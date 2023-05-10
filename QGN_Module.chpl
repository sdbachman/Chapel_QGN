use IO;
use FFTW;
use C_FFTW;
use CTypes;

use Math;
use LinearAlgebra;
use Random;

use parameters;
use domains;
use arrays;
use FFT_utils;
use ARK43;
use IO_Module;

use compare_fortran;
use Time;

////////////////////////////////////////////////////////////////
//  Initialize: Set up the initial PV field, modes, and FFTs  //
////////////////////////////////////////////////////////////////

proc Initialize() {

  if restart {
    read_initial_state(q);
  }
  else {
    create_initial_state(q);
  }

  read_background_state(background_file);

  writeln("---------------------------------------------------------------");
  writeln("dx = ", dx,"                dy = ", dy);

  /* Horizontal locations of grid points (physical space) */
    forall (i,j) in D {
        x[i,j] = -(Lx+dx)/2 + dx*(j+1);
    }

    forall (i,j) in D {
        y[i,j] = -(Ly+dy)/2 + dy*(i+1);
    }

    /* THIS IS NOT NEEDED SINCE WE DO THIS WHEN MAKING THE INPUT FILE
    /* Normalize so that the depth average of Modes[..,k]**2 is 1. */
      for k in 0..#nz {
        var tmp = + reduce ((H/Htot)*Modes[..,k]**2);
        Modes[..,k] = Modes[..,k] / sqrt(tmp);
      }
    */

      writeln("---------------------------------------------------------------");
      writeln(" The deformation radii are ");
      writeln( 1.0/sqrt(-EVals) );
      writeln("---------------------------------------------------------------");

    /* Initialize wavenumbers. We are going to do all spectral operations on the
       transposed domain, so kx will be constant in the first dimension and ky
       will be constant in the second dimension. */

      // These loops will only affect parts of kx, ky that
      // are in the local subdomain
      for i in 0..#nx2p {
        kx[..,i] = (2*pi/Lx) * i;
      }
      for j in 0..#ny2p {
        ky[j,..] = (2*pi/Ly) * j;
      }
      for j in ny2p..(ny-1) {
        ky[j,..] = (2*pi/Ly) * (j-ny);
      }

      k2 = kx**2+ky**2;

    /* Set up Fourier Transforms */
      set_up_forward_FFTs();
      set_up_backward_FFTs();

    /* Transform the initial condition */
      writeln(" Transforming the initial condition ");
      writeln("---------------------------------------------------------------");

      execute_forward_FFTs(q, q_hat);

      DeAlias(q_hat);

    /* Set the arrays for the ARK43 timestepping scheme on each Locale */
      writeln(" Setting ARK43 arrays ");
      writeln("---------------------------------------------------------------");
      set_ARK43_vars();


      writeln(" Initialization complete ");
      writeln("---------------------------------------------------------------");

}


////////////////////////////////////////////////////////////////
//          DeAlias: zeros out top 1/3 of coefficients        //
////////////////////////////////////////////////////////////////

proc DeAlias(ref field : [] complex(cp)) {

    forall (i,j,k) in D3_hat_sp1 {
      field[i,j,k] = 0;
    }
    forall (i,j,k) in D3_hat_sp2 {
      field[i,j,k] = 0;
    }
    forall (i,j,k) in D3_hat[zl.dim(0),0..0,0..0] {
      field[i,j,k] = 0;
    }
}

////////////////////////////////////////////////////////////////
//                          Jacobian                          //
//         Computes 2/3-rule dealiased jacobian, returns      //
//           Fourier coefficients thereof in jaco_spec        //
//       Also computes spectral coefficients of v for beta    //
//                  and of psi for viscosity                  //
//                   Assumes dealiased input                  //
////////////////////////////////////////////////////////////////

proc Jacobian(ref q_in : [] complex(cp), ref jaco_hat : [] complex(cp)) {

  /* Get psi_hat, u_hat, v_hat */
    GetPsi(q_in);

    forall (i,j,k) in D3_hat {
      u_hat[i,j,k] = -1i*ky[j,k]*psi_hat[i,j,k];
      v_hat[i,j,k] = 1i*kx[j,k]*psi_hat[i,j,k];
    }

  /* Get u, v, q */
    execute_backward_FFTs(q_in, q_phys);
    normalize(q_phys);

    execute_backward_FFTs(u_hat, u_phys);
    normalize(u_phys);

    execute_backward_FFTs(v_hat, v_phys);
    normalize(v_phys);

    forall (i,j,k) in D3 {
      uq_phys[i,j,k] = u_phys[i,j,k]*q_phys[i,j,k];
      vq_phys[i,j,k] = v_phys[i,j,k]*q_phys[i,j,k];
    }

  /* Get uq_hat, vq_hat */
    execute_forward_FFTs(uq_phys, uq_hat);
    execute_forward_FFTs(vq_phys, vq_hat);

  /* Compute jacobian_spec */
  /* The RHS term is the negative of the Jacobian, so I will put a
     minus sign here to avoid a needless array copy in the calling function. */
    forall (i,j,k) in D3_hat {
      jaco_hat[i,j,k] = -1i*(kx[j,k]*uq_hat[i,j,k]+ky[j,k]*vq_hat[i,j,k]);
    }

}


////////////////////////////////////////////////////////////////
//                  GetRHS: computes FFT of                   //
//     -J[psi,q]-uBar*qx-(beta+qyBar)*v-nu*del[q] - Ekman     //
////////////////////////////////////////////////////////////////

proc GetRHS(ref q_in : [] complex(cp), ref RHS : [] complex(cp)) {

  /* Advection */
  Jacobian(q_in,RHS);

  /* Mean advection, beta and viscosity */
    forall (i,j,k) in D3_hat {
      RHS[i,j,k] = RHS[i,j,k] - uBar[i]*1i*kx[j,k]*q_in[i,j,k]
                 - (beta + qyBar[i])*v_hat[i,j,k] - A8*(k2[j,k]**4)*q_in[i,j,k];
    }

  /* Ekman */
    forall (j,k) in D_hat {
      RHS[nz1m,j,k] = RHS[nz1m,j,k] + (r0*Htot/H[nz1m]) * k2[j,k] * psi_hat[nz1m,j,k];
    }

  /* Quadratic drag */
    forall (j,k) in D {
      Uu_drag[j,k] = sqrt(u_phys[nz1m,j,k]**2+v_phys[nz1m,j,k]**2)*u_phys[nz1m,j,k];
      Uv_drag[j,k] = sqrt(u_phys[nz1m,j,k]**2+v_phys[nz1m,j,k]**2)*v_phys[nz1m,j,k];
    }

    execute_forward_FFTs_single_level(Uu_drag, drag_tmp);
    forall (j,k) in D_hat {
      drag_hat[j,k] = 1i*ky[j,k]*drag_tmp[j,k];
    }

    execute_forward_FFTs_single_level(Uv_drag, drag_tmp);
    forall (j,k) in D_hat {
      drag_hat[j,k] = drag_hat[j,k] - 1i*kx[j,k]*drag_tmp[j,k];
    }

    forall (j,k) in D_hat {
      RHS[nz1m,j,k] = RHS[nz1m,j,k] + (C_d*Htot/H[nz1m])*drag_hat[j,k];
    }

  /* Dealias */
    DeAlias(RHS);

}


////////////////////////////////////////////////////////////////
//              GetPsi: Gets psi_hat from q_hat               //
////////////////////////////////////////////////////////////////

proc GetPsi(ref in_arr : [] complex(cp)) {

  q_hat_mode = 0;

  /* Get q_hat_mode and psi_hat_mode */
  forall (i,j,k) in D3_hat {
    for ii in 0..#nz {
      q_hat_mode[i,j,k] = q_hat_mode[i,j,k] + H[ii]*Modes[ii,i]*in_arr[ii,j,k];
    }
    q_hat_mode[i,j,k] = q_hat_mode[i,j,k]/Htot;
    psi_hat_mode[i,j,k] = q_hat_mode[i,j,k]/(-k2[j,k]+EVals[i]);
  }

  forall (i,j,k) in D3_hat[0..#nz, 0..0, 0..0] {
    psi_hat_mode[i,0,0] = 0;
  }

  psi_hat = 0;
  /* Get psi_hat */
    forall (i,j,k) in D3_hat {
        for ii in 0..#nz {
          psi_hat[i,j,k] = psi_hat[i,j,k] + psi_hat_mode[ii,j,k]*Modes[i,ii];
        }
    }

  /* Get b_hat = f0*dpsi/dz
    do k=1,nz-1
      b_hat(:,:,k) = (2._dp*f0/(H(k)+H(k+1)))*(psi_hat(:,:,k)-psi_hat(:,:,k+1))
    end do
  */
  //}
}


////////////////////////////////////////////////////////////////
//   create_initial_state: Fill initial q with random values  //
////////////////////////////////////////////////////////////////

proc create_initial_state(ref in_arr : [?dom] real(rp)) {

  /*
  var seed=17+(here.id : int);
  var D = dom.localSubdomain();
  var tmp = in_arr.localSlice(D);
  fillRandom(tmp, seed);

  tmp = 1e-6*(tmp - (+ reduce tmp)/(D.size));
  for i in 1..nx {
    tmp[..,..,i..i] = tmp[..,..,i..i] + 1e-5*sin(16.0*pi*(i-1) / nx);
  }
  */

  var D = dom.localSubdomain();
  var tmp = in_arr.localSlice(D);

  for (i,j,k) in D {
    tmp[i,j,k] = 1e-5*sin(48.0*pi*(k-1) / nx);
    tmp[i,j,k] = tmp[i,j,k] + 2e-5*cos(23.0*pi*(j-1) / ny);
  }

  in_arr[D] = tmp;
}
