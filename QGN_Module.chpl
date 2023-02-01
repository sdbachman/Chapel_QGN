use IO;
use FFTW;
use C_FFTW;
use CTypes;

use LAPACK;
use ClassicLAPACK;

use LinearAlgebra;
use Random;

use parameters;
use domains;
use arrays;
use FFT_utils;
use Stratification;
use Shear;

use compare_fortran;
use Time;

////////////////////////////////////////////////////////////////
//  Initialize: Set up the initial PV field, modes, and FFTs  //
////////////////////////////////////////////////////////////////

proc Initialize() {

create_initial_state(q);
//print_array_3D(q);

writeln("---------------------------------------------------------------");
writeln("dx = ", dx,"                dy = ", dy);

/* Horizontal locations of grid points (physical space) */
    forall (i,j) in D {
        x[i,j] = -(Lx+dx)/2 + dx*j;
    }

    forall (i,j) in D {
        y[i,j] = -(Ly+dy)/2 + dy*i;
    }

/* Initialize vertical modes */

    /* Set elements of L ([Del+L]psi=q) */
      L[1,2] = 2*S[1]/(H[1]*(H[1]+H[2]));
      L[1,1] =-L[1,2];
      if (nz > 2) {
        for k in 2..(nz-1) {
          L[k,k-1] = 2*S[k-1]/(H[k]*(H[k]+H[k-1]));
          L[k,k+1] = 2*S[k]/(H[k]*(H[k]+H[k+1]));
          L[k,k]   =-(L[k,k-1]+L[k,k+1]);
        }
      }
      L[nz,nz-1] = 2*S[nz-1]/(H[nz]*(H[nz]+H[nz-1]));
      L[nz,nz]   = -L[nz,nz-1];

    /* At this point need to set up mean zonal velocity and associated mean PV
    gradient before L gets over-written. */
      var tmp = -uBar;
      qyBar = dot(L,tmp);

    /* Eigenvalue decomposition; DGEEV over-writes L */
      LAPACKE_dgeev(lapack_memory_order.row_major,'N','V',nz : c_int, L, nz : c_int,
                    EVals, EValsI, ModesL, nz : c_int, Modes,nz : c_int);

    /* Now normalize so that depth average of Modes[..,k]**2 is 1. */
      for k in 1..nz {
        tmp = + reduce ((H/Htot)*Modes[..,k]**2);
        Modes[..,k] = Modes[..,k] / sqrt(tmp);
      }

    writeln("---------------------------------------------------------------");
    writeln(" The deformation radii are ");
    writeln( 1.0/sqrt(-EVals) );
    writeln("---------------------------------------------------------------");

    /* Initialize wavenumbers. We are going to do all spectral operations on the
       transposed domain, so kx will be constant in the first dimension and ky
       will be constant in the second dimension. */

      for i in 1..nx2p {
        kx[i,..] = (2*pi/Lx) * (i-1);
      }
      for j in 1..ny2p {
        ky[..,j] = (2*pi/Ly) * (j-1);
      }
      for j in (ny2p+1)..ny {
        ky[..,j] = (2*pi/Ly) * (j-ny-1);
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

      writeln(" Initialization complete ");
      writeln("---------------------------------------------------------------");


}


////////////////////////////////////////////////////////////////
//          DeAlias: zeros out top 1/3 of coefficients        //
////////////////////////////////////////////////////////////////

proc DeAlias(ref field : [] complex) {

    forall (i,j,k) in D3_hat_sp1 {
      field[i,j,k] = 0;
    }
    forall (i,j,k) in D3_hat_sp2 {
      field[i,j,k] = 0;
    }
    forall i in zl {
      field[i,1,1] = 0;
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

proc Jacobian(ref q_in : [] complex, ref jaco_hat : [] complex) {

  var uq_hat : [D3_hat] complex;
  var vq_hat : [D3_hat] complex;

  var uq_phys : [D3] real;
  var vq_phys : [D3] real;

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

proc GetRHS(ref q_in : [] complex, ref RHS : [] complex) {

  var jaco_hat : [D3_hat] complex;
  var drag_tmp : [D_hat] complex;
  var drag_hat : [D_hat] complex;
  var Uu_drag : [D] real;
  var Uv_drag : [D] real;

  /* Advection */
    Jacobian(q_in,RHS);

  /* Mean advection, beta and viscosity */
    forall (i,j,k) in D3_hat {
      RHS[i,j,k] = RHS[i,j,k] - uBar[i]*1i*kx[j,k]*q_in[i,j,k]
                 - (beta + qyBar[i])*v_hat[i,j,k] - A8*(k2[j,k]**4)*q_in[i,j,k];
    }

  /* Ekman */
    forall (j,k) in D_hat {
      RHS[nz,j,k] = RHS[nz,j,k] + (r0*Htot/H[nz]) * k2[j,k] * psi_hat[nz,j,k];
    }

  /* Quadratic drag */
    forall (j,k) in D {
      Uu_drag[j,k] = sqrt(u_phys[nz,j,k]**2+v_phys[nz,j,k]**2)*u_phys[nz,j,k];
      Uv_drag[j,k] = sqrt(u_phys[nz,j,k]**2+v_phys[nz,j,k]**2)*v_phys[nz,j,k];
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
      RHS[nz,j,k] = RHS[nz,j,k] + (C_d*Htot/H[nz])*drag_hat[j,k];
    }

  /* Dealias */
    DeAlias(RHS);

}


////////////////////////////////////////////////////////////////
//              GetPsi: Gets psi_hat from q_hat               //
////////////////////////////////////////////////////////////////

proc GetPsi(ref in_arr : [] complex) {

  var q_hat_mode : [D3_hat] complex;
  var psi_hat_mode : [D3_hat] complex;

  var t : [D_hat] complex;

  q_hat_mode = 0;
  /* Get q_hat_mode and psi_hat_mode */
    forall (i,j,k) in D3_hat {
        for ii in 1..nz {
          q_hat_mode[i,j,k] = q_hat_mode[i,j,k] + H[ii]*Modes[ii,i]*in_arr[ii,j,k];
        }
        q_hat_mode[i,j,k] = q_hat_mode[i,j,k]/Htot;
        psi_hat_mode[i,j,k] = q_hat_mode[i,j,k]/(-k2[j,k]+EVals[i]);
        psi_hat_mode[i,1,1] = 0;
    }

  psi_hat = 0;
  /* Get psi_hat */
    forall (i,j,k) in D3_hat {
        for ii in 1..nz {
          psi_hat[i,j,k] = psi_hat[i,j,k] + psi_hat_mode[ii,j,k]*Modes[i,ii];
        }
    }

  /* Get b_hat = f0*dpsi/dz
    do k=1,nz-1
      b_hat(:,:,k) = (2._dp*f0/(H(k)+H(k+1)))*(psi_hat(:,:,k)-psi_hat(:,:,k+1))
    end do
  */

}


////////////////////////////////////////////////////////////////
//   create_initial_state: Fill initial q with random values  //
////////////////////////////////////////////////////////////////

proc create_initial_state(ref in_arr : [?dom] real) {

  var seed=17;
  fillRandom(in_arr, seed);

  for k in 1..nz {
    var tmp = in_arr[k,..,..];
    in_arr[k,..,..] = 1e-6*(tmp - (+ reduce tmp)/(nx*ny));
  }

}
