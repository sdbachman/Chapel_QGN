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
use IO_Module;

use AllLocalesBarriers;
use compare_fortran;
use Time;

////////////////////////////////////////////////////////////////
//  Initialize: Set up the initial PV field, modes, and FFTs  //
////////////////////////////////////////////////////////////////

proc Initialize() {

/*
  if restart {
    read_initial_state(q);
  }
  else {
    create_initial_state(q);
  }
*/

load_fortran_grid("test_grid", q);

//print_array_3D(q);

writeln("---------------------------------------------------------------");
writeln("dx = ", dx,"                dy = ", dy);

  /* Horizontal locations of grid points (physical space) */
    forall (i,j) in _D.localSubdomain() {
        x[i,j] = -(Lx+dx)/2 + dx*(j+1);
    }

    forall (i,j) in _D.localSubdomain() {
        y[i,j] = -(Ly+dy)/2 + dy*(i+1);
    }

    //writeln("x on ", _D.localSubdomain(), ":");
    //print_array_2D(x);

    //writeln("y on ", _D.localSubdomain(), ":");
    //print_array_2D(y);

/* Initialize vertical modes */

    /* Set elements of L ([Del+L]psi=q) */
      L[0,1] = 2*S[0]/(H[0]*(H[0]+H[1]));
      L[0,0] =-L[0,1];
      if (nz > 2) {
        for k in 1..(nz-2) {
          L[k,k-1] = 2*S[k-1]/(H[k]*(H[k]+H[k-1]));
          L[k,k+1] = 2*S[k]/(H[k]*(H[k]+H[k+1]));
          L[k,k]   =-(L[k,k-1]+L[k,k+1]);
        }
      }
      L[nz-1,nz-2] = 2*S[nz-2]/(H[nz-1]*(H[nz-1]+H[nz-2]));
      L[nz-1,nz-1]   = -L[nz-1,nz-2];

      //writeln("L on ", _D.localSubdomain(), ":");
      //writeln(L);

    /* At this point need to set up mean zonal velocity and associated mean PV
    gradient before L gets over-written. */
    /* Cannot pass replicated arrays into LAPACK or BLAS functions, so I need
    to make local copies for now. */
      var tmpuBar : [zl] real = -uBar;
      var tmpL : [zl2] real = L;

      qyBar = dot(tmpL,tmpuBar);

    /* Eigenvalue decomposition; DGEEV over-writes L */
    /* Cannot pass replicated arrays into LAPACK or BLAS functions, so I need
    to make local copies for now. */
      var tmpEVals : [zl] real;
      var tmpEValsI : [zl] real;
      var tmpModesL : [zl2] real;
      var tmpModes : [zl2] real;
      LAPACKE_dgeev(lapack_memory_order.row_major,'N','V',nz : c_int, tmpL, nz : c_int,
                    tmpEVals, tmpEValsI, tmpModesL, nz : c_int, tmpModes,nz : c_int);

      EVals = tmpEVals;
      EValsI = tmpEValsI;
      ModesL = tmpModesL;
      Modes = tmpModes;

    /* Now normalize so that depth average of Modes[..,k]**2 is 1. */
      for k in 0..#nz {
        tmp = + reduce ((H/Htot)*Modes[..,k]**2);
        Modes[..,k] = Modes[..,k] / sqrt(tmp);
      }

      //writeln("EVals on ", _D.localSubdomain(), ":");
      //writeln(EVals);
      //writeln("Modes on ", _D.localSubdomain(), ":");
      //writeln(Modes);

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
        kx[i,..] = (2*pi/Lx) * i;
      }
      for j in 0..#ny2p {
        ky[..,j] = (2*pi/Ly) * j;
      }
      for j in ny2p..(ny-1) {
        ky[..,j] = (2*pi/Ly) * (j-ny);
      }

      k2 = kx**2+ky**2;

      //if (here.id == 0) {
      //writeln("k2 on ", _D_hat.localSubdomain(), ":");
      //print_array_2D(k2);
      //}

    /* Set up Fourier Transforms */
      set_up_forward_FFTs();
      set_up_backward_FFTs();

    /* Transform the initial condition */
      writeln(" Transforming the initial condition ");
      writeln("---------------------------------------------------------------");

      execute_forward_FFTs(q, q_hat);

      allLocalesBarrier.barrier();
      DeAlias(q_hat);

      writeln(" Initialization complete ");
      writeln("---------------------------------------------------------------");

}


////////////////////////////////////////////////////////////////
//          DeAlias: zeros out top 1/3 of coefficients        //
////////////////////////////////////////////////////////////////

proc DeAlias(ref field : [] complex) {

    forall (i,j,k) in _D3_hat_sp1 {
      field[i,j,k] = 0;
    }
    forall (i,j,k) in _D3_hat_sp2 {
      field[i,j,k] = 0;
    }
    forall i in _zl {
      field[i,0,0] = 0;
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

  //var uq_hat : [_D3_hat] complex;
  //var vq_hat : [_D3_hat] complex;

  //var uq_phys : [_D3] real;
  //var vq_phys : [_D3] real;

  /* Get psi_hat, u_hat, v_hat */
    GetPsi(q_in);

    forall (i,j,k) in _D3_hat.localSubdomain() {
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

    forall (i,j,k) in _D3.localSubdomain() {
      uq_phys[i,j,k] = u_phys[i,j,k]*q_phys[i,j,k];
      vq_phys[i,j,k] = v_phys[i,j,k]*q_phys[i,j,k];
    }

  /* Get uq_hat, vq_hat */
    execute_forward_FFTs(uq_phys, uq_hat);
    execute_forward_FFTs(vq_phys, vq_hat);

  /* Compute jacobian_spec */
  /* The RHS term is the negative of the Jacobian, so I will put a
     minus sign here to avoid a needless array copy in the calling function. */
    forall (i,j,k) in _D3_hat.localSubdomain() {
      jaco_hat[i,j,k] = -1i*(kx[j,k]*uq_hat[i,j,k]+ky[j,k]*vq_hat[i,j,k]);
    }

}


////////////////////////////////////////////////////////////////
//                  GetRHS: computes FFT of                   //
//     -J[psi,q]-uBar*qx-(beta+qyBar)*v-nu*del[q] - Ekman     //
////////////////////////////////////////////////////////////////

proc GetRHS(ref q_in : [] complex, ref RHS : [] complex) {

  //var drag_tmp : [_D_hat] complex;
  //var drag_hat : [_D_hat] complex;
  var Uu_drag : [_D] real;
  var Uv_drag : [_D] real;

  /* Advection */
    Jacobian(q_in,RHS);

  /* Mean advection, beta and viscosity */
    forall (i,j,k) in _D3_hat.localSubdomain() {
      RHS[i,j,k] = RHS[i,j,k] - uBar[i]*1i*kx[j,k]*q_in[i,j,k]
                 - (beta + qyBar[i])*v_hat[i,j,k] - A8*(k2[j,k]**4)*q_in[i,j,k];
    }

  /* Ekman */
    forall (j,k) in _D_hat.localSubdomain() {
      RHS[nz1m,j,k] = RHS[nz1m,j,k] + (r0*Htot/H[nz1m]) * k2[j,k] * psi_hat[nz1m,j,k];
    }

  /* Quadratic drag */
    forall (j,k) in _D.localSubdomain() {
      Uu_drag[j,k] = sqrt(u_phys[nz1m,j,k]**2+v_phys[nz1m,j,k]**2)*u_phys[nz1m,j,k];
      Uv_drag[j,k] = sqrt(u_phys[nz1m,j,k]**2+v_phys[nz1m,j,k]**2)*v_phys[nz1m,j,k];
    }

    execute_forward_FFTs_single_level(Uu_drag, drag_tmp);

    forall (j,k) in _D_hat.localSubdomain() {
      drag_hat[j,k] = 1i*ky[j,k]*drag_tmp[j,k];
    }

    execute_forward_FFTs_single_level(Uv_drag, drag_tmp);
    forall (j,k) in _D_hat.localSubdomain() {
      drag_hat[j,k] = drag_hat[j,k] - 1i*kx[j,k]*drag_tmp[j,k];
    }

    forall (j,k) in _D_hat.localSubdomain() {
      RHS[nz1m,j,k] = RHS[nz1m,j,k] + (C_d*Htot/H[nz1m])*drag_hat[j,k];
    }

  /* Dealias */
    allLocalesBarrier.barrier();
    DeAlias(RHS);

}


////////////////////////////////////////////////////////////////
//              GetPsi: Gets psi_hat from q_hat               //
////////////////////////////////////////////////////////////////

proc GetPsi(ref in_arr : [] complex) {

  var q_hat_mode : [_D3_hat] complex;
  var psi_hat_mode : [_D3_hat] complex;

  q_hat_mode.localSlice(_D3_hat.localSubdomain()) = 0;
  /* Get q_hat_mode and psi_hat_mode */
    forall (i,j,k) in _D3_hat.localSubdomain() {
        for ii in 0..#nz {
          q_hat_mode[i,j,k] = q_hat_mode[i,j,k] + H[ii]*Modes[ii,i]*in_arr[ii,j,k];
        }
        q_hat_mode[i,j,k] = q_hat_mode[i,j,k]/Htot;
        psi_hat_mode[i,j,k] = q_hat_mode[i,j,k]/(-k2[j,k]+EVals[i]);
        psi_hat_mode[i,0,0] = 0;
    }

  psi_hat.localSlice(_D3_hat.localSubdomain()) = 0;
  /* Get psi_hat */
    forall (i,j,k) in _D3_hat.localSubdomain() {
        for ii in 0..#nz {
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

  var seed=17+(here.id : int);
  var D = dom.localSubdomain();
  var tmp = in_arr.localSlice(D);
  fillRandom(tmp, seed);

  tmp = 1e-6*(tmp - (+ reduce tmp)/(D.size));
  in_arr[D] = tmp;
}
