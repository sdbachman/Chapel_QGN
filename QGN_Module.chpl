use IO;
use FFTW;
use C_FFTW;
use CTypes;

use LAPACK;
use ClassicLAPACK;

use LinearAlgebra;

use parameters;
use domains;
use arrays;
use FFT_utils;
use Stratification;
use Shear;

use compare_fortran;

////////////////////////////////////////////////////////////////
//  Initialize: Set up the initial PV field, modes, and FFTs  //
////////////////////////////////////////////////////////////////

proc Initialize() {

load_fortran_grid(q);
//load_1layer_test(q);
//print_array_3D(q);


writeln("---------------------------------------------------------------");
writeln("dx = ", dx,"                dy = ", dy);

/* Horizontal locations of grid points (physical space) */
    for j in 1..nx {
        x[..,j] = -(Lx+dx)/2 + dx*j;
    }

    for i in 1..ny {
        y[i,..] = -(Ly+dy)/2 + dy*i;
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

    for j in nx3p..nx2p {
      field[..,j,..] = 0;
    }
    for k in ny3p..ny3p2 {
      field[..,..,k] = 0;
    }
    field[..,1,1] = 0;

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

    for k in 1..nz {
      u_hat[k,..,..] =-1i*ky*psi_hat[k,..,..];
      v_hat[k,..,..] = 1i*kx*psi_hat[k,..,..];
    }
    //forall (i,j,k) in D_hat {
    //  u_hat[i,j,k] = -1i*ky[j,k]*psi_hat[i,j,k];
    //  v_hat[i,j,k] = 1i*kx[j,k]*psi_hat[i,j,k];
    //}

  /* Get u, v, q */
    execute_backward_FFTs(q_in, q_phys);
    normalize(q_phys, q_phys);

    execute_backward_FFTs(u_hat, u_phys);
    normalize(u_phys, u_phys);

    execute_backward_FFTs(v_hat, v_phys);
    normalize(v_phys, v_phys);

    uq_phys = u_phys*q_phys;
    vq_phys = v_phys*q_phys;

   // q_in matches
   // u_hat and v_hat match
   // q_phys matches initial grid when DeAlias is off
   // v_phys matches with Fortran when DeAlias is on, but u and q don't?

  /* Get uq_hat, vq_hat */
    execute_forward_FFTs(uq_phys, uq_hat);
    execute_forward_FFTs(vq_phys, vq_hat);

  /* Compute jacobian_spec */
    for k in 1..nz {
      jaco_hat[k,..,..] = 1i*(kx*uq_hat[k,..,..]+ky*vq_hat[k,..,..]);
    }

}


////////////////////////////////////////////////////////////////
//                  GetRHS: computes FFT of                   //
//     -J[psi,q]-uBar*qx-(beta+qyBar)*v-nu*del[q] - Ekman     //
////////////////////////////////////////////////////////////////

proc GetRHS(ref q_in : [] complex, ref RHS : [] complex) {

  var jaco_hat : [D3_hat] complex;
  var drag_tmp : [D3_hat] complex;
  var drag_hat : [D_hat] complex;
  var URMS : [D3] real;
  var Uu_drag : [D3] real;
  var Uv_drag : [D3] real;

  /* Advection */
    Jacobian(q_in,jaco_hat);
    RHS = -jaco_hat;

  /* Mean advection, beta and viscosity */
    for k in 1..nz {
      RHS[k,..,..] = RHS[k,..,..] - uBar[k]*1i*kx*q_in[k,..,..]
                 - (beta + qyBar[k])*v_hat[k,..,..] - A8*(k2**4)*q_in[k,..,..];
    }

  /* Ekman */
    RHS[nz,..,..] = RHS[nz,..,..] + (r0*Htot/H[nz]) * k2 * psi_hat[nz,..,..];

  /* Quadratic drag */
    URMS = sqrt(u_phys**2+v_phys**2);
    Uu_drag = URMS*u_phys;
    Uv_drag = URMS*v_phys;

    execute_forward_FFTs(Uu_drag, drag_tmp);
    drag_hat = 1i*ky*drag_tmp[nz,..,..];

    execute_forward_FFTs(Uv_drag, drag_tmp);
    drag_hat = drag_hat - 1i*kx*drag_tmp[nz,..,..];

    RHS[nz,..,..] = RHS[nz,..,..] + (C_d*Htot/H[nz])*drag_hat;

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
    for k in 1..nz {
      for kk in 1..nz {
        q_hat_mode[k,..,..] = q_hat_mode[k,..,..] + H[kk]*Modes[kk,k]*in_arr[kk,..,..];
      }
      q_hat_mode[k,..,..] = q_hat_mode[k,..,..]/Htot;
      psi_hat_mode[k,..,..] = q_hat_mode[k,..,..]/(-k2+EVals[k]);
      psi_hat_mode[k,1,1] = 0;
    }

  psi_hat = 0;
  /* Get psi_hat */
    for k in 1..nz {
      for kk in 1..nz {
        psi_hat[k,..,..] = psi_hat[k,..,..] + psi_hat_mode[kk,..,..]*Modes[k,kk];
      }
    }

  /* Get b_hat = f0*dpsi/dz
    do k=1,nz-1
      b_hat(:,:,k) = (2._dp*f0/(H(k)+H(k+1)))*(psi_hat(:,:,k)-psi_hat(:,:,k+1))
    end do
  */

}
