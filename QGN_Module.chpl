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

//load_fortran_grid();
load_1layer_test(q);
print_array_3D(q);


writeln("---------------------------------------------------------------");
writeln("dx = ", dx,"                dy = ", dy);

/* Horizontal locations of grid points (physical space) */
    for j in 1..ny {
        y[..,j] = -(Ly+dy)/2 + dy*j;
    }
    for i in 1..nx {
        x[i,..] = -(Lx+dx)/2 + dx*i;
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

      for i in 1..ny2p {
        ky[i,..] = (2*pi/Ly) * (i-1);
      }
      for j in 1..nx2p {
        kx[..,j] = (2*pi/Lx) * (j-1);
      }
      for j in (nx2p+1)..nx {
        kx[..,j] = (2*pi/Lx) * (j-nx-1);
      }

      k2 = kx**2+ky**2;

    /* Set up Fourier Transforms */
      set_up_forward_FFTs();
      set_up_backward_FFTs();

    /* Transform the initial condition */
      writeln(" Transforming the initial condition ");
      execute_forward_FFTs(q, q_hat);

      DeAlias(q_hat);


//print_array_3D(q_hat);


}


////////////////////////////////////////////////////////////////
//          DeAlias: zeros out top 1/3 of coefficients        //
////////////////////////////////////////////////////////////////

proc DeAlias(ref field : [] complex) {

    for j in ny3p..ny2p {
      field[..,j,..] = 0;
    }
    for k in nx3p..nx3p2 {
      field[..,..,k] = 0;
    }
    field[..,1,1] = 0;

}

