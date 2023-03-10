use IO;
use FFTW;
use C_FFTW;
use CTypes;

use parameters;
use domains;

use compare_fortran;


/* Forward FFT plans */
var plan_f1 : fftw_plan;
var plan_f2 : fftw_plan;

/* Backward FFT plans */
var plan_b1 : fftw_plan;
var plan_b2 : fftw_plan;

/* Dummy arrays for setting up plans */
var arr_nx : [D_nx] real;
var arr_nx2p : [D_nx2p] complex;
var arr_ny : [D_ny] complex;

/* 2D versions for testing */
var plan_2df : fftw_plan;
var plan_2db : fftw_plan;
var arr_2D_in : [D] real;
var arr_2D_out : [D_hat] complex;



proc set_up_forward_FFTs() {

  // First forward FFT is along the x dimension, which takes an nx-vector as input and returns an nx2p-vector as output
    plan_f1 = fftw_plan_dft_r2c_1d(nx : c_int, c_ptrTo(arr_nx), c_ptrTo(arr_nx2p), FFTW_ESTIMATE);

  // Second forward FFT is along the y dimension, for which both input and output are length ny
    plan_f2 = fftw_plan_dft_1d(ny : c_int, c_ptrTo(arr_ny), c_ptrTo(arr_ny), -1, FFTW_ESTIMATE);

  // Forward 2d (for testing only)
    plan_2df = fftw_plan_dft_r2c_2d(ny : c_int, nx : c_int, c_ptrTo(arr_2D_in), c_ptrTo(arr_2D_out), FFTW_ESTIMATE);
}

proc set_up_backward_FFTs() {

  // First backward FFT is along the y dimension, for which both input and output are length ny
    plan_b1 = fftw_plan_dft_1d(ny : c_int, c_ptrTo(arr_ny), c_ptrTo(arr_ny), 1, FFTW_ESTIMATE);

  // Second backward FFT is along the x dimension, which takes an nx2p-vector as input and returns an nx-vector as output
    plan_b2 = fftw_plan_dft_c2r_1d(nx : c_int, c_ptrTo(arr_nx2p), c_ptrTo(arr_nx), FFTW_ESTIMATE);

  // Backward 2D (for testing only)
    plan_2db = fftw_plan_dft_c2r_2d(ny : c_int, nx : c_int, c_ptrTo(arr_2D_out), c_ptrTo(arr_2D_in), FFTW_ESTIMATE);

}

proc execute_forward_FFTs(ref in_arr: [] real, ref out_arr : [] complex) {

  /* Temporary arrays to hold 1D transforms */
    var tmp_f1 : [D3_hatT] complex;
    var tmp_f1T : [D3_hat] complex;

  /* First forward FFT */
    for (i,j) in D_zy {
      fftw_execute_dft_r2c(plan_f1, c_ptrTo(in_arr[i,j,1]), c_ptrTo(tmp_f1[i,j,1]));
    }

  /* Transpose */
    transpose_3D(tmp_f1, tmp_f1T);

  /* Second forward FFT */
    for (i,j) in D_zxhat {
      fftw_execute_dft(plan_f2, c_ptrTo(tmp_f1T[i,j,1]), c_ptrTo(out_arr[i,j,1]));
    }

}

proc execute_forward_FFTs_single_level(ref in_arr: [] real, ref out_arr : [] complex) {

  /* Temporary arrays to hold 1D transforms */
    var tmp_f1 : [D_hatT] complex;
    var tmp_f1T : [D_hat] complex;

  /* First forward FFT */
    for i in D_ny {
      fftw_execute_dft_r2c(plan_f1, c_ptrTo(in_arr[i,1]), c_ptrTo(tmp_f1[i,1]));
    }

  /* Transpose */
    transpose_2D(tmp_f1, tmp_f1T);

  /* Second forward FFT */
    for i in D_nx2p {
      fftw_execute_dft(plan_f2, c_ptrTo(tmp_f1T[i,1]), c_ptrTo(out_arr[i,1]));
    }

}

proc execute_backward_FFTs(ref in_arr: [] complex, ref out_arr : [] real) {

    /* Temporary arrays to hold 1D transforms */
      var tmp_b1 : [D3_hat] complex;
      var tmp_b1T : [D3_hatT] complex;

  /* First backward FFT */
    for (i,j) in D_zxhat {
      fftw_execute_dft(plan_b1, c_ptrTo(in_arr[i,j,1]), c_ptrTo(tmp_b1[i,j,1]));
    }

  /* Transpose */
    transpose_3D(tmp_b1, tmp_b1T);

  /* Second backward FFT */
    for (i,j) in D_zy {
      fftw_execute_dft_c2r(plan_b2, c_ptrTo(tmp_b1T[i,j,1]), c_ptrTo(out_arr[i,j,1]));
    }

}

proc execute_backward_FFTs_single_level(ref in_arr: [] complex, ref out_arr : [] real) {

    /* Temporary arrays to hold 1D transforms */
      var tmp_b1 : [D_hat] complex;
      var tmp_b1T : [D_hatT] complex;

  /* First backward FFT */
    for i in D_nx2p {
      fftw_execute_dft(plan_b1, c_ptrTo(in_arr[i,1]), c_ptrTo(tmp_b1[i,1]));
    }

  /* Transpose */
    transpose_2D(tmp_b1, tmp_b1T);

  /* Second backward FFT */
    for i in D_ny {
      fftw_execute_dft_c2r(plan_b2, c_ptrTo(tmp_b1T[i,1]), c_ptrTo(out_arr[i,1]));
    }

}

proc transpose_2D(ref in_arr: [?D] complex, ref out_arr: [] complex) {

    for (i,j) in D {
      out_arr[j,i] = in_arr[i,j];
    }
}

proc transpose_3D(ref in_arr: [?D] complex, ref out_arr: [] complex) {

    for (i,j,k) in D {
      out_arr[i,k,j] = in_arr[i,j,k];
    }
}

proc normalize(ref in_arr: [] real, ref out_arr: [] real) {

    out_arr = in_arr / (nx*ny);

}

proc execute_forward_FFTs_2D(ref in_arr: [] real, out out_arr : [] complex) {

    var tmp : [D_hatT] complex;
    var tmpT : [D_hat] complex;
    fftw_execute_dft_r2c(plan_2df, c_ptrTo(in_arr), c_ptrTo(tmp));

    /* Transpose */
    transpose_2D(tmp, tmpT);

    out_arr = tmpT;
}

proc execute_backward_FFTs_2D(ref in_arr: [] complex, out out_arr : [] real) {

    var out_tmp : [D] real;
    var tmp : [D_hatT] complex;

    /* Transpose */
    transpose_2D(in_arr, tmp);

    fftw_execute_dft_c2r(plan_2db, c_ptrTo(tmp), c_ptrTo(out_tmp));
    out_arr = out_tmp;

}
