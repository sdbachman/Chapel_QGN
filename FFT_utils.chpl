use IO;
use FFTW;
use C_FFTW;
use CTypes;

use parameters;
use domains;

/* Forward FFT plans */
var plan_f1 : fftw_plan;
var plan_f2 : fftw_plan;

var plan_2d : fftw_plan;

/* Backward FFT plans */
var plan_b1 : fftw_plan;
var plan_b2 : fftw_plan;

/* Dummy arrays for setting up plans */
var arr_ny : [D_ny] real;
var arr_ny2p : [D_ny2p] complex;
var arr_nx : [D_nx] complex;

/* TO DO: REMOVE */
var arr_2D_in : [D] real;
var arr_2D_out : [D_hat] complex;

proc set_up_forward_FFTs() {

  // First forward FFT is along the y dimension, which takes an ny-vector as input and returns an ny2p-vector as output
    plan_f1 = fftw_plan_dft_r2c_1d(ny : c_int, c_ptrTo(arr_ny), c_ptrTo(arr_ny2p), FFTW_ESTIMATE);

  // Second forward FFT is along the x dimension, for which both input and output are length nx
    plan_f2 = fftw_plan_dft_1d(nx : c_int, c_ptrTo(arr_nx), c_ptrTo(arr_nx), -1, FFTW_ESTIMATE);

    /* TO DO: REMOVE */
    plan_2d = fftw_plan_dft_r2c_2d(nx : c_int, ny : c_int, c_ptrTo(arr_2D_in), c_ptrTo(arr_2D_out), FFTW_ESTIMATE);

}

proc set_up_backward_FFTs() {

  // First backward FFT is along the x dimension, for which both input and output are length nx
    plan_b1 = fftw_plan_dft_1d(nx : c_int, c_ptrTo(arr_nx), c_ptrTo(arr_nx), 1, FFTW_ESTIMATE);

  // Second backward FFT is along the y dimension, which takes an ny2p-vector as input and returns an ny-vector as output
    plan_b2 = fftw_plan_dft_c2r_1d(ny : c_int, c_ptrTo(arr_ny2p), c_ptrTo(arr_ny), FFTW_ESTIMATE);

}

/* TO DO: REMOVE */
proc execute_forward_FFT_2D(in in_arr: [] real, out out_arr : [] complex) {
    var tmp : [D_hat] complex;
    fftw_execute_dft_r2c(plan_2d, c_ptrTo(in_arr), c_ptrTo(tmp));
    out_arr = tmp;
}

proc execute_forward_FFTs(ref in_arr: [] real, ref out_arr : [] complex) {

  /* Temporary arrays to hold 1D transforms */
    var tmp_f1 : [D3_hatT] complex;
    var tmp_f1T : [D3_hat] complex;
    var tmp_vec : [D_ny] real;

  /* First forward FFT */
    // forall broken?
    for (i,j) in D_zx {
      tmp_vec = in_arr[i,j,..];
      writeln((i,j), ": ", tmp_vec);
      fftw_execute_dft_r2c(plan_f1, c_ptrTo(in_arr[i,j,1]), c_ptrTo(tmp_f1[i,j,1]));
    }

  /* Transpose */
    transpose_3D(tmp_f1, tmp_f1T);

  /* Second forward FFT */
    for (i,j) in D_zyhat {
      fftw_execute_dft(plan_f2, c_ptrTo(tmp_f1T[i,j,1]), c_ptrTo(out_arr[i,j,1]));
    }

}

proc execute_backward_FFTs2(ref in_arr: [] complex, ref out_arr : [] real) {

    /* Temporary arrays to hold 1D transforms */
      var tmp_b1 : [D3_hat] complex;
      var tmp_b1T : [D3_hatT] complex;

  /* First backward FFT */
    for (i,j) in D_zyhat {
      fftw_execute_dft(plan_b1, c_ptrTo(in_arr[i,j,1]), c_ptrTo(tmp_b1[i,j,1]));
    }

  /* Transpose */
    transpose_3D(tmp_b1, tmp_b1T);

  /* Second backward FFT */
    for (i,j) in D_zx {
      fftw_execute_dft_c2r(plan_b2, c_ptrTo(tmp_b1T[i,j,1]), c_ptrTo(out_arr[i,j,1]));
    }

}



proc execute_backward_FFTs(in in_arr: [] complex, out out_arr : [] real) {

    /* Temporary arrays to hold 1D transforms */
      var tmp_b1 : [D_hat] complex;
      var tmp_b1T : [D_hatT] complex;
      var tmp_out : [D] real;

  /* First backward FFT */
    for j in 1..ny2p {
      fftw_execute_dft(plan_b1, c_ptrTo(in_arr[j,..]), c_ptrTo(tmp_b1T[j,..]));
    }

  /* Transpose */
    transpose(tmp_b1T, tmp_b1);

  /* Second backward FFT */
    for i in 1..nx {
      fftw_execute_dft_c2r(plan_b2, c_ptrTo(tmp_b1[i,..]), c_ptrTo(tmp_out[i,..]));
    }

    out_arr = tmp_out;

}

proc transpose(ref in_arr: [?D] complex, ref out_arr: [] complex) {

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

