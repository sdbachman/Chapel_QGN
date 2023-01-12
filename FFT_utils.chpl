use IO;
use FFTW;
use C_FFTW;
use CTypes;

use parameters;
use domains;

/* Forward FFT plans */
var plan_f1 : fftw_plan;
var plan_f2 : fftw_plan;

/* Backward FFT plans */
var plan_b1 : fftw_plan;
var plan_b2 : fftw_plan;

/* Dummy arrays for setting up plans */
var arr_ny : [D_ny] real;
var arr_nyp : [D_nyp] complex;
var arr_nx : [D_nx] complex;



proc set_up_forward_FFTs() {

  // First forward FFT is along the y dimension, which takes an ny-vector as input and returns an nyp-vector as output
    plan_f1 = fftw_plan_dft_r2c_1d(ny : c_int, c_ptrTo(arr_ny), c_ptrTo(arr_nyp), FFTW_ESTIMATE);

  // Second forward FFT is along the x dimension, for which both input and output are length nx
    plan_f2 = fftw_plan_dft_1d(nx : c_int, c_ptrTo(arr_nx), c_ptrTo(arr_nx), -1, FFTW_ESTIMATE);

}

proc set_up_backward_FFTs() {

  // First backward FFT is along the x dimension, for which both input and output are length nx
    plan_b1 = fftw_plan_dft_1d(nx : c_int, c_ptrTo(arr_nx), c_ptrTo(arr_nx), 1, FFTW_ESTIMATE);

  // Second backward FFT is along the y dimension, which takes an nyp-vector as input and returns an ny-vector as output
    plan_b2 = fftw_plan_dft_c2r_1d(ny : c_int, c_ptrTo(arr_nyp), c_ptrTo(arr_ny), FFTW_ESTIMATE);

}

proc execute_forward_FFTs(ref in_arr: [] real, ref out_arr : [] complex) {

  /* Temporary arrays to hold 1D transforms */
    var tmp_f1 : [D_out] complex;
    var tmp_f1T : [D_outT] complex;

  /* First forward FFT */
    for i in 1..nx {
      fftw_execute_dft_r2c(plan_f1, c_ptrTo(in_arr[i,..]), c_ptrTo(tmp_f1[i,..]));
    }

  /* Transpose */
    transpose(tmp_f1, tmp_f1T);

  /* Second forward FFT */
    for j in 1..nyp {
      fftw_execute_dft(plan_f2, c_ptrTo(tmp_f1T[j,..]), c_ptrTo(out_arr[j,..]));
    }

}

proc execute_backward_FFTs(ref in_arr: [] complex, ref out_arr : [] real) {

    /* Temporary arrays to hold 1D transforms */
      var tmp_b1 : [D_out] complex;
      var tmp_b1T : [D_outT] complex;

  /* First backward FFT */
    for j in 1..nyp {
      fftw_execute_dft(plan_b1, c_ptrTo(in_arr[j,..]), c_ptrTo(tmp_b1T[j,..]));
    }

  /* Transpose */
    transpose(tmp_b1T, tmp_b1);

  /* Second backward FFT */
    for i in 1..nx {
      fftw_execute_dft_c2r(plan_b2, c_ptrTo(tmp_b1[i,..]), c_ptrTo(out_arr[i,..]));
    }

}

proc transpose(ref in_arr: [?D] complex, ref out_arr: [] complex) {

    for (i,j) in D {
      out_arr[j,i] = in_arr[i,j];
    }

}

proc normalize(ref in_arr: [] real, ref out_arr: [] real) {

    out_arr = in_arr / (nx*ny);

}

