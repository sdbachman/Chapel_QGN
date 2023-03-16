use IO;
use FFTW;
use C_FFTW;
use CTypes;
use AllLocalesBarriers;
use CopyAggregation;

use parameters;
use domains;
use arrays;

use compare_fortran;


/* Forward FFT plans */
var plan_f1 : fftw_plan;
var plan_f2 : fftw_plan;

/* Backward FFT plans */
var plan_b1 : fftw_plan;
var plan_b2 : fftw_plan;

/* Dummy arrays for setting up plans */
//var arr_nx : [D_nx] real;
//var arr_nx2p : [D_nx2p] complex;
//var arr_ny : [D_ny] complex;


proc set_up_forward_FFTs() {

  var arr_nx : [D_nx] real;
  var arr_nx2p : [D_nx2p] complex;
  var arr_ny : [D_ny] complex;

  // First forward FFT is along the x dimension, which takes an nx-vector as input and returns an nx2p-vector as output
    plan_f1 = fftw_plan_dft_r2c_1d(nx : c_int, c_ptrTo(arr_nx), c_ptrTo(arr_nx2p), FFTW_ESTIMATE);

  // Second forward FFT is along the y dimension, for which both input and output are length ny
    plan_f2 = fftw_plan_dft_1d(ny : c_int, c_ptrTo(arr_ny), c_ptrTo(arr_ny), -1, FFTW_ESTIMATE);

  //  writeln("On ", here.id, ":");
  //  writeln("setting up FFTs");
}

proc set_up_backward_FFTs() {

  var arr_nx : [D_nx] real;
  var arr_nx2p : [D_nx2p] complex;
  var arr_ny : [D_ny] complex;

  // First backward FFT is along the y dimension, for which both input and output are length ny
    plan_b1 = fftw_plan_dft_1d(ny : c_int, c_ptrTo(arr_ny), c_ptrTo(arr_ny), 1, FFTW_ESTIMATE);

  // Second backward FFT is along the x dimension, which takes an nx2p-vector as input and returns an nx-vector as output
    plan_b2 = fftw_plan_dft_c2r_1d(nx : c_int, c_ptrTo(arr_nx2p), c_ptrTo(arr_nx), FFTW_ESTIMATE);

  //  writeln("On ", here.id, ":");
  //  writeln("setting up backwards");
}

proc execute_forward_FFTs(ref in_arr: [] real, ref out_arr : [] complex) {

  var arr_nx : [D_nx] real;
  var arr_nx2p : [D_nx2p] complex;
  var arr_ny : [D_ny] complex;

  // First forward FFT is along the x dimension, which takes an nx-vector as input and returns an nx2p-vector as output
    var plan_f1 = fftw_plan_dft_r2c_1d(nx : c_int, c_ptrTo(arr_nx), c_ptrTo(arr_nx2p), FFTW_ESTIMATE);

  // Second forward FFT is along the y dimension, for which both input and output are length ny
    var plan_f2 = fftw_plan_dft_1d(ny : c_int, c_ptrTo(arr_ny), c_ptrTo(arr_ny), -1, FFTW_ESTIMATE);

  var loc_D3 = _D3.localSubdomain();
  var loc_D3_hat = _D3_hat.localSubdomain();
  var loc_D3_hatT = _D3_hatT.localSubdomain();

  /* First forward FFT */

    forall (i,j) in _D_zy.localSubdomain() {
      fftw_execute_dft_r2c(plan_f1, c_ptrTo(in_arr.localSlice(loc_D3)[i,j,0]), c_ptrTo(tmp_f1.localSlice(loc_D3_hatT)[i,j,0]));
    }

  /* Transpose */
    allLocalesBarrier.barrier();
    transpose_3D(tmp_f1, tmp_f1T);
    allLocalesBarrier.barrier();

  /* Second forward FFT */
    forall (i,j) in _D_zxhat.localSubdomain() {
      fftw_execute_dft(plan_f2, c_ptrTo(tmp_f1T.localSlice(loc_D3_hat)[i,j,0]), c_ptrTo(out_arr.localSlice(loc_D3_hat)[i,j,0]));
    }
    allLocalesBarrier.barrier();

}

proc execute_forward_FFTs_single_level(ref in_arr: [] real, ref out_arr : [] complex) {

  var arr_nx : [D_nx] real;
  var arr_nx2p : [D_nx2p] complex;
  var arr_ny : [D_ny] complex;

  // First forward FFT is along the x dimension, which takes an nx-vector as input and returns an nx2p-vector as output
    var plan_f1 = fftw_plan_dft_r2c_1d(nx : c_int, c_ptrTo(arr_nx), c_ptrTo(arr_nx2p), FFTW_ESTIMATE);

  // Second forward FFT is along the y dimension, for which both input and output are length ny
    var plan_f2 = fftw_plan_dft_1d(ny : c_int, c_ptrTo(arr_ny), c_ptrTo(arr_ny), -1, FFTW_ESTIMATE);

  var loc_D = _D.localSubdomain();
  var loc_D_hat = _D_hat.localSubdomain();
  var loc_D_hatT = _D_hatT.localSubdomain();

  /* First forward FFT */
    forall i in _D_ny.localSubdomain() {
      fftw_execute_dft_r2c(plan_f1, c_ptrTo(in_arr.localSlice(loc_D)[i,0]), c_ptrTo(tmp_f1_2D.localSlice(loc_D_hatT)[i,0]));
    }

  /* Transpose */
    allLocalesBarrier.barrier();
    transpose_2D(tmp_f1_2D, tmp_f1T_2D);
    allLocalesBarrier.barrier();

  /* Second forward FFT */
    forall i in _D_nx2p.localSubdomain() {
      fftw_execute_dft(plan_f2, c_ptrTo(tmp_f1T_2D.localSlice(loc_D_hat)[i,0]), c_ptrTo(out_arr.localSlice(loc_D_hat)[i,0]));
    }
    allLocalesBarrier.barrier();

}

proc execute_backward_FFTs(ref in_arr: [] complex, ref out_arr : [] real) {

  var arr_nx : [D_nx] real;
  var arr_nx2p : [D_nx2p] complex;
  var arr_ny : [D_ny] complex;

  // First backward FFT is along the y dimension, for which both input and output are length ny
    var plan_b1 = fftw_plan_dft_1d(ny : c_int, c_ptrTo(arr_ny), c_ptrTo(arr_ny), 1, FFTW_ESTIMATE);

  // Second backward FFT is along the x dimension, which takes an nx2p-vector as input and returns an nx-vector as output
    var plan_b2 = fftw_plan_dft_c2r_1d(nx : c_int, c_ptrTo(arr_nx2p), c_ptrTo(arr_nx), FFTW_ESTIMATE);

  var loc_D3 = _D3.localSubdomain();
  var loc_D3_hat = _D3_hat.localSubdomain();
  var loc_D3_hatT = _D3_hatT.localSubdomain();

  /* First backward FFT */
    forall (i,j) in _D_zxhat.localSubdomain() {
      fftw_execute_dft(plan_b1, c_ptrTo(in_arr.localSlice(loc_D3_hat)[i,j,0]), c_ptrTo(tmp_b1.localSlice(loc_D3_hat)[i,j,0]));
    }

  /* Transpose */
    allLocalesBarrier.barrier();
    transpose_3D(tmp_b1, tmp_b1T);
    allLocalesBarrier.barrier();

  /* Second backward FFT */
    forall (i,j) in _D_zy.localSubdomain() {
      fftw_execute_dft_c2r(plan_b2, c_ptrTo(tmp_b1T.localSlice(loc_D3_hatT)[i,j,0]), c_ptrTo(out_arr.localSlice(loc_D3)[i,j,0]));
    }
    allLocalesBarrier.barrier();

}

proc execute_backward_FFTs_single_level(ref in_arr: [] complex, ref out_arr : [] real) {

    /* Temporary arrays to hold 1D transforms */
      var tmp_b1 : [D_hat] complex;
      var tmp_b1T : [D_hatT] complex;

  /* First backward FFT */
    forall i in D_nx2p {
      fftw_execute_dft(plan_b1, c_ptrTo(in_arr[i,0]), c_ptrTo(tmp_b1[i,0]));
    }

  /* Transpose */
    transpose_2D(tmp_b1, tmp_b1T);

  /* Second backward FFT */
    forall i in D_ny {
      fftw_execute_dft_c2r(plan_b2, c_ptrTo(tmp_b1T[i,0]), c_ptrTo(out_arr[i,0]));
    }

}

proc transpose_2D(ref in_arr: [?D_in] complex, ref out_arr: [?D_out] complex) {

    forall (i,j) in D_in.localSubdomain() {
      out_arr[j,i] = in_arr[i,j];
    }

}

proc transpose_3D(ref in_arr: [?D_in] complex, ref out_arr: [?D_out] complex) {

//    forall (i,j,k) in D_in.localSubdomain() with (var agg = new DstAggregator(complex)) do {
//      agg.copy(out_arr[i,k,j], in_arr[i,j,k]);
//    }

//    forall (i,j,k) in D_in with (var agg = new DstAggregator(complex)) do {
//      agg.copy(out_arr[i,k,j], in_arr[i,j,k]);
//    }

    forall (i,j,k) in D_in.localSubdomain() {
      out_arr[i,k,j] = in_arr[i,j,k];
    }

}

proc normalize(ref in_arr: [?dom] real) {

    var norm = nx*ny;
    forall (i,j,k) in dom.localSubdomain() {
      in_arr[i,j,k] = in_arr[i,j,k] / norm;
    }

}

