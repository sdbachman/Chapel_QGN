use IO;
use FFTW;
use C_FFTW;
use CTypes;
use CopyAggregation;

//use Time;
use parameters;
use domains;
use arrays;

//use compare_fortran;


/* Forward FFT plans */
var plan_f : fftw_plan;

/* Backward FFT plans */
var plan_b : fftw_plan;


proc set_up_forward_FFTs() {

    var arr2_in : [D] real(rp);
    var arr2_out : [D_hat] complex(cp);

    if (rp == 64) {
      extern proc fftw_plan_dft_r2c_2d(n0: c_int, n1: c_int, in_arg: c_ptr(real(rp)), out_arg: c_ptr(complex(cp)), flags: c_uint) : fftw_plan;
      plan_f = fftw_plan_dft_r2c_2d(ny : c_int, nx : c_int, c_ptrTo(arr2_in), c_ptrTo(arr2_out), FFTW_ESTIMATE);
    } else {
      extern proc fftwf_plan_dft_r2c_2d(n0: c_int, n1: c_int, in_arg: c_ptr(real(rp)), out_arg: c_ptr(complex(cp)), flags: c_uint) : fftw_plan;
      plan_f = fftwf_plan_dft_r2c_2d(ny : c_int, nx : c_int, c_ptrTo(arr2_in), c_ptrTo(arr2_out), FFTW_ESTIMATE);
    }

}

proc set_up_backward_FFTs() {

    var arr2_in : [D_hat] complex(cp);
    var arr2_out : [D] real(rp);

    if (rp == 64) {
      extern proc fftw_plan_dft_c2r_2d(n0: c_int, n1: c_int, in_arg: c_ptr(complex(cp)), out_arg: c_ptr(real(rp)), flags: c_uint) : fftw_plan;
      plan_b = fftw_plan_dft_c2r_2d(ny : c_int, nx : c_int, c_ptrTo(arr2_in), c_ptrTo(arr2_out), FFTW_ESTIMATE);
    } else {
      extern proc fftwf_plan_dft_c2r_2d(n0: c_int, n1: c_int, in_arg: c_ptr(complex(cp)), out_arg: c_ptr(real(rp)), flags: c_uint) : fftw_plan;
      plan_b = fftwf_plan_dft_c2r_2d(ny : c_int, nx : c_int, c_ptrTo(arr2_in), c_ptrTo(arr2_out), FFTW_ESTIMATE);
    }

}

proc execute_forward_FFTs(ref in_arr: [] real(rp), ref out_arr : [] complex(cp)) {

    if (rp == 64) {
      extern proc fftw_execute_dft_r2c(p: fftw_plan, in_arg: c_ptr(real(rp)), out_arg: c_ptr(complex(cp)));

      forall k in zl {
        fftw_execute_dft_r2c(plan_f, c_ptrTo(in_arr[k,0,0]), c_ptrTo(out_arr[k,0,0]));
      }
    } else {
      extern proc fftwf_execute_dft_r2c(p: fftw_plan, in_arg: c_ptr(real(rp)), out_arg: c_ptr(complex(cp)));

      forall k in zl {
        fftwf_execute_dft_r2c(plan_f, c_ptrTo(in_arr[k,0,0]), c_ptrTo(out_arr[k,0,0]));
      }
    }

}

proc execute_forward_FFTs_single_level(ref in_arr: [] real(rp), ref out_arr : [] complex(cp)) {

    if (rp == 64) {
      extern proc fftw_execute_dft_r2c(p: fftw_plan, in_arg: c_ptr(real(rp)), out_arg: c_ptr(complex(cp)));

      fftw_execute_dft_r2c(plan_f, c_ptrTo(in_arr), c_ptrTo(out_arr));
    } else {
      extern proc fftwf_execute_dft_r2c(p: fftw_plan, in_arg: c_ptr(real(rp)), out_arg: c_ptr(complex(cp)));

      fftwf_execute_dft_r2c(plan_f, c_ptrTo(in_arr), c_ptrTo(out_arr));
    }

}

proc execute_backward_FFTs(in in_arr: [] complex(cp), ref out_arr : [] real(rp)) {

    if (rp == 64) {
      extern proc fftw_execute_dft_c2r(p: fftw_plan, in_arg: c_ptr(complex(cp)), out_arg: c_ptr(real(rp)));

      forall k in zl {
        fftw_execute_dft_c2r(plan_b, c_ptrTo(in_arr[k,0,0]), c_ptrTo(out_arr[k,0,0]));
      }
    } else {
      extern proc fftwf_execute_dft_c2r(p: fftw_plan, in_arg: c_ptr(complex(cp)), out_arg: c_ptr(real(rp)));

      forall k in zl {
        fftwf_execute_dft_c2r(plan_b, c_ptrTo(in_arr[k,0,0]), c_ptrTo(out_arr[k,0,0]));
      }
    }

}

proc transpose_2D(ref in_arr: [?D_in] complex(cp), ref out_arr: [?D_out] complex(cp)) {

    forall (i,j) in D_in.localSubdomain() {
      out_arr[j,i] = in_arr[i,j];
    }

}

proc transpose_3D(ref in_arr: [?D_in] complex(cp), ref out_arr: [?D_out] complex(cp)) {

    forall (i,j,k) in D_in.localSubdomain() with (var agg = new DstAggregator(complex(cp))) do {
      agg.copy(out_arr[i,k,j], in_arr[i,j,k]);
    }

}

proc normalize(ref in_arr: [?dom] real(rp)) {

    var norm = nx*ny;
    forall (i,j,k) in dom {
      in_arr[i,j,k] = in_arr[i,j,k] / norm;
    }

}

