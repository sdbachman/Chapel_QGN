use IO;
use FFTW;
use C_FFTW;
use CTypes;

use parameters;
use domains;
use FFT_utils;

proc Initialize() {

var arr : [D] real;


var arr_out_f1 : [D_out] complex;
var arr_f : [D_outT] complex;
var arr_b : [D] real;
var arr_b_norm : [D] real;

var f = open("../test_grid.dat", iomode.r);
var r = f.reader(kind=ionative);

for j in 1..ny {
  for i in 1..nx {
    var tmp : real;
    r.readBinary(tmp);
    //writef("%r\n", tmp);
    arr[i,j] = tmp;
  }
}
r.close();

writeln("grid:");
writeln(arr);
writeln();


set_up_forward_FFTs();
set_up_backward_FFTs();


execute_forward_FFTs(arr, arr_f);

execute_backward_FFTs(arr_f, arr_b);


normalize(arr_b, arr_b_norm);

writeln("Diff between grid before and after FFTs:");
writeln(arr - arr_b_norm);


}
