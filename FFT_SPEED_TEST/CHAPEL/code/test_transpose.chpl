use parameters;
use domains;
use arrays;
use FFT_utils;
//use QGN_Module;

use compare_fortran;
use Time;
use IO;

var q_hatT : [D3_hatT] complex;
var q_hatT2 : [D3_hatT] complex;
var q2 : [D3] real;
var q3 : [D3] real;

config const filename : string = "test_grid";

var t1 : stopwatch;
var t2 : stopwatch;
var t3 : stopwatch;
var t4 : stopwatch;
var t5 : stopwatch;
var t6 : stopwatch;
var t7 : stopwatch;
var t8 : stopwatch;
var t9 : stopwatch;

t1.start();

load_fortran_grid(filename, q);
writeln(q.shape);
set_up_forward_FFTs();
set_up_backward_FFTs();
t1.stop();

t2.start();
execute_forward_FFTs(q, q_hat);
t2.stop();

t3.start();
for i in 1..10000 {
transpose_3D(q_hat, q_hatT);
}
t3.stop();

/*
t4.start();
execute_backward_FFTs(q_hat, q2);
t4.stop();

t5.start();
for k in 1..nz {
  q2[k,..,..] = 2*q2[k,..,..]+1;
}
t5.stop();

create_initial_state(q3);

t6.start();
for i in 1..10000 {
execute_forward_FFTs_2D(q3, q_hatT2);
}
t6.stop();

t7.start();
execute_backward_FFTs_2D(q_hatT2, q2);
t7.stop();
normalize(q2);

t8.start();
execute_forward_FFTs_3D(q, q_hatT2);
t8.stop();
*/

//writeln("Setup: ", t1.elapsed());
//writeln("Forward FFT: ", t2.elapsed() - t3.elapsed());

openwriter("timings_tranpose.txt").writeln(t3.elapsed());

/*
writeln("Transpose: ", t3.elapsed());
writeln("Backward FFT: ", t4.elapsed() - t3.elapsed());
writeln("Rank reduction slice: ", t5.elapsed());
writeln("2D FFT: ", t6.elapsed());
writeln("2D FFT backwards: ", t7.elapsed());
writeln("3D FFT: ", t8.elapsed());
*/
