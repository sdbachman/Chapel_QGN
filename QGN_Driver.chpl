use parameters;
use domains;
use arrays;
use QGN_Module;
use ARK43;
use compare_fortran;
use FFT_utils;

use IO_Module;
use Time;
use IO;

proc main() {

  var t0 : stopwatch;
  t0.start();

  Initialize();

  for i in 1..Nt {

 //   var t2 : stopwatch;
 //   t2.start();

    TimeStep();
    DeAlias(q_hat);
    if ((i % diag_freq)==0) {
      var q2 : [D3] real;
      execute_backward_FFTs(q_hat, q2);
      normalize(q2);
      WriteOutput(q2, i);
    }

//    var q2 : [D3] real;
//    execute_backward_FFTs(q_hat, q2);
//    normalize(q2);

//difference("out", q2);

//    print_array_3D(q2);

//    writeln();
//    writeln(t2.elapsed());
//    writeln();
  }

  openwriter("timings.txt").writeln(t0.elapsed());

}
