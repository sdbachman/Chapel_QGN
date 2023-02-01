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

    TimeStep();
    DeAlias(q_hat);

    /* Write diagnostics */
      if ((i % diag_freq)==0) {
        var q2 : [D3] real;
        execute_backward_FFTs(q_hat, q2);
        normalize(q2);
        WriteOutput(q2, i);
      }

  }

  openwriter("timings.txt").writeln(t0.elapsed());

}
