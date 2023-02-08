use parameters;
use domains;
use arrays;
use QGN_Module;
use ARK43;
use IO_Module;

use compare_fortran;
use FFT_utils;
use Time;
use IO;

proc main() {

  var t0 : stopwatch;
  t0.start();

  Initialize();

/*
  for i in (Nt_start+1)..(Nt_start+Nt) {

    TimeStep();
    DeAlias(q_hat);
    Diagnostics(i);

  }

  openwriter("timings.txt").writeln(t0.elapsed());
*/
}
