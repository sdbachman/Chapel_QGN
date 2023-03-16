use parameters;
use domains;
use arrays;
use QGN_Module;
use ARK43;
use IO_Module;
use AllLocalesBarriers;

use compare_fortran;
use FFT_utils;
use Time;
use IO;

proc main() {

  var t0 : stopwatch;
  t0.start();

  coforall loc in Locales do on loc {

    //var storage = new myClass;

    // Put these in a separate subroutine?
    H = H.replicand(Locales[0]);
    S = S.replicand(Locales[0]);

    Initialize();

    for i in (Nt_start+1)..(Nt_start+Nt) {

      TimeStep();

      allLocalesBarrier.barrier();
      DeAlias(q_hat);

      Diagnostics(i);

    }

    openwriter("timings.txt").writeln(t0.elapsed());
  }
}
