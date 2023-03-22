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

  coforall loc in Locales with (ref Nt_start) do on loc {

    var t0 : stopwatch;
    t0.start();

    Initialize();
 writeln("In driver Nt_start is: ", Nt_start);

    for i in (Nt_start+1)..(Nt_start+Nt) {

      TimeStep();

      DeAlias(q_hat);

      Diagnostics(i);

      //if (here.id == 0) {
      //  writeln(t0.elapsed());
      //}

    } // Timestepping loop

  } // coforall loop

} // main
