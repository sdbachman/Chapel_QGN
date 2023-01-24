use parameters;
use arrays;
use QGN_Module;
use ARK43;
use compare_fortran;

use Time;
use IO;

proc main() {

  var t : Timer;
  t.start();

  Initialize();

  for i in 1..Nt {

    var t2 : Timer;
    t2.start();

    TimeStep();

    writeln(t2.elapsed());
  }

  openwriter("timings.txt").writeln(t.elapsed());

}
