use BlockDist;
use AllLocalesBarriers;
use CopyAggregation;
use Time;

///////////////////////////////////////
//   Set up domain size
///////////////////////////////////////

config const trials = 5;
config const nz = 3;
config const nx = 10;
const ny = nx; // Must be a square domain

if (nx%numLocales) {
  writeln("ERROR: Number of locales must evenly divide nx. Aborting.");
  exit();
}

///////////////////////////////////////
//   Set up blocks
///////////////////////////////////////

const block_size = nx / numLocales;
const D_block_ind = {0..0, 0..<block_size, 0..<block_size};


///////////////////////////////////////
//   Set up domains and arrays
///////////////////////////////////////

const myTargetLocales3D = reshape(Locales, {1..1, 1..Locales.size, 1..1});
var D = {0..<nz, 0..<ny, 0..<nx};
var _D = D dmapped Block(D, targetLocales=myTargetLocales3D);
var A : [_D] complex;
var AT : [_D] complex;

const nz_blocks = {0..<nz, 0..<numLocales};

///////////////////////////////////////
//   Set up timer
///////////////////////////////////////
var timings_tile : [0..trials] real;
var mean_timing_tile : real;

var timings_agg : [0..trials] real;
var mean_timing_agg : real;

///////////////////////////////////////
//   Parallel loop
///////////////////////////////////////

coforall loc in Locales do on loc {

  /////////////////////////////////////
  //   Initialize array
  /////////////////////////////////////
  for (i,j,k) in _D.localSubdomain() {
    A[i,j,k] = i*ny*nx + j*nx + k;
  }

  allLocalesBarrier.barrier();

  for trial in 0..trials {
    var t_tile : stopwatch;
    t_tile.start();

      forall block in 0..<numLocales {
          var D_block = {0..<nz,(here.id*block_size)..((here.id+1)*block_size-1), (block*block_size)..((block+1)*block_size-1)};
        if here.id == block {
          AT[D_block] = A[D_block];
        } else {
          const blockIdx = (here.id, block);
          const dest = Locales[block];
          var dest_block = {0..<nz, D_block.dim(2), D_block.dim(1)};
          AT[dest_block] = A[D_block];
        }
      }
      allLocalesBarrier.barrier();

      /////////////////////////////////////
      //   Transpose block by block
      /////////////////////////////////////
      // horizontal strips of blocks on the same locale
      forall (k,block) in nz_blocks {
        var D_block = {k..k,(here.id*block_size)..((here.id+1)*block_size-1), (block*block_size)..((block+1)*block_size-1)};
        const first = D_block.first;

        for i in 0..<block_size {
          for j in (i+1)..<block_size {
            ref x = AT[k,first[1]+i,first[2]+j];
            ref y = AT[k,first[1]+j,first[2]+i];
            var temp = x;
            x = y;
            y = temp;
          }
        }

      } // forall

   // } // k

    t_tile.stop();
    timings_tile[trial] = t_tile.elapsed();
    var t_agg : stopwatch;
    t_agg.start();
    forall (i,j,k) in _D.localSubdomain() with (var agg = new DstAggregator(complex)) do {
      agg.copy(AT[i,k,j], A[i,j,k]);
    }

    t_agg.stop();
    timings_agg[trial] = t_agg.elapsed();
  } // trials
} // coforall

mean_timing_tile = (+ reduce timings_tile[1..]) / trials;
mean_timing_agg = (+ reduce timings_agg[1..]) / trials;

writeln("Timings for tiled transpose are:");
writeln(timings_tile[1..]);
writeln("Mean timing for tiled transpose is:");
writeln(mean_timing_tile);
writeln();

writeln("Timings for aggregated transpose are:");
writeln(timings_agg[1..]);
writeln("Mean timing for aggregated transpose is:");
writeln(mean_timing_agg);

/*
for k in 0..<nz {
  writeln("A at layer ", k, " is:");
  writeln(A[k, .., ..]);
  writeln();
  writeln("AT at layer ", k, " is:");
  writeln(AT[k, .., ..]);
  writeln();
}
*/
