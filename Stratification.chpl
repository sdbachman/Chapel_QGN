use parameters;
use arrays;

/* Even grid spacing */
  H = Htot / nz;

writeln("-------------------------------------------------");
writeln(" Layer depths H/Htot ");
writeln(H/Htot);
writeln("-------------------------------------------------");

/* N(z) = constant */
  S[0] = 1.0 / 225;
  for k in 1..nz {
    S[k] = S[0];
  }

writeln("-------------------------------------------------");
writeln(" f^2/N^2 ");
writeln(S);
writeln("-------------------------------------------------");
