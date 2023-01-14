use parameters;
use arrays;

var tmp : real;

/* Set mean shear */
  uBar[1] = 0.1*cos(pi*H[1]/(2*Htot));

  for k in 2..nz {
    tmp = + reduce H[1..k-1];
    uBar[k] = 0.1*cos(pi*(tmp + H[k]/2)/Htot);
  }

  tmp = + reduce (H*uBar);
  uBar = uBar - tmp/Htot;
