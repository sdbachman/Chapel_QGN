use parameters;
use arrays;

/*
var UBarBC : [1..50] real;
var UBarOC : [1..50] real;

UBarBC[1..10] =  [ 217, 217, 216, 215, 213, 211, 208, 204, 200, 195 ];
UBarBC[11..20] = [ 190, 184, 178, 171, 163, 155, 147, 137, 128, 120 ];
UBarBC[21..30] = [ 111, 102,  92,  83,  74,  64,  53,  43,  32,  23 ];
UBarBC[31..40] = [  14,   5,  -4, -12, -19, -26, -32, -38, -43, -47 ];
UBarBC[41..50] = [ -51, -53, -55, -56, -57, -58, -58, -59, -59, -59 ];
UBarBC = (1.0/276)*UBarBC;

UBarOC[1..10]  = [ 446, 442, 437, 431, 423, 415, 406, 396, 385, 374 ];
UBarOC[11..20] = [ 363, 350, 337, 325, 311, 297, 282, 266, 252, 238 ];
UBarOC[21..30] = [ 225, 211, 198, 184, 171, 157, 142, 128, 115, 102 ];
UBarOC[31..40] = [  91,  79,  68,  58,  49,  40,  32,  26,  20,  14 ];
UBarOC[41..50] = [  10,   7,   5,   4,   2,   2,   1,   0,   0,   0 ];
UBarOC = (1.0/446)*UBarOC;

var total_amp = 5.0;
var bc_amp = 0.1;
var oc_amp = 0.5;
uBar = total_amp * ((bc_amp)*UBarBC + (oc_amp)*UBarOC);
*/


var tmp : real;

/* Set mean shear */
  uBar[0] = 0.1*cos(pi*H[0]/(2*Htot));

  for k in 1..(nz-1) {
    tmp = + reduce H[0..(k-1)];
    uBar[k] = 0.1*cos(pi*(tmp + H[k]/2)/Htot);
  }

  tmp = + reduce (H*uBar);
  uBar = uBar - tmp/Htot;

