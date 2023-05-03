use parameters;

////////////////////////////////////////////////////////////
//   NOTE: All spectral arrays will be transposed with    //
//   respect to their physical space counterparts. This   //
//   is to enable the code to use serial 1D FFTs.         //
////////////////////////////////////////////////////////////

/* Vertical domain at layer centers */
var zl = {0..#nz};
var zl2 = {0..#nz, 0..#nz};

/* Vertical domain at layer interfaces */
var zi = {0..#(nz+1)};

/* Horizontal domain (physical space) */
var D = {0..#ny,0..#nx};

/* Horizontal domain (spectral space) */
var D_hat = {0..#ny, 0..#nx2p};

/* 1D slices of the horizontal domain */
const D_nx = {0..#nx};
const D_ny = {0..#ny};
const D_nx2p = {0..#nx2p};
var D_ny2p = {0..#ny2p};

/* 2D vertical slices */
var D_zx = {0..#nz,0..#nx};
var D_zy = {0..#nz,0..#ny};
var D_zxhat = {0..#nz,0..#nx2p};

/* 3D domain */
var D3 = {0..#nz,0..#ny,0..#nx};
var D3_hat = {0..#nz,0..#ny,0..#nx2p};

/* Special domains for DeAliasing */
var D3_hat_sp1 = {0..#nz, 0..#ny, (nx3p-1)..(nx2p-1)};
var D3_hat_sp2 = {0..#nz, (ny3p-1)..(ny3p2-1), 0..#nx2p};

/* Special domains for ARK43 */
var ark1d = {1..6};
var ark2d = {1..6, 1..6};


