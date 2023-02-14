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
var D_hat = {0..#nx2p, 0..#ny};
var D_hatT = {0..#ny, 0..#nx2p};

/* 1D slices of the horizontal domain */
var D_nx = {0..#nx};
var D_ny = {0..#ny};
var D_nx2p = {0..#nx2p};
var D_ny2p = {0..#ny2p};

/* 2D vertical slices */
var D_zx = {0..#nz,0..#nx};
var D_zy = {0..#nz,0..#ny};
var D_zxhat = {0..#nz,0..#nx2p};

/* 3D domain */
var D3 = {0..#nz,0..#ny,0..#nx};
var D3_hat = {0..#nz,0..#nx2p,0..#ny};
var D3_hatT = {0..#nz,0..#ny,0..#nx2p};

/* Special domains for DeAliasing */
var D3_hat_sp1 = {0..#nz, (nx3p-1)..(nx2p-1), 0..#ny};
var D3_hat_sp2 = {0..#nz, 0..#nx2p, (ny3p-1)..(ny3p2-1)};
