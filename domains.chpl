use parameters;

////////////////////////////////////////////////////////////
//   NOTE: All spectral arrays will be transposed with    //
//   respect to their physical space counterparts. This   //
//   is to enable the code to use serial 1D FFTs.         //
////////////////////////////////////////////////////////////

/* Vertical domain at layer centers */
var zl = {1..nz};
var zl2 = {1..nz, 1..nz};

/* Vertical domain at layer interfaces */
var zi = {0..nz};

/* Horizontal domain (physical space) */
var D = {1..ny,1..nx};

/* Horizontal domain (spectral space) */
var D_hat = {1..nx2p, 1..ny};
var D_hatT = {1..ny, 1..nx2p};

/* 1D slices of the horizontal domain */
var D_nx = {1..nx};
var D_ny = {1..ny};
var D_nx2p = {1..nx2p};
var D_ny2p = {1..ny2p};

/* 2D vertical slices */
var D_zx = {1..nz,1..nx};
var D_zy = {1..nz,1..ny};
var D_zxhat = {1..nz,1..nx2p};

/* 3D domain */
var D3 = {1..nz,1..ny,1..nx};
var D3_hat = {1..nz,1..nx2p,1..ny};
var D3_hatT = {1..nz,1..ny,1..nx2p};

/* Special domains for DeAliasing */
var D3_hat_sp1 = {1..nz, nx3p..nx2p, 1..ny};
var D3_hat_sp2 = {1..nz, 1..nx2p, ny3p..ny3p2};
