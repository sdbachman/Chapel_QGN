use parameters;
use BlockDist;
use ReplicatedDist;

////////////////////////////////////////////////////////////
//   NOTE: All spectral arrays will be transposed with    //
//   respect to their physical space counterparts. This   //
//   is to enable the code to use serial 1D FFTs.         //
////////////////////////////////////////////////////////////

const myTargetLocales2D = reshape(Locales, {1..Locales.size, 1..1});
const myTargetLocales2D_hat = reshape(Locales, {1..1, 1..Locales.size});
const myTargetLocales3D = reshape(Locales, {1..1, 1..Locales.size, 1..1});
const myTargetLocales3D_hat = reshape(Locales, {1..1, 1..1, 1..Locales.size});

/* Vertical domain at layer centers */
var zl = {0..#nz};
var _zl = zl dmapped Replicated();
var zl2 = {0..#nz, 0..#nz};
var _zl2 = zl2 dmapped Replicated();

/* Vertical domain at layer interfaces */
var zi = {0..#(nz+1)};
var _zi = zi dmapped Replicated();

/* Horizontal domain (physical space) */
var D = {0..#ny,0..#nx};
var _D = D dmapped Block(D, targetLocales=myTargetLocales2D);

/* Horizontal domain (spectral space) */
var D_hat = {0..#nx2p, 0..#ny};
var _D_hat = D_hat dmapped Block(D_hat, targetLocales=myTargetLocales2D);
var D_hatT = {0..#ny, 0..#nx2p};
var _D_hatT = D_hatT dmapped Block(D_hatT, targetLocales=myTargetLocales2D);

/* 1D slices of the horizontal domain */
var D_nx = {0..#nx};
var _D_nx = D_nx dmapped Replicated();
var D_ny = {0..#ny};
//var _D_ny = D_ny dmapped Replicated();
var _D_ny = D_ny dmapped Block(D_ny);
var D_nx2p = {0..#nx2p};
//var _D_nx2p = D_nx2p dmapped Replicated();
var _D_nx2p = D_nx2p dmapped Block(D_nx2p);
var D_ny2p = {0..#ny2p};

/* 2D vertical slices */
var D_zx = {0..#nz,0..#nx};
var D_zy = {0..#nz,0..#ny};
var _D_zy = D_zy dmapped Block(D_zy, targetLocales=myTargetLocales2D_hat);
var D_zxhat = {0..#nz,0..#nx2p};
var _D_zxhat = D_zxhat dmapped Block(D_zxhat, targetLocales=myTargetLocales2D_hat);

/* 3D domain */
var D3 = {0..#nz,0..#ny,0..#nx};
var _D3 = D3 dmapped Block(D3, targetLocales=myTargetLocales3D);
var D3_hat = {0..#nz,0..#nx2p,0..#ny};
var _D3_hat = D3_hat dmapped Block(D3_hat, targetLocales=myTargetLocales3D);
var D3_hatT = {0..#nz,0..#ny,0..#nx2p};
var _D3_hatT = D3_hatT dmapped Block(D3_hatT, targetLocales=myTargetLocales3D);

/* Special domains for DeAliasing */
var D3_hat_sp1 = {0..#nz, (nx3p-1)..(nx2p-1), 0..#ny};
var _D3_hat_sp1 = D3_hat_sp1 dmapped Replicated();
var D3_hat_sp2 = {0..#nz, 0..#nx2p, (ny3p-1)..(ny3p2-1)};
var _D3_hat_sp2 = D3_hat_sp2 dmapped Replicated();
