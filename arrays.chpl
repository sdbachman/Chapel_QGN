use parameters;
use domains;

/* Horizontal grid spacing */
  const dx : real = Lx / nx;
  const dy : real = Ly / ny;

/* Vertical layer depths: H(nz) is the bottom layer; S=f^2/N^2(z) */
  var H : [zl] real;
  var S : [zi] real;

/* Zonal mean velocity profile and associated meridional PV gradient */
  var uBar : [zl] real;
  var qyBar : [zl] real;

/* Vertical modes, EVals = -kd^2 */
  var Modes : [zl2] real;
  var EVals : [zl] real;

/* L matrix, LeftEVs */
  var L : [zl2] real;
  var ModesL : [zl2] real;

/* Imaginary part of eigenvalues */
  var EValsI : [zl] real;

/* Wavenumbers */
  var kx : [D_hat] real;
  var ky : [D_hat] real;
  var k2 : [D_hat] real;

/* Spectral potential vorticity */
  var q_hat : [D3_hat] complex;
