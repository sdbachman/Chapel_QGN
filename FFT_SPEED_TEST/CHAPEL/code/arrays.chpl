use parameters;
use domains;

/* Location of horizontal grid */
  var x : [D] real;
  var y : [D] real;

/* Horizontal grid spacing */
  const dx : real = Lx / nx;
  const dy : real = Ly / ny;

/* Vertical grid */
  var z : [zl] real;

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

/* Potential vorticity */
  var q : [D3] real;

/* Spectral potential vorticity */
  var q_hat : [D3_hat] complex;

/* Spectral streamfunction */
  var psi_hat : [D3_hat] complex;

/* Spectral buoyancy */
  var b_hat : [D3_hat] complex;

/* Spectral velocities */
  var u_hat : [D3_hat] complex;
  var v_hat : [D3_hat] complex;

/* Physical space arrays */
  var q_phys : [D3] real;
  var u_phys : [D3] real;
  var v_phys : [D3] real;

