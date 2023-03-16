use parameters;
use domains;

/* Location of horizontal grid */
  var x : [_D] real;
  var y : [_D] real;

/* Horizontal grid spacing */
  const dx : real = Lx / nx;
  const dy : real = Ly / ny;

/* Vertical grid */
  var z : [_zl] real;

/* Vertical layer depths: H(nz) is the bottom layer; S=f^2/N^2(z) */
  var H : [_zl] real;
  var S : [_zi] real;

/* Zonal mean velocity profile and associated meridional PV gradient */
  var uBar : [_zl] real;
  var qyBar : [_zl] real;

/* Vertical modes, EVals = -kd^2 */
  var Modes : [_zl2] real;
  var EVals : [_zl] real;

/* L matrix, LeftEVs */
  var L : [_zl2] real;
  var ModesL : [_zl2] real;

/* Imaginary part of eigenvalues */
  var EValsI : [_zl] real;

/* Wavenumbers */
  var kx : [_D_hat] real;
  var ky : [_D_hat] real;
  var k2 : [_D_hat] real;

/* Potential vorticity */
  var q : [_D3] real;

/* Spectral potential vorticity */
  var q_hat : [_D3_hat] complex;

/* Spectral streamfunction */
  var psi_hat : [_D3_hat] complex;

/* Spectral buoyancy */
  var b_hat : [_D3_hat] complex;

/* Spectral velocities */
  var u_hat : [_D3_hat] complex;
  var v_hat : [_D3_hat] complex;

/* Physical space arrays */
  var q_phys : [_D3] real;
  var u_phys : [_D3] real;
  var v_phys : [_D3] real;

/* For the Jacobian */
  var uq_hat : [_D3_hat] complex;
  var vq_hat : [_D3_hat] complex;
  var uq_phys : [_D3] real;
  var vq_phys : [_D3] real;

/* For the drag term */
  var drag_tmp : [_D_hat] complex;
  var drag_hat : [_D_hat] complex;

/* Temporary arrays to hold 1D transforms */
  var tmp_f1 : [_D3_hatT] complex;
  var tmp_f1T : [_D3_hat] complex;
  var tmp_f1_2D : [_D_hatT] complex;
  var tmp_f1T_2D : [_D_hat] complex;
  var tmp_b1 : [_D3_hat] complex;
  var tmp_b1T : [_D3_hatT] complex;
