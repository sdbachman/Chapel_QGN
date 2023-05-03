use parameters;
use domains;

/* Location of horizontal grid */
  var x : [D] real(rp);
  var y : [D] real(rp);

/* Horizontal grid spacing */
  const dx : real(rp) = Lx / nx;
  const dy : real(rp) = Ly / ny;

/* Vertical grid */
  var z : [zl] real(rp);

/* Vertical layer depths: H(nz) is the bottom layer; S=f^2/N^2(z) */
  var H : [zl] real(rp);
  var S : [zi] real(rp);

/* Zonal mean velocity profile and associated meridional PV gradient */
  var uBar : [zl] real(rp);
  var qyBar : [zl] real(rp);

/* Vertical modes, EVals = -kd^2 */
  var Modes : [zl2] real(rp);
  var EVals : [zl] real(rp);

/* L matrix, LeftEVs */
  var L : [zl2] real(rp);
  var ModesL : [zl2] real(rp);

/* Imaginary part of eigenvalues */
  var EValsI : [zl] real(rp);

/* Wavenumbers */
  var kx : [D_hat] real(rp);
  var ky : [D_hat] real(rp);
  var k2 : [D_hat] real(rp);

/* Potential vorticity */
  var q : [D3] real(rp);

/* Spectral potential vorticity */
  var q_hat : [D3_hat] complex(cp);

/* Spectral streamfunction */
  var psi_hat : [D3_hat] complex(cp);

/* Modal decomposition */
  var q_hat_mode : [D3_hat] complex(cp);
  var psi_hat_mode : [D3_hat] complex(cp);

/* Spectral buoyancy */
  var b_hat : [D3_hat] complex(cp);

/* Spectral velocities */
  var u_hat : [D3_hat] complex(cp);
  var v_hat : [D3_hat] complex(cp);

/* Physical space arrays */
  var q_phys : [D3] real(rp);
  var u_phys : [D3] real(rp);
  var v_phys : [D3] real(rp);

/* For the Jacobian */
  var uq_hat : [D3_hat] complex(cp);
  var vq_hat : [D3_hat] complex(cp);
  var uq_phys : [D3] real(rp);
  var vq_phys : [D3] real(rp);

/* For the drag term */
  var drag_tmp : [D_hat] complex(cp);
  var drag_hat : [D_hat] complex(cp);
  var Uu_drag : [D] real(rp);
  var Uv_drag : [D] real(rp);

/* Arrays and variables for the ARK43 timestepping */
  var N1, N2, N3, N4, N5, N6 : [D3_hat] complex(cp);
  var L1, L2, L3, L4, L5, L6 : [D3_hat] complex(cp);
  var q_tmp : [D3_hat] complex(cp);
  var Mq : [D3_hat] complex(cp);
  var k8 : [D_hat] real(rp);

  var err : [D3] real(rp);
  var err_hat : [D3_hat] complex(cp);

  var ae, ai : [ark2d] real(rp);
  var b, be : [ark1d] real(rp);

  // These will be only on Locale 0
  var err0, err1 : real(rp);
  var reject : bool = false;

/* For QG Leith */
  var A8L : [zl] real(rp);
  var q8_tmp : [D3_hat] real(rp);
  var k4 : [D_hat] real(rp);
