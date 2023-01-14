/* Timestepping */
var t : real = 0;
var dt : real = 900;
var Nt : int = 50000;
var diag_freq : int = 1000;
var chkpt_freq : int = 1000;
var out_freq : int = 1000;

/* Domain, meters */
param Lx : real = 1.024e6;
param Ly : real = 1.024e6;
param Htot : real = 4e3;

/* Grid: Horizontal part must be divisible by 3. */
param nx : int = 4;
param ny : int = 8;
//param nx : int  = (3*128)/2;
//param ny : int  = (3*128)/2;
param nz : int  = 2;

param nx3p : int = (nx/3)+1;
param ny3p : int = (ny/3)+1;
param nx2p : int = (nx/2)+1;
param ny2p : int = (ny/2)+1;
param nx3p2 : int = (2*(nx/3)+1);

/* Coriolis coefficients */
param f0 : real = 1e-4;
param beta : real = 2e-11;

/* Gravity, reference density */
param g : real = 9.81;
param rho0 : real = 1028.8;

/* Viscosities & Ekman friction: r = r0*Htot/H(nz) fixes barotropic drag;r0=fdE/(2Htot) */
param A2 : real = 25;
param A8 : real = 3e22;
param r0 : real = 5e-8;

/* Quadratic drag parameter. C_d = c_d/h_BL where h_BL is the bottom
   boundary layer thickness. This is the drag felt by the barotropic mode,
   not by the bottom layer. C_d*L_d is a crude measure of the strength of the drag.
   > 1 is strong drag, < 1 is weak drag. */
param C_d : real = 1e-5;

/* Time step error tolerance; max-norm vorticity */
param TOL : real = 1e-8;

/* CFL for the fastest Rossby wave = 1 / (10 * frequency) */
param dt_max : real = 0.2 * pi / (beta * Lx);
