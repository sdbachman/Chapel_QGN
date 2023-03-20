/* Restart? */
config const restart : bool = false;

/* Background settings */
config const background_file = "background_eady.nc";

/* Timestepping */
config var Nt_start : int = 0;
config var Nt : int = 10000;
config var t : real = 0;
config var dt : real = 900;

/* Domain, meters */
config const Lx : real = 3.84e5;
config const Ly : real = 3.84e5;
config const Htot : real = 5e3;

/* Grid: Horizontal part must be divisible by 3. */
config const nx : int  = (3*128)/2;
config const ny : int  = (3*128)/2;
config const nz : int  = 8;  // 50

var nx3p : int = (nx/3)+1;
var ny3p : int = (ny/3)+1;
var nx2p : int = (nx/2)+1;
var ny2p : int = (ny/2)+1;
var ny3p2 : int = (2*(ny/3)+1);
var nz1m : int = nz-1;

/* Coriolis coefficients */
config const f0 : real = 1e-4;
config const beta : real = 2e-11;

/* Gravity, reference density */
config const g : real = 9.81;
config const rho0 : real = 1028.8;

/* Viscosities & Ekman friction: r = r0*Htot/H(nz) fixes barotropic drag;r0=fdE/(2Htot) */
config const A2 : real = 25;
config const A8 : real = 3e22;
config const r0 : real = 5e-8;

/* Quadratic drag parameter. C_d = c_d/h_BL where h_BL is the bottom
   boundary layer thickness. This is the drag felt by the barotropic mode,
   not by the bottom layer. C_d*L_d is a crude measure of the strength of the drag.
   > 1 is strong drag, < 1 is weak drag. */
config const C_d : real = 5e-5;

/* Time step error tolerance; max-norm vorticity */
param TOL : real = 1e-8;

/* CFL for the fastest Rossby wave = 1 / (10 * frequency) */
var dt_max : real = 0.2 * pi / (beta * Lx);
