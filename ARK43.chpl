use parameters;
use domains;
use arrays;
use QGN_Module;
use FFT_utils;

use compare_fortran;
use Time;

var N1, N2, N3, N4, N5, N6 : [D3_hat] complex;
var L1, L2, L3, L4, L5, L6 : [D3_hat] complex;
var q_tmp : [D3_hat] complex;
var Mq : [D3_hat] complex;
var k8 : [D_hat] real;

var err : [D3] real;
var err_hat : [D3_hat] complex;

/* Current and previous errors for adaptive time step */
  var err0, err1 : real;

/* Explicit & implicit RK coefficients */
  var ae, ai : [{1..6,1..6}] real;

/* Error coefficients */
  var b, be : [{1..6}] real;

  var reject : bool;

/* TODO: Check this */
  k8 = k2**4;

  err0 = TOL;
  ae[1,1] = 0;
  ae[2,1] = 0.5;
  ae[3,1] = 13861.0/62500.0;
  ae[3,2] = 6889.0/62500.0;
  ae[4,1] = -116923316275.0/2393684061468.0;
  ae[4,2] = -2731218467317.0/15368042101831.0;
  ae[4,3] = 9408046702089.0/11113171139209.0;
  ae[5,1] = -451086348788.0/2902428689909.0;
  ae[5,2] = -2682348792572.0/7519795681897.0;
  ae[5,3] = 12662868775082.0/11960479115383.0;
  ae[5,4] = 3355817975965.0/11060851509271.0;
  ae[6,1] = 647845179188.0/3216320057751.0;
  ae[6,2] = 73281519250.0/8382639484533.0;
  ae[6,3] = 552539513391.0/3454668386233.0;
  ae[6,4] = 3354512671639.0/8306763924573.0;
  ae[6,5] = 4040.0/17871.0;

  ai[1,1] = 0;
  ai[2,1] = 0.25;
  ai[2,2] = 0.25;
  ai[3,1] = 8611.0/62500.0;
  ai[3,2] = -1743.0/31250.0;
  ai[3,3] = 0.25;
  ai[4,1] = 5012029.0/34652500.0;
  ai[4,2] = -654441.0/2922500.0;
  ai[4,3] = 174375.0/388108.0;
  ai[4,4] = 0.25;
  ai[5,1] = 15267082809.0/155376265600.0;
  ai[5,2] = -71443401.0/120774400.0;
  ai[5,3] = 730878875.0/902184768.0;
  ai[5,4] = 2285395.0/8070912.0;
  ai[5,5] = 0.25;
  ai[6,1] = 82889.0/524892.0;
  ai[6,2] = 0;
  ai[6,3] = 15625.0/83664.0;
  ai[6,4] = 69875.0/102672.0;
  ai[6,5] = -2260.0/8211.0;
  ai[6,6] = 0.25;

  b[1] = 82889.0/524892.0;
  b[2] = 0;
  b[3] = 15625.0/83664.0;
  b[4] = 69875.0/102672.0;
  b[5] = -2260.0/8211.0;
  b[6] = 0.25;

  be[1] = 31666707.0/9881966720.0;
  be[2] = 0;
  be[3] = -256875.0/105007616.0;
  be[4] = -2768025.0/128864768.0;
  be[5] = 169839.0/3864644.0;
  be[6] = -5247.0/225920.0;

proc TimeStep() {

  /* First RK stage, t=0 */
    GetRHS(q_hat, N1);
    forall (i,j,k) in D3_hat {
      L1[i,j,k] = -(A2*k2[j,k] + A8*k8[j,k])*q_hat[i,j,k];
    }

  do {
    forall (i,j,k) in D3_hat {
      Mq[i,j,k] = 1.0/(1.0 + 0.25*dt*(A2*k2[j,k]+A8*k8[j,k]));
    }

    var t1 : Timer;
    var t2 : Timer;
    var t3 : Timer;
    var t4 : Timer;

    t1.start();
    /* Second RK stage */
      t2.start();
      q_tmp = Mq*(q_hat + dt*(ae[2,1]*N1+ai[2,1]*L1));
      t2.stop();

      t3.start();
      GetRHS(q_tmp,N2);
      t3.stop();

      t4.start();
      forall (i,j,k) in D3_hat {
        L2[i,j,k] = -(A2*k2[j,k] + A8*k8[j,k])*q_tmp[i,j,k];
      }
      t4.stop();
    t1.stop();

    writeln();
    writeln("Total RK stage: ", t1.elapsed());
    writeln("q_tmp: ", t2.elapsed());
    writeln("GetRHS: ", t3.elapsed());
    writeln("L2: ", t4.elapsed());
    writeln();

    /* Third RK stage */
      q_tmp = Mq*(q_hat + dt*(ae[3,1]*N1+ae[3,2]*N2
                             +ai[3,1]*L1+ai[3,2]*L2));
      GetRHS(q_tmp,N3);
      forall (i,j,k) in D3_hat {
        L3[i,j,k] = -(A2*k2[j,k] + A8*k8[j,k])*q_tmp[i,j,k];
      }

    /* Fourth RK stage */
      q_tmp = Mq*(q_hat + dt*(ae[4,1]*N1+ae[4,2]*N2+ae[4,3]*N3
                             +ai[4,1]*L1+ai[4,2]*L2+ai[4,3]*L3));
      GetRHS(q_tmp,N4);
      forall (i,j,k) in D3_hat {
        L4[i,j,k] = -(A2*k2[j,k] + A8*k8[j,k])*q_tmp[i,j,k];
      }

    /* Fifth RK stage */
      q_tmp = Mq*(q_hat+dt*(ae[5,1]*N1+ae[5,2]*N2+ae[5,3]*N3+ae[5,4]*N4
                           +ai[5,1]*L1+ai[5,2]*L2+ai[5,3]*L3+ai[5,4]*L4));
      GetRHS(q_tmp,N5);
      forall (i,j,k) in D3_hat {
        L5[i,j,k] = -(A2*k2[j,k] + A8*k8[j,k])*q_tmp[i,j,k];
      }

    /* Sixth RK stage */
      q_tmp = Mq*(q_hat + dt*(ae[6,1]*N1+ae[6,2]*N2+ae[6,3]*N3
                             +ae[6,4]*N4+ae[6,5]*N5+ai[6,1]*L1
                             +ai[6,2]*L2+ai[6,3]*L3+ai[6,4]*L4+ai[6,5]*L5));
      GetRHS(q_tmp,N6);
      forall (i,j,k) in D3_hat {
        L6[i,j,k] = -(A2*k2[j,k] + A8*k8[j,k])*q_tmp[i,j,k];
      }

    /* Error control */
      err_hat = be[1]*(N1+L1)+be[3]*(N3+L3)
               +be[4]*(N4+L4)+be[5]*(N5+L5)+be[6]*(N6+L6);
      execute_backward_FFTs(err_hat, err);
      normalize(err);
      err1 = dt*(max reduce (abs(err)));

      if (err1 > TOL) {
        dt = 0.75*dt;
        reject = true;
      }
      else {
        /* Compute update */
        q_hat = q_hat + dt*(b[1]*(N1+L1)+b[3]*(N3+L3)
                           +b[4]*(N4+L4)+b[5]*(N5+L5)
                           +b[6]*(N6+L6));
        forall i in zl {
          q_hat[i,1,1] = 0;
        }
        t = t + dt;

        /* Stepsize adjustment PI.3.4, divide by 4 for 4th order method with 3rd embedded */
        dt = min(dt*((0.75*TOL/err1)**0.075)*((err0/err1)**0.1),dt_max);
        err0 = err1;
        reject = false;
      }
  } while reject;

}
