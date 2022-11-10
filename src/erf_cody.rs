//!
//! Original Fortran code taken from http://www.netlib.org/specfun/erf, compiled with f2c, and adapted by hand.
//!
//! Created with command line f2c -C++ -c -a -krd -r8 cody_erf.f
//!
//! Translated by f2c (version 20100827).
//!
//!

use crate::definitions::*;

/// SUBROUTINE CALERF(ARG,RESULT,JINT)
#[allow(clippy::excessive_precision)]
fn calerf(x: f64, jint: i64) -> f64 {
  let a = [3.1611237438705656, 113.864154151050156, 377.485237685302021, 3209.37758913846947, 0.185777706184603153];
  let b = [23.6012909523441209, 244.024637934444173, 1282.61652607737228, 2844.23683343917062];
  let c__ = [
    0.564188496988670089,
    8.88314979438837594,
    66.1191906371416295,
    298.635138197400131,
    881.95222124176909,
    1712.04761263407058,
    2051.07837782607147,
    1230.33935479799725,
    2.15311535474403846e-8,
  ];
  let d__ = [
    15.7449261107098347,
    117.693950891312499,
    537.181101862009858,
    1621.38957456669019,
    3290.79923573345963,
    4362.61909014324716,
    3439.36767414372164,
    1230.33935480374942,
  ];
  let p = [
    0.305326634961232344,
    0.360344899949804439,
    0.125781726111229246,
    0.0160837851487422766,
    6.58749161529837803e-4,
    0.0163153871373020978,
  ];
  let q = [
    2.56852019228982242,
    1.87295284992346047,
    0.527905102951428412,
    0.0605183413124413191,
    0.00233520497626869185,
  ];

  const ZERO: f64 = 0.0;
  const HALF: f64 = 0.5;
  const ONE: f64 = 1.0;
  const TWO: f64 = 2.0;
  const FOUR: f64 = 4.0;
  const SQRPI: f64 = 0.56418958354775628695;
  const THRESH: f64 = 0.46875;
  const SIXTEN: f64 = 16.0;

  const XINF: f64 = 1.79e308;
  const XNEG: f64 = -26.628;
  const XSMALL: f64 = 1.11e-16;
  const XBIG: f64 = 26.543;
  const XHUGE: f64 = 6.71e7;
  const XMAX: f64 = 2.53e307;

  let fix_negative_result = |mut result: f64| {
    if jint == 0 {
      result = (HALF - result) + HALF;
      if x < ZERO {
        result = -result;
      }
    } else if jint == 1 {
      if x < ZERO {
        result = TWO - result;
      }
    } else if x < ZERO {
      if x < XNEG {
        result = XINF;
      } else {
        let d_1 = x * SIXTEN;
        let ysq = d_int(d_1) / SIXTEN;
        let del = (x - ysq) * (x + ysq);
        let y = exp(ysq * ysq) * exp(del);
        result = y + y - result;
      }
    }
    result
  };

  let y = fabs(x);
  if y <= THRESH {
    let mut ysq = ZERO;
    if y > XSMALL {
      ysq = y * y;
    }
    let mut xnum = a[4] * ysq;
    let mut xden = ysq;
    for i in 1..=3 {
      xnum = (xnum + a[i - 1]) * ysq;
      xden = (xden + b[i - 1]) * ysq;
    }
    let mut result = x * (xnum + a[3]) / (xden + b[3]);
    if jint != 0 {
      result = ONE - result;
    }
    if jint == 2 {
      result *= exp(ysq);
    }
    result
  } else if y <= FOUR {
    let mut xnum = c__[8] * y;
    let mut xden = y;
    for i in 1..=7 {
      xnum = (xnum + c__[i - 1]) * y;
      xden = (xden + d__[i - 1]) * y;
    }
    let mut result = (xnum + c__[7]) / (xden + d__[7]);
    if jint != 2 {
      let mut d_1 = y * SIXTEN;
      let ysq = d_int(d_1) / SIXTEN;
      let del = (y - ysq) * (y + ysq);
      d_1 = exp(-ysq * ysq) * exp(-del);
      result *= d_1;
    }
    fix_negative_result(result)
  } else {
    let mut result = ZERO;
    if y >= XBIG {
      if jint != 2 || y >= XMAX {
        return fix_negative_result(result);
      }
      if y >= XHUGE {
        result = SQRPI / y;
        return fix_negative_result(result);
      }
    }
    let ysq = ONE / (y * y);
    let mut xnum = p[5] * ysq;
    let mut xden = ysq;
    for i in 1..=4 {
      xnum = (xnum + p[i - 1]) * ysq;
      xden = (xden + q[i - 1]) * ysq;
    }
    result = ysq * (xnum + p[4]) / (xden + q[4]);
    result = (SQRPI - result) / y;
    if jint != 2 {
      let mut d_1 = y * SIXTEN;
      let ysq = d_int(d_1) / SIXTEN;
      let del = (y - ysq) * (y + ysq);
      d_1 = exp(-ysq * ysq) * exp(-del);
      result *= d_1;
    }
    fix_negative_result(result)
  }
}

pub fn erf_cody(x: f64) -> f64 {
  calerf(x, 0)
}

pub fn erfc_cody(x: f64) -> f64 {
  calerf(x, 1)
}

pub fn erfcx_cody(x: f64) -> f64 {
  calerf(x, 2)
}
