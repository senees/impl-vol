use lazy_static::lazy_static;

pub const VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_BELOW_INTRINSIC: f64 = f64::MIN;
pub const VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_ABOVE_MAXIMUM: f64 = f64::MAX;
pub const DBL_MIN: f64 = f64::MIN_POSITIVE;
pub const DBL_MAX: f64 = f64::MAX;
pub const DBL_EPSILON: f64 = f64::EPSILON;

lazy_static! {
  pub static ref SQRT_DBL_EPSILON: f64 = DBL_EPSILON.sqrt();
  pub static ref FOURTH_ROOT_DBL_EPSILON: f64 = SQRT_DBL_EPSILON.sqrt();
  pub static ref EIGHTH_ROOT_DBL_EPSILON: f64 = FOURTH_ROOT_DBL_EPSILON.sqrt();
  pub static ref SIXTEENTH_ROOT_DBL_EPSILON: f64 = EIGHTH_ROOT_DBL_EPSILON.sqrt();
  pub static ref SQRT_DBL_MIN: f64 = DBL_MIN.sqrt();
  pub static ref SQRT_DBL_MAX: f64 = DBL_MAX.sqrt();
  pub static ref ASYMPTOTIC_EXPANSION_ACCURACY_THRESHOLD: f64 = -10.0;
  pub static ref SMALL_T_EXPANSION_OF_NORMALISED_BLACK_THRESHOLD: f64 = 2.0 * *SIXTEENTH_ROOT_DBL_EPSILON;
}

#[inline(always)]
pub fn is_zero(x: f64) -> bool {
  x.abs() < DBL_MIN
}

#[inline(always)]
pub fn max(x: f64, y: f64) -> f64 {
  if x >= y {
    x
  } else {
    y
  }
}

#[inline(always)]
pub fn exp(x: f64) -> f64 {
  x.exp()
}

#[inline(always)]
pub fn pow(x: f64, n: f64) -> f64 {
  x.powf(n)
}

#[inline(always)]
pub fn log(x: f64) -> f64 {
  x.ln()
}

#[inline(always)]
pub fn sqrt(x: f64) -> f64 {
  x.sqrt()
}

#[inline(always)]
pub fn square(x: f64) -> f64 {
  x * x
}

#[inline(always)]
pub fn fabs(x: f64) -> f64 {
  x.abs()
}

#[inline(always)]
pub fn sel(c: bool, x: f64, y: f64) -> f64 {
  if c {
    x
  } else {
    y
  }
}

#[inline(always)]
pub fn d_int(x: f64) -> f64 {
  if x > 0.0 {
    x.floor()
  } else {
    -((-x).floor())
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn test_max() {
    assert_eq!(10.0, max(10.0, 9.9999));
    assert_eq!(10.0, max(10.0, 10.0));
    assert_eq!(10.0001, max(10.0, 10.0001));
  }

  #[test]
  fn test_exp() {
    assert_eq!((0.0_f64).exp(), exp(0.0));
    assert_eq!((1.0_f64).exp(), exp(1.0));
    assert_eq!((-1.0_f64).exp(), exp(-1.0));
  }

  #[test]
  fn test_sqrt() {
    assert_eq!((0.0_f64).sqrt(), sqrt(0.0));
    assert_eq!((1.0_f64).sqrt(), sqrt(1.0));
    assert_eq!((12.54_f64).sqrt(), sqrt(12.54));
  }

  #[test]
  fn test_fabs() {
    assert_eq!(0.0, fabs(0.0));
    assert_eq!(1.0, fabs(1.0));
    assert_eq!(1.0, fabs(-1.0));
  }

  #[test]
  fn test_sel() {
    assert_eq!(1.0, sel(true, 1.0, 2.0));
    assert_eq!(2.0, sel(false, 1.0, 2.0));
  }

  #[test]
  fn test_d_int() {
    assert_eq!(2.0, d_int(2.0));
    assert_eq!(2.0, d_int(2.1));
    assert_eq!(2.0, d_int(2.4));
    assert_eq!(2.0, d_int(2.5));
    assert_eq!(2.0, d_int(2.8));
    assert_eq!(2.0, d_int(2.9));
    assert_eq!(-2.0, d_int(-2.0));
    assert_eq!(-2.0, d_int(-2.1));
    assert_eq!(-2.0, d_int(-2.4));
    assert_eq!(-2.0, d_int(-2.5));
    assert_eq!(-2.0, d_int(-2.8));
    assert_eq!(-2.0, d_int(-2.9));
  }
}
