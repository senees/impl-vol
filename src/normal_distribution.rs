use crate::definitions::*;
use crate::erf_cody::*;
use lazy_static::lazy_static;

lazy_static! {
   ///``` text
   /// The asymptotic expansion  Φ(z) = φ(z)/|z|·[1-1/z^2+...],  Abramowitz & Stegun (26.2.12), suffices for Φ(z) to have
   /// relative accuracy of 1.64E-16 for z<=-10 with 17 terms inside the square brackets (not counting the leading 1).
   /// This translates to a maximum of about 9 iterations below, which is competitive with a call to erfc() and never
   /// less accurate when z<=-10. Note that, as mentioned in section 4 (and discussion of figures 2 and 3) of George
   /// Marsaglia's article "Evaluating the Normal Distribution" (available at http://www.jstatsoft.org/v11/a05/paper),
   /// for values of x approaching -8 and below, the error of any cumulative normal function is actually dominated by
   /// the hardware (or compiler implementation) accuracy of exp(-x²/2) which is not reliably more than 14 digits when
   /// x becomes large. Still, we should switch to the asymptotic only when it is beneficial to do so.
   /// ```
   static ref NORM_CDF_ASYMPTOTIC_EXPANSION_FIRST_THRESHOLD: f64 = -10.0;
   static ref NORM_CDF_ASYMPTOTIC_EXPANSION_SECOND_THRESHOLD: f64 = -1.0/sqrt(DBL_EPSILON);
}

pub const ONE_OVER_SQRT_TWO: f64 = std::f64::consts::FRAC_1_SQRT_2;
#[allow(clippy::excessive_precision)]
pub const ONE_OVER_SQRT_TWO_PI: f64 = 0.3989422804014326779399460599343818684758586311649;
#[allow(clippy::excessive_precision)]
pub const SQRT_TWO_PI: f64 = 2.506628274631000502415765284811045253006986740610;

pub fn norm_cdf(z: f64) -> f64 {
  if z <= *NORM_CDF_ASYMPTOTIC_EXPANSION_FIRST_THRESHOLD {
    // Asymptotic expansion for very negative z following (26.2.12) on page 408
    // in M. Abramowitz and A. Stegun, Pocketbook of Mathematical Functions, ISBN 3-87144818-4.
    let mut sum = 1.0;
    if z >= *NORM_CDF_ASYMPTOTIC_EXPANSION_SECOND_THRESHOLD {
      let z_square = z * z;
      let mut i = 1;
      let mut g = 1.0;
      let mut a = DBL_MAX;
      let mut last_a;
      loop {
        last_a = a;
        let x = (4 * i - 3) as f64 / z_square;
        let y = x * ((4 * i - 1) as f64 / z_square);
        a = g * (x - y);
        sum -= a;
        g *= y;
        i += 1;
        a = fabs(a);
        if !(last_a > a && a >= fabs(sum * DBL_EPSILON)) {
          break;
        };
      }
      return -norm_pdf(z) * sum / z;
    }
  }
  0.5 * erfc_cody(-z * ONE_OVER_SQRT_TWO)
}

#[allow(clippy::excessive_precision)]
pub fn inverse_norm_cdf(u: f64) -> f64 {
  //
  // ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
  //
  // Produces the normal deviate Z corresponding to a given lower
  // tail area of u; Z is accurate to about 1 part in 10**16.
  // see http://lib.stat.cmu.edu/apstat/241
  //
  const SPLIT1: f64 = 0.425;
  const SPLIT2: f64 = 5.0;
  const CONST1: f64 = 0.180625;
  const CONST2: f64 = 1.6;

  // Coefficients for P close to 0.5
  const A0: f64 = 3.3871328727963666080E0;
  const A1: f64 = 1.3314166789178437745E+2;
  const A2: f64 = 1.9715909503065514427E+3;
  const A3: f64 = 1.3731693765509461125E+4;
  const A4: f64 = 4.5921953931549871457E+4;
  const A5: f64 = 6.7265770927008700853E+4;
  const A6: f64 = 3.3430575583588128105E+4;
  const A7: f64 = 2.5090809287301226727E+3;
  const B1: f64 = 4.2313330701600911252E+1;
  const B2: f64 = 6.8718700749205790830E+2;
  const B3: f64 = 5.3941960214247511077E+3;
  const B4: f64 = 2.1213794301586595867E+4;
  const B5: f64 = 3.9307895800092710610E+4;
  const B6: f64 = 2.8729085735721942674E+4;
  const B7: f64 = 5.2264952788528545610E+3;
  // Coefficients for P not close to 0, 0.5 or 1.
  const C0: f64 = 1.42343711074968357734E0;
  const C1: f64 = 4.63033784615654529590E0;
  const C2: f64 = 5.76949722146069140550E0;
  const C3: f64 = 3.64784832476320460504E0;
  const C4: f64 = 1.27045825245236838258E0;
  const C5: f64 = 2.41780725177450611770E-1;
  const C6: f64 = 2.27238449892691845833E-2;
  const C7: f64 = 7.74545014278341407640E-4;
  const D1: f64 = 2.05319162663775882187E0;
  const D2: f64 = 1.67638483018380384940E0;
  const D3: f64 = 6.89767334985100004550E-1;
  const D4: f64 = 1.48103976427480074590E-1;
  const D5: f64 = 1.51986665636164571966E-2;
  const D6: f64 = 5.47593808499534494600E-4;
  const D7: f64 = 1.05075007164441684324E-9;
  // Coefficients for P very close to 0 or 1
  const E0: f64 = 6.65790464350110377720E0;
  const E1: f64 = 5.46378491116411436990E0;
  const E2: f64 = 1.78482653991729133580E0;
  const E3: f64 = 2.96560571828504891230E-1;
  const E4: f64 = 2.65321895265761230930E-2;
  const E5: f64 = 1.24266094738807843860E-3;
  const E6: f64 = 2.71155556874348757815E-5;
  const E7: f64 = 2.01033439929228813265E-7;
  const F1: f64 = 5.99832206555887937690E-1;
  const F2: f64 = 1.36929880922735805310E-1;
  const F3: f64 = 1.48753612908506148525E-2;
  const F4: f64 = 7.86869131145613259100E-4;
  const F5: f64 = 1.84631831751005468180E-5;
  const F6: f64 = 1.42151175831644588870E-7;
  const F7: f64 = 2.044_263_103_389_939_785_64E-15;

  if u <= 0.0 {
    return log(u);
  }
  if u >= 1.0 {
    return log(1.0 - u);
  }

  let q = u - 0.5;
  if fabs(q) <= SPLIT1 {
    let r = CONST1 - q * q;
    q * (((((((A7 * r + A6) * r + A5) * r + A4) * r + A3) * r + A2) * r + A1) * r + A0) / (((((((B7 * r + B6) * r + B5) * r + B4) * r + B3) * r + B2) * r + B1) * r + 1.0)
  } else {
    let mut r = sel(q < 0.0, u, 1.0 - u);
    r = sqrt(-log(r));

    let ret = if r < SPLIT2 {
      r -= CONST2;
      (((((((C7 * r + C6) * r + C5) * r + C4) * r + C3) * r + C2) * r + C1) * r + C0) / (((((((D7 * r + D6) * r + D5) * r + D4) * r + D3) * r + D2) * r + D1) * r + 1.0)
    } else {
      r -= SPLIT2;
      (((((((E7 * r + E6) * r + E5) * r + E4) * r + E3) * r + E2) * r + E1) * r + E0) / (((((((F7 * r + F6) * r + F5) * r + F4) * r + F3) * r + F2) * r + F1) * r + 1.0)
    };
    sel(q < 0.0, -ret, ret)
  }
}

pub fn norm_pdf(x: f64) -> f64 {
  ONE_OVER_SQRT_TWO_PI * exp(-0.5 * x * x)
}
