use crate::definitions::*;
use crate::erf_cody::*;
use crate::normal_distribution::*;
use crate::rational_cubic::*;
use std::sync::atomic::{AtomicUsize, Ordering};

#[allow(clippy::excessive_precision)]
const TWO_PI: f64 = 6.283185307179586476925286766559005768394338798750;
#[allow(clippy::excessive_precision)]
const SQRT_THREE: f64 = 1.732050807568877293527446341505872366942805253810;
#[allow(clippy::excessive_precision)]
const SQRT_PI_OVER_TWO: f64 = 1.253314137315500251207882642405522626503493370305; // sqrt(pi/2) to avoid misinterpretation.
#[allow(clippy::excessive_precision)]
const SQRT_ONE_OVER_THREE: f64 = 0.577350269189625764509148780501957455647601751270;
#[allow(clippy::excessive_precision)]
const TWO_PI_OVER_SQRT_TWENTY_SEVEN: f64 = 1.209199576156145233729385505094770488189377498728; // 2*pi/sqrt(27)
const PI_OVER_SIX: f64 = std::f64::consts::FRAC_PI_6;

/// Set this to 0 if you want positive results for (positive) denormalised inputs, else to DBL_MIN.
/// Note that you cannot achieve full machine accuracy from denormalised inputs!
const DENORMALISATION_CUTOFF: f64 = 0.0;

/// (DBL_DIG*20)/3 ≈ 100.
/// Only needed when the iteration effectively alternates Householder/Halley/Newton steps
/// and binary nesting due to roundoff truncation.
static IMPLIED_VOLATILITY_MAXIMUM_ITERATIONS: AtomicUsize = AtomicUsize::new(2);

fn get_implied_volatility_maximum_iterations() -> usize {
  IMPLIED_VOLATILITY_MAXIMUM_ITERATIONS.load(Ordering::Relaxed)
}

#[cfg(feature = "ENABLE_SWITCHING_THE_OUTPUT_TO_ITERATION_COUNT")]
static IMPLIED_VOLATILITY_OUTPUT_TYPE: AtomicUsize = AtomicUsize::new(0);

#[cfg(feature = "ENABLE_SWITCHING_THE_OUTPUT_TO_ITERATION_COUNT")]
fn get_implied_volatility_output_type() -> usize {
  IMPLIED_VOLATILITY_OUTPUT_TYPE.load(Ordering::Relaxed)
}

#[cfg(feature = "ENABLE_SWITCHING_THE_OUTPUT_TO_ITERATION_COUNT")]
fn implied_volatility_output(count: usize, volatility: f64) -> f64 {
  sel(get_implied_volatility_output_type() > 0, count as f64, volatility)
}
#[cfg(not(feature = "ENABLE_SWITCHING_THE_OUTPUT_TO_ITERATION_COUNT"))]
fn implied_volatility_output(count: usize, volatility: f64) -> f64 {
  volatility
}

#[cfg(feature = "DO_NOT_OPTIMISE_NORMALISED_BLACK_IN_REGIONS_3_AND_4_FOR_CODYS_FUNCTIONS")]
const DO_NOT_OPTIMISE_NORMALISED_BLACK_IN_REGIONS_3_AND_4_FOR_CODYS_FUNCTIONS: bool = true;
#[cfg(not(feature = "DO_NOT_OPTIMISE_NORMALISED_BLACK_IN_REGIONS_3_AND_4_FOR_CODYS_FUNCTIONS"))]
const DO_NOT_OPTIMISE_NORMALISED_BLACK_IN_REGIONS_3_AND_4_FOR_CODYS_FUNCTIONS: bool = false;

#[cfg(feature = "ENABLE_CHANGING_THE_HOUSEHOLDER_METHOD_ORDER")]
static IMPLIED_VOLATILITY_HOUSEHOLDER_METHOD_ORDER: AtomicUsize = AtomicUsize::new(4);

#[cfg(feature = "ENABLE_CHANGING_THE_HOUSEHOLDER_METHOD_ORDER")]
fn get_implied_volatility_householder_method_order() -> usize {
  IMPLIED_VOLATILITY_HOUSEHOLDER_METHOD_ORDER.load(Ordering::Relaxed)
}

#[cfg(feature = "ENABLE_CHANGING_THE_HOUSEHOLDER_METHOD_ORDER")]
fn householder_factor(newton: f64, halley: f64, hh3: f64) -> f64 {
  sel(
    get_implied_volatility_householder_method_order() > 3,
    (1.0 + 0.5 * halley * newton) / (1.0 + newton * (halley + hh3 * newton / 6.0)),
    sel(get_implied_volatility_householder_method_order() > 2, 1.0 / (1.0 + 0.5 * halley * newton), 1.0),
  )
}

#[cfg(not(feature = "ENABLE_CHANGING_THE_HOUSEHOLDER_METHOD_ORDER"))]
fn householder_factor(newton: f64, halley: f64, hh3: f64) -> f64 {
  (1 + 0.5 * halley * newton) / (1 + newton * (halley + hh3 * newton / 6))
}

///```text
/// Asymptotic expansion of
///
///              b  =  Φ(h+t)·exp(x/2) - Φ(h-t)·exp(-x/2)
/// with
///              h  =  x/s   and   t  =  s/2
/// which makes
///              b  =  Φ(h+t)·exp(h·t) - Φ(h-t)·exp(-h·t)
///
///                    exp(-(h²+t²)/2)
///                 =  ---------------  ·  [ Y(h+t) - Y(h-t) ]
///                        √(2π)
/// with
///           Y(z) := Φ(z)/φ(z)
///
/// for large negative (t-|h|) by the aid of Abramowitz & Stegun (26.2.12) where Φ(z) = φ(z)/|z|·[1-1/z^2+...].
/// We define
///                     r
///         A(h,t) :=  --- · [ Y(h+t) - Y(h-t) ]
///                     t
///
/// with r := (h+t)·(h-t) and give an expansion for A(h,t) in q:=(h/r)² expressed in terms of e:=(t/h)² .
/// ```
#[rustfmt::skip]
fn asymptotic_expansion_of_normalised_black_call(h: f64, t: f64) -> f64 {
  let e = (t / h) * (t / h);
  let r = (h + t) * (h - t);
  let q = (h / r) * (h / r);
  // 17th order asymptotic expansion of A(h,t) in q, sufficient for Φ(h) [and thus y(h)] to have relative accuracy of 1.64E-16 for h <= η  with  η:=-10.
  let asymptotic_expansion_sum = 2.0+q*(-6.0E0-2.0*e+3.0*q*(1.0E1+e*(2.0E1+2.0*e)+5.0*q*(-1.4E1+e*(-7.0E1+e*(-4.2E1-2.0*e))+7.0*q*(1.8E1+e*(1.68E2+e*(2.52E2+e*(7.2E1+2.0*e)))+9.0*q*(-2.2E1+e*(-3.3E2+e*(-9.24E2+e*(-6.6E2+e*(-1.1E2-2.0*e))))+1.1E1*q*(2.6E1+e*(5.72E2+e*(2.574E3+e*(3.432E3+e*(1.43E3+e*(1.56E2+2.0*e)))))+1.3E1*q*(-3.0E1+e*(-9.1E2+e*(-6.006E3+e*(-1.287E4+e*(-1.001E4+e*(-2.73E3+e*(-2.1E2-2.0*e))))))+1.5E1*q*(3.4E1+e*(1.36E3+e*(1.2376E4+e*(3.8896E4+e*(4.862E4+e*(2.4752E4+e*(4.76E3+e*(2.72E2+2.0*e)))))))+1.7E1*q*(-3.8E1+e*(-1.938E3+e*(-2.3256E4+e*(-1.00776E5+e*(-1.84756E5+e*(-1.51164E5+e*(-5.4264E4+e*(-7.752E3+e*(-3.42E2-2.0*e))))))))+1.9E1*q*(4.2E1+e*(2.66E3+e*(4.0698E4+e*(2.3256E5+e*(5.8786E5+e*(7.05432E5+e*(4.0698E5+e*(1.08528E5+e*(1.197E4+e*(4.2E2+2.0*e)))))))))+2.1E1*q*(-4.6E1+e*(-3.542E3+e*(-6.7298E4+e*(-4.90314E5+e*(-1.63438E6+e*(-2.704156E6+e*(-2.288132E6+e*(-9.80628E5+e*(-2.01894E5+e*(-1.771E4+e*(-5.06E2-2.0*e))))))))))+2.3E1*q*(5.0E1+e*(4.6E3+e*(1.0626E5+e*(9.614E5+e*(4.08595E6+e*(8.9148E6+e*(1.04006E7+e*(6.53752E6+e*(2.16315E6+e*(3.542E5+e*(2.53E4+e*(6.0E2+2.0*e)))))))))))+2.5E1*q*(-5.4E1+e*(-5.85E3+e*(-1.6146E5+e*(-1.77606E6+e*(-9.37365E6+e*(-2.607579E7+e*(-4.01166E7+e*(-3.476772E7+e*(-1.687257E7+e*(-4.44015E6+e*(-5.9202E5+e*(-3.51E4+e*(-7.02E2-2.0*e))))))))))))+2.7E1*q*(5.8E1+e*(7.308E3+e*(2.3751E5+e*(3.12156E6+e*(2.003001E7+e*(6.919458E7+e*(1.3572783E8+e*(1.5511752E8+e*(1.0379187E8+e*(4.006002E7+e*(8.58429E6+e*(9.5004E5+e*(4.7502E4+e*(8.12E2+2.0*e)))))))))))))+2.9E1*q*(-6.2E1+e*(-8.99E3+e*(-3.39822E5+e*(-5.25915E6+e*(-4.032015E7+e*(-1.6934463E8+e*(-4.1250615E8+e*(-6.0108039E8+e*(-5.3036505E8+e*(-2.8224105E8+e*(-8.870433E7+e*(-1.577745E7+e*(-1.472562E6+e*(-6.293E4+e*(-9.3E2-2.0*e))))))))))))))+3.1E1*q*(6.6E1+e*(1.0912E4+e*(4.74672E5+e*(8.544096E6+e*(7.71342E7+e*(3.8707344E8+e*(1.14633288E9+e*(2.07431664E9+e*(2.33360622E9+e*(1.6376184E9+e*(7.0963464E8+e*(1.8512208E8+e*(2.7768312E7+e*(2.215136E6+e*(8.184E4+e*(1.056E3+2.0*e)))))))))))))))+3.3E1*(-7.0E1+e*(-1.309E4+e*(-6.49264E5+e*(-1.344904E7+e*(-1.4121492E8+e*(-8.344518E8+e*(-2.9526756E9+e*(-6.49588632E9+e*(-9.0751353E9+e*(-8.1198579E9+e*(-4.6399188E9+e*(-1.6689036E9+e*(-3.67158792E8+e*(-4.707164E7+e*(-3.24632E6+e*(-1.0472E5+e*(-1.19E3-2.0*e)))))))))))))))))*q))))))))))))))));
  let b = ONE_OVER_SQRT_TWO_PI * exp(-0.5 * (h * h + t * t)) * (t / r) * asymptotic_expansion_sum;
  fabs(max(b, 0.0))
}

///```text
/// Calculation of
///
///              b  =  Φ(h+t)·exp(h·t) - Φ(h-t)·exp(-h·t)
///
///                    exp(-(h²+t²)/2)
///                 =  --------------- ·  [ Y(h+t) - Y(h-t) ]
///                        √(2π)
/// with
///           Y(z) := Φ(z)/φ(z)
///
/// using an expansion of Y(h+t)-Y(h-t) for small t to twelvth order in t.
/// Theoretically accurate to (better than) precision  ε = 2.23E-16  when  h<=0  and  t < τ  with  τ := 2·ε^(1/16) ≈ 0.21.
/// The main bottleneck for precision is the coefficient a:=1+h·Y(h) when |h|>1 .
/// ```
#[rustfmt::skip]
fn small_t_expansion_of_normalised_black_call(h: f64, t: f64) -> f64 {
  // Y(h) := Φ(h)/φ(h) = √(π/2)·erfcx(-h/√2)
  // a := 1+h·Y(h)  --- Note that due to h<0, and h·Y(h) -> -1 (from above) as h -> -∞, we also have that a>0 and a -> 0 as h -> -∞
  // w := t² , h2 := h²
  let a = 1.0 + h * (0.5 * SQRT_TWO_PI) * erfcx_cody(-ONE_OVER_SQRT_TWO * h);
  let w = t * t;
  let h2 = h * h;
  let expansion = 2.0*t*(a+w*((-1.0+3.0*a+a*h2)/6.0+w*((-7.0+15.0*a+h2*(-1.0+10.0*a+a*h2))/120.0+w*((-57.0+105.0*a+h2*(-18.0+105.0*a+h2*(-1.0+21.0*a+a*h2)))/5040.0+w*((-561.0+945.0*a+h2*(-285.0+1260.0*a+h2*(-33.0+378.0*a+h2*(-1.0+36.0*a+a*h2))))/362880.0+w*((-6555.0+10395.0*a+h2*(-4680.0+17325.0*a+h2*(-840.0+6930.0*a+h2*(-52.0+990.0*a+h2*(-1.0+55.0*a+a*h2)))))/39916800.0+((-89055.0+135135.0*a+h2*(-82845.0+270270.0*a+h2*(-20370.0+135135.0*a+h2*(-1926.0+25740.0*a+h2*(-75.0+2145.0*a+h2*(-1.0+78.0*a+a*h2))))))*w)/6227020800.0))))));
  let b = ONE_OVER_SQRT_TWO_PI * exp(-0.5 * (h * h + t * t)) * expansion;
  fabs(max(b, 0.0))
}

///
fn normalised_intrinsic(x: f64, q: f64 /* q=±1 */) -> f64 {
  if q * x <= 0.0 {
    return 0.0;
  }
  let x2 = x * x;
  // The factor 98 is computed from last coefficient: √√92897280 = 98.1749
  if x2 < 98.0 * *FOURTH_ROOT_DBL_EPSILON {
    return fabs(max(
      sel(q < 0.0, -1.0, 1.0) * x * (1.0 + x2 * ((1.0 / 24.0) + x2 * ((1.0 / 1920.0) + x2 * ((1.0 / 322560.0) + (1.0 / 92897280.0) * x2)))),
      0.0,
    ));
  }
  let b_max = (0.5 * x).exp();
  let one_over_b_max = 1.0 / b_max;
  fabs(max(sel(q < 0.0, -1.0, 1.0) * (b_max - one_over_b_max), 0.0))
}

///
fn normalised_intrinsic_call(x: f64) -> f64 {
  normalised_intrinsic(x, 1.0)
}

///```text
///     b(x,s)  =  Φ(x/s+s/2)·exp(x/2)  -   Φ(x/s-s/2)·exp(-x/2)
///             =  Φ(h+t)·exp(x/2)      -   Φ(h-t)·exp(-x/2)
/// with
///             h  =  x/s   and   t  =  s/2
/// ```
fn normalised_black_call_using_norm_cdf(x: f64, s: f64) -> f64 {
  let h = x / s;
  let t = 0.5 * s;
  let b_max = exp(0.5 * x);
  let b = norm_cdf(h + t) * b_max - norm_cdf(h - t) / b_max;
  fabs(max(b, 0.0))
}

///```text
/// Given h = x/s and t = s/2, the normalised Black function can be written as
///
///     b(x,s)  =  Φ(x/s+s/2)·exp(x/2)  -   Φ(x/s-s/2)·exp(-x/2)
///             =  Φ(h+t)·exp(h·t)      -   Φ(h-t)·exp(-h·t) .                     (*)
///
/// It is mentioned in section 4 (and discussion of figures 2 and 3) of George Marsaglia's article "Evaluating the
/// Normal Distribution" (available at www.jstatsoft.org/v11/a05/paper) that the error of any cumulative normal
/// function Φ(z) is dominated by the hardware (or compiler implementation) accuracy of exp(-z²/2) which is not
/// reliably more than 14 digits when z is large. The accuracy of Φ(z) typically starts coming down to 14 digits when
/// z is around -8. For the (normalised) Black function, as above in (*), this means that we are subtracting two terms
/// that are each products of terms with about 14 digits of accuracy. The net result, in each of the products, is even
/// less accuracy, and then we are taking the difference of these terms, resulting in even less accuracy. When we are
/// using the asymptotic expansion asymptotic_expansion_of_normalised_black_call() invoked in the second branch at the
/// beginning of this function, we are using only *one* exponential instead of 4, and this improves accuracy. It
/// actually improves it a bit more than you would expect from the above logic, namely, almost the full two missing
/// digits (in 64 bit IEEE floating point).  Unfortunately, going higher order in the asymptotic expansion will not
/// enable us to gain more accuracy (by extending the range in which we could use the expansion) since the asymptotic
/// expansion, being a divergent series, can never gain 16 digits of accuracy for z=-8 or just below. The best you can
/// get is about 15 digits (just), for about 35 terms in the series (26.2.12), which would result in an prohibitively
/// long expression in function asymptotic expansion asymptotic_expansion_of_normalised_black_call(). In this last branch,
/// here, we therefore take a different tack as follows.
///     The "scaled complementary error function" is defined as erfcx(z) = exp(z²)·erfc(z). Cody's implementation of this
/// function as published in "Rational Chebyshev approximations for the error function", W. J. Cody, Math. Comp., 1969, pp.
/// 631-638, uses rational functions that theoretically approximates erfcx(x) to at least 18 significant decimal digits,
/// *without* the use of the exponential function when x>4, which translates to about z<-5.66 in Φ(z). To make use of it,
/// we write
///             Φ(z) = exp(-z²/2)·erfcx(-z/√2)/2
///
/// to transform the normalised black function to
///
///   b   =  ½ · exp(-½(h²+t²)) · [ erfcx(-(h+t)/√2) -  erfcx(-(h-t)/√2) ]
///
/// which now involves only one exponential, instead of three, when |h|+|t| > 5.66 , and the difference inside the
/// square bracket is between the evaluation of two rational functions, which, typically, according to Marsaglia,
/// retains the full 16 digits of accuracy (or just a little less than that).
///```
fn normalised_black_call_using_erfcx(h: f64, t: f64) -> f64 {
  let b = 0.5 * exp(-0.5 * (h * h + t * t)) * (erfcx_cody(-ONE_OVER_SQRT_TWO * (h + t)) - erfcx_cody(-ONE_OVER_SQRT_TWO * (h - t)));
  fabs(max(b, 0.0))
}

///```text
/// Introduced on 2017-02-18
///
///     b(x,s)  =  Φ(x/s+s/2)·exp(x/2)  -   Φ(x/s-s/2)·exp(-x/2)
///             =  Φ(h+t)·exp(x/2)      -   Φ(h-t)·exp(-x/2)
///             =  ½ · exp(-u²-v²) · [ erfcx(u-v) -  erfcx(u+v) ]
///             =  ½ · [ exp(x/2)·erfc(u-v)     -  exp(-x/2)·erfc(u+v)    ]
///             =  ½ · [ exp(x/2)·erfc(u-v)     -  exp(-u²-v²)·erfcx(u+v) ]
///             =  ½ · [ exp(-u²-v²)·erfcx(u-v) -  exp(-x/2)·erfc(u+v)    ]
/// with
///              h  =  x/s ,       t  =  s/2 ,
/// and
///              u  = -h/√2  and   v  =  t/√2 .
///
/// Cody's erfc() and erfcx() functions each, for some values of their argument, involve the evaluation
/// of the exponential function exp(). The normalised Black function requires additional evaluation(s)
/// of the exponential function irrespective of which of the above formulations is used. However, the total
/// number of exponential function evaluations can be minimised by a judicious choice of one of the above
/// formulations depending on the input values and the branch logic in Cody's erfc() and erfcx().
///```
fn normalised_black_call_with_optimal_use_of_codys_functions(x: f64, s: f64) -> f64 {
  let codys_threshold = 0.46875;
  let h = x / s;
  let t = 0.5 * s;
  let q1 = -ONE_OVER_SQRT_TWO * (h + t);
  let q2 = -ONE_OVER_SQRT_TWO * (h - t);
  let two_b;
  if q1 < codys_threshold {
    if q2 < codys_threshold {
      two_b = exp(0.5 * x) * erfc_cody(q1) - exp(-0.5 * x) * erfc_cody(q2);
    } else {
      two_b = exp(0.5 * x) * erfc_cody(q1) - exp(-0.5 * (h * h + t * t)) * erfcx_cody(q2);
    }
  } else if q2 < codys_threshold {
    two_b = exp(-0.5 * (h * h + t * t)) * erfcx_cody(q1) - exp(-0.5 * x) * erfc_cody(q2);
  } else {
    two_b = exp(-0.5 * (h * h + t * t)) * (erfcx_cody(q1) - erfcx_cody(q2));
  }
  fabs(max(0.5 * two_b, 0.0))
}

///
pub fn black(f: f64, k: f64, sigma: f64, t: f64, q: f64 /* q=±1 */) -> f64 {
  let intrinsic = fabs(max(sel(q < 0.0, k - f, f - k), 0.0));
  // Map in-the-money to out-of-the-money
  if q * (f - k) > 0.0 {
    return intrinsic + black(f, k, sigma, t, -q);
  }
  max(intrinsic, (sqrt(f) * sqrt(k)) * normalised_black(log(f / k), sigma * sqrt(t), q))
}

///
pub fn normalised_black_call(x: f64, s: f64) -> f64 {
  if x > 0.0 {
    return normalised_intrinsic_call(x) + normalised_black_call(-x, s); // In the money.
  }
  if s <= fabs(x) * DENORMALISATION_CUTOFF {
    return normalised_intrinsic_call(x); // sigma=0 -> intrinsic value.
  }
  // Denote h := x/s and t := s/2.
  // We evaluate the condition |h|>|η|, i.e., h<η  &&  t < τ+|h|-|η| avoiding any divisions by s,
  // where η = ASYMPTOTIC_EXPANSION_ACCURACY_THRESHOLD  and τ = SMALL_T_EXPANSION_OF_NORMALISED_BLACK_THRESHOLD.
  let eta = *ASYMPTOTIC_EXPANSION_ACCURACY_THRESHOLD;
  let tau = *SMALL_T_EXPANSION_OF_NORMALISED_BLACK_THRESHOLD;
  if x < s * eta && 0.5 * s * s + x < s * (tau + eta) {
    return asymptotic_expansion_of_normalised_black_call(x / s, 0.5 * s);
  }
  if 0.5 * s < tau {
    return small_t_expansion_of_normalised_black_call(x / s, 0.5 * s);
  }
  if DO_NOT_OPTIMISE_NORMALISED_BLACK_IN_REGIONS_3_AND_4_FOR_CODYS_FUNCTIONS {
    // When b is more than, say, about 85% of b_max=exp(x/2), then b is dominated by the first of the two terms in the Black formula, and we retain more accuracy by not attempting to combine the two terms in any way.
    // We evaluate the condition h+t>0.85  avoiding any divisions by s.
    if x + 0.5 * s * s > s * 0.85 {
      normalised_black_call_using_norm_cdf(x, s)
    } else {
      normalised_black_call_using_erfcx(x / s, 0.5 * s)
    }
  } else {
    normalised_black_call_with_optimal_use_of_codys_functions(x, s)
  }
}

///
pub fn normalised_black(x: f64, s: f64, q: f64 /* q=±1 */) -> f64 {
  /* Reciprocal-strike call-put equivalence */
  normalised_black_call(sel(q < 0.0, -x, x), s)
}

///
pub fn normalised_vega(x: f64, s: f64) -> f64 {
  let ax = fabs(x);
  if ax <= 0.0 {
    ONE_OVER_SQRT_TWO_PI * exp(-0.125 * s * s)
  } else if s <= 0.0 || s <= ax * *SQRT_DBL_MIN {
    0.0
  } else {
    ONE_OVER_SQRT_TWO_PI * exp(-0.5 * (square(x / s) + square(0.5 * s)))
  }
}

/// This weeds out denormalised (a.k.a. 'subnormal') numbers.
fn is_below_horizon(x: f64) -> bool {
  fabs(x) < DENORMALISATION_CUTOFF
}

///
fn compute_f_lower_map_and_first_two_derivatives(x: f64, s: f64) -> (f64, f64, f64) {
  let ax = fabs(x);
  let z = SQRT_ONE_OVER_THREE * ax / s;
  let y = z * z;
  let s2 = s * s;
  let phi_minus = norm_cdf(-z);
  let phi_plus = norm_pdf(z);
  let fpp = PI_OVER_SIX * y / (s2 * s) * phi_minus * (8.0 * SQRT_THREE * s * ax + (3.0 * s2 * (s2 - 8.0) - 8.0 * x * x) * phi_minus / phi_plus) * exp(2.0 * y + 0.25 * s2);
  let (f, fp) = if is_below_horizon(s) {
    (0.0, 1.0)
  } else {
    let phi2 = phi_minus * phi_minus;
    let fp = TWO_PI * y * phi2 * exp(y + 0.125 * s * s);
    let f = if is_below_horizon(x) {
      0.0
    } else {
      TWO_PI_OVER_SQRT_TWENTY_SEVEN * ax * (phi2 * phi_minus)
    };
    (f, fp)
  };
  (f, fp, fpp)
}

///
fn compute_f_upper_map_and_first_two_derivatives(x: f64, s: f64) -> (f64, f64, f64) {
  let f = norm_cdf(-0.5 * s);
  let (fp, fpp) = if is_below_horizon(x) {
    let fp = -0.5;
    let fpp = 0.0;
    (fp, fpp)
  } else {
    let w = square(x / s);
    let fp = -0.5 * exp(0.5 * w);
    let fpp = SQRT_PI_OVER_TWO * exp(w + 0.125 * s * s) * w / s;
    (fp, fpp)
  };
  (f, fp, fpp)
}

///
fn inverse_f_lower_map(x: f64, f: f64) -> f64 {
  sel(
    is_below_horizon(f),
    0.0,
    fabs(x / (SQRT_THREE * inverse_norm_cdf(pow(f / (TWO_PI_OVER_SQRT_TWENTY_SEVEN * fabs(x)), 1.0 / 3.0)))),
  )
}

///
fn inverse_f_upper_map(f: f64) -> f64 {
  -2.0 * inverse_norm_cdf(f)
}

///```text
/// See www.en.wikipedia.org/wiki/Householder%27s_method for a detailed explanation of the third order Householder iteration.
///
/// Given the objective function g(s) whose root x such that 0 = g(s) we seek, iterate
///
///     s_n+1  =  s_n  -  (g/g') · [ 1 - (g''/g')·(g/g') ] / [ 1 - (g/g')·( (g''/g') - (g'''/g')·(g/g')/6 ) ]
///
/// Denoting  newton:=-(g/g'), halley:=(g''/g'), and hh3:=(g'''/g'), this reads
///
///     s_n+1  =  s_n  +  newton · [ 1 + halley·newton/2 ] / [ 1 + newton·( halley + hh3·newton/6 ) ]
///
/// NOTE that this function returns 0 when beta<intrinsic without any safety checks.
///```
fn unchecked_normalised_implied_volatility_from_a_transformed_rational_guess_with_limited_iterations(mut beta: f64, mut x: f64, mut q: f64 /* q=±1 */, n: usize) -> f64 {
  // Subtract intrinsic.
  if q * x > 0.0 {
    beta = fabs(max(beta - normalised_intrinsic(x, q), 0.0));
    q = -q;
  }
  // Map puts to calls
  if q < 0.0 {
    x = -x;
  }
  // For negative or zero prices we return 0.
  if beta <= 0.0 {
    return implied_volatility_output(0, 0.0);
  }
  // For positive but denormalised (a.k.a. 'subnormal') prices, we return 0 since it would be impossible to converge to full machine accuracy anyway.
  if beta < DENORMALISATION_CUTOFF {
    return implied_volatility_output(0, 0.0);
  }
  let b_max = exp(0.5 * x);
  if beta >= b_max {
    return implied_volatility_output(0, VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_ABOVE_MAXIMUM);
  }
  let mut iterations = 0_usize;
  let mut direction_reversal_count = 0_usize;

  let mut f = -DBL_MAX;
  let mut s = -DBL_MAX;
  let mut ds = s;
  let mut ds_previous = 0_f64;
  let mut s_left = DBL_MIN;
  let mut s_right = DBL_MAX;
  // The temptation is great to use the optimised form b_c = exp(x/2)/2-exp(-x/2)·Phi(sqrt(-2·x)) but that would require implementing all of the above types of round-off and over/underflow handling for this expression, too.
  let s_c = sqrt(fabs(2.0 * x));
  let b_c = normalised_black_call(x, s_c);
  let v_c = normalised_vega(x, s_c);
  if beta < b_c {
    let s_l = s_c - b_c / v_c;
    let b_l = normalised_black_call(x, s_l);
    if beta < b_l {
      let (f_lower_map_l, d_f_lower_map_l_d_beta, d2_f_lower_map_l_d_beta2) = compute_f_lower_map_and_first_two_derivatives(x, s_l);
      let r_ll =
        convex_rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(0., b_l, 0., f_lower_map_l, 1., d_f_lower_map_l_d_beta, d2_f_lower_map_l_d_beta2, true);
      f = rational_cubic_interpolation(beta, 0., b_l, 0., f_lower_map_l, 1., d_f_lower_map_l_d_beta, r_ll);
      if f.le(&0.0) {
        // This can happen due to roundoff truncation for extreme values such as |x|>500.
        // We switch to quadratic interpolation using f(0)≡0, f(b_l), and f'(0)≡1 to specify the quadratic.
        let t = beta / b_l;
        f = (f_lower_map_l * t + b_l * (1.0 - t)) * t;
      }
      s = inverse_f_lower_map(x, f);
      s_right = s_l;
      //
      // In this branch, which comprises the lowest segment, the objective function is
      //     g(s) = 1/ln(b(x,s)) - 1/ln(beta)
      //          ≡ 1/ln(b(s)) - 1/ln(beta)
      // This makes
      //              g'               =   -b'/(b·ln(b)²)
      //              newton = -g/g'   =   (ln(beta)-ln(b))·ln(b)/ln(beta)·b/b'
      //              halley = g''/g'  =   b''/b'  -  b'/b·(1+2/ln(b))
      //              hh3    = g'''/g' =   b'''/b' +  2(b'/b)²·(1+3/ln(b)·(1+1/ln(b)))  -  3(b''/b)·(1+2/ln(b))
      //
      // The Householder(3) iteration is
      //     s_n+1  =  s_n  +  newton · [ 1 + halley·newton/2 ] / [ 1 + newton·( halley + hh3·newton/6 ) ]
      //
      while iterations < n && fabs(ds) > DBL_EPSILON * s {
        if ds * ds_previous < 0.0 {
          direction_reversal_count += 1;
        }
        if iterations > 0 && (3 == direction_reversal_count || !(s > s_left && s < s_right)) {
          // If looping inefficiently, or the forecast step takes us outside the bracket, or onto its edges, switch to binary nesting.
          // NOTE that this can only really happen for very extreme values of |x|, such as |x| = |ln(F/K)| > 500.
          s = 0.5 * (s_left + s_right);
          if s_right - s_left <= DBL_EPSILON * s {
            break;
          };
          direction_reversal_count = 0;
          ds = 0.0;
        }
        ds_previous = ds;
        let b = normalised_black_call(x, s);
        let bp = normalised_vega(x, s);
        if b > beta && s < s_right {
          s_right = s;
        } else if b < beta && s > s_left {
          s_left = s;
        } // Tighten the bracket if applicable.
        if b <= 0.0 || bp <= 0.0 {
          // Numerical underflow. Switch to binary nesting for this iteration.
          ds = 0.5 * (s_left + s_right) - s;
        } else {
          let ln_b = log(b);
          let ln_beta = log(beta);
          let bpob = bp / b;
          let h = x / s;
          let b_halley = h * h / s - s / 4.0;
          let newton = (ln_beta - ln_b) * ln_b / ln_beta / bpob;
          let halley = b_halley - bpob * (1.0 + 2.0 / ln_b);
          let b_hh3 = b_halley * b_halley - 3.0 * square(h / s) - 0.25;
          let hh3 = b_hh3 + 2.0 * square(bpob) * (1.0 + 3.0 / ln_b * (1.0 + 1.0 / ln_b)) - 3.0 * b_halley * bpob * (1.0 + 2.0 / ln_b);
          ds = newton * householder_factor(newton, halley, hh3);
        }
        ds = max(-0.5 * s, ds);
        s += ds;
        iterations += 1;
      }
      return implied_volatility_output(iterations, s);
    } else {
      let v_l = normalised_vega(x, s_l);
      let r_lm = convex_rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(b_l, b_c, s_l, s_c, 1.0 / v_l, 1.0 / v_c, 0.0, false);
      s = rational_cubic_interpolation(beta, b_l, b_c, s_l, s_c, 1.0 / v_l, 1.0 / v_c, r_lm);
      s_left = s_l;
      s_right = s_c;
    }
  } else {
    let s_h = sel(v_c > DBL_MIN, s_c + (b_max - b_c) / v_c, s_c);
    let b_h = normalised_black_call(x, s_h);
    if beta <= b_h {
      let v_h = normalised_vega(x, s_h);
      let r_hm = convex_rational_cubic_control_parameter_to_fit_second_derivative_at_left_side(b_c, b_h, s_c, s_h, 1.0 / v_c, 1.0 / v_h, 0.0, false);
      s = rational_cubic_interpolation(beta, b_c, b_h, s_c, s_h, 1.0 / v_c, 1.0 / v_h, r_hm);
      s_left = s_c;
      s_right = s_h;
    } else {
      let (f_upper_map_h, d_f_upper_map_h_d_beta, d2_f_upper_map_h_d_beta2) = compute_f_upper_map_and_first_two_derivatives(x, s_h);
      if d2_f_upper_map_h_d_beta2 > -*SQRT_DBL_MAX && d2_f_upper_map_h_d_beta2 < *SQRT_DBL_MAX {
        let r_hh = convex_rational_cubic_control_parameter_to_fit_second_derivative_at_left_side(
          b_h,
          b_max,
          f_upper_map_h,
          0.,
          d_f_upper_map_h_d_beta,
          -0.5,
          d2_f_upper_map_h_d_beta2,
          true,
        );
        f = rational_cubic_interpolation(beta, b_h, b_max, f_upper_map_h, 0., d_f_upper_map_h_d_beta, -0.5, r_hh);
      }
      if f <= 0.0 {
        let h = b_max - b_h;
        let t = (beta - b_h) / h;
        f = (f_upper_map_h * (1.0 - t) + 0.5 * h * t) * (1.0 - t); // We switch to quadratic interpolation using f(b_h), f(b_max)≡0, and f'(b_max)≡-1/2 to specify the quadratic.
      }
      s = inverse_f_upper_map(f);
      s_left = s_h;
      if beta > 0.5 * b_max {
        // Else we better drop through and let the objective function be g(s) = b(x,s)-beta.
        //
        // In this branch, which comprises the upper segment, the objective function is
        //     g(s) = ln(b_max-beta)-ln(b_max-b(x,s))
        //          ≡ ln((b_max-beta)/(b_max-b(s)))
        // This makes
        //              g'               =   b'/(b_max-b)
        //              newton = -g/g'   =   ln((b_max-b)/(b_max-beta))·(b_max-b)/b'
        //              halley = g''/g'  =   b''/b'  +  b'/(b_max-b)
        //              hh3    = g'''/g' =   b'''/b' +  g'·(2g'+3b''/b')
        // and the iteration is
        //     s_n+1  =  s_n  +  newton · [ 1 + halley·newton/2 ] / [ 1 + newton·( halley + hh3·newton/6 ) ].
        //

        while iterations < n && fabs(ds) > DBL_EPSILON * s {
          if ds * ds_previous < 0.0 {
            direction_reversal_count += 1;
          }
          if iterations > 0 && (3 == direction_reversal_count || !(s > s_left && s < s_right)) {
            // If looping inefficiently, or the forecast step takes us outside the bracket, or onto its edges, switch to binary nesting.
            // NOTE that this can only really happen for very extreme values of |x|, such as |x| = |ln(F/K)| > 500.
            s = 0.5 * (s_left + s_right);
            if s_right - s_left <= DBL_EPSILON * s {
              break;
            };
            direction_reversal_count = 0;
            ds = 0.0;
          }
          ds_previous = ds;
          let b = normalised_black_call(x, s);
          let bp = normalised_vega(x, s);
          if b > beta && s < s_right {
            // Tighten the bracket if applicable.
            s_right = s;
          } else if b < beta && s > s_left {
            s_left = s;
          }
          if b >= b_max || bp <= DBL_MIN {
            // Numerical underflow. Switch to binary nesting for this iteration.
            ds = 0.5 * (s_left + s_right) - s;
          } else {
            let b_max_minus_b = b_max - b;
            let g = log((b_max - beta) / b_max_minus_b);
            let gp = bp / b_max_minus_b;
            let b_halley = square(x / s) / s - s / 4.0;
            let b_hh3 = b_halley * b_halley - 3.0 * square(x / (s * s)) - 0.25;
            let newton = -g / gp;
            let halley = b_halley + gp;
            let hh3 = b_hh3 + gp * (2.0 * gp + 3.0 * b_halley);
            ds = newton * householder_factor(newton, halley, hh3);
          }
          ds = max(-0.5 * s, ds);
          s += ds;
          iterations += 1;
        }
        return implied_volatility_output(iterations, s);
      }
    }
  }
  // In this branch, which comprises the two middle segments, the objective function is g(s) = b(x,s)-beta, or g(s) = b(s) - beta, for short.
  // This makes
  //              newton = -g/g'   =  -(b-beta)/b'
  //              halley = g''/g'  =    b''/b'    =  x²/s³-s/4
  //              hh3    = g'''/g' =    b'''/b'   =  halley² - 3·(x/s²)² - 1/4
  // and the iteration is
  //     s_n+1  =  s_n  +  newton · [ 1 + halley·newton/2 ] / [ 1 + newton·( halley + hh3·newton/6 ) ].
  //
  while iterations < n && fabs(ds) > DBL_EPSILON * s {
    if ds * ds_previous < 0.0 {
      direction_reversal_count += 1;
    }
    if iterations > 0 && (3 == direction_reversal_count || !(s > s_left && s < s_right)) {
      // If looping inefficiently, or the forecast step takes us outside the bracket, or onto its edges, switch to binary nesting.
      // NOTE that this can only really happen for very extreme values of |x|, such as |x| = |ln(F/K)| > 500.
      s = 0.5 * (s_left + s_right);
      if s_right - s_left <= DBL_EPSILON * s {
        break;
      };
      direction_reversal_count = 0;
      ds = 0.0;
    }
    ds_previous = ds;
    let b = normalised_black_call(x, s);
    let bp = normalised_vega(x, s);

    if b > beta && s < s_right {
      // Tighten the bracket if applicable.
      s_right = s;
    } else if b < beta && s > s_left {
      s_left = s;
    }
    let newton = (beta - b) / bp;
    let halley = square(x / s) / s - s / 4.0;
    let hh3 = halley * halley - 3.0 * square(x / (s * s)) - 0.25;
    ds = max(-0.5 * s, newton * householder_factor(newton, halley, hh3));
    s += ds;
    iterations += 1;
  }
  implied_volatility_output(iterations, s)
}

///
pub fn implied_volatility_from_a_transformed_rational_guess(price: f64, f: f64, k: f64, t: f64, q: f64 /* q=±1 */) -> f64 {
  implied_volatility_from_a_transformed_rational_guess_with_limited_iterations(price, f, k, t, q, get_implied_volatility_maximum_iterations())
}

///
pub fn implied_volatility_from_a_transformed_rational_guess_with_limited_iterations(mut price: f64, f: f64, k: f64, t: f64, mut q: f64 /* q=±1 */, n: usize) -> f64 {
  let intrinsic = fabs(max(sel(q < 0.0, k - f, f - k), 0.0));
  if price < intrinsic {
    return implied_volatility_output(0, VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_BELOW_INTRINSIC);
  }
  let max_price = sel(q < 0.0, k, f);
  if price >= max_price {
    return implied_volatility_output(0, VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_ABOVE_MAXIMUM);
  }
  let x = (f / k).ln();
  // Map in-the-money to out-of-the-money
  if q * x > 0.0 {
    price = fabs(max(price - intrinsic, 0.0));
    q = -q;
  }
  unchecked_normalised_implied_volatility_from_a_transformed_rational_guess_with_limited_iterations(price / (f.sqrt() * k.sqrt()), x, q, n) / t.sqrt()
}

pub fn normalised_implied_volatility_from_a_transformed_rational_guess_with_limited_iterations(mut beta: f64, x: f64, mut q: f64 /* q=±1 */, n: usize) -> f64 {
  // Map in-the-money to out-of-the-money
  if q * x > 0.0 {
    beta -= normalised_intrinsic(x, q);
    q = -q;
  }
  if beta < 0.0 {
    return implied_volatility_output(0, VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_BELOW_INTRINSIC);
  }
  unchecked_normalised_implied_volatility_from_a_transformed_rational_guess_with_limited_iterations(beta, x, q, n)
}

pub fn normalised_implied_volatility_from_a_transformed_rational_guess(beta: f64, x: f64, q: f64 /* q=±1 */) -> f64 {
  normalised_implied_volatility_from_a_transformed_rational_guess_with_limited_iterations(beta, x, q, get_implied_volatility_maximum_iterations())
}
