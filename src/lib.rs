extern crate lazy_static;

mod definitions;
mod erf_cody;
mod lets_be_rational;
mod normal_distribution;
mod rational_cubic;

pub use erf_cody::{erf_cody, erfc_cody, erfcx_cody};
pub use lets_be_rational::{
  black, implied_volatility_from_a_transformed_rational_guess, implied_volatility_from_a_transformed_rational_guess_with_limited_iterations, normalised_black,
  normalised_black_call, normalised_implied_volatility_from_a_transformed_rational_guess, normalised_implied_volatility_from_a_transformed_rational_guess_with_limited_iterations,
  normalised_vega,
};
