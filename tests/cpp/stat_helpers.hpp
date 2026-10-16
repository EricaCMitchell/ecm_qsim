#ifndef ECM_QSIM_TESTS_STAT_HELPERS_HPP_
#define ECM_QSIM_TESTS_STAT_HELPERS_HPP_

// Catch2 helpers for statistical checks on sampled counts.

#include <catch2/catch_test_macros.hpp>
#include <cmath>

// |count - trials * p| within 5 binomial standard deviations.
inline void require_binomial(int count, int trials, double p) {
  const double mean = trials * p;
  const double sigma = std::sqrt(trials * p * (1 - p));
  CAPTURE(count, trials, p, mean, sigma);
  REQUIRE(std::abs(count - mean) <= 5 * sigma);
}

#endif  // ECM_QSIM_TESTS_STAT_HELPERS_HPP_
