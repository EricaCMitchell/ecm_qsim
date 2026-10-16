#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "ecm_qsim/noise.hpp"

using namespace ecm_qsim;
using Catch::Matchers::ContainsSubstring;

TEST_CASE("NoiseModel defaults to noiseless") {
  NoiseModel noise;
  REQUIRE(noise.p_gate == 0.0);
  REQUIRE(noise.p_idle == 0.0);
  REQUIRE(noise.p_readout == 0.0);
}

TEST_CASE("validate_noise_model accepts probabilities in [0, 1]") {
  CHECK_NOTHROW(validate_noise_model(NoiseModel{}));
  CHECK_NOTHROW(validate_noise_model(NoiseModel{1.0, 1.0, 1.0}));
  CHECK_NOTHROW(validate_noise_model(NoiseModel{0.001, 0.5, 0.01}));
}

TEST_CASE("validate_noise_model rejects values outside [0, 1] and NaN") {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  const double bad = GENERATE_COPY(-0.1, 1.5, nan, inf, -inf);
  CAPTURE(bad);
  CHECK_THROWS_WITH(validate_noise_model(NoiseModel{bad, 0, 0}),
                    ContainsSubstring("p_gate must be in [0, 1], got "));
  CHECK_THROWS_WITH(validate_noise_model(NoiseModel{0, bad, 0}),
                    ContainsSubstring("p_idle must be in [0, 1], got "));
  CHECK_THROWS_WITH(validate_noise_model(NoiseModel{0, 0, bad}),
                    ContainsSubstring("p_readout must be in [0, 1], got "));
  CHECK_THROWS_AS(validate_noise_model(NoiseModel{bad, 0, 0}),
                  std::invalid_argument);
}

TEST_CASE("validate_noise_model reports the value as given") {
  CHECK_THROWS_WITH(validate_noise_model(NoiseModel{0, 1.5, 0}),
                    "p_idle must be in [0, 1], got 1.5");
  // Default ostream formatting would print "got 1" here.
  CHECK_THROWS_WITH(validate_noise_model(NoiseModel{1.0000001, 0, 0}),
                    "p_gate must be in [0, 1], got 1.0000001");
  CHECK_THROWS_WITH(validate_noise_model(NoiseModel{
                        0, 0, -std::numeric_limits<double>::infinity()}),
                    "p_readout must be in [0, 1], got -inf");
}

TEST_CASE("validate_noise_model reports the first bad field") {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  CHECK_THROWS_WITH(validate_noise_model(NoiseModel{-0.1, 1.5, nan}),
                    "p_gate must be in [0, 1], got -0.1");
  CHECK_THROWS_WITH(validate_noise_model(NoiseModel{0, 1.5, nan}),
                    "p_idle must be in [0, 1], got 1.5");
}
