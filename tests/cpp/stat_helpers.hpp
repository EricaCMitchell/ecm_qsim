// Copyright 2026 Erica C. Mitchell
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// SPDX-License-Identifier: Apache-2.0

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
