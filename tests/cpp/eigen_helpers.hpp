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

#ifndef ECM_QSIM_TESTS_EIGEN_HELPERS_HPP_
#define ECM_QSIM_TESTS_EIGEN_HELPERS_HPP_

// Catch2 helpers for comparing Eigen matrices and vectors.

#include <Eigen/Core>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

// Absolute tolerance only: epsilon(0) disables Approx's default relative
// tolerance, which would otherwise also be accepted.
inline Catch::Approx near(double expected) {
  return Catch::Approx(expected).epsilon(0).margin(1e-12);
}

// Entrywise comparison. Eigen's isApprox is relative, which is meaningless
// when the expected matrix is mostly zeros.
inline void require_close(const Eigen::MatrixXcd& actual,
                          const Eigen::MatrixXcd& expected) {
  INFO("actual:\n" << actual << "\nexpected:\n" << expected);
  REQUIRE(actual.rows() == expected.rows());
  REQUIRE(actual.cols() == expected.cols());
  for (Eigen::Index row = 0; row < actual.rows(); ++row) {
    for (Eigen::Index col = 0; col < actual.cols(); ++col) {
      CAPTURE(row, col);
      REQUIRE(actual(row, col).real() == near(expected(row, col).real()));
      REQUIRE(actual(row, col).imag() == near(expected(row, col).imag()));
    }
  }
}

#endif  // ECM_QSIM_TESTS_EIGEN_HELPERS_HPP_
