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

#ifndef ECM_QSIM_MEASUREMENT_HPP_
#define ECM_QSIM_MEASUREMENT_HPP_

// Measurement of all qubits. Sample a probability vector to draw from one
// distribution many times; sample a StateVector to draw from a state once
// without a 2^n probabilities copy.

#include <Eigen/Core>
#include <cstdint>

#include "ecm_qsim/state_vector.hpp"

namespace ecm_qsim {

// Maps a uniform draw u in [0, 1) to the index whose interval in the
// cumulative distribution contains u * (sum of probabilities). The
// probabilities need not sum to 1. Never returns an index with probability 0.
// Throws std::invalid_argument if u is not in [0, 1), if an entry is negative
// or NaN, or if the sum is not finite. Throws std::domain_error if all
// entries are zero (including an empty vector).
std::uint64_t sample_index(const Eigen::VectorXd& probabilities, double u);

// sample_index(state.probabilities(), u) without the copy: each probability
// is computed from its amplitude as it is read. Same result, same exceptions.
std::uint64_t sample_index(const StateVector& state, double u);

// A uniform draw in [0, 1) from 64 random bits: the top 53 bits scaled by
// 2^-53. Use it instead of std::uniform_real_distribution, which can return
// 1.0 on libc++ (LWG 2524) and gives different values on other standard
// libraries.
constexpr double uniform_from_bits(std::uint64_t bits) noexcept {
  return static_cast<double>(bits >> 11) * 0x1p-53;
}

}  // namespace ecm_qsim

#endif  // ECM_QSIM_MEASUREMENT_HPP_
