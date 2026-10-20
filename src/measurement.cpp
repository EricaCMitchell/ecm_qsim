#include "ecm_qsim/measurement.hpp"

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>

namespace ecm_qsim {

namespace {

// The algorithm behind both overloads. probabilities is any dense Eigen
// expression of double; each entry is read once in each pass.
template <typename Derived>
std::uint64_t sample_index_impl(const Eigen::DenseBase<Derived>& probabilities,
                                double u) {
  if (!(u >= 0.0 && u < 1.0)) {  // also rejects NaN
    throw std::invalid_argument("sample_index: u must be in [0, 1), got " +
                                std::to_string(u));
  }
  double total = 0.0;
  for (Eigen::Index i = 0; i < probabilities.size(); ++i) {
    const double p = probabilities(i);
    if (!(p >= 0.0)) {  // also rejects NaN
      throw std::invalid_argument("sample_index: probability " +
                                  std::to_string(i) + " is " +
                                  std::to_string(p));
    }
    total += p;
  }
  if (!std::isfinite(total)) {
    throw std::invalid_argument(
        "sample_index: total probability is not finite");
  }
  if (total == 0.0) {
    throw std::domain_error("sample_index: all probabilities are zero");
  }
  const double target = u * total;
  double cumulative = 0.0;
  Eigen::Index last_nonzero = 0;
  for (Eigen::Index i = 0; i < probabilities.size(); ++i) {
    const double p = probabilities(i);
    if (p == 0.0) {
      continue;
    }
    cumulative += p;
    last_nonzero = i;
    if (cumulative > target) {
      return static_cast<std::uint64_t>(i);
    }
  }
  // The running sum ends at exactly total, which is greater than u * total
  // unless total is subnormal and u * total rounds up to it.
  return static_cast<std::uint64_t>(last_nonzero);
}

}  // namespace

std::uint64_t sample_index(const Eigen::VectorXd& probabilities, double u) {
  return sample_index_impl(probabilities, u);
}

std::uint64_t sample_index(const StateVector& state, double u) {
  return sample_index_impl(state.amplitudes_.cwiseAbs2(), u);
}

}  // namespace ecm_qsim
