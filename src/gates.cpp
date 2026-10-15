#include "ecm_qsim/gates.hpp"

#include <cmath>
#include <complex>
#include <numbers>
#include <stdexcept>

namespace ecm_qsim {
namespace {

using cd = std::complex<double>;

constexpr cd kI{0.0, 1.0};

Eigen::Matrix2cd make_2x2(cd a, cd b, cd c, cd d) {
  Eigen::Matrix2cd m;
  m << a, b, c, d;
  return m;
}

}  // namespace

Eigen::Matrix2cd matrix_1q(GateType gate, double theta) {
  const double c = std::cos(theta / 2);
  const double s = std::sin(theta / 2);
  const double r = 1 / std::numbers::sqrt2;
  switch (gate) {
    case GateType::X:
      return make_2x2(0, 1, 1, 0);
    case GateType::Y:
      return make_2x2(0, -kI, kI, 0);
    case GateType::Z:
      return make_2x2(1, 0, 0, -1);
    case GateType::H:
      return make_2x2(r, r, r, -r);
    case GateType::S:
      return make_2x2(1, 0, 0, kI);
    case GateType::T:
      return make_2x2(1, 0, 0, std::polar(1.0, std::numbers::pi / 4));
    case GateType::RX:
      return make_2x2(c, -kI * s, -kI * s, c);
    case GateType::RY:
      return make_2x2(c, -s, s, c);
    case GateType::RZ:
      return make_2x2(cd(c, -s), 0, 0, cd(c, s));
    case GateType::CNOT:
    case GateType::CZ:
    case GateType::CY:
      break;
  }
  throw std::invalid_argument("matrix_1q: not a one-qubit gate");
}

}  // namespace ecm_qsim
