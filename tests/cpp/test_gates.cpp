#include <Eigen/Core>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <cmath>
#include <complex>
#include <numbers>
#include <stdexcept>
#include <unsupported/Eigen/MatrixFunctions>

#include "ecm_qsim/circuit.hpp"
#include "ecm_qsim/gates.hpp"

using namespace ecm_qsim;
using Catch::Approx;
using std::numbers::pi;

namespace {

using cd = std::complex<double>;

constexpr cd kI{0.0, 1.0};

// Absolute tolerance only: epsilon(0) disables Approx's default relative
// tolerance, which would otherwise also be accepted.
Approx near(double expected) {
  return Approx(expected).epsilon(0).margin(1e-12);
}

// Entrywise comparison. Eigen's isApprox is relative, which is meaningless
// when the expected matrix is mostly zeros.
void require_close(const Eigen::MatrixXcd& actual,
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

Eigen::Matrix2cd m2(cd a, cd b, cd c, cd d) {
  Eigen::Matrix2cd m;
  m << a, b, c, d;
  return m;
}

const Eigen::Matrix2cd kId2 = Eigen::Matrix2cd::Identity();

Eigen::Matrix2cd X() {
  return matrix_1q(GateType::X);
}
Eigen::Matrix2cd Y() {
  return matrix_1q(GateType::Y);
}
Eigen::Matrix2cd Z() {
  return matrix_1q(GateType::Z);
}
Eigen::Matrix2cd H() {
  return matrix_1q(GateType::H);
}
Eigen::Matrix2cd S() {
  return matrix_1q(GateType::S);
}
Eigen::Matrix2cd T() {
  return matrix_1q(GateType::T);
}

// A value outside the GateType enumerators (valid for a scoped enum).
constexpr GateType kInvalidGate = static_cast<GateType>(99);

}  // namespace

TEST_CASE("Fixed one-qubit gates match their definitions") {
  const double r = 1 / std::sqrt(2.0);
  require_close(X(), m2(0, 1, 1, 0));
  require_close(Y(), m2(0, -kI, kI, 0));
  require_close(Z(), m2(1, 0, 0, -1));
  require_close(H(), m2(r, r, r, -r));
  require_close(S(), m2(1, 0, 0, kI));
  require_close(T(), m2(1, 0, 0, cd(r, r)));
}

TEST_CASE("Rotations equal exp(-i theta P / 2)") {
  const double theta = GENERATE(0.0, 0.7, pi / 2, pi, -1.3, 2 * pi, 10.0);
  CAPTURE(theta);
  const cd k = -kI * theta / 2.0;
  require_close(matrix_1q(GateType::RX, theta), (k * X()).exp());
  require_close(matrix_1q(GateType::RY, theta), (k * Y()).exp());
  require_close(matrix_1q(GateType::RZ, theta), (k * Z()).exp());
}

TEST_CASE("Rotations at special angles") {
  require_close(matrix_1q(GateType::RX, 0), kId2);
  require_close(matrix_1q(GateType::RY, 0), kId2);
  require_close(matrix_1q(GateType::RZ, 0), kId2);
  require_close(matrix_1q(GateType::RX, pi), -kI * X());
  require_close(matrix_1q(GateType::RY, pi), -kI * Y());
  require_close(matrix_1q(GateType::RZ, pi), -kI * Z());
  require_close(matrix_1q(GateType::RX, 2 * pi), -kId2);
  require_close(matrix_1q(GateType::RY, pi / 2), H() * Z());
  require_close(matrix_1q(GateType::RZ, pi / 2),
                std::polar(1.0, -pi / 4) * S());
  require_close(matrix_1q(GateType::RZ, pi / 4),
                std::polar(1.0, -pi / 8) * T());
}

TEST_CASE("One-qubit gates are unitary") {
  const GateType gate =
      GENERATE(GateType::X, GateType::Y, GateType::Z, GateType::H, GateType::S,
               GateType::T, GateType::RX, GateType::RY, GateType::RZ);
  const double theta = GENERATE(0.0, 0.7, -2.5);
  CAPTURE(gate, theta);
  const Eigen::Matrix2cd u = matrix_1q(gate, theta);
  require_close(u.adjoint() * u, kId2);
  require_close(u * u.adjoint(), kId2);
}

TEST_CASE("One-qubit gate identities") {
  // Nielsen & Chuang Eq. 4.18.
  require_close(H() * X() * H(), Z());
  require_close(H() * Y() * H(), -Y());
  require_close(H() * Z() * H(), X());
  require_close(X() * Y(), kI * Z());
  require_close(S() * S(), Z());
  require_close(T() * T(), S());
}

TEST_CASE("Rotation by -theta undoes rotation by theta") {
  const GateType gate = GENERATE(GateType::RX, GateType::RY, GateType::RZ);
  const double theta = GENERATE(0.7, -2.5, 10.0);
  CAPTURE(gate, theta);
  const Eigen::Matrix2cd u = matrix_1q(gate, theta);
  require_close(matrix_1q(gate, -theta), u.adjoint());
  require_close(matrix_1q(gate, -theta) * u, kId2);
}

TEST_CASE("Fixed gates ignore theta") {
  const GateType gate = GENERATE(GateType::X, GateType::Y, GateType::Z,
                                 GateType::H, GateType::S, GateType::T);
  CAPTURE(gate);
  require_close(matrix_1q(gate, 1.23), matrix_1q(gate));
}

TEST_CASE("matrix_1q rejects gates that are not one-qubit gates") {
  const GateType gate =
      GENERATE(GateType::CNOT, GateType::CZ, GateType::CY, kInvalidGate);
  CAPTURE(gate);
  REQUIRE_THROWS_AS(matrix_1q(gate), std::invalid_argument);
  REQUIRE_THROWS_WITH(matrix_1q(gate), "matrix_1q: not a one-qubit gate");
}
