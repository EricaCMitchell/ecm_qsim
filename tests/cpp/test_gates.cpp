#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <cmath>
#include <complex>
#include <initializer_list>
#include <numbers>
#include <stdexcept>
#include <unsupported/Eigen/KroneckerProduct>
#include <unsupported/Eigen/MatrixFunctions>

#include "ecm_qsim/circuit.hpp"
#include "ecm_qsim/gates.hpp"
#include "eigen_helpers.hpp"

using namespace ecm_qsim;
using std::numbers::pi;

namespace {

using cd = std::complex<double>;

constexpr cd kI{0.0, 1.0};

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

Eigen::Matrix4cd m4(std::initializer_list<cd> entries) {
  Eigen::Matrix4cd m;
  auto it = entries.begin();
  for (int row = 0; row < 4; ++row) {
    for (int col = 0; col < 4; ++col) {
      m(row, col) = *it++;
    }
  }
  return m;
}

// u on the control (first-listed qubit, the high bit of the 4x4 row index)
// or on the target. kroneckerProduct(A, B) puts A on the high bit.
Eigen::Matrix4cd on_control(const Eigen::Matrix2cd& u) {
  return Eigen::kroneckerProduct(u, kId2);
}
Eigen::Matrix4cd on_target(const Eigen::Matrix2cd& u) {
  return Eigen::kroneckerProduct(kId2, u);
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

TEST_CASE("Two-qubit gates match their definitions") {
  require_close(matrix_2q(GateType::CNOT), m4({1, 0, 0, 0,  //
                                               0, 1, 0, 0,  //
                                               0, 0, 0, 1,  //
                                               0, 0, 1, 0}));
  require_close(matrix_2q(GateType::CY), m4({1, 0, 0, 0,    //
                                             0, 1, 0, 0,    //
                                             0, 0, 0, -kI,  //
                                             0, 0, kI, 0}));
  require_close(matrix_2q(GateType::CZ), m4({1, 0, 0, 0,  //
                                             0, 1, 0, 0,  //
                                             0, 0, 1, 0,  //
                                             0, 0, 0, -1}));
}

TEST_CASE("Controlled gates are P0 (x) I + P1 (x) U, control as high bit") {
  // kroneckerProduct(A, B): A acts on the high bit of the row index.
  const Eigen::Matrix2cd p0 = m2(1, 0, 0, 0);
  const Eigen::Matrix2cd p1 = m2(0, 0, 0, 1);
  const auto controlled = [&](const Eigen::Matrix2cd& u) -> Eigen::Matrix4cd {
    return Eigen::kroneckerProduct(p0, kId2).eval() +
           Eigen::kroneckerProduct(p1, u).eval();
  };
  require_close(matrix_2q(GateType::CNOT), controlled(X()));
  require_close(matrix_2q(GateType::CY), controlled(Y()));
  require_close(matrix_2q(GateType::CZ), controlled(Z()));
}

TEST_CASE("CNOT maps basis state |c t> to |c, t xor c>") {
  const Eigen::Matrix4cd cnot = matrix_2q(GateType::CNOT);
  for (int c = 0; c < 2; ++c) {
    for (int t = 0; t < 2; ++t) {
      CAPTURE(c, t);
      const int in = 2 * c + t;
      const int out = 2 * c + (t ^ c);
      require_close(cnot.col(in), Eigen::Vector4cd::Unit(out));
    }
  }
}

TEST_CASE("Two-qubit gates are unitary and Hermitian") {
  const GateType gate = GENERATE(GateType::CNOT, GateType::CZ, GateType::CY);
  CAPTURE(gate);
  const Eigen::Matrix4cd m = matrix_2q(gate);
  require_close(m.adjoint() * m, Eigen::Matrix4cd::Identity());
  require_close(m.adjoint(), m);
}

TEST_CASE("CNOT propagates Paulis between control and target") {
  // Nielsen & Chuang Eqs. 4.32-4.37, with subscript 0 = control and
  // 1 = target. Those subscripts are not the simulator's qubit numbers.
  const Eigen::Matrix4cd c = matrix_2q(GateType::CNOT);
  require_close(c * on_control(X()) * c, on_control(X()) * on_target(X()));
  require_close(c * on_control(Y()) * c, on_control(Y()) * on_target(X()));
  require_close(c * on_control(Z()) * c, on_control(Z()));
  require_close(c * on_target(X()) * c, on_target(X()));
  require_close(c * on_target(Y()) * c, on_control(Z()) * on_target(Y()));
  require_close(c * on_target(Z()) * c, on_control(Z()) * on_target(Z()));
}

TEST_CASE("CNOT commutes with RZ on the control and RX on the target") {
  // Nielsen & Chuang Eqs. 4.38-4.39.
  const double theta = GENERATE(0.7, -2.5);
  CAPTURE(theta);
  const Eigen::Matrix4cd c = matrix_2q(GateType::CNOT);
  const Eigen::Matrix4cd rz = on_control(matrix_1q(GateType::RZ, theta));
  const Eigen::Matrix4cd rx = on_target(matrix_1q(GateType::RX, theta));
  require_close(rz * c, c * rz);
  require_close(rx * c, c * rx);
}

TEST_CASE("CZ is symmetric in its qubits; CNOT is not") {
  const Eigen::Matrix4cd swap = m4({1, 0, 0, 0,  //
                                    0, 0, 1, 0,  //
                                    0, 1, 0, 0,  //
                                    0, 0, 0, 1});
  const Eigen::Matrix4cd cz = matrix_2q(GateType::CZ);
  const Eigen::Matrix4cd cnot = matrix_2q(GateType::CNOT);
  require_close(swap * cz * swap, cz);
  REQUIRE((swap * cnot * swap - cnot).cwiseAbs().maxCoeff() > 0.5);
}

TEST_CASE("matrix_2q rejects gates that are not two-qubit gates") {
  const GateType gate = GENERATE(
      GateType::X, GateType::Y, GateType::Z, GateType::H, GateType::S,
      GateType::T, GateType::RX, GateType::RY, GateType::RZ, kInvalidGate);
  CAPTURE(gate);
  REQUIRE_THROWS_AS(matrix_2q(gate), std::invalid_argument);
  REQUIRE_THROWS_WITH(matrix_2q(gate), "matrix_2q: not a two-qubit gate");
}

TEST_CASE("is_two_qubit agrees with matrix_1q and matrix_2q") {
  // Every enumerator, from X to CY (the last one).
  for (int i = 0; i <= static_cast<int>(GateType::CY); ++i) {
    const auto gate = static_cast<GateType>(i);
    CAPTURE(i);
    if (is_two_qubit(gate)) {
      CHECK_THROWS_AS(matrix_1q(gate), std::invalid_argument);
      CHECK_NOTHROW(matrix_2q(gate));
    } else {
      CHECK_NOTHROW(matrix_1q(gate));
      CHECK_THROWS_AS(matrix_2q(gate), std::invalid_argument);
    }
  }
}
