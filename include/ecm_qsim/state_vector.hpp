#ifndef ECM_QSIM_STATE_VECTOR_HPP_
#define ECM_QSIM_STATE_VECTOR_HPP_

// StateVector: owns the 2^n amplitudes; apply_1q / apply_2q kernels.
// Amplitude index i holds the basis state in which qubit k equals bit k of i,
// so qubit 0 is the least significant bit.

#include <Eigen/Core>
#include <complex>
#include <cstdint>
#include <iosfwd>
#include <string>

namespace ecm_qsim {

class StateVector {
public:
  // |0...0> on n_qubits qubits. Throws std::invalid_argument unless
  // 1 <= n_qubits <= kMaxQubits, and std::length_error if
  // n_qubits > max_qubits() or allocation otherwise fails.
  explicit StateVector(int n_qubits);

  // Takes the amplitudes as given, without normalising them. Throws
  // std::invalid_argument unless the size is 2^n for some n >= 1.
  explicit StateVector(Eigen::VectorXcd amplitudes);

  // Largest qubit count the first constructor accepts on this machine: the
  // amplitudes (16 bytes each) plus one probabilities() copy (8 bytes each)
  // must fit in half of physical memory. At most kMaxQubits.
  static int max_qubits();

  int n_qubits() const noexcept { return n_qubits_; }
  std::uint64_t size() const noexcept { return std::uint64_t{1} << n_qubits_; }

  // Throws std::out_of_range if index >= size().
  std::complex<double> amplitude(std::uint64_t index) const;

  // A copy, so callers do not depend on the internal layout.
  Eigen::VectorXcd to_vector() const;

  // |amplitude|^2 for each basis state.
  Eigen::VectorXd probabilities() const;

  // Applies u to qubit k in place: for each i with bit k clear,
  // (psi[i], psi[i + 2^k]) <- u * (psi[i], psi[i + 2^k]).
  // Throws std::out_of_range unless 0 <= k < n_qubits().
  void apply_1q(const Eigen::Matrix2cd& u, int k);

  // Applies m to qubits j and k in place. j is the high bit of m's row index
  // (for CNOT c t, j = c): for each i with bits j and k clear, applies m to
  // (psi[i], psi[i + 2^k], psi[i + 2^j], psi[i + 2^j + 2^k]).
  // Throws std::out_of_range if j or k is not a qubit, std::invalid_argument
  // if j == k.
  void apply_2q(const Eigen::Matrix4cd& m, int j, int k);

  // One line "|<bitstring>>  p=<4 decimals>" per basis state, in index
  // order. States with negligible probability are omitted.
  std::string to_string() const;

private:
  int n_qubits_;  // declared before amplitudes_: the constructors rely on it
  Eigen::VectorXcd amplitudes_;
};

std::ostream& operator<<(std::ostream& os, const StateVector& state);

// index as an n_qubits-character bitstring, most significant bit first, so
// qubit 0 is the rightmost character. Shared by StateVector::to_string() and
// the counts output. Throws std::invalid_argument unless
// 1 <= n_qubits <= kMaxQubits, std::out_of_range if index >= 2^n_qubits.
std::string to_bitstring(std::uint64_t index, int n_qubits);

}  // namespace ecm_qsim

#endif  // ECM_QSIM_STATE_VECTOR_HPP_
