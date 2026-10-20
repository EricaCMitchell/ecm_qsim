#include "ecm_qsim/fusion.hpp"

#include <algorithm>
#include <array>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "ecm_qsim/gates.hpp"

namespace ecm_qsim {

namespace {

// True if the gate's matrix is real, so that its cluster can stay real. No
// default case, so -Wswitch flags a gate added later.
bool is_real_gate(GateType gate) {
  switch (gate) {
    case GateType::X:
    case GateType::Z:
    case GateType::H:
    case GateType::RY:
    case GateType::CNOT:
    case GateType::CZ:
      return true;
    case GateType::Y:
    case GateType::S:
    case GateType::T:
    case GateType::RX:
    case GateType::RZ:
    case GateType::CY:
      return false;
  }
  return false;
}

bool is_real(const FusedGate& gate) {
  return std::holds_alternative<Eigen::Matrix2d>(gate.matrix) ||
         std::holds_alternative<Eigen::Matrix4d>(gate.matrix);
}

int arity(const FusedGate& gate) {
  return std::holds_alternative<Eigen::Matrix2d>(gate.matrix) ||
                 std::holds_alternative<Eigen::Matrix2cd>(gate.matrix)
             ? 1
             : 2;
}

// The same check and message as the StateVector kernel named by where, so
// that simulate() reports a bad qubit as it did before fusion.
void check_qubit(int qubit, int n_qubits, const char* where) {
  if (qubit < 0 || qubit >= n_qubits) {
    throw std::out_of_range(std::string(where) + ": qubit " +
                            std::to_string(qubit) + " out of range for " +
                            std::to_string(n_qubits) + " qubits");
  }
}

// P m P, where P swaps the two bits of the matrix index: m with its
// operands in the other order.
Eigen::Matrix4cd swap_operands(const Eigen::Matrix4cd& m) {
  constexpr std::array<int, 4> kSwapped = {0, 2, 1, 3};
  Eigen::Matrix4cd result;
  for (int r = 0; r < 4; ++r) {
    for (int c = 0; c < 4; ++c) {
      result(r, c) = m(kSwapped[r], kSwapped[c]);
    }
  }
  return result;
}

// An operation as a fused gate of its own. Checks its qubits first, with
// the kernel that simulate() used for it before fusion: apply_ry for RY,
// apply_cnot for CNOT, apply_2q for the other two-qubit gates and apply_1q
// for the rest.
FusedGate gate_of(const Op& op, int n_qubits) {
  if (is_two_qubit(op.gate)) {
    const char* where = op.gate == GateType::CNOT ? "apply_cnot" : "apply_2q";
    check_qubit(op.qubits[0], n_qubits, where);
    check_qubit(op.qubits[1], n_qubits, where);
    if (op.qubits[0] == op.qubits[1]) {
      throw std::invalid_argument(std::string(where) + ": both qubits are " +
                                  std::to_string(op.qubits[0]));
    }
    Eigen::Matrix4cd m = matrix_2q(op.gate);
    if (op.qubits[0] > op.qubits[1]) {
      m = swap_operands(m);
    }
    const std::array<int, 2> qubits = {std::min(op.qubits[0], op.qubits[1]),
                                       std::max(op.qubits[0], op.qubits[1])};
    if (is_real_gate(op.gate)) {
      return FusedGate{qubits, Eigen::Matrix4d(m.real())};
    }
    return FusedGate{qubits, m};
  }
  check_qubit(op.qubits[0], n_qubits,
              op.gate == GateType::RY ? "apply_ry" : "apply_1q");
  const Eigen::Matrix2cd u = matrix_1q(op.gate, op.theta);
  if (is_real_gate(op.gate)) {
    return FusedGate{{op.qubits[0], 0}, Eigen::Matrix2d(u.real())};
  }
  return FusedGate{{op.qubits[0], 0}, u};
}

FusedGate gate_of(const PauliError& error, int n_qubits) {
  return gate_of(Op{error.pauli, {error.qubit, 0}, 0.0}, n_qubits);
}

// matrix's alternative of type Matrix. A real matrix is promoted when Matrix
// is complex; a complex one is never made real.
template <typename Matrix>
Matrix get_as(const FusedMatrix& matrix) {
  using Real = Eigen::Matrix<double, Matrix::RowsAtCompileTime,
                             Matrix::ColsAtCompileTime>;
  if constexpr (!std::is_same_v<Matrix, Real>) {
    if (const auto* real = std::get_if<Real>(&matrix)) {
      return real->template cast<std::complex<double>>();
    }
  }
  return std::get<Matrix>(matrix);
}

// u (x) I or I (x) u: u on the high or the low bit of a 4x4 index.
template <typename Matrix2>
Eigen::Matrix<typename Matrix2::Scalar, 4, 4> embed(const Matrix2& u,
                                                    bool high_bit) {
  Eigen::Matrix<typename Matrix2::Scalar, 4, 4> result;
  result.setZero();
  for (int a = 0; a < 2; ++a) {
    for (int b = 0; b < 2; ++b) {
      if (high_bit) {
        result(2 * a, 2 * b) = u(a, b);
        result(2 * a + 1, 2 * b + 1) = u(a, b);
      } else {
        result(a, b) = u(a, b);
        result(2 + a, 2 + b) = u(a, b);
      }
    }
  }
  return result;
}

// The product of parts, the last part leftmost, as one gate on qubits
// (sorted, qubits[1] unused if n_qubits is 1). Every part acts on a subset
// of qubits. Scalar is double only if every part is real.
template <typename Scalar>
FusedGate product(const std::array<int, 2>& qubits, int n_qubits,
                  const std::array<const FusedGate*, 3>& parts,
                  std::size_t n_parts) {
  using Matrix2 = Eigen::Matrix<Scalar, 2, 2>;
  using Matrix4 = Eigen::Matrix<Scalar, 4, 4>;
  if (n_qubits == 1) {
    Matrix2 result = Matrix2::Identity();
    for (std::size_t i = 0; i < n_parts; ++i) {
      result = get_as<Matrix2>(parts[i]->matrix) * result;
    }
    return FusedGate{{qubits[0], 0}, result};
  }
  Matrix4 result = Matrix4::Identity();
  for (std::size_t i = 0; i < n_parts; ++i) {
    const FusedGate& part = *parts[i];
    // A two-qubit part already acts on (qubits[0], qubits[1]) in that order.
    const Matrix4 m = arity(part) == 2 ? get_as<Matrix4>(part.matrix)
                                       : embed(get_as<Matrix2>(part.matrix),
                                               part.qubits[0] == qubits[0]);
    result = m * result;
  }
  return FusedGate{qubits, result};
}

// Greedy fusion to width two over the operations in execution order.
class Fuser {
public:
  explicit Fuser(std::size_t n_operations) { done_.reserve(n_operations); }

  void add(const FusedGate& operation) {
    // The open clusters are disjoint, so at most two touch the operation.
    std::array<std::size_t, 2> touching{};
    std::size_t n_touching = 0;
    // The operation's qubits and those of up to two touching clusters.
    std::array<int, 6> qubits{};
    std::size_t n_qubits = 0;
    const auto add_qubits = [&qubits, &n_qubits](const FusedGate& gate) {
      for (int i = 0; i < arity(gate); ++i) {
        qubits[n_qubits++] = gate.qubits[static_cast<std::size_t>(i)];
      }
    };
    add_qubits(operation);
    for (std::size_t c = 0; c < open_.size(); ++c) {
      if (shares_qubit(open_[c].gate, operation)) {
        touching[n_touching++] = c;
        add_qubits(open_[c].gate);
      }
    }
    std::sort(qubits.begin(), qubits.begin() + n_qubits);
    n_qubits = static_cast<std::size_t>(
        std::unique(qubits.begin(), qubits.begin() + n_qubits) -
        qubits.begin());
    // Earliest first, so that both branches below keep creation order.
    if (n_touching == 2 &&
        open_[touching[0]].created > open_[touching[1]].created) {
      std::swap(touching[0], touching[1]);
    }

    if (n_qubits > 2) {
      for (std::size_t i = 0; i < n_touching; ++i) {
        done_.push_back(open_[touching[i]].gate);
      }
      remove(touching, n_touching);
      open_.push_back(Open{next_created_++, operation});
      return;
    }

    // The touching clusters are disjoint, so they commute, and the
    // operation follows both.
    std::array<const FusedGate*, 3> parts{};
    std::size_t n_parts = 0;
    bool real = is_real(operation);
    for (std::size_t i = 0; i < n_touching; ++i) {
      parts[n_parts++] = &open_[touching[i]].gate;
      real = real && is_real(open_[touching[i]].gate);
    }
    parts[n_parts++] = &operation;
    const std::array<int, 2> union_qubits = {qubits[0],
                                             n_qubits == 2 ? qubits[1] : 0};
    const int width = static_cast<int>(n_qubits);
    Open merged{n_touching == 0 ? next_created_++ : open_[touching[0]].created,
                real ? product<double>(union_qubits, width, parts, n_parts)
                     : product<std::complex<double>>(union_qubits, width, parts,
                                                     n_parts)};
    remove(touching, n_touching);
    open_.push_back(std::move(merged));
  }

  // The emitted clusters, then the open ones in creation order.
  std::vector<FusedGate> finish() {
    std::sort(open_.begin(), open_.end(), [](const Open& a, const Open& b) {
      return a.created < b.created;
    });
    for (Open& cluster : open_) {
      done_.push_back(std::move(cluster.gate));
    }
    open_.clear();
    return std::move(done_);
  }

private:
  struct Open {
    std::size_t created;  // position of the cluster's first operation
    FusedGate gate;
  };

  static bool shares_qubit(const FusedGate& a, const FusedGate& b) {
    for (int i = 0; i < arity(a); ++i) {
      for (int j = 0; j < arity(b); ++j) {
        if (a.qubits[static_cast<std::size_t>(i)] ==
            b.qubits[static_cast<std::size_t>(j)]) {
          return true;
        }
      }
    }
    return false;
  }

  // Erases open_[indices[i]] for i < count, highest index first.
  void remove(const std::array<std::size_t, 2>& indices, std::size_t count) {
    std::array<std::size_t, 2> sorted = indices;
    if (count == 2 && sorted[0] < sorted[1]) {
      std::swap(sorted[0], sorted[1]);
    }
    for (std::size_t i = 0; i < count; ++i) {
      open_.erase(open_.begin() + static_cast<std::ptrdiff_t>(sorted[i]));
    }
  }

  std::vector<Open> open_;
  std::vector<FusedGate> done_;
  std::size_t next_created_ = 0;
};

}  // namespace

void validate_pauli_errors(const Circuit& circuit,
                           const std::vector<PauliError>& errors) {
  const std::size_t n_steps = circuit.timesteps.size();
  for (std::size_t i = 0; i < errors.size(); ++i) {
    const PauliError& error = errors[i];
    if (error.timestep >= n_steps) {
      throw std::out_of_range(
          "error timestep " + std::to_string(error.timestep) +
          " out of range for " + std::to_string(n_steps) + " time steps");
    }
    if (i > 0 && error.timestep < errors[i - 1].timestep) {
      throw std::invalid_argument("errors must be sorted by timestep");
    }
    if (error.pauli != GateType::X && error.pauli != GateType::Y &&
        error.pauli != GateType::Z) {
      throw std::invalid_argument("error pauli must be X, Y or Z");
    }
  }
}

std::vector<FusedGate> fuse_trajectory(const Circuit& circuit,
                                       const std::vector<PauliError>& errors) {
  validate_pauli_errors(circuit, errors);
  std::size_t n_operations = errors.size();
  for (const TimeStep& step : circuit.timesteps) {
    n_operations += step.ops.size();
  }
  Fuser fuser(n_operations);
  std::size_t next = 0;
  for (std::size_t t = 0; t < circuit.timesteps.size(); ++t) {
    for (const Op& op : circuit.timesteps[t].ops) {
      fuser.add(gate_of(op, circuit.n_qubits));
    }
    for (; next < errors.size() && errors[next].timestep == t; ++next) {
      fuser.add(gate_of(errors[next], circuit.n_qubits));
    }
  }
  return fuser.finish();
}

void apply_fused(StateVector& state, const FusedGate& gate) {
  const int j = gate.qubits[0];
  const int k = gate.qubits[1];
  std::visit(
      [&state, j, k](const auto& m) {
        using Matrix = std::decay_t<decltype(m)>;
        if constexpr (std::is_same_v<Matrix, Eigen::Matrix2d>) {
          state.apply_real_1q(m, j);
        } else if constexpr (std::is_same_v<Matrix, Eigen::Matrix4d>) {
          state.apply_real_2q(m, j, k);
        } else if constexpr (std::is_same_v<Matrix, Eigen::Matrix2cd>) {
          state.apply_1q(m, j);
        } else {
          state.apply_2q(m, j, k);
        }
      },
      gate.matrix);
}

}  // namespace ecm_qsim
