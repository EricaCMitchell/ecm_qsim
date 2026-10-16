#include "ecm_qsim/noise.hpp"

#include <algorithm>
#include <charconv>
#include <limits>
#include <stdexcept>
#include <string>

#include "ecm_qsim/measurement.hpp"

namespace ecm_qsim {

namespace {

void check_probability(const char* name, double p) {
  // Written so that NaN fails the check.
  if (!(p >= 0.0 && p <= 1.0)) {
    // Shortest round-trip form, so a rejected 1.0000001 is not shown as 1.
    char buf[32];
    const auto result = std::to_chars(buf, buf + sizeof buf, p);
    throw std::invalid_argument(std::string(name) + " must be in [0, 1], got " +
                                std::string(buf, result.ptr));
  }
}

void check_qubit(int qubit, int n_qubits) {
  if (qubit < 0 || qubit >= n_qubits) {
    throw std::out_of_range("qubit " + std::to_string(qubit) +
                            " out of range for " + std::to_string(n_qubits) +
                            " qubits");
  }
}

// One occurrence draw; no draw at all when p is 0. u is in [0, 1), so p = 1
// always fires.
bool occurs(double p, Rng& rng) {
  return p != 0.0 && uniform_from_bits(rng()) < p;
}

// Reject UINT64_MAX, then map bits % 3 values 0, 1, 2 to X, Y, Z. The
// accepted range has 2^64 - 1 values, a multiple of 3, so the choice is
// exactly uniform.
GateType random_pauli(Rng& rng) {
  std::uint64_t bits = rng();
  while (bits == std::numeric_limits<std::uint64_t>::max()) {
    bits = rng();
  }
  switch (bits % 3) {
    case 0:
      return GateType::X;
    case 1:
      return GateType::Y;
    default:
      return GateType::Z;
  }
}

}  // namespace

void validate_noise_model(const NoiseModel& noise) {
  check_probability("p_gate", noise.p_gate);
  check_probability("p_idle", noise.p_idle);
  check_probability("p_readout", noise.p_readout);
}

std::vector<PauliError> sample_pauli_errors(const Circuit& circuit,
                                            const NoiseModel& noise, Rng& rng) {
  // Check everything before the first draw, so a failure leaves rng
  // unchanged.
  validate_noise_model(noise);
  const int n = circuit.n_qubits;
  if (n < 1 || n > kMaxQubits) {
    throw std::invalid_argument("circuit.n_qubits must be in [1, " +
                                std::to_string(kMaxQubits) + "], got " +
                                std::to_string(n));
  }
  for (const TimeStep& step : circuit.timesteps) {
    for (const Op& op : step.ops) {
      check_qubit(op.qubits[0], n);
      if (is_two_qubit(op.gate)) {
        check_qubit(op.qubits[1], n);
      }
    }
  }

  std::vector<PauliError> errors;
  std::vector<char> acted(static_cast<std::size_t>(n), 0);
  for (std::size_t t = 0; t < circuit.timesteps.size(); ++t) {
    std::fill(acted.begin(), acted.end(), 0);
    for (const Op& op : circuit.timesteps[t].ops) {
      acted[static_cast<std::size_t>(op.qubits[0])] = 1;
      // qubits[1] of a one-qubit gate is 0 and unused. Marking it would give
      // qubit 0 gate noise instead of idle noise.
      if (is_two_qubit(op.gate)) {
        acted[static_cast<std::size_t>(op.qubits[1])] = 1;
      }
    }
    for (int q = 0; q < n; ++q) {
      const double p =
          acted[static_cast<std::size_t>(q)] ? noise.p_gate : noise.p_idle;
      if (occurs(p, rng)) {
        errors.push_back(PauliError{t, q, random_pauli(rng)});
      }
    }
  }
  return errors;
}

}  // namespace ecm_qsim
