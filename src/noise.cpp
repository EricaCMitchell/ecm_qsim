#include "ecm_qsim/noise.hpp"

#include <charconv>
#include <stdexcept>
#include <string>

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

}  // namespace

void validate_noise_model(const NoiseModel& noise) {
  check_probability("p_gate", noise.p_gate);
  check_probability("p_idle", noise.p_idle);
  check_probability("p_readout", noise.p_readout);
}

}  // namespace ecm_qsim
