#ifndef ECM_QSIM_NOISE_HPP_
#define ECM_QSIM_NOISE_HPP_

// Gate, idle and readout noise: trajectory sampling and readout flips.

namespace ecm_qsim {

struct NoiseModel {
  double p_gate = 0.0;
  double p_idle = 0.0;
  double p_readout = 0.0;
};

// Throws std::invalid_argument("<field> must be in [0, 1], got <value>") for
// the first field, in the order p_gate, p_idle, p_readout, that is not in
// [0, 1] (NaN included).
void validate_noise_model(const NoiseModel& noise);

}  // namespace ecm_qsim

#endif  // ECM_QSIM_NOISE_HPP_
