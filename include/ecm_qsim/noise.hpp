#ifndef ECM_QSIM_NOISE_HPP_
#define ECM_QSIM_NOISE_HPP_

// Gate, idle and readout noise: trajectory sampling and readout flips.

namespace ecm_qsim {

struct NoiseModel {
  double p_gate = 0.0;
  double p_idle = 0.0;
  double p_readout = 0.0;
};

}  // namespace ecm_qsim

#endif  // ECM_QSIM_NOISE_HPP_
