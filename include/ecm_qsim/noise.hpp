#pragma once

// Gate, idle and readout noise: trajectory sampling and readout flips.

namespace ecm_qsim {

struct NoiseModel {
  double p_gate = 0.0;
  double p_idle = 0.0;
  double p_readout = 0.0;
};

}  // namespace ecm_qsim
