#include <catch2/catch_test_macros.hpp>

#include "ecm_qsim/circuit.hpp"
#include "ecm_qsim/noise.hpp"

using namespace ecm_qsim;

TEST_CASE("Circuit holds empty timesteps") {
  Circuit c{3, {TimeStep{}, TimeStep{}}};
  REQUIRE(c.n_qubits == 3);
  REQUIRE(c.timesteps.size() == 2);
  REQUIRE(c.timesteps[0].ops.empty());
}

TEST_CASE("NoiseModel defaults to noiseless") {
  NoiseModel noise;
  REQUIRE(noise.p_gate == 0.0);
  REQUIRE(noise.p_idle == 0.0);
  REQUIRE(noise.p_readout == 0.0);
}
