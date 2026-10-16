#include <catch2/catch_test_macros.hpp>

#include "ecm_qsim/circuit.hpp"

using namespace ecm_qsim;

TEST_CASE("Circuit holds empty timesteps") {
  Circuit c{3, {TimeStep{}, TimeStep{}}};
  REQUIRE(c.n_qubits == 3);
  REQUIRE(c.timesteps.size() == 2);
  REQUIRE(c.timesteps[0].ops.empty());
}
