#include <catch2/catch_test_macros.hpp>
#include <initializer_list>

#include "ecm_qsim/circuit.hpp"

using namespace ecm_qsim;

TEST_CASE("Circuit holds empty timesteps") {
  Circuit c{3, {TimeStep{}, TimeStep{}}};
  REQUIRE(c.n_qubits == 3);
  REQUIRE(c.timesteps.size() == 2);
  REQUIRE(c.timesteps[0].ops.empty());
}

TEST_CASE("is_two_qubit is true only for controlled gates") {
  for (GateType g : {GateType::CNOT, GateType::CZ, GateType::CY}) {
    CHECK(is_two_qubit(g));
  }
  for (GateType g :
       {GateType::X, GateType::Y, GateType::Z, GateType::H, GateType::S,
        GateType::T, GateType::RX, GateType::RY, GateType::RZ}) {
    CHECK_FALSE(is_two_qubit(g));
  }
}
