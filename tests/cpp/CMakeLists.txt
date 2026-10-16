include(FetchContent)
FetchContent_Declare(Catch2
  URL https://github.com/catchorg/Catch2/archive/refs/tags/v3.16.1.tar.gz
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  FIND_PACKAGE_ARGS 3
)
FetchContent_MakeAvailable(Catch2)

add_executable(unit_tests
  test_circuit.cpp
  test_cli.cpp
  test_gates.cpp
  test_measurement.cpp
  test_noise.cpp
  test_parser.cpp
  test_simulator.cpp
  test_state_vector.cpp
)
target_link_libraries(unit_tests PRIVATE ecm_qsim_cli Catch2::Catch2WithMain)
target_compile_definitions(unit_tests PRIVATE
  ECM_QSIM_EXAMPLES_DIR="${PROJECT_SOURCE_DIR}/examples"
)
target_compile_options(unit_tests PRIVATE -Wall -Wextra -Wpedantic)

include(Catch)
catch_discover_tests(unit_tests)

add_test(NAME sim_example_command
  COMMAND sim ${PROJECT_SOURCE_DIR}/examples/ghz3.txt --shots 1000
          --p-gate 0.001 --p-idle 0.001 --p-readout 0.01 --seed 1)
set_tests_properties(sim_example_command PROPERTIES
  PASS_REGULAR_EXPRESSION "^000 [0-9]+\n")
add_test(NAME sim_usage_error COMMAND sim)
set_tests_properties(sim_usage_error PROPERTIES WILL_FAIL TRUE)
