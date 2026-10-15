include(FetchContent)
FetchContent_Declare(Catch2
  URL https://github.com/catchorg/Catch2/archive/refs/tags/v3.16.1.tar.gz
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  FIND_PACKAGE_ARGS 3
)
FetchContent_MakeAvailable(Catch2)

add_executable(unit_tests
  test_circuit.cpp
  test_gates.cpp
  test_parser.cpp
)
target_link_libraries(unit_tests PRIVATE ecm_qsim_core Catch2::Catch2WithMain)
target_compile_definitions(unit_tests PRIVATE
  ECM_QSIM_EXAMPLES_DIR="${PROJECT_SOURCE_DIR}/examples"
)
target_compile_options(unit_tests PRIVATE -Wall -Wextra -Wpedantic)

include(Catch)
catch_discover_tests(unit_tests)
