cmake_minimum_required(VERSION 3.24)
project(ecm_qsim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

find_package(Eigen3 REQUIRED NO_MODULE)

add_library(ecm_qsim_core
  src/parser.cpp
  src/gates.cpp
  src/state_vector.cpp
  src/noise.cpp
  src/simulator.cpp
)
target_include_directories(ecm_qsim_core PUBLIC include)
target_link_libraries(ecm_qsim_core PUBLIC Eigen3::Eigen)
target_compile_options(ecm_qsim_core PRIVATE -Wall -Wextra -Wpedantic)

add_executable(sim src/main.cpp)
target_link_libraries(sim PRIVATE ecm_qsim_core)
target_compile_options(sim PRIVATE -Wall -Wextra -Wpedantic)

option(ECM_QSIM_BUILD_TESTS "Build the C++ unit tests" ON)
if(ECM_QSIM_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests/cpp)
endif()
