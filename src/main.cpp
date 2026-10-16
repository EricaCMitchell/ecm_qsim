// sim CLI: parse arguments, call run(), print counts.

#include <iostream>
#include <string>
#include <vector>

#include "ecm_qsim/cli.hpp"

int main(int argc, char* argv[]) {
  return ecm_qsim::sim_main(std::vector<std::string>(argv + 1, argv + argc),
                            std::cout, std::cerr);
}
