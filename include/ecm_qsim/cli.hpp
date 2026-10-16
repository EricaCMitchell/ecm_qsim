#ifndef ECM_QSIM_CLI_HPP_
#define ECM_QSIM_CLI_HPP_

// The sim command line as a function, so tests can call it in-process.

#include <ostream>
#include <string>
#include <vector>

namespace ecm_qsim {

// The sim command line. args excludes the program name. Writes the counts
// to out and messages to err, prefixed "sim: ". Returns the exit status:
// 0 on success, 1 if the circuit cannot be read, parsed or run, 2 for a
// usage error. Catches every std::exception.
int sim_main(const std::vector<std::string>& args, std::ostream& out,
             std::ostream& err);

}  // namespace ecm_qsim

#endif  // ECM_QSIM_CLI_HPP_
