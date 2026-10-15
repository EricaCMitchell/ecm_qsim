#ifndef ECM_QSIM_PARSER_HPP_
#define ECM_QSIM_PARSER_HPP_

// Reads the circuit text format into a Circuit.

#include <filesystem>
#include <istream>
#include <stdexcept>
#include <string>

#include "ecm_qsim/circuit.hpp"

namespace ecm_qsim {

// Malformed circuit text. what() is "line N: message"; line() is 1-based and
// counts every physical line, including comments and blank lines.
class ParseError : public std::runtime_error {
public:
  ParseError(int line, const std::string& message);
  int line() const noexcept { return line_; }

private:
  int line_;
};

// Parses a whole circuit. Throws ParseError on malformed input and
// std::runtime_error if the input stream fails before normal EOF.
Circuit parse_circuit(std::istream& in);

// Opens and parses a circuit file. Throws std::runtime_error if the file
// cannot be opened or read, ParseError if it is malformed.
Circuit parse_circuit_file(const std::filesystem::path& path);

}  // namespace ecm_qsim

#endif  // ECM_QSIM_PARSER_HPP_
