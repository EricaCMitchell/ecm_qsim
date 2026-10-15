#include "ecm_qsim/parser.hpp"

#include <cctype>
#include <charconv>
#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

namespace ecm_qsim {

ParseError::ParseError(int line, const std::string& message)
    : std::runtime_error("line " + std::to_string(line) + ": " + message),
      line_(line) {}

namespace {

using Tokens = std::vector<std::string_view>;

// Splits on any whitespace, so leading/trailing blanks, tabs and a CR from a
// CRLF line ending never reach a token.
Tokens split_whitespace(std::string_view line) {
  Tokens tokens;
  std::size_t i = 0;
  while (i < line.size()) {
    while (i < line.size() &&
           std::isspace(static_cast<unsigned char>(line[i]))) {
      ++i;
    }
    const std::size_t start = i;
    while (i < line.size() &&
           !std::isspace(static_cast<unsigned char>(line[i]))) {
      ++i;
    }
    if (i > start) {
      tokens.push_back(line.substr(start, i - start));
    }
  }
  return tokens;
}

std::string quoted(std::string_view token) {
  return "'" + std::string(token) + "'";
}

// Parses a base-10 integer that fills the whole token. Values too large for
// long long saturate, so callers report them as out of range rather than as
// malformed.
std::optional<long long> parse_integer(std::string_view token) {
  long long value = 0;
  const char* end = token.data() + token.size();
  const auto [ptr, ec] = std::from_chars(token.data(), end, value);
  if (ptr != end) {
    return std::nullopt;
  }
  if (ec == std::errc::result_out_of_range) {
    return token.front() == '-' ? std::numeric_limits<long long>::min()
                                : std::numeric_limits<long long>::max();
  }
  if (ec != std::errc{}) {
    return std::nullopt;
  }
  return value;
}

class Parser {
public:
  explicit Parser(std::istream& in) : in_(in) {}

  Circuit parse() {
    std::string text;
    while (std::getline(in_, text)) {
      ++line_;
      const Tokens tokens = split_whitespace(text);
      if (tokens.empty() || tokens[0].front() == '#') {
        continue;
      }
      if (header_line_ == 0) {
        parse_header(tokens);
      } else {
        parse_statement(tokens);
      }
    }
    if (in_.bad() || !in_.eof()) {
      throw std::runtime_error("error reading circuit input");
    }
    if (header_line_ == 0) {
      throw ParseError(line_ + 1, "missing 'QUBITS n' declaration");
    }
    // A trailing TICK already closed the last step; only an open step with
    // gates is closed by end of file, so no phantom all-idle step is added.
    if (!current_.ops.empty()) {
      circuit_.timesteps.push_back(std::move(current_));
    }
    return std::move(circuit_);
  }

private:
  [[noreturn]] void fail(const std::string& message) const {
    throw ParseError(line_, message);
  }

  void expect_arguments(const Tokens& tokens, std::size_t count,
                        const std::string& usage) const {
    const std::size_t got = tokens.size() - 1;
    if (got != count) {
      fail(std::string(tokens[0]) + " expects " + std::to_string(count) +
           (count == 1 ? " argument (" : " arguments (") + usage + "), got " +
           std::to_string(got));
    }
  }

  void parse_header(const Tokens& tokens) {
    if (tokens[0] != "QUBITS") {
      fail("expected 'QUBITS n' before any gate or TICK, got " +
           quoted(tokens[0]));
    }
    expect_arguments(tokens, 1, "QUBITS n");
    const auto n = parse_integer(tokens[1]);
    if (!n || *n <= 0) {
      fail("qubit count must be a positive integer, got " + quoted(tokens[1]));
    }
    if (*n > kMaxQubits) {
      fail("qubit count " + quoted(tokens[1]) + " exceeds the maximum of " +
           std::to_string(kMaxQubits));
    }
    circuit_.n_qubits = static_cast<int>(*n);
    header_line_ = line_;
  }

  void parse_statement(const Tokens& tokens) {
    if (tokens[0] == "QUBITS") {
      fail("QUBITS already declared on line " + std::to_string(header_line_));
    }
    if (tokens[0] == "TICK") {
      expect_arguments(tokens, 0, "TICK");
      circuit_.timesteps.push_back(std::move(current_));
      current_ = TimeStep{};
      return;
    }
    fail("unknown gate " + quoted(tokens[0]));
  }

  std::istream& in_;
  int line_ = 0;
  int header_line_ = 0;  // 0 until QUBITS is seen
  Circuit circuit_{};
  TimeStep current_;
};

}  // namespace

Circuit parse_circuit(std::istream& in) {
  return Parser(in).parse();
}

}  // namespace ecm_qsim
