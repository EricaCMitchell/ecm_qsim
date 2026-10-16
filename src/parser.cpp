#include "ecm_qsim/parser.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
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

// Decimal syntax accepted by the circuit format. Hexadecimal floats and
// spellings such as nan and inf are not part of the format.
bool has_decimal_syntax(std::string_view token) {
  std::size_t i = 0;
  if (i < token.size() && (token[i] == '+' || token[i] == '-')) {
    ++i;
  }

  bool has_digits = false;
  while (i < token.size() && token[i] >= '0' && token[i] <= '9') {
    has_digits = true;
    ++i;
  }
  if (i < token.size() && token[i] == '.') {
    ++i;
    while (i < token.size() && token[i] >= '0' && token[i] <= '9') {
      has_digits = true;
      ++i;
    }
  }
  if (!has_digits) {
    return false;
  }

  if (i < token.size() && (token[i] == 'e' || token[i] == 'E')) {
    ++i;
    if (i < token.size() && (token[i] == '+' || token[i] == '-')) {
      ++i;
    }
    const std::size_t exponent_start = i;
    while (i < token.size() && token[i] >= '0' && token[i] <= '9') {
      ++i;
    }
    if (i == exponent_start) {
      return false;
    }
  }
  return i == token.size();
}

// Converts a validated decimal token and rejects overflow to infinity.
// strtod is used because libc++ floating-point from_chars is recent.
std::optional<double> parse_angle(std::string_view token) {
  if (!has_decimal_syntax(token)) {
    return std::nullopt;
  }
  const std::string text(token);
  char* end = nullptr;
  const double value = std::strtod(text.c_str(), &end);
  if (end != text.c_str() + text.size() || !std::isfinite(value)) {
    return std::nullopt;
  }
  return value;
}

struct GateSpec {
  std::string_view name;
  GateType type;
  bool has_angle;
  int n_qubits;
};

constexpr std::array<GateSpec, 12> kGates{{{"X", GateType::X, false, 1},
                                           {"Y", GateType::Y, false, 1},
                                           {"Z", GateType::Z, false, 1},
                                           {"H", GateType::H, false, 1},
                                           {"S", GateType::S, false, 1},
                                           {"T", GateType::T, false, 1},
                                           {"RX", GateType::RX, true, 1},
                                           {"RY", GateType::RY, true, 1},
                                           {"RZ", GateType::RZ, true, 1},
                                           {"CNOT", GateType::CNOT, false, 2},
                                           {"CZ", GateType::CZ, false, 2},
                                           {"CY", GateType::CY, false, 2}}};

const GateSpec* find_gate(std::string_view name) {
  for (const GateSpec& spec : kGates) {
    if (spec.name == name) {
      return &spec;
    }
  }
  return nullptr;
}

// Usage shown in arity errors, e.g. "RX theta q" or "CNOT control target".
std::string usage(const GateSpec& spec) {
  std::string text(spec.name);
  if (spec.has_angle) {
    text += " theta";
  }
  text += spec.n_qubits == 1 ? " q" : " control target";
  return text;
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
    // A trailing TICK already closed the last step. Otherwise EOF closes the
    // current step, including an empty initial step after QUBITS.
    if (close_step_at_eof_) {
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
    used_on_line_.assign(circuit_.n_qubits, 0);
    header_line_ = line_;
    close_step_at_eof_ = true;
  }

  void parse_statement(const Tokens& tokens) {
    if (tokens[0] == "QUBITS") {
      fail("QUBITS already declared on line " + std::to_string(header_line_));
    }
    if (tokens[0] == "TICK") {
      expect_arguments(tokens, 0, "TICK");
      circuit_.timesteps.push_back(std::move(current_));
      current_ = TimeStep{};
      std::fill(used_on_line_.begin(), used_on_line_.end(), 0);
      close_step_at_eof_ = false;
      return;
    }
    parse_gate(tokens);
  }

  void parse_gate(const Tokens& tokens) {
    const GateSpec* spec = find_gate(tokens[0]);
    if (spec == nullptr) {
      fail("unknown gate " + quoted(tokens[0]));
    }
    const std::size_t n_args = (spec->has_angle ? 1 : 0) + spec->n_qubits;
    expect_arguments(tokens, n_args, usage(*spec));

    Op op{spec->type, {0, 0}, 0.0};
    std::size_t next = 1;
    if (spec->has_angle) {
      op.theta = parse_theta(tokens[next++]);
    }
    for (int i = 0; i < spec->n_qubits; ++i) {
      op.qubits[i] = parse_qubit(tokens[next++]);
    }
    if (spec->n_qubits == 2 && op.qubits[0] == op.qubits[1]) {
      fail(std::string(spec->name) + " control and target must differ, got " +
           std::to_string(op.qubits[0]) + " and " +
           std::to_string(op.qubits[1]));
    }
    for (int i = 0; i < spec->n_qubits; ++i) {
      claim(op.qubits[i]);
    }
    current_.ops.push_back(op);
    close_step_at_eof_ = true;
  }

  double parse_theta(std::string_view token) const {
    const auto theta = parse_angle(token);
    if (!theta) {
      fail("angle must be a finite decimal number, got " + quoted(token));
    }
    return *theta;
  }

  int parse_qubit(std::string_view token) const {
    const auto q = parse_integer(token);
    if (!q || *q < 0) {
      fail("qubit index must be a non-negative integer, got " + quoted(token));
    }
    if (*q >= circuit_.n_qubits) {
      fail("qubit index " + quoted(token) + " out of range for " +
           std::to_string(circuit_.n_qubits) + " qubits");
    }
    return static_cast<int>(*q);
  }

  // Marks qubit as acted on in the current time step.
  void claim(int qubit) {
    int& used = used_on_line_[qubit];
    if (used != 0) {
      fail("qubit " + std::to_string(qubit) +
           " already used in this time step (line " + std::to_string(used) +
           ")");
    }
    used = line_;
  }

  std::istream& in_;
  int line_ = 0;
  int header_line_ = 0;  // 0 until QUBITS is seen
  bool close_step_at_eof_ = false;
  Circuit circuit_{};
  TimeStep current_;
  std::vector<int> used_on_line_;  // per qubit: line of its gate this step
};

}  // namespace

Circuit parse_circuit(std::istream& in) {
  return Parser(in).parse();
}

Circuit parse_circuit_file(const std::filesystem::path& path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("cannot open circuit file '" + path.string() +
                             "'");
  }
  return parse_circuit(in);
}

}  // namespace ecm_qsim
