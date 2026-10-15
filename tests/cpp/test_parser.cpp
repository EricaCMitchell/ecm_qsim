#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <optional>
#include <sstream>
#include <streambuf>
#include <string>
#include <string_view>

#include "ecm_qsim/circuit.hpp"
#include "ecm_qsim/parser.hpp"

using namespace ecm_qsim;
using Catch::Matchers::ContainsSubstring;

namespace {

class FailingStreambuf : public std::streambuf {
public:
  explicit FailingStreambuf(std::string text) : text_(text) {
    setg(text_.data(), text_.data(), text_.data() + text_.size());
  }

protected:
  int_type underflow() override {
    throw std::runtime_error("simulated read failure");
  }

private:
  std::string text_;
};

Circuit parse(std::string_view text) {
  std::istringstream in{std::string(text)};
  return parse_circuit(in);
}

// Returns the ParseError thrown for text, or nullopt if it parses.
std::optional<ParseError> parse_error(std::string_view text) {
  try {
    parse(text);
  } catch (const ParseError& e) {
    return e;
  }
  return std::nullopt;
}

void require_error(std::string_view text, int line, std::string_view message) {
  CAPTURE(text);
  const auto err = parse_error(text);
  REQUIRE(err.has_value());
  CHECK(err->line() == line);
  CHECK_THAT(std::string(err->what()), ContainsSubstring(std::string(message)));
}

}  // namespace

TEST_CASE("QUBITS declares the qubit count") {
  const Circuit c = parse("QUBITS 3\n");
  CHECK(c.n_qubits == 3);
  CHECK(c.timesteps.empty());
}

TEST_CASE("QUBITS accepts the maximum qubit count") {
  CHECK(parse("QUBITS 63\n").n_qubits == kMaxQubits);
}

TEST_CASE("Comments and blank lines are skipped, including before QUBITS") {
  const Circuit c = parse(
      "\n"
      "# header comment\n"
      "   # indented comment\n"
      "\t\n"
      "QUBITS 2\n"
      "\n"
      "#no space after hash\n"
      "TICK\n");
  CHECK(c.n_qubits == 2);
  REQUIRE(c.timesteps.size() == 1);
  CHECK(c.timesteps[0].ops.empty());
}

TEST_CASE("Consecutive TICKs create empty time steps") {
  const Circuit c = parse("QUBITS 1\nTICK\nTICK\nTICK\n");
  REQUIRE(c.timesteps.size() == 3);
  for (const TimeStep& step : c.timesteps) {
    CHECK(step.ops.empty());
  }
}

TEST_CASE("CRLF line endings, tabs and missing final newline are accepted") {
  const Circuit c = parse("QUBITS\t2  \r\n  TICK\r\nTICK");
  CHECK(c.n_qubits == 2);
  CHECK(c.timesteps.size() == 2);
}

TEST_CASE("A stream read failure does not return a partial circuit") {
  FailingStreambuf buffer("QUBITS 1\nTICK\n");
  std::istream in(&buffer);
  CHECK_THROWS_WITH(parse_circuit(in), "error reading circuit input");
}

TEST_CASE("ParseError message is prefixed with the line number") {
  const auto err = parse_error("QUBITS 2\nQUBITS 2\n");
  REQUIRE(err.has_value());
  CHECK(err->line() == 2);
  CHECK(std::string(err->what()) ==
        "line 2: QUBITS already declared on line 1");
}

TEST_CASE("Missing QUBITS is rejected at end of input") {
  require_error("", 1, "missing 'QUBITS n' declaration");
  require_error("# only a comment\n\n", 3, "missing 'QUBITS n' declaration");
}

TEST_CASE("QUBITS must come before any other statement") {
  require_error("TICK\nQUBITS 2\n", 1, "expected 'QUBITS n' before any gate");
  require_error("# c\n\nH 0\n", 3, "got 'H'");
  require_error("qubits 2\n", 1, "expected 'QUBITS n'");
}

TEST_CASE("QUBITS repeated later in the file is rejected") {
  require_error("QUBITS 2\nTICK\n# c\nQUBITS 3\n", 4,
                "QUBITS already declared on line 1");
}

TEST_CASE("QUBITS with the wrong number of arguments is rejected") {
  require_error("QUBITS\n", 1, "QUBITS expects 1 argument (QUBITS n), got 0");
  require_error("QUBITS 2 3\n", 1,
                "QUBITS expects 1 argument (QUBITS n), got 2");
}

TEST_CASE("Qubit count must be a positive integer") {
  const std::string count =
      GENERATE(as<std::string>{}, "0", "-1", "2.5", "two", "+2", "0x2", "2q");
  require_error("QUBITS " + count + "\n", 1,
                "qubit count must be a positive integer, got '" + count + "'");
}

TEST_CASE("Qubit count above the maximum is rejected") {
  const std::string count =
      GENERATE(as<std::string>{}, "64", "100", "99999999999999999999");
  require_error("QUBITS " + count + "\n", 1,
                "qubit count '" + count + "' exceeds the maximum of 63");
}

TEST_CASE("TICK takes no arguments") {
  require_error("QUBITS 2\nTICK 1\n", 2,
                "TICK expects 0 arguments (TICK), got 1");
}
