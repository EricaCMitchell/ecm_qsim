#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <optional>
#include <sstream>
#include <streambuf>
#include <string>
#include <string_view>
#include <vector>

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

namespace Catch {
template <>
struct StringMaker<ecm_qsim::Op> {
  static std::string convert(const ecm_qsim::Op& op) {
    return "Op{gate=" + std::to_string(static_cast<int>(op.gate)) +
           ", qubits={" + std::to_string(op.qubits[0]) + ", " +
           std::to_string(op.qubits[1]) +
           "}, theta=" + StringMaker<double>::convert(op.theta) + "}";
  }
};
}  // namespace Catch

TEST_CASE("QUBITS declares the qubit count") {
  const Circuit c = parse("QUBITS 3\n");
  CHECK(c.n_qubits == 3);
}

TEST_CASE("End of file closes an empty initial time step") {
  const Circuit c = parse("QUBITS 3\n");
  REQUIRE(c.timesteps.size() == 1);
  CHECK(c.timesteps[0].ops.empty());
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

TEST_CASE("Every gate parses to its type, angle and qubits") {
  const Circuit c = parse(
      "QUBITS 3\n"
      "X 0\n"
      "Y 1\n"
      "Z 2\n"
      "TICK\n"
      "H 0\n"
      "S 1\n"
      "T 2\n"
      "TICK\n"
      "RX 0.5 0\n"
      "RY -1.25 1\n"
      "RZ 1e-3 2\n"
      "TICK\n"
      "CNOT 0 1\n"
      "TICK\n"
      "CZ 2 1\n"
      "TICK\n"
      "CY 1 0\n");
  REQUIRE(c.timesteps.size() == 6);
  CHECK(c.timesteps[0].ops == std::vector<Op>{{GateType::X, {0, 0}, 0.0},
                                              {GateType::Y, {1, 0}, 0.0},
                                              {GateType::Z, {2, 0}, 0.0}});
  CHECK(c.timesteps[1].ops == std::vector<Op>{{GateType::H, {0, 0}, 0.0},
                                              {GateType::S, {1, 0}, 0.0},
                                              {GateType::T, {2, 0}, 0.0}});
  CHECK(c.timesteps[2].ops == std::vector<Op>{{GateType::RX, {0, 0}, 0.5},
                                              {GateType::RY, {1, 0}, -1.25},
                                              {GateType::RZ, {2, 0}, 1e-3}});
  CHECK(c.timesteps[3].ops == std::vector<Op>{{GateType::CNOT, {0, 1}, 0.0}});
  CHECK(c.timesteps[4].ops == std::vector<Op>{{GateType::CZ, {2, 1}, 0.0}});
  CHECK(c.timesteps[5].ops == std::vector<Op>{{GateType::CY, {1, 0}, 0.0}});
}

TEST_CASE("End of file closes the open time step; a trailing TICK adds none") {
  CHECK(parse("QUBITS 1\nH 0\n").timesteps.size() == 1);
  CHECK(parse("QUBITS 1\nH 0\nTICK\n").timesteps.size() == 1);
  CHECK(parse("QUBITS 1\nH 0\nTICK\n# trailing comment\n\n").timesteps.size() ==
        1);

  const Circuit c = parse("QUBITS 1\nH 0\nTICK\nTICK\n");
  REQUIRE(c.timesteps.size() == 2);
  CHECK(c.timesteps[1].ops.empty());
}

TEST_CASE("A qubit may be reused after TICK") {
  const Circuit c = parse("QUBITS 2\nCNOT 0 1\nTICK\nH 0\nX 1\nTICK\nX 0\n");
  REQUIRE(c.timesteps.size() == 3);
  CHECK(c.timesteps[1].ops.size() == 2);
  CHECK(c.timesteps[2].ops == std::vector<Op>{{GateType::X, {0, 0}, 0.0}});
}

TEST_CASE("Extra whitespace between tokens is accepted") {
  const Circuit c = parse("QUBITS 2\n  CNOT\t0    1  \n");
  REQUIRE(c.timesteps.size() == 1);
  CHECK(c.timesteps[0].ops == std::vector<Op>{{GateType::CNOT, {0, 1}, 0.0}});
}

TEST_CASE("Unknown and lowercase gate names are rejected") {
  require_error("QUBITS 2\nFOO 0\n", 2, "unknown gate 'FOO'");
  require_error("QUBITS 2\nh 0\n", 2, "unknown gate 'h'");
  require_error("QUBITS 2\nCX 0 1\n", 2, "unknown gate 'CX'");
}

TEST_CASE("Gates with the wrong number of arguments are rejected") {
  require_error("QUBITS 3\nH\n", 2, "H expects 1 argument (H q), got 0");
  require_error("QUBITS 3\nH 0 1\n", 2, "H expects 1 argument (H q), got 2");
  require_error("QUBITS 3\nRX 0\n", 2,
                "RX expects 2 arguments (RX theta q), got 1");
  require_error("QUBITS 3\nCNOT 0\n", 2,
                "CNOT expects 2 arguments (CNOT control target), got 1");
  require_error("QUBITS 3\nCZ 0 1 2\n", 2,
                "CZ expects 2 arguments (CZ control target), got 3");
}

TEST_CASE("A trailing comment on a gate line is an argument-count error") {
  require_error("QUBITS 1\nH 0  # hadamard\n", 2,
                "H expects 1 argument (H q), got 3");
}

TEST_CASE("Angles must be finite decimal numbers") {
  const std::string angle =
      GENERATE(as<std::string>{}, "pi", "1.0rad", "0x1p0", "0x1.8p+1", "nan",
               "inf", "-inf", "1e999", "--1");
  require_error("QUBITS 1\nRY " + angle + " 0\n", 2,
                "angle must be a finite decimal number, got '" + angle + "'");
}

TEST_CASE("Qubit indices must be non-negative integers") {
  const std::string index =
      GENERATE(as<std::string>{}, "-1", "1.0", "q0", "+1", "0x1");
  require_error(
      "QUBITS 3\nH " + index + "\n", 2,
      "qubit index must be a non-negative integer, got '" + index + "'");
}

TEST_CASE("Qubit indices must be below the qubit count") {
  require_error("QUBITS 3\nH 3\n", 2,
                "qubit index '3' out of range for 3 qubits");
  require_error("QUBITS 3\nCNOT 0 7\n", 2,
                "qubit index '7' out of range for 3 qubits");
  require_error("QUBITS 3\nX 4294967296\n", 2,
                "qubit index '4294967296' out of range for 3 qubits");
  require_error("QUBITS 3\nX 99999999999999999999\n", 2,
                "qubit index '99999999999999999999' out of range for 3 qubits");
}

TEST_CASE("Two-qubit gates need distinct control and target") {
  const std::string gate = GENERATE(as<std::string>{}, "CNOT", "CZ", "CY");
  require_error("QUBITS 3\n" + gate + " 1 1\n", 2,
                gate + " control and target must differ, got 1 and 1");
}

TEST_CASE("Two gates on one qubit in a time step are rejected") {
  const auto err = parse_error("QUBITS 3\nH 0\n# comment\nCNOT 1 0\n");
  REQUIRE(err.has_value());
  CHECK(err->line() == 4);
  CHECK(std::string(err->what()) ==
        "line 4: qubit 0 already used in this time step (line 2)");

  require_error("QUBITS 1\nTICK\nX 0\nX 0\n", 4,
                "qubit 0 already used in this time step (line 3)");
}

TEST_CASE("The GHZ example file parses") {
  const Circuit c = parse_circuit_file(ECM_QSIM_EXAMPLES_DIR "/ghz3.txt");
  CHECK(c.n_qubits == 3);
  REQUIRE(c.timesteps.size() == 4);
  CHECK(c.timesteps[0].ops == std::vector<Op>{{GateType::H, {0, 0}, 0.0}});
  CHECK(c.timesteps[1].ops == std::vector<Op>{{GateType::CNOT, {0, 1}, 0.0}});
  CHECK(c.timesteps[2].ops == std::vector<Op>{{GateType::CNOT, {1, 2}, 0.0}});
  CHECK(c.timesteps[3].ops == std::vector<Op>{{GateType::RZ, {0, 0}, 0.785}});
}

TEST_CASE("A missing circuit file is reported with its path") {
  const std::string path = ECM_QSIM_EXAMPLES_DIR "/does_not_exist.txt";
  CHECK_THROWS_WITH(parse_circuit_file(path),
                    "cannot open circuit file '" + path + "'");
}
