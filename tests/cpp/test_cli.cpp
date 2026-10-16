#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "ecm_qsim/circuit.hpp"
#include "ecm_qsim/cli.hpp"
#include "ecm_qsim/noise.hpp"
#include "ecm_qsim/parser.hpp"
#include "ecm_qsim/simulator.hpp"
#include "ecm_qsim/state_vector.hpp"

using namespace ecm_qsim;
using Catch::Matchers::ContainsSubstring;

namespace {

struct CliResult {
  int status;
  std::string out;
  std::string err;
};

// Calls sim_main with std::ostringstream for out and err.
CliResult sim(const std::vector<std::string>& args) {
  std::ostringstream out;
  std::ostringstream err;
  const int status = sim_main(args, out, err);
  return CliResult{status, out.str(), err.str()};
}

// Writes text to temp_directory_path() / "ecm_qsim_test_cli" / name and
// returns the path (creates the directory).
std::filesystem::path write_circuit(const std::string& name,
                                    const std::string& text) {
  const std::filesystem::path dir =
      std::filesystem::temp_directory_path() / "ecm_qsim_test_cli";
  std::filesystem::create_directories(dir);
  const std::filesystem::path path = dir / name;
  std::ofstream(path) << text;
  return path;
}

// Parses "<bitstring> <count>\n" lines into Counts. REQUIREs that every line
// has this form and that out ends with '\n'.
Counts parse_table(const std::string& out) {
  REQUIRE(!out.empty());
  REQUIRE(out.back() == '\n');
  Counts counts;
  std::istringstream lines(out);
  std::string line;
  while (std::getline(lines, line)) {
    CAPTURE(line);
    const std::size_t space = line.find(' ');
    REQUIRE(space != std::string::npos);
    const std::string bits = line.substr(0, space);
    const std::string count = line.substr(space + 1);
    REQUIRE(!bits.empty());
    REQUIRE(bits.find_first_not_of("01") == std::string::npos);
    REQUIRE(!count.empty());
    REQUIRE(count.find_first_not_of("0123456789") == std::string::npos);
    REQUIRE(counts.count(bits) == 0);
    counts[bits] = std::stoi(count);
  }
  return counts;
}

const std::string kGhz = ECM_QSIM_EXAMPLES_DIR "/ghz3.txt";

}  // namespace

TEST_CASE("sim accepts the example command line") {
  const CliResult result =
      sim({kGhz, "--shots", "1000", "--p-gate", "0.001", "--p-idle", "0.001",
           "--p-readout", "0.01", "--seed", "1"});
  REQUIRE(result.status == 0);
  CHECK(result.err.empty());
  CHECK(parse_table(result.out) ==
        run(parse_circuit_file(kGhz), NoiseModel{0.001, 0.001, 0.01}, 1000, 1));
}

TEST_CASE("sim passes each noise flag to its own field") {
  const auto [flag, noise] = GENERATE(table<std::string, NoiseModel>({
      {"--p-gate", NoiseModel{0.2, 0, 0}},
      {"--p-idle", NoiseModel{0, 0.2, 0}},
      {"--p-readout", NoiseModel{0, 0, 0.2}},
  }));
  CAPTURE(flag);
  const CliResult result =
      sim({kGhz, "--shots", "2000", "--seed", "3", flag, "0.2"});
  REQUIRE(result.status == 0);
  CHECK(parse_table(result.out) ==
        run(parse_circuit_file(kGhz), noise, 2000, 3));
}

TEST_CASE("sim defaults the noise flags to 0") {
  const CliResult result = sim({kGhz, "--shots", "500", "--seed", "1"});
  REQUIRE(result.status == 0);
  CHECK(parse_table(result.out) ==
        run(parse_circuit_file(kGhz), NoiseModel{}, 500, 1));
}

TEST_CASE("sim prints one line per outcome, in bitstring order") {
  const auto x0 = write_circuit("x0.txt", "QUBITS 3\nX 0\n");
  CHECK(sim({x0.string(), "--shots", "5", "--seed", "1"}).out == "001 5\n");
  // Noiseless GHZ has exactly two outcomes; 000 must come first.
  const std::string out = sim({kGhz, "--shots", "100", "--seed", "1"}).out;
  CHECK(out.starts_with("000 "));
  CHECK(out.find("\n111 ") != std::string::npos);
}

TEST_CASE("sim --json prints the counts as one JSON object") {
  const auto x0 = write_circuit("x0.txt", "QUBITS 3\nX 0\n");
  CHECK(sim({x0.string(), "--shots", "5", "--seed", "1", "--json"}).out ==
        "{\"counts\": {\"001\": 5}}\n");
  const auto h0 = write_circuit("h0.txt", "QUBITS 1\nH 0\n");
  const Counts counts = run(parse_circuit_file(h0), NoiseModel{}, 100, 1);
  REQUIRE(counts.size() == 2);
  CHECK(sim({h0.string(), "--shots", "100", "--seed", "1", "--json"}).out ==
        "{\"counts\": {\"0\": " + std::to_string(counts.at("0")) +
            ", \"1\": " + std::to_string(counts.at("1")) + "}}\n");
}

TEST_CASE("sim accepts options before the path") {
  const auto x0 = write_circuit("x0.txt", "QUBITS 3\nX 0\n");
  CHECK(sim({"--seed", "1", "--shots", "5", x0.string()}).out == "001 5\n");
}

TEST_CASE("sim --help prints the usage") {
  const CliResult result = sim({"--help"});
  CHECK(result.status == 0);
  CHECK_THAT(result.out, ContainsSubstring("usage: sim"));
  CHECK(result.err.empty());
}

TEST_CASE("sim reports usage errors with status 2") {
  using Args = std::vector<std::string>;
  const auto [args, message] = GENERATE(table<Args, std::string>({
      {Args{}, "missing circuit file"},
      {Args{kGhz, "--seed", "1"}, "missing --shots"},
      {Args{kGhz, "--shots", "5"}, "missing --seed"},
      {Args{kGhz, "--shots", "0", "--seed", "1"},
       "--shots must be a positive integer, got '0'"},
      {Args{kGhz, "--shots", "abc", "--seed", "1"},
       "--shots must be a positive integer, got 'abc'"},
      {Args{kGhz, "--shots", "1e3", "--seed", "1"},
       "--shots must be a positive integer, got '1e3'"},
      {Args{kGhz, "--shots", "+5", "--seed", "1"},
       "--shots must be a positive integer, got '+5'"},
      {Args{kGhz, "--shots", "99999999999", "--seed", "1"},
       "--shots must be a positive integer, got '99999999999'"},
      {Args{kGhz, "--shots", "5", "--seed", "-1"},
       "--seed must be a non-negative integer, got '-1'"},
      {Args{kGhz, "--shots", "5", "--seed", "1", "--p-gate", "abc"},
       "--p-gate must be a number, got 'abc'"},
      {Args{kGhz, "--shots", "5", "--seed", "1", "--p-gate", "0.1x"},
       "--p-gate must be a number, got '0.1x'"},
      {Args{kGhz, "--shots", "5", "--seed", "1", "--p-gate", ""},
       "--p-gate must be a number, got ''"},
      {Args{kGhz, "--shots", "5", "--seed", "1", "--p-gate", "1.5"},
       "p_gate must be in [0, 1], got 1.5"},
      {Args{kGhz, "--shots", "5", "--seed", "1", "--p-idle", "nan"},
       "p_idle must be in [0, 1]"},
      {Args{kGhz, "--shots", "5", "--seed", "1", "--foo", "1"},
       "unknown option '--foo'"},
      {Args{kGhz, "--seed", "1", "--shots"}, "--shots needs a value"},
      {Args{kGhz, "--shots", "5", "--seed", "1", "--shots", "5"},
       "--shots given twice"},
      {Args{kGhz, kGhz, "--shots", "5", "--seed", "1"}, "unexpected argument"},
      // Arguments are checked before the file is opened.
      {Args{"no_such_file.txt", "--shots", "0", "--seed", "1"},
       "--shots must be a positive integer"},
  }));
  CAPTURE(args);
  const CliResult result = sim(args);
  CHECK(result.status == 2);
  CHECK(result.out.empty());
  CHECK(result.err.starts_with("sim: "));
  CHECK_THAT(result.err, ContainsSubstring(message));
  CHECK_THAT(result.err, ContainsSubstring("usage: sim"));
}

TEST_CASE("sim reports circuit errors with status 1") {
  std::string path;
  std::string message;
  SECTION("missing file") {
    path = "no_such_file.txt";
    message = "cannot open circuit file";
  }
  SECTION("malformed file") {
    path = write_circuit("bad.txt", "QUBITS 2\nX 5\n").string();
    message = "sim: line 2:";
  }
  SECTION("register too large for memory") {
    const int n = StateVector::max_qubits() + 1;
    if (n > kMaxQubits) {
      SKIP("this machine fits every register the IR can represent");
    }
    path =
        write_circuit("big.txt", "QUBITS " + std::to_string(n) + "\n").string();
    message = "sim: ";
  }
  const CliResult result = sim({path, "--shots", "1", "--seed", "1"});
  CHECK(result.status == 1);
  CHECK(result.out.empty());
  CHECK(result.err.starts_with("sim: "));
  CHECK_THAT(result.err, ContainsSubstring(message));
  CHECK_THAT(result.err, !ContainsSubstring("usage:"));
}
