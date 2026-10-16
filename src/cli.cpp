#include "ecm_qsim/cli.hpp"

#include <charconv>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

#include "ecm_qsim/circuit.hpp"
#include "ecm_qsim/noise.hpp"
#include "ecm_qsim/parser.hpp"
#include "ecm_qsim/simulator.hpp"

namespace ecm_qsim {

namespace {

constexpr const char* kUsage =
    "usage: sim CIRCUIT --shots N --seed S [--p-gate P] [--p-idle P] "
    "[--p-readout P] [--json]";

// Bad command line arguments: exit status 2, message followed by the usage.
class UsageError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

struct Options {
  bool help = false;
  std::string path;
  int shots = 0;
  std::uint64_t seed = 0;
  NoiseModel noise;
  bool json = false;
};

// Parses all of token as an integer of type T with std::from_chars, which
// rejects a sign other than '-', whitespace, exponents and overflow.
template <typename T>
std::optional<T> parse_integer(const std::string& token) {
  T value{};
  const char* first = token.data();
  const char* last = first + token.size();
  const auto [end, ec] = std::from_chars(first, last, value);
  if (ec != std::errc{} || end != last) {
    return std::nullopt;
  }
  return value;
}

// Parses all of token as a double with std::strtod. Overflow gives inf,
// which validate_noise_model rejects.
std::optional<double> parse_double(const std::string& token) {
  if (token.empty()) {
    return std::nullopt;
  }
  char* end = nullptr;
  const double value = std::strtod(token.c_str(), &end);
  if (end != token.c_str() + token.size()) {
    return std::nullopt;
  }
  return value;
}

Options parse_args(const std::vector<std::string>& args) {
  Options options;
  std::optional<std::string> path;
  std::optional<int> shots;
  std::optional<std::uint64_t> seed;
  std::optional<double> p_gate;
  std::optional<double> p_idle;
  std::optional<double> p_readout;
  bool json = false;

  for (std::size_t i = 0; i < args.size(); ++i) {
    const std::string& arg = args[i];
    if (arg == "--help") {
      options.help = true;
      return options;
    }
    if (arg == "--json") {
      if (json) {
        throw UsageError(arg + " given twice");
      }
      json = true;
      continue;
    }
    if (!arg.starts_with("--")) {
      if (path) {
        throw UsageError("unexpected argument '" + arg + "'");
      }
      path = arg;
      continue;
    }
    if (arg != "--shots" && arg != "--seed" && arg != "--p-gate" &&
        arg != "--p-idle" && arg != "--p-readout") {
      throw UsageError("unknown option '" + arg + "'");
    }
    if (i + 1 == args.size()) {
      throw UsageError(arg + " needs a value");
    }
    const std::string& value = args[++i];
    if (arg == "--shots") {
      if (shots) {
        throw UsageError(arg + " given twice");
      }
      shots = parse_integer<int>(value);
      if (!shots || *shots < 1) {
        throw UsageError("--shots must be a positive integer, got '" + value +
                         "'");
      }
    } else if (arg == "--seed") {
      if (seed) {
        throw UsageError(arg + " given twice");
      }
      seed = parse_integer<std::uint64_t>(value);
      if (!seed) {
        throw UsageError("--seed must be a non-negative integer, got '" +
                         value + "'");
      }
    } else {
      std::optional<double>& p = arg == "--p-gate"   ? p_gate
                                 : arg == "--p-idle" ? p_idle
                                                     : p_readout;
      if (p) {
        throw UsageError(arg + " given twice");
      }
      p = parse_double(value);
      if (!p) {
        throw UsageError(arg + " must be a number, got '" + value + "'");
      }
    }
  }

  if (!path) {
    throw UsageError("missing circuit file");
  }
  if (!shots) {
    throw UsageError("missing --shots");
  }
  if (!seed) {
    throw UsageError("missing --seed");
  }
  options.path = *path;
  options.shots = *shots;
  options.seed = *seed;
  options.noise = NoiseModel{p_gate.value_or(0.0), p_idle.value_or(0.0),
                             p_readout.value_or(0.0)};
  options.json = json;
  return options;
}

void write_table(const Counts& counts, std::ostream& out) {
  for (const auto& [bits, count] : counts) {
    out << bits << ' ' << count << '\n';
  }
}

void write_json(const Counts& counts, std::ostream& out) {
  out << "{\"counts\": {";
  const char* separator = "";
  for (const auto& [bits, count] : counts) {
    out << separator << '"' << bits << "\": " << count;
    separator = ", ";
  }
  out << "}}\n";
}

}  // namespace

int sim_main(const std::vector<std::string>& args, std::ostream& out,
             std::ostream& err) {
  try {
    const Options options = parse_args(args);
    if (options.help) {
      out << kUsage << '\n';
      return 0;
    }
    try {
      validate_noise_model(options.noise);
    } catch (const std::invalid_argument& e) {
      throw UsageError(e.what());
    }
    const Counts counts = run(parse_circuit_file(options.path), options.noise,
                              options.shots, options.seed);
    if (options.json) {
      write_json(counts, out);
    } else {
      write_table(counts, out);
    }
    return 0;
  } catch (const UsageError& e) {
    err << "sim: " << e.what() << '\n' << kUsage << '\n';
    return 2;
  } catch (const std::exception& e) {
    err << "sim: " << e.what() << '\n';
    return 1;
  }
}

}  // namespace ecm_qsim
