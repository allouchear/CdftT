#ifndef CDFTT_TESTS_HELPERS_HPP
#define CDFTT_TESTS_HELPERS_HPP

#include <string>
#include <vector>

namespace cdftt_tests {

// Gathers the outcome of a single cdftt process invocation.
//  - exit_code: the process exit status (0 = success, non-zero = failure).
//  - stdout_text: full standard output captured from the run.
//  - stderr_text: full standard error captured from the run.
struct RunResult {
    int exit_code = -1;
    std::string stdout_text;
    std::string stderr_text;
};

struct Test {
    const char* name;
    void (*fn)();
};

// Invoke the cdftt executable with input_file_path as its input file.
// Captures exit code, stdout, and stderr into the returned RunResult.
// Must not throw on non-zero exit; callers use assert_zero_exit / assert_nonzero_exit.
RunResult run_cdftt(const std::string& input_file_path);

// Write content to the file at path, creating it if absent.
// Used to produce temporary input files used by run_cdftt.
// Throws on I/O failure (ie. for open or write issues).
void write_input_file(const std::string& path, const std::string& content);

// Assert that a regular file exists at path.
// Throws error if the file is absent; used by output-existence checks.
void assert_file_exists(const std::string& path);

// Assert that a run succeeded (exit code is zero).
void assert_zero_exit(const RunResult& result);

// Assert that a run failed (exit code is non-zero).
void assert_nonzero_exit(const RunResult& result);

// Assert that stdout contains a text fragment.
void assert_stdout_contains(const RunResult& result, const std::string& fragment);

// Assert that stdout does not contain a text fragment.
void assert_stdout_does_not_contain(const RunResult& result, const std::string& fragment);

// Assert that stderr contains a text fragment.
void assert_stderr_contains(const RunResult& result, const std::string& fragment);

// Assert that a cube file has a parseable header (core geometry/grid lines).
void assert_cube_header_parseable(const std::string& path);

// Assert that all parsed voxel values in a cube file are finite (no NaN/Inf).
void assert_cube_has_finite_values(const std::string& path);

// Asserts a given file's size is within expected bounds for a cube of given dimensions and data type.
void assert_file_size_less_than(const std::string& path, std::size_t max_bytes);

// Count how many lines in `text` match the given regular expression pattern.
// Returns the number of matching lines.
int count_lines_matching_regex(const std::string& text, const std::string& regex_pattern);

// Extract the first floating-point number that appears after the given label
// in `text`. Throws if the label is not found or a number cannot be parsed.
double extract_double_after_label(const std::string& text, const std::string& label);

// Read the first N numeric voxel values from a cube file (after header).
// Returns a vector with up to N values (fewer if file shorter). Throws on parse errors.
std::vector<double> read_first_n_voxels(const std::string& path, std::size_t n);

// Get the repository root directory, either from the CDFTT_REPO_ROOT environment variable or defaulting to ".".
std::string repo_root();

} // namespace cdftt_tests

#endif // CDFTT_TESTS_HELPERS_HPP
