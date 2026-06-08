#include "../helpers.hpp"

#include <iostream>
#include <stdexcept>
#include <string>

using namespace cdftt_tests;

/*
  About this job type:
    - It is a special case that only prints help information and does not require
      any other input parameters.
    - It is the only job type that should be runnable with just "RunType=Help" in
      the input file, without any additional parameters.
    - It is not expected to perform any complex computations, so a "smoke test"
      approach is sufficient for the happy path.
    - The main negative path to test is providing an unrecognized RunType value,
      which should be handled gracefully with a clear error message.

    Detailed header command expectations:
    - For the happy path, the input file should contain exactly one line:
          RunType=Help
      This is the minimal valid input that should trigger the help output.
    - For the negative path, we can use a similar input file but with an invalid RunType, e.g.:
          RunType=NotARealJob
      This should cause the program to exit with an error and print a message about the unknown job type.

    Details about the expected output:
    - For the happy path, we expect:
        1. The program to exit with code 0 (success).
        2. The standard output to contain a header line that lists available jobs, e.g. "Available jobs (runType=) :".
        3. The standard output to include the names of all public jobs defined in Job::_jobsList, which are:
           "Help",
           "computeDescriptors",
           "computePartialCharges",
           "computeIntegrals",
           "computeGridDifference",
           "MakeDensityCube",
           "MakeOrbitalsCube",
           "MakeELFCube",
           "ConvertOrbitals".
    - For the negative path, we expect:
        1. The program to exit with a non-zero code (indicating failure).
        2. The standard error to contain a message in the format:
           "Error: Run type \"<value>\" unknown."
           where <value> is the invalid RunType provided in the input file.
*/

/*************************************************************************
    Local assertion helpers
    These are thin wrappers that produce clear failure messages.
    They throw on failure so the test runner can catch and report them.
**************************************************************************/

// Check that the process exited with the exact expected code.
static void assert_exit_code(int actual, int expected) {
    if (actual != expected)
        throw std::runtime_error("Expected exit code " + std::to_string(expected)
                                 + " but got " + std::to_string(actual));
}


/*************************************************************************
    Actual tests
**************************************************************************/

// _____________________ TEST 1 _____________________
// Mandatory [P0], fast for pull requests [PR-fast] 
//
// What we test:
//   Running cdftt with a one-line input file "RunType=Help" should:
//     1. Exit successfully (exit code 0).
//     2. Print a header line listing available jobs.
//     3. Include every public job name in that output.
//
// Why these job names specifically:
//   They are the exact strings in Job::_jobsList (src/lib/JobControl/Job.cpp),
//   which Job::printListOfRunTypes() iterates and prints.
//   Casing matters - the source uses mixed case (e.g. "computeDescriptors").
// -----------------------------------------------------------------------------
static void test_help_smoke() {
    // Write a minimal input file that asks for the Help job.
    const std::string input_path = "/tmp/cdftt_test_help_smoke.txt";
    write_input_file(input_path, "RunType=Help\n");

    const RunResult r = run_cdftt(input_path);

    // The program should exit cleanly.
    assert_zero_exit(r);

    // The header line is always printed first.
    assert_stdout_contains(r, "Available jobs (runType=) :");

    // All 9 public job names must appear (hidden jobs are intentionally excluded here).
    assert_stdout_contains(r, "Help");
    assert_stdout_contains(r, "computeDescriptors");
    assert_stdout_contains(r, "computePartialCharges");
    assert_stdout_contains(r, "computeIntegrals");
    assert_stdout_contains(r, "computeGridDifference");
    assert_stdout_contains(r, "MakeDensityCube");
    assert_stdout_contains(r, "MakeOrbitalsCube");
    assert_stdout_contains(r, "MakeELFCube");
    assert_stdout_contains(r, "ConvertOrbitals");

    // Clean up the temporary input file.
    std::remove(input_path.c_str());
}

// -----------------------------------------------------------------------------
// [P0 | PR-fast] Negative path: unrecognised RunType value
//
// What we test:
//   Providing a RunType that cdftt does not know should:
//     1. Exit with code 1 (failure).
//     2. Print an explicit, user-readable error message to stderr.
//
// Why stderr and not stdout:
//   cdftt routes errors through print_error() (src/lib/Utils/Utils.cpp),
//   which writes to std::cerr. The exact message format is:
//       Error: Run type "<value>" unknown.
//   This is assembled in Job::readRunType() (src/lib/JobControl/Job.cpp).
// -----------------------------------------------------------------------------
static void test_help_invalid_runtype() {
    // Write an input file with a nonsense RunType value.
    const std::string input_path = "/tmp/cdftt_test_help_invalid.txt";
    write_input_file(input_path, "RunType=NotARealJob\n");

    const RunResult r = run_cdftt(input_path);

    // The program must signal failure to the caller (e.g. the CI system).
    assert_exit_code(r.exit_code, 1);

    // The error message must name the bad value so the user knows what to fix.
    assert_stderr_contains(r, "Error: Run type \"NotARealJob\" unknown.");

    // Clean up the temporary input file.
    std::remove(input_path.c_str());
}

// =============================================================================
// Minimal test runner
//
// No external framework is used — each test is a plain function.
// The runner calls each in turn, catches any exception as a failure,
// and returns exit code 1 if any test failed (CI-friendly).
// =============================================================================

int main() {
    const Test tests[] = {
        { "test_help_smoke",           test_help_smoke           },
        { "test_help_invalid_runtype", test_help_invalid_runtype },
    };

    int passed = 0, failed = 0;
    for (const auto& t : tests) {
        std::cout << "[ RUN ] " << t.name << std::endl;
        try {
            t.fn();
            std::cout << "[ OK  ] " << t.name << std::endl;
            ++passed;
        } catch (const std::exception& e) {
            std::cerr << "[FAIL ] " << t.name << ": " << e.what() << std::endl;
            ++failed;
        }
    }
    std::cout << "\n" << passed << " passed, " << failed << " failed." << std::endl;
    return failed > 0 ? 1 : 0;
}
