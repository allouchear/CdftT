#include "../helpers.hpp"

#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
#include <iterator>

using namespace cdftt_tests;

/*
    About this job type (RunLambdaDiagnostic):
        - It computes the lambda diagnostic from orbital data and excited-state
            transition information.
        - This workflow requires both a valid analytic source and a transitions
            definition source (file-backed in this test file).
        - Missing transitions input should fail early with a clear file-open error.

        
    Detailed header command expectations:
        - Mandatory header line:
                    RunType=RunLambdaDiagnostic
        - Required companion lines for this smoke scenario:
                    AnalyticFiles=<inputOrbitalFile>
                    Size=<Small|Medium|Large|...>
                    TransitionsFile=<path>

    Expected behavior to anchor assertions:
        - Happy path:
                1. Exit code 0.
                2. stdout includes domain-build marker, state count, and lambda value.
        - Negative path:
                1. Non-zero exit code.
                2. stderr reports inability to open transitions file.
*/

// [P0 | PR-fast] Valid lambda diagnostic run with explicit transitions file.
static void test_lambda_diagnostic_smoke() {
    // As input files are required, this test is fixture-backed with checked-in analytic
    // and transitions files from existing test assets. The input content explicitly defines
    // the path to the transitions file, which is required for this workflow and should be
    // successfully read. The analytic file is also required but the test focuses on the
    // transitions input as the critical dependency for this workflow.
    const std::string root = repo_root();
    const std::string analytic = root + "/src/applications/cdftt/tests/CC-pvtz/CO_TDDFT.fchk";
    const std::string transitions = root + "/src/applications/cdftt/tests/transitionsFiles/transitions_CO_twoFirstStates.txt";
    const std::string input_path = "/tmp/cdftt_test_lambda_diagnostic_smoke.txt";

    write_input_file(input_path,
                     std::string("RunType=RunLambdaDiagnostic\n")
                   + "AnalyticFiles=" + analytic + "\n"
                   + "Size=Medium\n"
                   + "TransitionsFile=" + transitions + "\n");

    // The run should succeed and print markers about the current job, domain build, number of
    // excited states read, and the computed lambda value.
    const RunResult r = run_cdftt(input_path);

    assert_zero_exit(r);
    assert_stdout_contains(r, "Building domain and grid, please wait...");
    assert_stdout_contains(r, "Number of excited states read:");
    assert_stdout_contains(r, "Lambda =");

    // Clean up the temporary input file.
    std::remove(input_path.c_str());
}

// [P0 | PR-fast] Missing required transitions source should fail clearly.
static void test_lambda_diagnostic_missing_transitions() {
    // The run should fail and print an explicit error mentioning inability to open the transitions file.
    const std::string root = repo_root();
    const std::string analytic = root + "/src/applications/cdftt/tests/CC-pvtz/H2CO_TDDFT.fchk";
    const std::string input_path = "/tmp/cdftt_test_lambda_diagnostic_missing_transitions.txt";

    write_input_file(input_path,
                     std::string("RunType=RunLambdaDiagnostic\n")
                   + "AnalyticFiles=" + analytic + "\n"
                   + "Size=Medium\n");

    const RunResult r = run_cdftt(input_path);

    std::cout << "STDOUT:" << r.stdout_text << std::endl;

    assert_nonzero_exit(r);
    assert_stdout_contains(r, "Note: the \"TransitionsFile\" parameter is not specified in the provided input file (" + input_path + ").");

    // TODO -> The current behaviour is odd  it reports the file is not found, tries to read another one as default, segfaults, thus returns nonzero but the stderr is empty.
    assert_stderr_contains(r, "Error: could not open transitions file");

    // Clean up the temporary input file.
    std::remove(input_path.c_str());
}

int main() {
    const Test tests[] = {
        { "test_lambda_diagnostic_smoke",               test_lambda_diagnostic_smoke               },
        { "test_lambda_diagnostic_missing_transitions", test_lambda_diagnostic_missing_transitions },
    };

    int passed = 0;
    int failed = 0;
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

    std::cout << passed << " passed, " << failed << " failed.\n" << std::endl;
    return failed > 0 ? 1 : 0;
}
