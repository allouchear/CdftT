#include "../helpers.hpp"

#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
#include <iterator>

using namespace cdftt_tests;

/*
    About this job type (LambdaDiagnostic):
        - It computes the lambda diagnostic from orbital data and excited-state
            transition information.
        - This workflow requires both a valid analytic source and a transitions
            definition source (file-backed in this test file).
        - Missing transitions input should fail early with a clear file-open error.

        
    Detailed header command expectations:
        - Mandatory header line:
                    RunType=LambdaDiagnostic
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
                     std::string("RunType=LambdaDiagnostic\n")
                   + "AnalyticFiles=" + analytic + "\n"
                   + "Size=Medium\n"
                   + "TransitionsFile=" + transitions + "\n");

    // The run should succeed and print markers about the current job, domain build, number of
    // excited states read, and the computed lambda value.
    const RunResult r = run_cdftt(input_path);

    assert_zero_exit(r);
    assert_stdout_contains(r, "Current job: LambdaDiagnostic");
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
    const std::string analytic = root + "/src/applications/cdftt/tests/CC-pvtz/CO_TDDFT.fchk";
    const std::string input_path = "/tmp/cdftt_test_lambda_diagnostic_missing_transitions.txt";

    write_input_file(input_path,
                     std::string("RunType=LambdaDiagnostic\n")
                   + "AnalyticFiles=" + analytic + "\n"
                   + "Size=Medium\n");

    const RunResult r = run_cdftt(input_path);

    assert_nonzero_exit(r);
    assert_stderr_contains(r, "could not open transitions file");

    // Clean up the temporary input file.
    std::remove(input_path.c_str());
}

// [P1 | nightly] Test transitions-from-file path vs transitions-from-analytic path
static void test_lambda_diagnostic_transitions_comparison() {
    // This test compares the two possible transitions sources:
    // 1. Explicit file-based transitions (already tested in smoke)
    // 2. Analytic file-derived transitions (when supported by the analytic format)
    
    const std::string root = repo_root();
    const std::string analytic = root + "/src/applications/cdftt/tests/CC-pvtz/CO_TDDFT.fchk";
    const std::string transitions_file = root + "/src/applications/cdftt/tests/transitionsFiles/transitions_CO_twoFirstStates.txt";
    const std::string input_path_file = "/tmp/cdftt_test_lambda_diagnostic_transitions_file.txt";
    const std::string input_path_analytic = "/tmp/cdftt_test_lambda_diagnostic_transitions_analytic.txt";

    // Test 1: transitions from file (baseline - should work)
    write_input_file(input_path_file,
                     std::string("RunType=LambdaDiagnostic\n")
                   + "AnalyticFiles=" + analytic + "\n"
                   + "Size=Medium\n"
                   + "TransitionsFile=" + transitions_file + "\n");

    const RunResult r_file = run_cdftt(input_path_file);

    assert_zero_exit(r_file);
    assert_stdout_contains(r_file, "Current job: LambdaDiagnostic");
    assert_stdout_contains(r_file, "Building domain and grid, please wait...");
    assert_stdout_contains(r_file, "Number of excited states read:");
    assert_stdout_contains(r_file, "Lambda =");

    // Test 2: Try transitions from analytic (if supported by this format)
    // Note: This may not work for all analytic formats, so we'll check if it succeeds
    write_input_file(input_path_analytic,
                     std::string("RunType=LambdaDiagnostic\n")
                   + "AnalyticFiles=" + analytic + "\n"
                   + "Size=Medium\n");
                   // No TransitionsFile specified - should use analytic-derived transitions if supported

    const RunResult r_analytic = run_cdftt(input_path_analytic);

    // TODO: analytic derivation should in principle be available but is not achieved so far.

    // Either both should succeed or both should fail gracefully
    // The important thing is that the system handles both approaches
    if (r_analytic.exit_code == 0) {
        // If analytic-derived transitions work, verify the basic output
        assert_stdout_contains(r_analytic, "Current job: LambdaDiagnostic");
        assert_stdout_contains(r_analytic, "Building domain and grid, please wait...");
        assert_stdout_contains(r_analytic, "Number of excited states read:");
        assert_stdout_contains(r_analytic, "Lambda =");
        // Print captured outputs for easier diagnosis when tests run

        // Try to extract numeric Lambda values and compare within tolerance
        try {
            const double lambda_file = extract_double_after_label(r_file.stdout_text, "Lambda =");
            const double lambda_analytic = extract_double_after_label(r_analytic.stdout_text, "Lambda =");
            const double diff = std::abs(lambda_file - lambda_analytic);
            if (diff > 1e-6) {
                throw std::runtime_error("Lambda values differ between file-derived and analytic-derived transitions: diff=" + std::to_string(diff));
            }
        } catch (const std::exception &e) {
            std::cerr << "Failed to compare Lambda values: " << e.what() << std::endl;
            throw;
        }
    } else {
        // Expect the current behavior: program reports it "could not open transitions file"
        // with an empty filename when analytic-derived transitions are not available.
        assert_stderr_contains(r_analytic, "could not open transitions file");

        // Quick probe: scan the analytic file to see if it contains any transition-like markers.
        std::ifstream af(analytic);
        if (!af.is_open()) {
            std::cout << "Cannot open analytic file for inspection: " << analytic << std::endl;
        } else {
            std::string content((std::istreambuf_iterator<char>(af)), std::istreambuf_iterator<char>());
            const bool has_transition_word = content.find("Transition") != std::string::npos ||
                                       content.find("TRANSITION") != std::string::npos ||
                                       content.find("Excited") != std::string::npos ||
                                       content.find("Alpha Orbital Energies") != std::string::npos;
            if (has_transition_word) {
                // Print baseline and analytic captured outputs to help debugging
                // std::cout << "--- baseline (file) stdout ---\n" << r_file.stdout_text << "\n--- baseline stderr ---\n" << r_file.stderr_text << std::endl;
                std::cout << "--- analytic-derived stdout ---\n" << r_analytic.stdout_text << "\n--- analytic-derived stderr ---\n" << r_analytic.stderr_text << std::endl;

                throw std::runtime_error("Analytic file " + analytic + " appears to contain transitions but analytic-derived transitions path failed");
            }
        }
    }

    // Clean up temporary files
    std::remove(input_path_file.c_str());
    std::remove(input_path_analytic.c_str());
}

int main() {
    const Test tests[] = {
        { "test_lambda_diagnostic_smoke",               test_lambda_diagnostic_smoke               },
        { "test_lambda_diagnostic_missing_transitions", test_lambda_diagnostic_missing_transitions },
        { "test_lambda_diagnostic_transitions_comparison", test_lambda_diagnostic_transitions_comparison },
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

    std::cout << "\n" << passed << " passed, " << failed << " failed." << std::endl;
    return failed > 0 ? 1 : 0;
}
