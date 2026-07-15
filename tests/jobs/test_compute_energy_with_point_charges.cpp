#include "../helpers.hpp"

#include <cstdlib>
#include <iostream>
#include <string>

using namespace cdftt_tests;

/*
    About this job type (ComputeEnergyWithPointCharges):
        - It computes excited-state/energy-related quantities in presence of one or
            more external point charges.
        - Inputs combine an analytic orbital source, charge values, their 3D
            positions, and optionally explicit transitions data.
        - Mismatch between number of charges and number of position triplets must be
            treated as a hard input error.

    Detailed header command expectations:
        - Mandatory header line:
                    RunType=ComputeEnergyWithPointCharges
        - Required companion lines used in this file:
                    Charges=<q1>[,<q2>,...]
                    AnalyticFiles=<inputOrbitalFile>
                    Positions=(<x1>,<y1>,<z1>)[,(<x2>,<y2>,<z2>),...]
        - Optional but exercised in smoke path:
                    TransitionsFile=<path>

    Expected behavior to anchor assertions:
        - Happy path:
                1. Exit code 0.
                2. stdout includes charge count and transitions/eigenvalue markers.
        - Negative path:
                1. Non-zero exit code.
                2. stderr reports incorrect number of point-charge positions.
*/

// [P0 | PR-fast] Minimal valid run with one analytic file, one point charge,
// explicit position and transitions source.
static void test_compute_energy_with_point_charges_smoke() {
    // As input files are required, this test is fixture-backed with checked-in analytic and
    // transitions files from existing test assets. The point charge and position are explicitly
    // defined in the input content.
    const std::string root = repo_root();
    const std::string analytic = root + "/src/applications/cdftt/tests/Guegan2020/acrolein.fchk";
    const std::string transitions = root + "/src/applications/cdftt/tests/transitionsFiles/transitions_acrolein.txt";
    const std::string input_path = "/tmp/cdftt_test_energy_point_charges_smoke.txt";

    // Keep Charges before RunType to match parser behavior used by existing project inputs.
    write_input_file(input_path,
                     std::string("Charges=-0.1\n")
                   + "RunType=ComputeEnergyWithPointCharges\n"
                   + "AnalyticFiles=" + analytic + "\n"
                   + "TransitionsFile=" + transitions + "\n"
                   + "Positions=(-1.752102,-0.142545,-0.000106)\n");

    // The run should succeed and print the number of point charges, read transitions, and computed eigenvalues.
    const RunResult r = run_cdftt(input_path);

    assert_zero_exit(r);
    assert_stdout_contains(r, "Number of runs: 1");
    assert_stdout_contains(r, "Reading transitions from file:");
    assert_stdout_contains(r, "Total number of states:");
    assert_stdout_contains(r, "Sorted Eigenvalues:");

    // Clean up the temporary input file.
    std::remove(input_path.c_str());
}

// [P0 | PR-fast] Charges/positions mismatch should fail clearly.
static void test_compute_energy_with_point_charges_positions_mismatch() {
    // The run should fail and print an explicit error mentioning incorrect number of point-charge positions.
    // The analytic file is required for the parser to reach the point of validating charges vs. positions count
    // so we use a fixture-backed path with a checked-in analytic file from existing test assets. The input content
    // deliberately includes one charge but only one position triplet, which should trigger the error.
    const std::string root = repo_root();
    const std::string analytic = root + "/src/applications/cdftt/tests/Guegan2020/acrolein.fchk";
    const std::string transitions = root + "/src/applications/cdftt/tests/transitionsFiles/transitions_acrolein.txt";
    const std::string input_path = "/tmp/cdftt_test_energy_point_charges_positions_mismatch.txt";

    write_input_file(input_path,
                     std::string("Charges=-0.1,-0.2\n")
                   + "RunType=ComputeEnergyWithPointCharges\n"
                   + "AnalyticFiles=" + analytic + "\n"
                   + "TransitionsFile=" + transitions + "\n"
                   + "ChargesPositionsBijections=True\n"
                   + "Positions=(-1.752102,-0.142545,-0.000106),(0.0,0.0,0.0),(1.0,1.0,1.0)\n"
                   + "Becke=default\n");

    // The run should fail and print an explicit error mentioning incorrect number of point-charge positions.
    const RunResult r = run_cdftt(input_path);

    assert_nonzero_exit(r);
    assert_stderr_contains(r, "the number of positions is not a multiple of the number of charges in the input file");

    // Clean up the temporary input file.
    std::remove(input_path.c_str());
}

// [P1 | nightly] Test Loop-on-atoms mode (no Positions provided)
static void test_compute_energy_with_point_charges_loop_on_atoms() {
    // This test exercises the loop-on-atoms mode where positions are derived from atomic coordinates
    // rather than explicitly provided. This requires an analytic file that supports atomic coordinate extraction.
    const std::string root = repo_root();
    const std::string analytic = root + "/src/applications/cdftt/tests/Guegan2020/acrolein.fchk";
    const std::string transitions = root + "/src/applications/cdftt/tests/transitionsFiles/transitions_acrolein.txt";
    const std::string input_path = "/tmp/cdftt_test_energy_point_charges_loop_on_atoms.txt";

    // Test the loop-on-atoms mode by omitting the Positions line
    // This should trigger automatic atom-based positioning
    write_input_file(input_path,
                     std::string("Charges=-0.1\n")
                   + "RunType=ComputeEnergyWithPointCharges\n"
                   + "AnalyticFiles=" + analytic + "\n"
                   + "TransitionsFile=" + transitions + "\n"
                   // Note: No Positions line - this triggers loop-on-atoms mode
                   + "Becke=default\n");

    // The run should succeed in loop-on-atoms mode and show appropriate messaging
    const RunResult r = run_cdftt(input_path);

    assert_zero_exit(r);
    assert_stdout_contains(r, "Number of runs: 8"); // Assuming 8 atoms in the analytic file for this test
    // In loop-on-atoms mode, we might see different output indicating automatic positioning
    assert_stdout_contains(r, "Note: the \"Positions\" parameter is not specified in the provided input file");
    assert_stdout_contains(r, "The program will use atom positions.");
    assert_stdout_contains(r, "Total number of states:");
    assert_stdout_contains(r, "Sorted Eigenvalues:");

    // Clean up the temporary input file
    std::remove(input_path.c_str());
}

// [P1 | nightly] Test Multiple charges and multiple states case
static void test_compute_energy_with_point_charges_multiple_charges_states() {
    // This test exercises a more complex scenario with multiple charges and multiple states
    // using existing test fixtures that support this configuration
    const std::string root = repo_root();
    const std::string analytic = root + "/src/applications/cdftt/tests/Guegan2020/acrolein.fchk";
    const std::string transitions = root + "/src/applications/cdftt/tests/transitionsFiles/transitions_acrolein.txt";
    const std::string input_path = "/tmp/cdftt_test_energy_point_charges_multiple_charges.txt";

    // Test with multiple charges and multiple positions (should work with same transitions file)
    write_input_file(input_path,
                     std::string("Charges=-0.1,-0.2,+0.1\n")
                   + "RunType=ComputeEnergyWithPointCharges\n"
                   + "AnalyticFiles=" + analytic + "\n"
                   + "TransitionsFile=" + transitions + "\n"
                   + "Positions=(-1.752102,-0.142545,-0.000106),(0.0,0.0,0.0),(1.0,1.0,1.0)\n"
                   + "Becke=default\n");

    // The run should succeed with multiple charges and states
    const RunResult r = run_cdftt(input_path);

    assert_zero_exit(r);
    assert_stdout_contains(r, "Number of runs: 9"); // 3 charges x 3 positions = 9 runs
    assert_stdout_contains(r, "Reading transitions from file:");
    assert_stdout_contains(r, "Total number of states:");
    assert_stdout_contains(r, "Sorted Eigenvalues:");
    // Should also show information about multiple charge positions
    assert_stdout_contains(r, "Point charge of -0.1 e at position");
    assert_stdout_contains(r, "Point charge of -0.2 e at position");
    assert_stdout_contains(r, "Point charge of 0.1 e at position");
    //! Note -> the input file provides distances in Bohr radii, but the produced output
    //!         currently shows them in Angstroms. While confusing for users, this is not
    //!         a parsing failure - we don't assert the exact position here for that reason.
    // Clean up the temporary input file
    std::remove(input_path.c_str());
}

// [P1 | nightly] Test Missing transitions source behavior
static void test_compute_energy_with_point_charges_missing_transitions() {
    // This test verifies the behavior when transitions file is missing or not derivable
    // It should fail with a clear error message about missing transitions
    const std::string root = repo_root();
    const std::string analytic = root + "/src/applications/cdftt/tests/Guegan2020/acrolein.fchk";
    const std::string non_existent_transitions = "/tmp/non_existent_transitions_file.txt";
    const std::string input_path = "/tmp/cdftt_test_energy_point_charges_missing_transitions.txt";

    // Test with a non-existent transitions file
    write_input_file(input_path,
                     std::string("Charges=-0.1\n")
                   + "RunType=ComputeEnergyWithPointCharges\n"
                   + "AnalyticFiles=" + analytic + "\n"
                   + "TransitionsFile=" + non_existent_transitions + "\n"
                   + "Positions=(-1.752102,-0.142545,-0.000106)\n"
                   + "Becke=default\n");

    // The run should fail because the transitions file doesn't exist
    const RunResult r = run_cdftt(input_path);

    assert_nonzero_exit(r);
    // Should report file not found or unable to read transitions
    assert_stderr_contains(r, "could not open transitions file " + non_existent_transitions);

    // Clean up the temporary input file
    std::remove(input_path.c_str());
}

int main() {
    const Test tests[] = {
        { "test_compute_energy_with_point_charges_smoke",              test_compute_energy_with_point_charges_smoke              },
        { "test_compute_energy_with_point_charges_positions_mismatch", test_compute_energy_with_point_charges_positions_mismatch },
        { "test_compute_energy_with_point_charges_loop_on_atoms",      test_compute_energy_with_point_charges_loop_on_atoms      },
        { "test_compute_energy_with_point_charges_multiple_charges_states", test_compute_energy_with_point_charges_multiple_charges_states },
        { "test_compute_energy_with_point_charges_missing_transitions", test_compute_energy_with_point_charges_missing_transitions },
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
