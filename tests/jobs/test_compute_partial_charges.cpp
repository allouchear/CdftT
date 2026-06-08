#include "../helpers.hpp"

#include <cstdlib>
#include <iostream>
#include <string>
#include <regex>
#include <sstream>

using namespace cdftt_tests;

/*
    About this job type (ComputePartialCharges):
        - It computes atomic partial charges from an input scalar grid and an
            accepted partition method.
        - The smoke path in this file uses one fixture-backed density cube with
            PartitionMethod=On-Grid.
        - Some methods are forbidden for this workflow and must fail explicitly,
            including BBS and B2S.

    Detailed header command expectations:
        - Mandatory header line:
                    RunType=ComputePartialCharges
        - Required companion lines:
                    Grids=<densityCube>
                    PartitionMethod=<method>

    Expected behavior to anchor assertions:
        - Happy path:
                1. Exit code 0.
                2. stdout includes partition echo and charge-table markers.
        - Negative path:
                1. Non-zero exit code.
                2. stderr reports invalid partition method for this job.
*/

// [P0 | PR-fast] Fixture-backed happy path using one cube grid and an allowed
// partition method.
static void test_compute_partial_charges_smoke() {
    // As input files are required, this test is fixture-backed with a checked-in cube grid from
    // existing test assets. The partition method and output paths are explicitly defined in the
    // input content.
    const std::string root = repo_root();
    const std::string grid = root + "/examples/Density/Density.cube";
    const std::string input_path = "/tmp/cdftt_test_partial_charges_smoke.txt";

    write_input_file(input_path,
                     std::string("RunType=ComputePartialCharges\n")
                   + "Grids=" + grid + "\n"
                   + "PartitionMethod=On-Grid\n");

    // The run should succeed and print markers about the chosen partition method and the partial
    // charges table.
    const RunResult r = run_cdftt(input_path);

    assert_zero_exit(r);
    assert_stdout_contains(r, "Volume partition method: On-Grid");
    assert_stdout_contains(r, "Partial charges (Pos in Angstrom)");
    assert_stdout_contains(r, "Total charge =");

    // Clean up the temporary input file.
    std::remove(input_path.c_str());
}

// [P0 | PR-fast] Negative path: BBS and B2S are forbidden for ComputePartialCharges.
static void test_compute_partial_charges_invalid_partition_bbs() {
    // The run should fail and print an explicit error mentioning invalid partitionMethod.
    const std::string input_path = "/tmp/cdftt_test_partial_charges_invalid_bbs.txt";
    write_input_file(input_path,
                     "RunType=ComputePartialCharges\n"
                     "Grids=dummy_density.cube\n"
                     "PartitionMethod=BBS\n");

    const RunResult r = run_cdftt(input_path);

    assert_nonzero_exit(r);
    assert_stderr_contains(r, "partitionMethod \"BBS\" invalid for this job");

    // Clean up the temporary input file.
    std::remove(input_path.c_str());
}

// [P0 | PR-fast] Negative path: same rule for B2S.
static void test_compute_partial_charges_invalid_partition_b2s() {
    // The run should fail and print an explicit error mentioning invalid partitionMethod.
    const std::string input_path = "/tmp/cdftt_test_partial_charges_invalid_b2s.txt";
    write_input_file(input_path,
                     "RunType=ComputePartialCharges\n"
                     "Grids=dummy_density.cube\n"
                     "PartitionMethod=B2S\n");

    const RunResult r = run_cdftt(input_path);

    assert_nonzero_exit(r);
    assert_stderr_contains(r, "partitionMethod \"B2S\" invalid for this job");

    // Clean up the temporary input file.
    std::remove(input_path.c_str());
}

// [P1 | nightly] Becke method path with analytic input-backed workflow.
static void test_compute_partial_charges_becke_analytic() {
    // Uses a small deterministic fixture (H2O) so we can assert atom count and
    // total partial charge sanity.
    const std::string root = repo_root();
    const std::string analytic = root + "/tests/fixtures/Orbitals/h2o.fchk";
    const std::string g1 = root + "/examples/Density/Density.cube";
    const std::string input_path = "/tmp/cdftt_test_partial_charges_becke_analytic.txt";

    write_input_file(input_path,
                     std::string("RunType=ComputePartialCharges\n")
                   + "AnalyticFiles=" + analytic + "\n"
                   + "PartitionMethod=Becke\n"
                   + "Grids=" + g1 + "\n"
                   + "Size=Medium\n");

    const RunResult r = run_cdftt(input_path);
    assert_zero_exit(r);

    // Confirm partition method applied and partial-charge block present.
    assert_stdout_contains(r, "Volume partition method: Becke");

    // Count per-atom rows. H2O has 3 atoms; expect 3 per-atom charge lines.
    // TODO: "value" should be replaced with "charge" or "partial charge" in the output for clearer parsing and user understanding.
    const int atom_lines = count_lines_matching_regex(r.stdout_text, R"(^ *Atom = [A-Za-z]+ *, *value = *[-+]?\d*\.\d+)");
    if (atom_lines != 3)
        throw std::runtime_error("Expected 3 atom partial-charge lines, got " + std::to_string(atom_lines));

    // TODO: for now, the job does not return the Total charge in an explicit way.
    //? As a placeholder, we can compute it by summing the per-atom charges
    double total_charge = 0.0;
    std::istringstream iss(r.stdout_text);
    std::string line;
    while (std::getline(iss, line)) {
        std::smatch match;
        if (std::regex_search(line, match, std::regex(R"(^ *Atom = [A-Za-z]+ *, *value = *([-+]?\d*\.\d+))"))) {
            total_charge += std::stod(match[1].str());
        }
    }
    // For this exact exemple, the charge is expected to be -0.077409
    // TODO -> One could expect the charge to be exacly 0 for a neutral molecule -> find a better example
    double expected_total_charge = -0.077409;
    if (std::abs(total_charge - expected_total_charge) > 1e-6)
        throw std::runtime_error("Total partial charge differs from expected "
                                    + std::to_string(expected_total_charge)
                                    + ": " + std::to_string(total_charge));

    // Clean up the temporary input file.
    std::remove(input_path.c_str());
}

int main() {
    const Test tests[] = {
        { "test_compute_partial_charges_smoke",                 test_compute_partial_charges_smoke                 },
        { "test_compute_partial_charges_invalid_partition_bbs", test_compute_partial_charges_invalid_partition_bbs },
        { "test_compute_partial_charges_invalid_partition_b2s", test_compute_partial_charges_invalid_partition_b2s },
        { "test_compute_partial_charges_becke_analytic",       test_compute_partial_charges_becke_analytic       },
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
