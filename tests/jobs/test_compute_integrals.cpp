#include "../helpers.hpp"

#include <cstdlib>
#include <iostream>
#include <string>

using namespace cdftt_tests;

/*
    About this job type (ComputeIntegrals):
        - It integrates one or more scalar fields over basins defined from a primary
            grid and a selected partition method.
        - The smoke path here uses three fixture-backed grids and an allowed
            partition method (On-Grid).
        - Several methods are intentionally forbidden for this job and should fail
            with clear diagnostics (Becke, FD, FMO in this test file).

    Detailed header command expectations:
        - Mandatory header line:
                    RunType=ComputeIntegrals
        - Required companion line:
                    GridFilesNames=<grid1>,<grid2>,<grid3>
            where the first grid is used to build basins.
        - Required partition line:
                    PartitionMethod=<method>

    Expected behavior to anchor assertions:
        - Happy path:
                1. Exit code 0.
                2. stdout contains partition echo and basin/integration progress markers.
        - Negative path:
                1. Non-zero exit code.
                2. stderr reports invalid partition method for this specific job.
*/

static std::string make_invalid_method_input(const std::string& method) {
    return std::string("RunType=ComputeIntegrals\n")
         + "GridFilesNames=grid1.cube,grid2.cube,grid3.cube\n"
         + "PartitionMethod=" + method + "\n";
}

// [P0 | PR-fast] Fixture-backed happy path using a basin-defining grid plus
// additional grids to integrate.
static void test_compute_integrals_smoke() {
    // As input files are required, this test is fixture-backed with checked-in cube
    // grids from existing test assets. The partition method and output paths are
    // explicitly defined in the input content.
    const std::string root = repo_root();
    const std::string g1 = root + "/examples/Density/Density.cube";
    const std::string g2 = root + "/examples/ELF/ELF.cube";
    const std::string g3 = root + "/examples/Density/Density.cube";
    const std::string input_path = "/tmp/cdftt_test_integrals_smoke.txt";

    write_input_file(input_path,
                     std::string("RunType=ComputeIntegrals\n")
                   + "GridFilesNames=" + g1 + "," + g2 + "," + g3 + "\n"
                   + "PartitionMethod=On-Grid\n");

    // The run should succeed and print markers about the chosen partition method, the
    // number of critical points found, and the progress of building basins and performing
    // integrations.
    const RunResult r = run_cdftt(input_path);

    assert_zero_exit(r);
    assert_stdout_contains(r, "Volume partition method: On-Grid");
    assert_stdout_contains(r, "to build basins.");
    assert_stdout_contains(r, "Number of critical points =");

    // Clean up the temporary input file.
    std::remove(input_path.c_str());
}

// [P0 | PR-fast] BECKE is forbidden for ComputeIntegrals.
static void test_compute_integrals_invalid_method_becke() {
    // The run should fail and print an explicit error mentioning invalid partitionMethod.
    const std::string input_path = "/tmp/cdftt_test_integrals_invalid_becke.txt";
    write_input_file(input_path, make_invalid_method_input("BECKE"));

    const RunResult r = run_cdftt(input_path);
    assert_nonzero_exit(r);
    assert_stderr_contains(r, "partitionMethod \"Becke\" invalid");

    // Clean up the temporary input file.
    std::remove(input_path.c_str());
}

// [P0 | PR-fast] FD is forbidden for ComputeIntegrals.
static void test_compute_integrals_invalid_method_fd() {
    // The run should fail and print an explicit error mentioning invalid partitionMethod.
    const std::string input_path = "/tmp/cdftt_test_integrals_invalid_fd.txt";
    write_input_file(input_path, make_invalid_method_input("FD"));

    const RunResult r = run_cdftt(input_path);
    assert_nonzero_exit(r);
    assert_stderr_contains(r, "partitionMethod \"FD\" invalid");

    // Clean up the temporary input file.
    std::remove(input_path.c_str());
}

// [P0 | PR-fast] FMO is forbidden for ComputeIntegrals.
static void test_compute_integrals_invalid_method_fmo() {
    // The run should fail and print an explicit error mentioning invalid partitionMethod.
    const std::string input_path = "/tmp/cdftt_test_integrals_invalid_fmo.txt";
    write_input_file(input_path, make_invalid_method_input("FMO"));

    const RunResult r = run_cdftt(input_path);
    assert_nonzero_exit(r);
    assert_stderr_contains(r, "partitionMethod \"FMO\" invalid");

    // Clean up the temporary input file.
    std::remove(input_path.c_str());
}

// [P1 | nightly] BBS partition path with cutoff parameter should succeed.
static void test_compute_integrals_bbs_cutoff() {
        // The run should succeed and print markers about the chosen partition method, building basins,
        // and the cutoff value used.
        // As input files are required, this test is fixture-backed with checked-in cube grids from existing test assets.
        // The partition method and cutoff parameter are explicitly defined in the input content.
        const std::string root = repo_root();
        const std::string g1 = root + "/examples/Density/Density.cube";
        const std::string g2 = root + "/examples/ELF/ELF.cube";
        const std::string g3 = root + "/examples/Density/Density.cube";
        const std::string input_path = "/tmp/cdftt_test_integrals_bbs_cutoff.txt";

        write_input_file(input_path,
                         std::string("RunType=ComputeIntegrals\n")
                         + "GridFilesNames=" + g1 + "," + g2 + "," + g3 + "\n"
                         + "PartitionMethod=BBS\n"
                         + "Cutoff=100.0\n");

        const RunResult r = run_cdftt(input_path);

        // Expect the BBS partition to be accepted and basins to be built.
        assert_zero_exit(r);
        assert_stdout_contains(r, "Volume partition method: BBS");
        assert_stdout_contains(r, "to build basins.");
        // assert_stdout_contains(r, "Cutoff"); //! Not reported in stdout yet
        // TODO: add the reporting of the cutoff value in the output for clearer diagnostics and test anchoring.

        std::remove(input_path.c_str());
}

// [P1 | nightly] B2S partition path with optional cutoff should succeed.
static void test_compute_integrals_b2s_cutoff() {
    // As input files are required, this test is fixture-backed with checked-in cube grids from existing test assets.
    // The partition method and cutoff parameter are explicitly defined in the input content.
    const std::string root = repo_root();
    const std::string g1 = root + "/examples/Density/Density.cube";
    const std::string g2 = root + "/examples/ELF/ELF.cube";
    const std::string g3 = root + "/examples/Density/Density.cube";
    const std::string input_path = "/tmp/cdftt_test_integrals_b2s_cutoff.txt";

    write_input_file(input_path,
                     std::string("RunType=ComputeIntegrals\n")
                   + "GridFilesNames=" + g1 + "," + g2 + "," + g3 + "\n"
                   + "PartitionMethod=B2S\n"
                   + "Cutoff=0.05\n");

    const RunResult r = run_cdftt(input_path);

    assert_zero_exit(r);
    assert_stdout_contains(r, "Volume partition method: B2S");
    assert_stdout_contains(r, "to build basins.");
    // assert_stdout_contains(r, "Cutoff"); //! Not reported in stdout yet
    // TODO: add the reporting of the cutoff value in the output for clearer diagnostics and test anchoring.

    // ? Note -> it seems the cutoff value, despite being accepted as input, and raising a warning
    //?          with default value when not provided, does not change the output in a noticeable way
    //?          for the tested case. This may be worth investigating further to confirm the cutoff is
    //?          properly applied in the B2S partitioning logic.

    std::remove(input_path.c_str());
}

int main() {
    const Test tests[] = {
        { "test_compute_integrals_smoke",                test_compute_integrals_smoke                },
        { "test_compute_integrals_invalid_method_becke", test_compute_integrals_invalid_method_becke },
        { "test_compute_integrals_invalid_method_fd",     test_compute_integrals_invalid_method_fd     },
        { "test_compute_integrals_invalid_method_fmo",    test_compute_integrals_invalid_method_fmo    },
        { "test_compute_integrals_bbs_cutoff",            test_compute_integrals_bbs_cutoff            },
        { "test_compute_integrals_b2s_cutoff",            test_compute_integrals_b2s_cutoff            },
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
