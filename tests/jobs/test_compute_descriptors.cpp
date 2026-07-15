#include "../helpers.hpp"

#include <cstdlib>
#include <iostream>
#include <string>

using namespace cdftt_tests;

/*
    About this job type (ComputeDescriptors):
        - It computes conceptual DFT descriptors from exactly three cube grids and,
            for this test shape, a pair of electronic energies.
        - It accepts partition strategies dedicated to descriptor workflows
            (the happy path here uses On-Grid).
        - It must reject partition strategies that are explicitly disallowed for this
            job, such as BBS and B2S.

    Detailed header command expectations:
        - Mandatory header line:
                    RunType=ComputeDescriptors
        - Required companion lines for this smoke scenario:
                    GridFilesNames=<grid1>,<grid2>,<grid3>
                    PartitionMethod=On-Grid
                    Energies=<IP>,<EA>

    Expected behavior to anchor assertions:
        - Happy path:
                1. Exit code 0.
                2. stdout includes partition method echo and descriptor-table markers.
        - Negative path:
                1. Non-zero exit code.
                2. stderr explicitly reports forbidden partition method for this job.
*/

// [P0 | PR-fast] Fixture-backed happy path with three grids + energies.
static void test_compute_descriptors_smoke_3grid_energies() {
    // Use repository fixtures: three checked-in cube grids and explicit energies.
    const std::string root = repo_root();
    const std::string g1 = root + "/examples/Density/Density.cube";
    const std::string g2 = root + "/examples/ELF/ELF.cube";
    const std::string g3 = root + "/examples/Density/Density.cube";
    const std::string input_path = "/tmp/cdftt_test_descriptors_smoke_3grid_energies.txt";

    write_input_file(input_path,
                     std::string("RunType=ComputeDescriptors\n")
                   + "GridFilesNames=" + g1 + "," + g2 + "," + g3 + "\n"
                   + "PartitionMethod=On-Grid\n"
                   + "Energies=0.30075,0.02092\n");

    // The run should succeed and echo the chosen partition method; the descriptor
    // table header is the canonical anchor used across descriptor tests.
    const RunResult r = run_cdftt(input_path);
    assert_zero_exit(r);
    assert_stdout_contains(r, "Volume partition method: On-Grid");
    assert_stdout_contains(r, "Reading Ionization potential I =");
    assert_stdout_contains(r, "Time in ms:");

    // Clean up the temporary input file.
    std::remove(input_path.c_str());
}

// [P0 | PR-fast] Negative path: BBS is forbidden for ComputeDescriptors.
static void test_compute_descriptors_invalid_partition_bbs() {
    // The run should fail and print an explicit error mentioning invalid partitionMethod.
    const std::string input_path = "/tmp/cdftt_test_descriptors_invalid_bbs.txt";
    write_input_file(input_path,
                     "RunType=ComputeDescriptors\n"
                     "PartitionMethod=BBS\n");

    // The run should fail and print an explicit error mentioning invalid partitionMethod.
    const RunResult r = run_cdftt(input_path);
    assert_nonzero_exit(r);
    assert_stderr_contains(r, "partitionMethod \"BBS\" invalid for this job");

    // Clean up the temporary input file.
    std::remove(input_path.c_str());
}

// [P0 | PR-fast] Negative path: B2S is forbidden for ComputeDescriptors.
static void test_compute_descriptors_invalid_partition_b2s() {
    // The run should fail and print an explicit error mentioning invalid partitionMethod.
    const std::string input_path = "/tmp/cdftt_test_descriptors_invalid_b2s.txt";
    write_input_file(input_path,
                     "RunType=ComputeDescriptors\n"
                     "PartitionMethod=B2S\n");
    
    // The run should fail and print an explicit error mentioning invalid partitionMethod.
    const RunResult r = run_cdftt(input_path);
    assert_nonzero_exit(r);
    assert_stderr_contains(r, "partitionMethod \"B2S\" invalid for this job");

    // Clean up the temporary input file.
    std::remove(input_path.c_str());
}

// [P1 | nightly] FD path with 3 analytic files should succeed and produce expected output structure.
static void test_fd_path_with_3_analytic_files() {

    // Use repository fixtures: three checked-in cube grids and three analytic
    // orbital files. FD partition reads energies from the provided analytic
    // files rather than from an explicit Energies= line.
    const std::string root = repo_root();
    const std::string g1 = root + "/examples/Density/Density.cube";
    const std::string g2 = root + "/examples/ELF/ELF.cube";
    const std::string g3 = root + "/examples/Density/Density.cube";
    const std::string a1 = root + "/tests/fixtures/Orbitals/h2o.molden";
    const std::string a2 = root + "/examples/ELF/h2o.molden";
    const std::string a3 = root + "/examples/Density/h2o.molden";

    const std::string input_path = "/tmp/cdftt_test_descriptors_fd_3_analytic.txt";

    write_input_file(input_path,
        "RunType=ComputeDescriptors\n"
        "Grid=" + g1 + "," + g2 + "," + g3 + "\n"
        "PartitionMethod=FD\n"
        "AnalyticFiles=" + a1 + "," + a2 + "," + a3 + "\n");

    const RunResult r = run_cdftt(input_path);

    // Successful run expected for this extended FD path.
    assert_zero_exit(r);

    // The run should echo the chosen partition method and indicate it read
    // energies from analytic files; the descriptor table header remains the
    // canonical anchor used across descriptor tests.
    assert_stdout_contains(r, "Volume partition method: FD");
    assert_stdout_contains(r, "Reading data from " + a1);
    assert_stdout_contains(r, "Reading data from " + a2);
    assert_stdout_contains(r, "Reading data from " + a3);
    assert_stdout_contains(r, "Time in ms:");


    // Clean up the temporary input file.
    std::remove(input_path.c_str());
}

// [P1 | nightly] FMO path with 1 analytic file should succeed and produce expected output structure.
static void test_fmo_path_with_1_analytic_file() {

    // Use repository fixtures: three checked-in cube grids and one analytic
    // orbital file. FMO partition reads energies from the provided analytic
    // file rather than from an explicit Energies= line.
    const std::string root = repo_root();
    const std::string g1 = root + "/examples/Density/Density.cube";
    const std::string g2 = root + "/examples/ELF/ELF.cube";
    const std::string g3 = root + "/examples/Density/Density.cube";

    //? Molden orbitals fail wihout clear error with FMO.
    //? WFX and FCHK orbitals succeed

    // const std::string a1 = root + "/examples/Orbitals/h2o.molden";
    const std::string a1 = root + "/tests/fixtures/Orbitals/h2o.fchk";

    const std::string input_path = "/tmp/cdftt_test_descriptors_fmo_1_analytic.txt";

    write_input_file(input_path,
        "RunType=ComputeDescriptors\n"
        "GridFilesNames=" + g1 + "," + g2 + "," + g3 + "\n"
        "PartitionMethod=FMO\n"
        "AnalyticFiles=" + a1 + "\n");

    const RunResult r = run_cdftt(input_path);

    // Successful run expected for this extended FMO path.
    assert_zero_exit(r);

    // The run should echo the chosen partition method and indicate it read
    // energies from the analytic file; the descriptor table header remains the
    // canonical anchor used across descriptor tests.
    assert_stdout_contains(r, "Volume partition method: FMO");
    assert_stdout_contains(r, "Reading data from " + a1);
    assert_stdout_contains(r, "Time in ms:");

    // Clean up the temporary input file.
    std::remove(input_path.c_str());
}

// [ P1 | nightly ] verify behavior for Energies list of size 2 and 3
static void test_compute_descriptors_various_energies_list_sizes() {
    // Use repository fixtures: three checked-in cube grids and explicit energies.
    const std::string root = repo_root();
    const std::string g1 = root + "/examples/Density/Density.cube";
    const std::string g2 = root + "/examples/ELF/ELF.cube";
    const std::string g3 = root + "/examples/Density/Density.cube";
    const std::string input_path = "/tmp/cdftt_test_descriptors_various_energies.txt";

    // Test with two energies (IP and EA).
    write_input_file(input_path,
        std::string("RunType=ComputeDescriptors\n")
        + "GridFilesNames=" + g1 + "," + g2 + "," + g3 + "\n"
        + "PartitionMethod=On-Grid\n"
        + "Energies=0.30075,0.02092\n");
        
    const RunResult r2 = run_cdftt(input_path);
    assert_zero_exit(r2);
    // TODO: refactor the output for better casing (capital letters are randomly at the start of words)
    assert_stdout_contains(r2, "Reading Ionization potential I = 0.30075");
    assert_stdout_contains(r2, "and electron Affinity A = 0.02092");
    assert_stdout_contains(r2, "Time in ms:");

    // Test with three energies (which will just be read as total energies).
    write_input_file(input_path,
                     std::string("RunType=ComputeDescriptors\n")
                   + "GridFilesNames=" + g1 + "," + g2 + "," + g3 + "\n"
                   + "PartitionMethod=On-Grid\n"
                   + "Energies=0.30075,0.02092,0.12345\n");

    const RunResult r3 = run_cdftt(input_path);
    assert_zero_exit(r3);
    assert_stdout_contains(r3, "Reading Total Energies: E1 = 0.30075");
    assert_stdout_contains(r3, "E2 = 0.02092");
    assert_stdout_contains(r3, "and E3 = 0.12345");
    assert_stdout_contains(r3, "Time in ms:");

    // Clean up the temporary input file.
    std::remove(input_path.c_str());
}

int main() {
    const Test tests[] = {
        { "test_compute_descriptors_smoke_3grid_energies",   test_compute_descriptors_smoke_3grid_energies   },
        { "test_compute_descriptors_invalid_partition_bbs", test_compute_descriptors_invalid_partition_bbs },
        { "test_compute_descriptors_invalid_partition_b2s", test_compute_descriptors_invalid_partition_b2s },
        { "test_fd_path_with_3_analytic_files",             test_fd_path_with_3_analytic_files             },
        { "test_fmo_path_with_1_analytic_file",             test_fmo_path_with_1_analytic_file             },
        { "test_compute_descriptors_various_energies_list_sizes", test_compute_descriptors_various_energies_list_sizes },
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
