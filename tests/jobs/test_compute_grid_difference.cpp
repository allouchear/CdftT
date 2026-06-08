#include "../helpers.hpp"

#include <cstdlib>
#include <iostream>
#include <string>
// #include <cmath>

using namespace cdftt_tests;

/*
    About this job type (ComputeGridDifference):
        - It subtracts two cube grids and writes the resulting difference into a
            third output cube file.
        - This job is strict on argument count: exactly three paths are expected in
            the Grids field (left operand, right operand, output file).

    Detailed header command expectations:
        - Mandatory header line:
                    RunType=ComputeGridDifference
        - Required companion line:
                    Grids=<gridA>,<gridB>,<outputCube>

    Expected behavior to anchor assertions:
        - Happy path:
                1. Exit code 0.
                2. stdout confirms difference computation and output save.
                3. Output file exists at requested path.
        - Negative path:
                1. Non-zero exit code.
                2. stderr reports incorrect number of grid file names.
*/

// [P0 | PR-fast] Fixture-backed happy path.
// Uses checked-in grids from examples and verifies output file creation.
static void test_compute_grid_difference_smoke() {
    // As input files are required, this test is fixture-backed with checked-in cube grids from
    // existing test assets. The output path is explicitly defined in the input content and checked
    // for existence after the run.
    const std::string root = repo_root();
    const std::string in1 = root + "/examples/Density/Density.cube";
    const std::string in2 = root + "/examples/ELF/ELF.cube";
    const std::string out = "/tmp/cdftt_test_grid_difference_output.cube";
    const std::string input_path = "/tmp/cdftt_test_grid_difference_smoke.txt";

    write_input_file(input_path,
                     std::string("RunType=ComputeGridDifference\n")
                   + "Grids=" + in1 + "," + in2 + "," + out + "\n");

    // The run should succeed and print markers about computing the difference and saving the output
    // file; the output file should be created at the requested path.
    const RunResult r = run_cdftt(input_path);

    assert_zero_exit(r);
    assert_stdout_contains(r, "Computing difference between");
    assert_stdout_contains(r, "Difference grid saved to file");
    assert_file_exists(out);

    // Clean up the temporary input and output files.
    std::remove(input_path.c_str());
    std::remove(out.c_str());
}

// [P0 | PR-fast] Negative path: exactly three grid paths are required.
static void test_compute_grid_difference_invalid_grid_count() {
    // The run should fail and print an explicit error mentioning incorrect number of grid file names.
    const std::string input_path = "/tmp/cdftt_test_grid_difference_invalid_count.txt";
    write_input_file(input_path,
                     "RunType=ComputeGridDifference\n"
                     "Grids=a.cube,b.cube\n");

    const RunResult r = run_cdftt(input_path);

    assert_nonzero_exit(r);
    assert_stderr_contains(r, "incorrect number of grid files names (three files expected)");

    // Clean up the temporary input file.
    std::remove(input_path.c_str());
}

// [P1 | nightly] Numeric spot-check: same input file twice so the expected difference is ~0.0.
static void test_compute_grid_difference_numeric_spotcheck() {
    const std::string root = repo_root();
    const std::string in1 = root + "/examples/Density/Density.cube";
    const std::string in2 = root + "/examples/Density/Density.cube"; // same file -> zero difference
    const std::string out = "/tmp/cdftt_test_grid_difference_output_spotcheck.cube";
    const std::string input_path = "/tmp/cdftt_test_grid_difference_spotcheck.txt";

    write_input_file(input_path,
                     std::string("RunType=ComputeGridDifference\n")
                   + "Grids=" + in1 + "," + in2 + "," + out + "\n");

    const RunResult r = run_cdftt(input_path);

    assert_zero_exit(r);
    assert_stdout_contains(r, "Difference grid saved to file");
    assert_file_exists(out);

    // Parse first N voxel values and ensure they are approximately zero.
    const auto vox = read_first_n_voxels(out, 32);
    if (vox.empty())
        throw std::runtime_error("Expected at least one voxel in output cube for spot-check");

    for (std::size_t i = 0; i < vox.size(); ++i) {
        if (std::abs(vox[i]) > 1e-12)
            throw std::runtime_error("Voxel " + std::to_string(i) + " differs from expected 0 by " + std::to_string(vox[i]));
    }

    // Clean up
    std::remove(out.c_str());
    std::remove(input_path.c_str());
}

// [P1 | nightly] Numeric verification: sampled voxels in output approx equal in1 - in2
static void test_compute_grid_difference_numeric_compare_inputs() {
    const std::string root = repo_root();
    const std::string in1 = root + "/examples/Density/Density.cube";
    const std::string in2 = root + "/examples/ELF/ELF.cube";
    const std::string out = "/tmp/cdftt_test_grid_difference_output_compare.cube";
    const std::string input_path = "/tmp/cdftt_test_grid_difference_compare.txt";

    write_input_file(input_path,
                     std::string("RunType=ComputeGridDifference\n")
                   + "Grids=" + in1 + "," + in2 + "," + out + "\n");

    const RunResult r = run_cdftt(input_path);

    assert_zero_exit(r);
    assert_stdout_contains(r, "Difference grid saved to file");
    assert_file_exists(out);

    // Read first N voxels from each file and compare
    const std::size_t N = 64;
    const auto v1 = read_first_n_voxels(in1, N);
    const auto v2 = read_first_n_voxels(in2, N);
    const auto vo = read_first_n_voxels(out, N);

    const std::size_t m = std::min(std::min(v1.size(), v2.size()), vo.size());
    if (m == 0)
        throw std::runtime_error("Insufficient voxel data for comparison");

    const double tol = 1e-12; // tolerant epsilon for numeric differences
    for (std::size_t i = 0; i < m; ++i) {
        const double expected = v1[i] - v2[i];
        const double actual = vo[i];
        if (std::abs(actual - expected) > tol)
            throw std::runtime_error("Voxel " + std::to_string(i) + " mismatch: expected "
                                     + std::to_string(expected) + " got " + std::to_string(actual));
    }

    // Clean up
    std::remove(out.c_str());
    std::remove(input_path.c_str());
}

int main() {
    const Test tests[] = {
        { "test_compute_grid_difference_smoke",              test_compute_grid_difference_smoke              },
        { "test_compute_grid_difference_invalid_grid_count", test_compute_grid_difference_invalid_grid_count },
        { "test_compute_grid_difference_numeric_spotcheck",  test_compute_grid_difference_numeric_spotcheck  },
        { "test_compute_grid_difference_numeric_compare",    test_compute_grid_difference_numeric_compare_inputs },
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

