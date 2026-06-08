#include "../helpers.hpp"

#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>

using namespace cdftt_tests;

/*
    About this job type (MakeDensityCube):
        - It builds a 3D grid domain from an analytic orbital source and writes an
            electron-density cube file.
        - The happy path here is fixture-backed and checks not only process success
            but also output artifact quality (parseable cube header and finite values).
        - Unsupported analytic file formats must fail with a clear parser error.

    Detailed header command expectations:
        - Mandatory header line:
                    RunType=MakeDensityCube
        - Required companion lines:
                    AnalyticFiles=<inputOrbitalFile>
                    Size=<Small|Medium|Large|...>
                    Grids=<outputCubePath>

    Expected behavior to anchor assertions:
        - Happy path:
                1. Exit code 0.
                2. stdout contains build and generation progress markers.
                3. Output cube exists, is parseable, and contains finite voxel values.
        - Negative path:
                1. Non-zero exit code.
                2. stderr reports invalid analytic file format.
*/

// [P0 | PR-fast] Fixture-backed happy path for medium-size density cube creation.
static void test_make_density_cube_smoke() {
    // As input files are required, this test is fixture-backed with a checked-in molden file from
    // existing test assets. The output path is explicitly defined in the input content and checked
    // for existence and content quality after the run.
    const std::string root = repo_root();
    const std::string analytic = root + "/examples/Density/h2o.molden";
    const std::string output = "/tmp/cdftt_test_make_density_cube_output.cube";
    const std::string input_path = "/tmp/cdftt_test_make_density_cube_smoke.txt";

    write_input_file(input_path,
                     std::string("RunType=MakeDensityCube\n")
                   + "AnalyticFiles=" + analytic + "\n"
                   + "Size=Medium\n"
                   + "Grids=" + output + "\n");

    // The run should succeed and print markers about building the domain and creating the density
    // cube; the output file should be created at the requested path, have a parseable header, and
    // contain finite values.
    const RunResult r = run_cdftt(input_path);

    assert_zero_exit(r);
    assert_stdout_contains(r, "Building domain, please wait...");
    assert_stdout_contains(r, "Creating density cube, please wait...");
    assert_stdout_contains(r, "Density cube saved to file");
    assert_file_exists(output);
    assert_cube_header_parseable(output);
    assert_cube_has_finite_values(output);

    // Clean up the temporary input and output files.
    std::remove(input_path.c_str());
    std::remove(output.c_str());
}

// [P0 | PR-fast] Negative path: unsupported analytic extension fails clearly.
static void test_make_density_cube_invalid_input_extension() {
    // The run should fail and print an explicit error mentioning invalid file format for
    // the unsupported extension.
    const std::string input_path = "/tmp/cdftt_test_make_density_cube_invalid_extension.txt";

    write_input_file(input_path,
                     "RunType=MakeDensityCube\n"
                     "AnalyticFiles=dummy.xyz\n"
                     "Size=Medium\n"
                     "Grids=/tmp/cdftt_should_not_exist_density.cube\n");

    const RunResult r = run_cdftt(input_path);

    assert_nonzero_exit(r);
    assert_stderr_contains(r, "invalid file format for file \"dummy.xyz\"");

    // Clean up the temporary input file.
    std::remove(input_path.c_str());
}

// [P1 | nightly] Size=Coarse produces expected dimensions
static void test_make_density_cube_size_coarse() {
    const std::string root = repo_root();
    const std::string analytic = root + "/examples/Density/h2o.molden";
    const std::string output = "/tmp/cdftt_test_make_density_cube_coarse.cube";
    const std::string input_path = "/tmp/cdftt_test_make_density_cube_size_coarse.txt";

    write_input_file(input_path,
                     std::string("RunType=MakeDensityCube\n")
                   + "AnalyticFiles=" + analytic + "\n"
                   + "Size=Coarse\n"
                   + "Grids=" + output + "\n");

    const RunResult r = run_cdftt(input_path);

    assert_zero_exit(r);
    assert_stdout_contains(r, "Building domain, please wait...");
    assert_stdout_contains(r, "Creating density cube, please wait...");
    assert_stdout_contains(r, "Density cube saved to file");
    assert_file_exists(output);
    assert_cube_header_parseable(output);
    assert_cube_has_finite_values(output);

    // Verify that coarse size produces reasonable small dimensions
    // Coarse -> 3 points per Bohr radius, so we expect <100Kb for a small molecule like H2O
    assert_file_size_less_than(output, 100 * 1024);
    
    // Clean up
    std::remove(input_path.c_str());
    std::remove(output.c_str());
}

// [P1 | nightly] Size=Medium produces expected dimensions
static void test_make_density_cube_size_medium() {
    const std::string root = repo_root();
    const std::string analytic = root + "/examples/Density/h2o.molden";
    const std::string output = "/tmp/cdftt_test_make_density_cube_medium.cube";
    const std::string input_path = "/tmp/cdftt_test_make_density_cube_size_medium.txt";

    write_input_file(input_path,
                     std::string("RunType=MakeDensityCube\n")
                   + "AnalyticFiles=" + analytic + "\n"
                   + "Size=Medium\n"
                   + "Grids=" + output + "\n");

    const RunResult r = run_cdftt(input_path);

    assert_zero_exit(r);
    assert_stdout_contains(r, "Building domain, please wait...");
    assert_stdout_contains(r, "Creating density cube, please wait...");
    assert_stdout_contains(r, "Density cube saved to file");
    assert_file_exists(output);
    assert_cube_header_parseable(output);
    assert_cube_has_finite_values(output);

    // Verify that medium size produces reasonable dimensions
    // Medium -> 6 points per Bohr radius, so we expect <1Mb for a small molecule like H2O
    assert_file_size_less_than(output, 1024 * 1024);
    
    // Clean up
    std::remove(input_path.c_str());
    std::remove(output.c_str());
}

// [P1 | nightly] Size=Fine produces expected dimensions
static void test_make_density_cube_size_fine() {
    const std::string root = repo_root();
    const std::string analytic = root + "/examples/Density/h2o.molden";
    const std::string output = "/tmp/cdftt_test_make_density_cube_fine.cube";
    const std::string input_path = "/tmp/cdftt_test_make_density_cube_size_fine.txt";

    write_input_file(input_path,
                     std::string("RunType=MakeDensityCube\n")
                   + "AnalyticFiles=" + analytic + "\n"
                   + "Size=Fine\n"
                   + "Grids=" + output + "\n");

    const RunResult r = run_cdftt(input_path);

    assert_zero_exit(r);
    assert_stdout_contains(r, "Building domain, please wait...");
    assert_stdout_contains(r, "Creating density cube, please wait...");
    assert_stdout_contains(r, "Density cube saved to file");
    assert_file_exists(output);
    assert_cube_header_parseable(output);
    assert_cube_has_finite_values(output);

    // Verify that fine size produces reasonable large dimensions
    // Fine -> 12 points per Bohr radius, so we expect <10Mb for a small molecule like H2O
    assert_file_size_less_than(output, 10 * 1024 * 1024);
    
    // Clean up
    std::remove(input_path.c_str());
    std::remove(output.c_str());
}


// [P1 | nightly] Size=Custom produces expected dimensions
static void test_make_density_cube_size_custom() {
    const std::string root = repo_root();
    const std::string analytic = root + "/examples/Density/h2o.molden";
    const std::string output = "/tmp/cdftt_test_make_density_cube_custom.cube";
    const std::string input_path = "/tmp/cdftt_test_make_density_cube_size_custom.txt";

    write_input_file(input_path,
                     std::string("RunType=MakeDensityCube\n")
                   + "AnalyticFiles=" + analytic + "\n"
                   + "Size=Custom\n"
                   + "CustomSizeData=80,80,80,5,5,5,0.15,0,0,0,0.15,0,0,0,0.15\n"
                   + "Grids=" + output + "\n");

    const RunResult r = run_cdftt(input_path);

    assert_zero_exit(r);
    assert_stdout_contains(r, "Building domain, please wait...");
    assert_stdout_contains(r, "Creating density cube, please wait...");
    assert_stdout_contains(r, "Density cube saved to file");
    assert_file_exists(output);
    assert_cube_header_parseable(output);
    assert_cube_has_finite_values(output);

    // Verify that custom size produces expected dimensions
    // For this exact grid and for water, we expect around 11Mb
    assert_file_size_less_than(output, 12 * 1024 * 1024);
    
    // Clean up
    std::remove(input_path.c_str());
    std::remove(output.c_str());
}

// [P1 | nightly] CustomSizeData malformed input fails clearly
static void test_make_density_cube_size_custom_malformed() {
    const std::string root = repo_root();
    const std::string analytic = root + "/examples/Density/h2o.molden";
    const std::string output = "/tmp/cdftt_test_make_density_cube_custom_malformed.cube";
    const std::string input_path = "/tmp/cdftt_test_make_density_cube_size_custom_malformed.txt";

    write_input_file(input_path,
                     std::string("RunType=MakeDensityCube\n")
                   + "AnalyticFiles=" + analytic + "\n"
                   + "Size=Custom\n"
                   + "CustomSizeData=80,80,80,5,5,\n"
                   + "Grids=" + output + "\n");

    const RunResult r = run_cdftt(input_path);

    // Assert failure with an explicit error about malformed CustomSizeData (incorrect number of values)
    assert_nonzero_exit(r);
    assert_stderr_contains(r, "Error: incorrect number of values for the \"CustomSizeData\" parameter (fifteen values expected)");

    // Clean up
    std::remove(input_path.c_str());
    std::remove(output.c_str());
}

int main() {
    const Test tests[] = {
        { "test_make_density_cube_smoke",                   test_make_density_cube_smoke                   },
        { "test_make_density_cube_invalid_input_extension", test_make_density_cube_invalid_input_extension },
        { "test_make_density_cube_size_coarse",            test_make_density_cube_size_coarse            },
        { "test_make_density_cube_size_medium",           test_make_density_cube_size_medium           },
        { "test_make_density_cube_size_fine",              test_make_density_cube_size_fine              },
        { "test_make_density_cube_size_custom",             test_make_density_cube_size_custom             },
        { "test_make_density_cube_size_custom_malformed",  test_make_density_cube_size_custom_malformed  },
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
