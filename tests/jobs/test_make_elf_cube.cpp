#include "../helpers.hpp"

#include <cstdlib>
#include <iostream>
#include <string>

using namespace cdftt_tests;

/*
    About this job type (MakeELFCube):
        - It computes ELF values from analytic orbitals and writes them to a cube.
        - The happy path uses a fixture-backed analytic source and a valid ELF method
            (Becke in this file).
        - Unknown or unsupported ELF methods must fail explicitly.

    Detailed header command expectations:
        - Mandatory header line:
                    RunType=MakeELFCube
        - Required companion lines:
                    AnalyticFiles=<inputOrbitalFile>
                    Size=<Small|Medium|Large|...>
                    ELFMethod=<method>
                    Grids=<outputCubePath>

    Expected behavior to anchor assertions:
        - Happy path:
                1. Exit code 0.
                2. stdout contains domain build and ELF generation markers.
                3. Output cube exists, header parses, and values are finite.
        - Negative path:
                1. Non-zero exit code.
                2. stderr reports unknown ELF method.
*/

// [P0 | PR-fast] Fixture-backed happy path for ELF cube generation.
static void test_make_elf_cube_smoke() {
    // As input files are required, this test is fixture-backed with a checked-in molden
    // file from existing test assets. The ELF method and output paths are explicitly
    // defined in the input content and checked for existence and content quality after the run.
    const std::string root = repo_root();
    const std::string analytic = root + "/examples/ELF/h2o.molden";
    const std::string output = "/tmp/cdftt_test_make_elf_cube_output.cube";
    const std::string input_path = "/tmp/cdftt_test_make_elf_cube_smoke.txt";

    write_input_file(input_path,
                     std::string("RunType=MakeELFCube\n")
                   + "AnalyticFiles=" + analytic + "\n"
                   + "Size=Medium\n"
                   + "ELFMethod=Becke\n"
                   + "Grids=" + output + "\n");

    // The run should succeed and print markers about building the domain and creating the ELF
    // cube; the output file should be created at the requested path, have a parseable header, and
    // contain finite values.
    const RunResult r = run_cdftt(input_path);

    assert_zero_exit(r);
    assert_stdout_contains(r, "Building domain, please wait...");
    assert_stdout_contains(r, "Creating ELF cube, please wait...");
    assert_stdout_contains(r, "ELF cube saved to file");
    assert_file_exists(output);
    assert_cube_header_parseable(output);
    assert_cube_has_finite_values(output);

    // Clean up the temporary input and output files.
    std::remove(input_path.c_str());
    std::remove(output.c_str());
}

// [P0 | PR-fast] Negative path: invalid ELF method fails explicitly.
static void test_make_elf_cube_invalid_method() {
    // The run should fail and print an explicit error mentioning unknown ELF method.
    const std::string root = repo_root();
    const std::string analytic = root + "/examples/ELF/h2o.molden";
    const std::string input_path = "/tmp/cdftt_test_make_elf_cube_invalid_method.txt";

    write_input_file(input_path,
                     std::string("RunType=MakeELFCube\n")
                   + "AnalyticFiles=" + analytic + "\n"
                   + "Size=Medium\n"
                   + "ELFMethod=BadMethod\n"
                   + "Grids=/tmp/cdftt_should_not_exist_elf.cube\n");

    const RunResult r = run_cdftt(input_path);

    assert_nonzero_exit(r);
    assert_stderr_contains(r, "ELF method \"BadMethod\" unknown");

    // Clean up the temporary input file.
    std::remove(input_path.c_str());
}

// [P1 | nightly] Test ELFMethod=Becke specifically
static void test_make_elf_cube_elf_method_becke() {
    const std::string root = repo_root();
    const std::string analytic = root + "/examples/ELF/h2o.molden";
    const std::string output = "/tmp/cdftt_test_make_elf_cube_becke.cube";
    const std::string input_path = "/tmp/cdftt_test_make_elf_cube_elf_method_becke.txt";

    write_input_file(input_path,
                     std::string("RunType=MakeELFCube\n")
                   + "AnalyticFiles=" + analytic + "\n"
                   + "Size=Medium\n"
                   + "ELFMethod=Becke\n"
                   + "Grids=" + output + "\n");

    const RunResult r = run_cdftt(input_path);

    assert_zero_exit(r);
    assert_stdout_contains(r, "Building domain, please wait...");
    assert_stdout_contains(r, "Creating ELF cube, please wait...");
    assert_stdout_contains(r, "ELF cube saved to file");
    assert_file_exists(output);
    assert_cube_header_parseable(output);
    assert_cube_has_finite_values(output);

    // Clean up
    std::remove(input_path.c_str());
    std::remove(output.c_str());
}

// [P1 | nightly] Test Size=Coarse produces expected dimensions for ELF cube
static void test_make_elf_cube_size_coarse() {
    const std::string root = repo_root();
    const std::string analytic = root + "/examples/ELF/h2o.molden";
    const std::string output = "/tmp/cdftt_test_make_elf_cube_size_coarse.cube";
    const std::string input_path = "/tmp/cdftt_test_make_elf_cube_size_coarse.txt";

    write_input_file(input_path,
                     std::string("RunType=MakeELFCube\n")
                   + "AnalyticFiles=" + analytic + "\n"
                   + "Size=Coarse\n"
                   + "ELFMethod=Becke\n"
                   + "Grids=" + output + "\n");

    const RunResult r = run_cdftt(input_path);

    assert_zero_exit(r);
    assert_stdout_contains(r, "Building domain, please wait...");
    assert_stdout_contains(r, "Creating ELF cube, please wait...");
    assert_stdout_contains(r, "ELF cube saved to file");
    assert_file_exists(output);
    assert_cube_header_parseable(output);
    assert_cube_has_finite_values(output);

    // Clean up
    std::remove(input_path.c_str());
    std::remove(output.c_str());
}

// [P1 | nightly] Test Size=Medium produces expected dimensions for ELF cube
static void test_make_elf_cube_size_medium() {
    const std::string root = repo_root();
    const std::string analytic = root + "/examples/ELF/h2o.molden";
    const std::string output = "/tmp/cdftt_test_make_elf_cube_size_medium.cube";
    const std::string input_path = "/tmp/cdftt_test_make_elf_cube_size_medium.txt";

    write_input_file(input_path,
                     std::string("RunType=MakeELFCube\n")
                   + "AnalyticFiles=" + analytic + "\n"
                   + "Size=Medium\n"
                   + "ELFMethod=Becke\n"
                   + "Grids=" + output + "\n");

    const RunResult r = run_cdftt(input_path);

    assert_zero_exit(r);
    assert_stdout_contains(r, "Building domain, please wait...");
    assert_stdout_contains(r, "Creating ELF cube, please wait...");
    assert_stdout_contains(r, "ELF cube saved to file");
    assert_file_exists(output);
    assert_cube_header_parseable(output);
    assert_cube_has_finite_values(output);

    // Clean up
    std::remove(input_path.c_str());
    std::remove(output.c_str());
}

// [P1 | nightly] Test Size=Fine produces expected dimensions for ELF cube
static void test_make_elf_cube_size_fine() {
    const std::string root = repo_root();
    const std::string analytic = root + "/examples/ELF/h2o.molden";
    const std::string output = "/tmp/cdftt_test_make_elf_cube_size_fine.cube";
    const std::string input_path = "/tmp/cdftt_test_make_elf_cube_size_fine.txt";

    write_input_file(input_path,
                     std::string("RunType=MakeELFCube\n")
                   + "AnalyticFiles=" + analytic + "\n"
                   + "Size=Fine\n"
                   + "ELFMethod=Becke\n"
                   + "Grids=" + output + "\n");

    const RunResult r = run_cdftt(input_path);

    assert_zero_exit(r);
    assert_stdout_contains(r, "Building domain, please wait...");
    assert_stdout_contains(r, "Creating ELF cube, please wait...");
    assert_stdout_contains(r, "ELF cube saved to file");
    assert_file_exists(output);
    assert_cube_header_parseable(output);
    assert_cube_has_finite_values(output);

    // Clean up
    std::remove(input_path.c_str());
    std::remove(output.c_str());
}

// [P1 | nightly] Test Size=Custom produces expected dimensions for ELF cube
static void test_make_elf_cube_size_custom() {
    const std::string root = repo_root();
    const std::string analytic = root + "/examples/ELF/h2o.molden";
    const std::string output = "/tmp/cdftt_test_make_elf_cube_size_custom.cube";
    const std::string input_path = "/tmp/cdftt_test_make_elf_cube_size_custom.txt";

    write_input_file(input_path,
                     std::string("RunType=MakeELFCube\n")
                   + "AnalyticFiles=" + analytic + "\n"
                   + "Size=Custom\n"
                   + "CustomSizeData=80,80,80,5,5,5,0.15,0,0,0,0.15,0,0,0,0.15\n"
                   + "ELFMethod=Becke\n"
                   + "Grids=" + output + "\n");

    const RunResult r = run_cdftt(input_path);

    assert_zero_exit(r);
    assert_stdout_contains(r, "Building domain, please wait...");
    assert_stdout_contains(r, "Creating ELF cube, please wait...");
    assert_stdout_contains(r, "ELF cube saved to file");
    assert_file_exists(output);
    assert_cube_header_parseable(output);
    assert_cube_has_finite_values(output);

    // Clean up
    std::remove(input_path.c_str());
    std::remove(output.c_str());
}

int main() {
    const Test tests[] = {
        { "test_make_elf_cube_smoke",          test_make_elf_cube_smoke          },
        { "test_make_elf_cube_invalid_method", test_make_elf_cube_invalid_method },
        { "test_make_elf_cube_elf_method_becke", test_make_elf_cube_elf_method_becke },
        { "test_make_elf_cube_size_coarse", test_make_elf_cube_size_coarse },
        { "test_make_elf_cube_size_medium", test_make_elf_cube_size_medium },
        { "test_make_elf_cube_size_fine", test_make_elf_cube_size_fine },
        { "test_make_elf_cube_size_custom", test_make_elf_cube_size_custom }
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
