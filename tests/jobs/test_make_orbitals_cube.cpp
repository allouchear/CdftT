#include "../helpers.hpp"

#include <cstdlib>
#include <iostream>
#include <string>

using namespace cdftt_tests;

/*
    About this job type (MakeOrbitalsCube):
        - It evaluates orbital-related scalar data on a generated domain and writes
            one or more output cube artifacts.
        - The happy path in this suite is currently anchored to a .molden fixture,
            which is stable in this environment.
        - Invalid orbital selection directives (OrbitalType) must fail explicitly.

    Detailed header command expectations:
        - Mandatory header line:
                    RunType=MakeOrbitalsCube
        - Required companion lines:
                    AnalyticFiles=<inputOrbitalFile>
                    Size=<Small|Medium|Large|...>
                    GridFilesNames=<outputCubePath>
        - Optional selector line (validated by negative test):
                    OrbitalType=<selection>

    Expected behavior to anchor assertions:
        - Happy path:
                1. Exit code 0.
                2. stdout contains domain build and cube generation markers.
                3. Output cube exists, header parses, and all checked values are finite.
        - Negative path:
                1. Non-zero exit code.
                2. stderr reports unknown orbital type for invalid selector value.
*/

// [P0 | PR-fast] Fixture-backed happy path for default/all orbitals cube generation.
static void test_make_orbitals_cube_smoke() {
    // As input files are required, this test is fixture-backed with a checked-in molden file
    // from existing test assets. The output path is explicitly defined in the input content
    // and checked for existence and content quality after the run.
    const std::string root = repo_root();

    const std::string analytic = root + "/tests/fixtures/Orbitals/h2o.fchk";
    // ? Manually verified working formats: `.molden`, `.fchk`, `.wfx`

    const std::string output = "/tmp/cdftt_test_make_orbitals_cube_output.cube";
    const std::string input_path = "/tmp/cdftt_test_make_orbitals_cube_smoke.txt";

    write_input_file(input_path,
                     std::string("RunType=MakeOrbitalsCube\n")
                   + "AnalyticFiles=" + analytic + "\n"
                   + "Size=Medium\n"
                   + "GridFilesNames=" + output + "\n");

    // The run should succeed and print markers about building the domain and creating the orbitals
    const RunResult r = run_cdftt(input_path);

    assert_zero_exit(r);
    assert_stdout_contains(r, "Reading data from " + analytic + "... Please wait.");
    assert_stdout_contains(r, "Writing cube file, please wait...");
    assert_stdout_contains(r, "Density cube saved to " + output);
    assert_file_exists(output);
    assert_cube_header_parseable(output);
    assert_cube_has_finite_values(output);

    // Clean up the temporary input and output files.
    std::remove(input_path.c_str());
    std::remove(output.c_str());
}

// [P0 | PR-fast] Negative path: invalid OrbitalType value fails explicitly.
static void test_make_orbitals_cube_invalid_selection() {
    // The run should fail and print an explicit error mentioning unknown orbital type.
    const std::string root = repo_root();
    const std::string analytic = root + "/tests/fixtures/Orbitals/h2o.molden";
    const std::string input_path = "/tmp/cdftt_test_make_orbitals_cube_invalid_selection.txt";

    write_input_file(input_path,
                     std::string("RunType=MakeOrbitalsCube\n")
                   + "AnalyticFiles=" + analytic + "\n"
                   + "Size=Medium\n"
                   + "GridFilesNames=/tmp/cdftt_should_not_exist_orbitals.cube\n"
                   + "OrbitalType=BadType\n");

    const RunResult r = run_cdftt(input_path);

    assert_nonzero_exit(r);
    assert_stderr_contains(r, "Orbital type \"BadType\" unknown");

    // Clean up the temporary input file.
    std::remove(input_path.c_str());
}

// [P1 | nightly] Test OrbitalType=All mode
static void test_make_orbitals_cube_orbital_type_all() {
    const std::string root = repo_root();
    const std::string analytic = root + "/tests/fixtures/Orbitals/h2o.molden";
    const std::string output = "/tmp/cdftt_test_make_orbitals_cube_orbital_type_all.cube";
    const std::string input_path = "/tmp/cdftt_test_make_orbitals_cube_orbital_type_all.txt";

    write_input_file(input_path,
                     std::string("RunType=MakeOrbitalsCube\n")
                   + "AnalyticFiles=" + analytic + "\n"
                   + "Size=Medium\n"
                   + "OrbitalType=All\n"
                   + "GridFilesNames=" + output + "\n");

    const RunResult r = run_cdftt(input_path);

    assert_zero_exit(r);
    assert_stdout_contains(r, "Reading data from " + analytic + "... Please wait.");
    assert_stdout_contains(r, "Writing cube file, please wait...");
    assert_stdout_contains(r, "Density cube saved to " + output);
    assert_file_exists(output);
    assert_cube_header_parseable(output);
    assert_cube_has_finite_values(output);

    // Clean up
    std::remove(input_path.c_str());
    std::remove(output.c_str());
}

// [P1 | nightly] Test OrbitalType=Occupied mode
static void test_make_orbitals_cube_orbital_type_occupied() {
    const std::string root = repo_root();
    const std::string analytic = root + "/tests/fixtures/Orbitals/h2o.molden";
    const std::string output = "/tmp/cdftt_test_make_orbitals_cube_orbital_type_occupied.cube";
    const std::string input_path = "/tmp/cdftt_test_make_orbitals_cube_orbital_type_occupied.txt";

    write_input_file(input_path,
                     std::string("RunType=MakeOrbitalsCube\n")
                   + "AnalyticFiles=" + analytic + "\n"
                   + "Size=Medium\n"
                   + "OrbitalType=Occupied\n"
                   + "SpinType=Alpha\n"
                   + "GridFilesNames=" + output + "\n");

    //! Using "Occ" as documentend returns `Error: Orbital type "Occ" unknown.`
    const RunResult r = run_cdftt(input_path);
    // TODO : rewrite documentation to clarify that "Occupied" is the correct value, not "Occ".

    assert_zero_exit(r);
    assert_stdout_contains(r, "Reading data from " + analytic + "... Please wait.");
    assert_stdout_contains(r, "Writing cube file, please wait...");
    assert_stdout_contains(r, "Density cube saved to " + output);
    assert_file_exists(output);
    assert_cube_header_parseable(output);
    assert_cube_has_finite_values(output);

    // Clean up
    std::remove(input_path.c_str());
    std::remove(output.c_str());
}

// [P1 | nightly] Test OrbitalType=Virtual mode
static void test_make_orbitals_cube_orbital_type_virtual() {
    const std::string root = repo_root();
    const std::string analytic = root + "/tests/fixtures/Orbitals/h2o.molden";
    const std::string output = "/tmp/cdftt_test_make_orbitals_cube_orbital_type_virtual.cube";
    const std::string input_path = "/tmp/cdftt_test_make_orbitals_cube_orbital_type_virtual.txt";

    write_input_file(input_path,
                     std::string("RunType=MakeOrbitalsCube\n")
                   + "AnalyticFiles=" + analytic + "\n"
                   + "Size=Medium\n"
                   + "OrbitalType=Virtual\n"
                   + "SpinType=Alpha\n"
                   + "GridFilesNames=" + output + "\n");

    const RunResult r = run_cdftt(input_path);
    // ! Exit code -1 with stderr :
    // !     free(): invalid next size (fast)

    assert_zero_exit(r);
    assert_stdout_contains(r, "Reading data from " + analytic + "... Please wait.");
    assert_stdout_contains(r, "Writing cube file, please wait...");
    assert_stdout_contains(r, "Density cube saved to " + output);
    assert_file_exists(output);
    assert_cube_header_parseable(output);
    assert_cube_has_finite_values(output);

    // Clean up
    std::remove(input_path.c_str());
    std::remove(output.c_str());
}

// [P1 | nightly] Test OrbitalType=Homo mode
static void test_make_orbitals_cube_orbital_type_homo() {
    const std::string root = repo_root();
    const std::string analytic = root + "/tests/fixtures/Orbitals/h2o.molden";
    const std::string output = "/tmp/cdftt_test_make_orbitals_cube_orbital_type_homo.cube";
    const std::string input_path = "/tmp/cdftt_test_make_orbitals_cube_orbital_type_homo.txt";

    write_input_file(input_path,
                     std::string("RunType=MakeOrbitalsCube\n")
                   + "AnalyticFiles=" + analytic + "\n"
                   + "Size=Medium\n"
                   + "OrbitalType=Homo\n"
                   + "GridFilesNames=" + output + "\n");

    const RunResult r = run_cdftt(input_path);

    assert_zero_exit(r);
    assert_stdout_contains(r, "Reading data from " + analytic + "... Please wait.");
    assert_stdout_contains(r, "Writing cube file, please wait...");
    assert_stdout_contains(r, "Density cube saved to " + output);
    assert_file_exists(output);
    assert_cube_header_parseable(output);
    assert_cube_has_finite_values(output);

    // Clean up
    std::remove(input_path.c_str());
    std::remove(output.c_str());
}

// [P1 | nightly] Test OrbitalType=Lumo mode
static void test_make_orbitals_cube_orbital_type_lumo() {
    const std::string root = repo_root();
    const std::string analytic = root + "/tests/fixtures/Orbitals/h2o.molden";
    const std::string output = "/tmp/cdftt_test_make_orbitals_cube_orbital_type_lumo.cube";
    const std::string input_path = "/tmp/cdftt_test_make_orbitals_cube_orbital_type_lumo.txt";

    write_input_file(input_path,
                     std::string("RunType=MakeOrbitalsCube\n")
                   + "AnalyticFiles=" + analytic + "\n"
                   + "Size=Medium\n"
                   + "OrbitalType=Lumo\n"
                   + "GridFilesNames=" + output + "\n");

    const RunResult r = run_cdftt(input_path);

    assert_zero_exit(r);
    assert_stdout_contains(r, "Reading data from " + analytic + "... Please wait.");
    assert_stdout_contains(r, "Writing cube file, please wait...");
    assert_stdout_contains(r, "Density cube saved to " + output);
    assert_file_exists(output);
    assert_cube_header_parseable(output);
    assert_cube_has_finite_values(output);

    // Clean up
    std::remove(input_path.c_str());
    std::remove(output.c_str());
}

// [P1 | nightly] Test OrbitalType=Homo-Lumo mode
static void test_make_orbitals_cube_orbital_type_homo_lumo() {
    const std::string root = repo_root();
    const std::string analytic = root + "/tests/fixtures/Orbitals/h2o.molden";
    const std::string output = "/tmp/cdftt_test_make_orbitals_cube_orbital_type_homo_lumo.cube";
    const std::string input_path = "/tmp/cdftt_test_make_orbitals_cube_orbital_type_homo_lumo.txt";

    write_input_file(input_path,
                     std::string("RunType=MakeOrbitalsCube\n")
                   + "AnalyticFiles=" + analytic + "\n"
                   + "Size=Medium\n"
                   + "OrbitalType=Homo-Lumo\n"
                   + "GridFilesNames=" + output + "\n");

    const RunResult r = run_cdftt(input_path);

    assert_zero_exit(r);
    assert_stdout_contains(r, "Reading data from " + analytic + "... Please wait.");
    assert_stdout_contains(r, "Writing cube file, please wait...");
    assert_stdout_contains(r, "Density cube saved to " + output);
    assert_file_exists(output);
    assert_cube_header_parseable(output);
    assert_cube_has_finite_values(output);

    // Clean up
    std::remove(input_path.c_str());
    std::remove(output.c_str());
}

// [P1 | nightly] Test OrbitalType=Custom mode
static void test_make_orbitals_cube_orbital_type_custom() {
    const std::string root = repo_root();
    const std::string analytic = root + "/examples/Orbitals/h2o.molden";
    const std::string output = "/tmp/cdftt_test_make_orbitals_cube_orbital_type_custom.cube";
    const std::string input_path = "/tmp/cdftt_test_make_orbitals_cube_orbital_type_custom.txt";

    write_input_file(input_path,
                     std::string("RunType=MakeOrbitalsCube\n")
                   + "AnalyticFiles=" + analytic + "\n"
                   + "Size=Medium\n"
                   + "OrbitalType=Custom\n"
                   + "OrbitalsNumbers=1,2,3\n"
                   + "SpinList=Alpha,Beta,Alpha\n"
                   + "GridFilesNames=" + output + "\n");

    const RunResult r = run_cdftt(input_path);

    std::cout << "STDOUT:\n" << r.stdout_text << std::endl;
    std::cout << "STDERR:\n" << r.stderr_text << std::endl;

    // ! `OrbitalType=Custom` crashes without any stdout or stderr output, returning -1 exit code.

    assert_zero_exit(r);
    assert_stdout_contains(r, "Reading data from " + analytic + "... Please wait.");
    assert_stdout_contains(r, "Writing cube file, please wait...");
    assert_stdout_contains(r, "Density cube saved to " + output);
    assert_file_exists(output);
    assert_cube_header_parseable(output);
    assert_cube_has_finite_values(output);

    // Clean up
    std::remove(input_path.c_str());
    std::remove(output.c_str());
}

// [P1 | nightly] Cover SpinType mode: Alpha which only works with OrbitalType=All
static void test_make_orbitals_cube_spin_type_alpha() {
    const std::string root = repo_root();
    const std::string analytic = root + "/examples/Orbitals/h2o.molden";
    const std::string output = "/tmp/cdftt_test_make_orbitals_cube_spin_type_alpha.cube";
    const std::string input_path = "/tmp/cdftt_test_make_orbitals_cube_spin_type_alpha.txt";

    write_input_file(input_path,
                     std::string("RunType=MakeOrbitalsCube\n")
                   + "AnalyticFiles=" + analytic + "\n"
                   + "Size=Medium\n"
                   + "OrbitalType=All\n"
                   + "SpinType=Alpha\n"
                   + "GridFilesNames=" + output + "\n");

    const RunResult r = run_cdftt(input_path);

    assert_zero_exit(r);
    assert_stdout_contains(r, "Reading data from " + analytic + "... Please wait.");
    assert_stdout_contains(r, "Writing cube file, please wait...");
    assert_stdout_contains(r, "Density cube saved to " + output);
    assert_file_exists(output);
    assert_cube_header_parseable(output);
    assert_cube_has_finite_values(output);

    // Clean up
    std::remove(input_path.c_str());
    std::remove(output.c_str());
}

// [P1 | nightly] Cover SpinType mode: Beta which only works with OrbitalType=All
static void test_make_orbitals_cube_spin_type_beta() {
    const std::string root = repo_root();
    const std::string analytic = root + "/examples/Orbitals/h2o.molden";
    const std::string output = "/tmp/cdftt_test_make_orbitals_cube_spin_type_beta.cube";
    const std::string input_path = "/tmp/cdftt_test_make_orbitals_cube_spin_type_beta.txt";

    write_input_file(input_path,
                     std::string("RunType=MakeOrbitalsCube\n")
                   + "AnalyticFiles=" + analytic + "\n"
                   + "Size=Medium\n"
                   + "OrbitalType=All\n"
                   + "SpinType=Beta\n"
                   + "GridFilesNames=" + output + "\n");
    
    const RunResult r = run_cdftt(input_path);

    assert_zero_exit(r);
    assert_stdout_contains(r, "Reading data from " + analytic + "... Please wait.");
    assert_stdout_contains(r, "Writing cube file, please wait...");
    assert_stdout_contains(r, "Density cube saved to " + output);
    assert_file_exists(output);
    assert_cube_header_parseable(output);
    assert_cube_has_finite_values(output);

    // Clean up
    std::remove(input_path.c_str());
    std::remove(output.c_str());
}

// [P1 | nightly] Cover SpinType mode: Alpha-Beta which only works with OrbitalType=All
static void test_make_orbitals_cube_spin_type_alpha_beta() {
    const std::string root = repo_root();
    const std::string analytic = root + "/examples/Orbitals/h2o.molden";
    const std::string output = "/tmp/cdftt_test_make_orbitals_cube_spin_type_alpha_beta.cube";
    const std::string input_path = "/tmp/cdftt_test_make_orbitals_cube_spin_type_alpha_beta.txt";

    write_input_file(input_path,
                     std::string("RunType=MakeOrbitalsCube\n")
                   + "AnalyticFiles=" + analytic + "\n"
                   + "Size=Medium\n"
                   + "OrbitalType=All\n"
                   + "SpinType=Alpha-Beta\n"
                   + "GridFilesNames=" + output + "\n");
    
    const RunResult r = run_cdftt(input_path);

    assert_zero_exit(r);
    assert_stdout_contains(r, "Reading data from " + analytic + "... Please wait.");
    assert_stdout_contains(r, "Writing cube file, please wait...");
    assert_stdout_contains(r, "Density cube saved to " + output);
    assert_file_exists(output);
    assert_cube_header_parseable(output);
    assert_cube_has_finite_values(output);

    // Clean up
    std::remove(input_path.c_str());
    std::remove(output.c_str());
}

// [P1 | nightly] Cover SpinList mode: Alpha which only works with OrbitalType=Custom
static void test_make_orbitals_cube_spin_list_alpha() {
    const std::string root = repo_root();
    const std::string analytic = root + "/examples/Orbitals/h2o.molden";
    const std::string output = "/tmp/cdftt_test_make_orbitals_cube_spin_list_alpha.cube";
    const std::string input_path = "/tmp/cdftt_test_make_orbitals_cube_spin_list_alpha.txt";

    write_input_file(input_path,
                     std::string("RunType=MakeOrbitalsCube\n")
                   + "AnalyticFiles=" + analytic + "\n"
                   + "Size=Medium\n"
                   + "OrbitalType=Custom\n"
                   + "OrbitalsNumbers=1,2,3\n"
                   + "SpinList=Alpha,Alpha,Alpha\n"
                   + "GridFilesNames=" + output + "\n");

    const RunResult r = run_cdftt(input_path);

    // ! `OrbitalType=Custom` crashes without any stdout or stderr output, returning -1 exit code.

    assert_zero_exit(r);
    assert_stdout_contains(r, "Reading data from " + analytic + "... Please wait.");
    assert_stdout_contains(r, "Writing cube file, please wait...");
    assert_stdout_contains(r, "Density cube saved to " + output);
    assert_file_exists(output);
    assert_cube_header_parseable(output);
    assert_cube_has_finite_values(output);

    // Clean up
    std::remove(input_path.c_str());
    std::remove(output.c_str());
}

// [P1 | nightly] Cover SpinList mode: Beta which only works with OrbitalType=Custom
static void test_make_orbitals_cube_spin_list_beta() {
    const std::string root = repo_root();
    const std::string analytic = root + "/examples/Orbitals/h2o.molden";
    const std::string output = "/tmp/cdftt_test_make_orbitals_cube_spin_list_beta.cube";
    const std::string input_path = "/tmp/cdftt_test_make_orbitals_cube_spin_list_beta.txt";

    write_input_file(input_path,
                     std::string("RunType=MakeOrbitalsCube\n")
                   + "AnalyticFiles=" + analytic + "\n"
                   + "Size=Medium\n"
                   + "OrbitalType=Custom\n"
                   + "OrbitalsNumbers=1,2,3\n"
                   + "SpinList=Beta,Beta,Beta\n"
                   + "GridFilesNames=" + output + "\n");

    const RunResult r = run_cdftt(input_path);

    // ! `OrbitalType=Custom` crashes without any stdout or stderr output, returning -1 exit code.

    assert_zero_exit(r);
    assert_stdout_contains(r, "Reading data from " + analytic + "... Please wait.");
    assert_stdout_contains(r, "Writing cube file, please wait...");
    assert_stdout_contains(r, "Density cube saved to " + output);
    assert_file_exists(output);
    assert_cube_header_parseable(output);
    assert_cube_has_finite_values(output);

    // Clean up
    std::remove(input_path.c_str());
    std::remove(output.c_str());
}

// [P1 | nightly] Cover SpinList mode: mixed which only works with OrbitalType=Custom
static void test_make_orbitals_cube_spin_list_mixed() {
    const std::string root = repo_root();
    const std::string analytic = root + "/examples/Orbitals/h2o.molden";
    const std::string output = "/tmp/cdftt_test_make_orbitals_cube_spin_list_mixed.cube";
    const std::string input_path = "/tmp/cdftt_test_make_orbitals_cube_spin_list_mixed.txt";

    write_input_file(input_path,
                     std::string("RunType=MakeOrbitalsCube\n")
                   + "AnalyticFiles=" + analytic + "\n"
                   + "Size=Medium\n"
                   + "OrbitalType=Custom\n"
                   + "OrbitalsNumbers=1,2,3\n"
                   + "SpinList=Alpha,Beta,Alpha\n"
                   + "GridFilesNames=" + output + "\n");

    const RunResult r = run_cdftt(input_path);

    // ! `OrbitalType=Custom` crashes without any stdout or stderr output, returning -1 exit code.

    assert_zero_exit(r);
    assert_stdout_contains(r, "Reading data from " + analytic + "... Please wait.");
    assert_stdout_contains(r, "Writing cube file, please wait...");
    assert_stdout_contains(r, "Density cube saved to " + output);
    assert_file_exists(output);
    assert_cube_header_parseable(output);
    assert_cube_has_finite_values(output);

    // Clean up
    std::remove(input_path.c_str());
    std::remove(output.c_str());
}

// [P1 | nightly] Cover SpinList shorter than OrbitalsNumbers
static void test_make_orbitals_cube_spin_list_shorter_than_orbitals_numbers() {
    // The expected behavior should be to raise an explicit error and produce non-zero exit code
    const std::string root = repo_root();
    const std::string analytic = root + "/examples/Orbitals/h2o.molden";
    const std::string output = "/tmp/cdftt_test_make_orbitals_cube_spin_list_shorter_than_orbitals_numbers.cube";
    const std::string input_path = "/tmp/cdftt_test_make_orbitals_cube_spin_list_shorter_than_orbitals_numbers.txt";

    write_input_file(input_path,
                     std::string("RunType=MakeOrbitalsCube\n")
                   + "AnalyticFiles=" + analytic + "\n"
                   + "Size=Medium\n"
                   + "OrbitalType=Custom\n"
                   + "OrbitalsNumbers=1,2,3,4,5\n"
                   + "SpinList=Alpha,Beta\n" // shorter than OrbitalsNumbers
                   + "GridFilesNames=" + output + "\n");
    
    const RunResult r = run_cdftt(input_path);

    assert_nonzero_exit(r);
    assert_stderr_contains(r, "Error: sizes of orbitals numbers list and spins list do not match");

    // Clean up
    std::remove(input_path.c_str());
    std::remove(output.c_str());  //? Should not be created ?
}

int main() {
    const Test tests[] = {
        { "test_make_orbitals_cube_smoke",             test_make_orbitals_cube_smoke             },
        { "test_make_orbitals_cube_invalid_selection", test_make_orbitals_cube_invalid_selection },
        { "test_make_orbitals_cube_orbital_type_all", test_make_orbitals_cube_orbital_type_all },
        { "test_make_orbitals_cube_orbital_type_occupied", test_make_orbitals_cube_orbital_type_occupied },
        { "test_make_orbitals_cube_orbital_type_virtual", test_make_orbitals_cube_orbital_type_virtual },
        { "test_make_orbitals_cube_orbital_type_homo", test_make_orbitals_cube_orbital_type_homo },
        { "test_make_orbitals_cube_orbital_type_lumo", test_make_orbitals_cube_orbital_type_lumo },
        { "test_make_orbitals_cube_orbital_type_homo_lumo", test_make_orbitals_cube_orbital_type_homo_lumo },
        { "test_make_orbitals_cube_orbital_type_custom", test_make_orbitals_cube_orbital_type_custom },
        { "test_make_orbitals_cube_spin_type_alpha", test_make_orbitals_cube_spin_type_alpha },
        { "test_make_orbitals_cube_spin_type_beta", test_make_orbitals_cube_spin_type_beta },
        { "test_make_orbitals_cube_spin_type_alpha_beta", test_make_orbitals_cube_spin_type_alpha_beta },
        { "test_make_orbitals_cube_spin_list_alpha", test_make_orbitals_cube_spin_list_alpha },
        { "test_make_orbitals_cube_spin_list_beta", test_make_orbitals_cube_spin_list_beta },
        {"test_make_orbitals_cube_spin_list_mixed", test_make_orbitals_cube_spin_list_mixed},
        { "test_make_orbitals_cube_spin_list_shorter_than_orbitals_numbers", test_make_orbitals_cube_spin_list_shorter_than_orbitals_numbers }
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
