#include "../helpers.hpp"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cdftt_tests;

/*
    About this job type (ConvertOrbitals):
        - It converts an orbital file from one supported format to another supported
            format and writes a new analytic file.

    Detailed header command expectations:
        - Mandatory header line:
                    RunType=ConvertOrbitals
        - Required companion line:
                    AnalyticFiles=<inputOrbitalFile>,<outputOrbitalFile>

    Expected behavior to anchor assertions:
        - Happy path:
                1. Exit code 0.
                2. stdout announces ConvertOrbitals job execution.
                3. Output file is created and contains non-zero content.
        - Negative path:
                1. Non-zero exit code.
                2. stderr reports invalid file format for unsupported extension.
*/

static void assert_file_nonempty(const std::string& path) {
    std::ifstream f(path, std::ios::binary | std::ios::ate);
    if (!f.is_open())
        throw std::runtime_error("Cannot open file to check non-empty content: " + path);
    if (f.tellg() <= 0)
        throw std::runtime_error("Expected non-empty file but size was zero: " + path);
}

// [P0 | PR-fast] Happy path: supported conversion from molden to wfx.
static void test_convert_orbitals_smoke() {
    const std::string root = repo_root();
    const std::string input_orbitals = root + "/examples/Orbitals/h2o.molden";
    const std::string output_orbitals = "/tmp/cdftt_test_convert_orbitals_output.wfx";
    const std::string input_path = "/tmp/cdftt_test_convert_orbitals_smoke.txt";

    write_input_file(input_path,
                     std::string("RunType=ConvertOrbitals\n")
                   + "AnalyticFiles=" + input_orbitals + "," + output_orbitals + "\n");

    const RunResult r = run_cdftt(input_path);

    assert_zero_exit(r);
    assert_stdout_contains(r, "has been converted");
    assert_file_exists(output_orbitals);
    assert_file_nonempty(output_orbitals);

    std::remove(input_path.c_str());
    std::remove(output_orbitals.c_str());
}

// [P0 | PR-fast] Negative path: unsupported input format fails clearly.
static void test_convert_orbitals_invalid_extension() {
    const std::string input_path = "/tmp/cdftt_test_convert_orbitals_invalid_extension.txt";

    write_input_file(input_path,
                     "RunType=ConvertOrbitals\n"
                     "AnalyticFiles=dummy.xyz,/tmp/cdftt_should_not_exist_convert_output.wfx\n");

    const RunResult r = run_cdftt(input_path);

    assert_nonzero_exit(r);
    assert_stderr_contains(r, "invalid file format for file \"dummy.xyz\"");

    std::remove(input_path.c_str());
}

// [P1 | nightly] Consolidated matrix-driven conversion test.
// Iterates over many format conversions, skipping missing fixtures, continues on failures,
// and reports a summary at the end so the test binary fails if any case failed.
static void test_convert_orbitals_matrix() {
    const std::string root = repo_root();

    // Matrix entries: input extension -> output extension
    // ? Note: As per documentation, for now {molden, wfx, fchk, log, gab} are all supported as input
    //?        Still as documented, only {molden, wfx, gab} are currently supported as output formats
    const std::vector<std::pair<std::string,std::string>> matrix = {
        // Molden <-> WFX
        {"molden","wfx"},  // Ok
        {"wfx","molden"},  //! "This option is nnot implemente." [UNDOCUMENTED]
        // Molden <-> FCHK
        // {"molden","fchk"},  //? "Format not recognize, please choose a valide format" [DOCUMENTED]
        {"fchk","molden"}, // Ok
        // Molden <-> Log
        // {"molden", "log"},  //? "Format not recognize, please choose a valide format" [DOCUMENTED]
        {"log","molden"},  // Ok
        // Molden <-> Gabedit
        {"molden","gab"},  //! "Gabedit Format can't read mixte basis." [UNDOCUMENTED]
        {"gab","molden"},  // Ok
        // WFX <-> FCHK
        // {"wfx","fchk"},     //? "Format not recognize, please choose a valide format" [DOCUMENTED]
        {"fchk","wfx"},    // Ok
        // WFX <-> Log
        // {"wfx", "log"},     //? "Format not recognize, please choose a valide format" [DOCUMENTED]
        {"log","wfx"},     // Ok
        // WFX <-> Gabedit
        {"wfx","gab"},     //! "This option is nnot implemente." [UNDOCUMENTED]
        {"gab","wfx"},     // Ok
        // FCHK <-> Log
        // {"fchk","log"},     //? "Format not recognize, please choose a valide format" [DOCUMENTED]
        // {"log","fchk"},   //? Format not recognize, please choose a valide format" [DOCUMENTED]
        // FCHK <-> Gabedit
        {"fchk","gab"},   // Ok
        // {"gab","fchk"},     //? "Format not recognize, please choose a valide format" [DOCUMENTED]
        // Log <-> Gabedit
        {"log","gab"},    // Ok
        // {"gab","log"},      //? "Format not recognize, please choose a valide format" [DOCUMENTED]
    };

    std::vector<std::string> failures;

    // Iterate over the defined format pairs
    for (const auto &p : matrix) {
        const std::string input_extension = p.first;
        const std::string output_extension = p.second;

        const std::string input_orbitals = root + "/tests/fixtures/Orbitals/h2o." + input_extension;
        std::ifstream inputOrbitalsStream(input_orbitals);
        if (!inputOrbitalsStream.good()) {
            std::cout << "[ SKIPPED ] Convert " << input_extension << " -> " << output_extension << " (fixture missing)" << std::endl;
            continue;
        }

        const std::string output_orbitals = "/tmp/cdftt_test_convert_orbitals_h2o_" + input_extension + "_to_" + output_extension + "." + output_extension;
        const std::string input_path = "/tmp/cdftt_test_convert_orbitals_matrix_" + input_extension + "_to_" + output_extension + ".txt";

        write_input_file(input_path,
                         std::string("RunType=ConvertOrbitals\n")
                       + "AnalyticFiles=" + input_orbitals + "," + output_orbitals + "\n");

        // Run the conversion and assert expected behavior
        try {
            const RunResult r = run_cdftt(input_path);
            // std::cout << "STDOUT:\n" << r.stdout_text << std::endl;
            assert_zero_exit(r);
            assert_stdout_contains(r, "has been converted");
            assert_file_exists(output_orbitals);
            assert_file_nonempty(output_orbitals);
            // Clean up for this iteration
            std::remove(output_orbitals.c_str());
            std::remove(input_path.c_str());
            std::cout << "[ OK ] Convert " << input_extension << " -> " << output_extension << std::endl;
        } catch (const std::exception &e) {
            // Explicit error reporting for any failing format pair without aborting the test
            std::string msg = "Convert " + input_extension + "->" + output_extension + " failed: " + e.what();
            std::cerr << "[FAIL] " << msg << std::endl;
            failures.push_back(msg);
            // Clean up input file -- no output to clean up on failure
            std::remove(input_path.c_str());
        }
    }

    // Aggregated reporting of any failures at the end of the matrix test
    if (!failures.empty()) {
        std::string summary = "Conversion matrix had failures:\n";
        for (const auto &s : failures) summary += " - " + s + "\n";
        throw std::runtime_error(summary);
    }
}

int main() {
    const Test tests[] = {
        { "test_convert_orbitals_smoke",             test_convert_orbitals_smoke             },
        { "test_convert_orbitals_invalid_extension", test_convert_orbitals_invalid_extension },
        { "test_convert_orbitals_matrix",            test_convert_orbitals_matrix            },
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
