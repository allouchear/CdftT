#!/usr/bin/env bash

# A simple test runner for C++ tests in the `tests/jobs` directory.
# Each test is a standalone C++ file with a main() function. This script can
# build and run them, passing the CDFTT_BINARY path via environment variable.

set -euo pipefail  # Strict mode: exit on error, undefined variable, or failed pipe. 

# Useful environment variables
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
TESTS_DIR="$REPO_DIR/tests"
JOBS_DIR="$TESTS_DIR/jobs"
BUILD_DIR="$REPO_DIR/build/tests"
HELPERS_CPP="$TESTS_DIR/helpers.cpp"

# Defaults; can be overridden by env var.
CDFTT_BINARY_DEFAULT="$REPO_DIR/src/applications/cdftt/cdftt"
CDFTT_BINARY="${CDFTT_BINARY:-$CDFTT_BINARY_DEFAULT}"

# Compiler settings; can be overridden by env var.
CXX="${CXX:-g++}"
CXXFLAGS=( -std=c++17 -O -I "$REPO_DIR/src/lib" )
LDFLAGS=( -lpthread )

# Print usage information.
usage() {
    echo "Usage: $0 [list|build|run] [test_name|all]"
    echo
    echo "Commands:"
    echo "  list                 List runnable tests (files containing int main())."
    echo "  build [name|all]     Build one test or all runnable tests."
    echo "  run [name|all]       Build then run one test or all runnable tests."
    echo
    echo "Examples:"
    echo "  $0 list"
    echo "  $0 build all"
    echo "  $0 run all"
    echo "  $0 run test_help"
    echo
    echo "Env vars:"
    echo "  CDFTT_BINARY=/absolute/path/to/cdftt"
    echo "  CXX=clang++"
}

# Check that the CDFTT_BINARY is set to an executable file before running tests.
ensure_cdftt_binary() {
    if [[ ! -f "$CDFTT_BINARY" || ! -x "$CDFTT_BINARY" ]]; then
        echo "Error: cdftt binary not found or not executable: $CDFTT_BINARY" >&2
        echo "Build cdftt first, or set CDFTT_BINARY to a valid executable path." >&2
        exit 2
    fi
}

# Return a list of runnable test names (without .cpp extension) in the JOBS_DIR.
# A test is considered runnable if it has a main() function.
runnable_tests() {
    shopt -s nullglob
    local file
    for file in "$JOBS_DIR"/*.cpp; do
        if grep -Eq '^[[:space:]]*int[[:space:]]+main[[:space:]]*\(' "$file"; then
            basename "$file" .cpp
        fi
    done
}

# Build a single test by name (without .cpp extension).
build_one() {
    local test_name="$1"
    local src="$JOBS_DIR/$test_name.cpp"
    local out="$BUILD_DIR/$test_name"

    if [[ ! -f "$src" ]]; then
        echo "Error: test source not found: $src" >&2
        return 1
    fi

    mkdir -p "$BUILD_DIR"
    echo "[BUILD] $test_name"
    "$CXX" "${CXXFLAGS[@]}" "$HELPERS_CPP" "$src" -o "$out" "${LDFLAGS[@]}"
}

# Build multiple tests by name.
build_many() {
    local names=("$@")
    local name
    for name in "${names[@]}"; do
        build_one "$name"
    done
}

# Run a single test by name (expects it to be built already).
run_one() {
    local test_name="$1"
    local exe="$BUILD_DIR/$test_name"

    if [[ ! -x "$exe" ]]; then
        echo "Error: executable not found: $exe" >&2
        echo "Run '$0 build $test_name' first, or use '$0 run $test_name'." >&2
        return 1
    fi

    echo "[RUN  ] $test_name"
    CDFTT_BINARY="$CDFTT_BINARY" "$exe"
}

# Run multiple tests by name.
run_many() {
    local names=("$@")
    local name
    for name in "${names[@]}"; do
        run_one "$name"
    done
}

# Parse test results from the last line of output
parse_test_results() {
    local output="$1"
    local passed=0
    local failed=0
    
    # Extract numbers from the last line of output
    if echo "$output" | grep -qE '[0-9]+ passed, [0-9]+ failed'; then
        passed=$(echo "$output" | tail -n1 | grep -oE '[0-9]+ passed' | cut -d' ' -f1)
        failed=$(echo "$output" | tail -n1 | grep -oE '[0-9]+ failed' | cut -d' ' -f1)
    fi
    
    echo "$passed $failed"
}

main() {
    local cmd="${1:-run}"  # Either "list", "build", or (default) "run"
    local target="${2:-all}" # Either (default) "all" or a specific test name

    # Discover runnable tests first, so "list" command can work without building.
    # Detail of how it works:
    #   - mapfile allows us to read the output of runnable_tests into an array variable.
    #   - The < <(...) syntax is process substitution, which runs the command and provides its output as a file-like input to mapfile.
    #   - discovered will be an array of test names (without .cpp) that have a main() function and are thus runnable.
    mapfile -t discovered < <(runnable_tests)

    # Handle the "list" command immediately, since it doesn't require building or running.
    if [[ "$cmd" == "list" ]]; then
        # If no runnable tests are found, print a message and exit with code 0 (not an error).
        if [[ ${#discovered[@]} -eq 0 ]]; then
            echo "No runnable tests found in $JOBS_DIR"
            exit 0
        fi
        # Just print the list of test names and exit.
        printf '%s\n' "${discovered[@]}"
        exit 0
    fi

    # If no runnable tests are found, print a message and exit with code 1 (error).
    if [[ ${#discovered[@]} -eq 0 ]]; then
        echo "No runnable tests found in $JOBS_DIR" >&2
        exit 1
    fi

    # Which tests are actually needed ?
    local selected=()
    if [[ "$target" == "all" ]]; then
        selected=("${discovered[@]}")
    else
        selected=("$target")
    fi

    # Validate that the selected tests are actually in the discovered list.
    case "$cmd" in
        build)
            build_many "${selected[@]}"
            ;;
        run)
            ensure_cdftt_binary
            build_many "${selected[@]}"
            
            # Track total results
            local total_passed=0
            local total_failed=0
            
            # Run each test and collect results
            for test_name in "${selected[@]}"; do
                local exe="$BUILD_DIR/$test_name"
                if [[ -x "$exe" ]]; then
                    echo "[RUN  ] $test_name"
                    local output
                    local status=0
                    if output=$(CDFTT_BINARY="$CDFTT_BINARY" "$exe" 2>&1); then
                        status=0
                    else
                        status=$?
                    fi

                    # Always display test output, even when the test binary failed.
                    echo "$output"

                    local results
                    results=$(parse_test_results "$output")
                    local passed=$(echo "$results" | cut -d' ' -f1)
                    local failed=$(echo "$results" | cut -d' ' -f2)

                    # If the binary failed before printing a summary, count it as one failed test binary.
                    if [[ "$status" -ne 0 && "$passed" -eq 0 && "$failed" -eq 0 ]]; then
                        failed=1
                    fi
                    
                    total_passed=$((total_passed + passed))
                    total_failed=$((total_failed + failed))
                else
                    echo "Error: executable not found: $exe" >&2
                    exit 1
                fi
            done
            
            # Display summary
            echo "=== SUMMARY ==="
            echo "Total: $total_passed passed, $total_failed failed"
            [[ "$total_failed" -eq 0 ]]
            ;;
        *)
            usage
            exit 1
            ;;
    esac
}

main "$@"