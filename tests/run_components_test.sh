#!/usr/bin/env bash

# A test runner for C++ components tests in the `tests/components` directory.
# This script can build and run the components test file, gathering outputs.

set -euo pipefail  # Strict mode: exit on error, undefined variable, or failed pipe. 

# Useful environment variables
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
TESTS_DIR="$REPO_DIR/tests"
COMPONENTS_DIR="$TESTS_DIR/components"
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
    echo "Usage: $0 [build|run] [test_name|all]"
    echo
    echo "Commands:"
    echo "  build                Build the components test."
    echo "  run                  Build and run the components test."
    echo "  run-only             Run the already built components test."
    echo
    echo "Examples:"
    echo "  $0 build"
    echo "  $0 run"
    echo "  $0 run-only"
    echo
    echo "Env vars:"
    echo "  CDFTT_BINARY=/absolute/path/to/cdftt"
    echo "  CXX=clang++"
}

# Check that the CDFTT_BINARY is set to an executable file before running tests.
ensure_cdftt_binary() {
    if [[ ! -x "$CDFTT_BINARY" ]]; then
        echo "Error: cdftt binary not found or not executable: $CDFTT_BINARY" >&2
        echo "Build cdftt first, or set CDFTT_BINARY to a valid executable path." >&2
        exit 2
    fi
}

# Build the components test
build_components_test() {
    local src="$COMPONENTS_DIR/test_job_components.cpp"
    local out="$BUILD_DIR/test_job_components"
    
    if [[ ! -f "$src" ]]; then
        echo "Error: components test source not found: $src" >&2
        return 1
    fi

    mkdir -p "$BUILD_DIR"
    echo "[BUILD] test_job_components"
    "$CXX" "${CXXFLAGS[@]}" "$HELPERS_CPP" "$src" -o "$out" "${LDFLAGS[@]}"
}

# Run the components test
run_components_test() {
    local exe="$BUILD_DIR/test_job_components"
    
    if [[ ! -x "$exe" ]]; then
        echo "Error: executable not found: $exe" >&2
        echo "Run '$0 build' first." >&2
        return 1
    fi

    echo "[RUN  ] test_job_components"
    echo "==========================="
    echo "Running components test..."
    echo ""
    
    # Capture both stdout and stderr
    CDFTT_BINARY="$CDFTT_BINARY" "$exe" 2>&1
    
    local exit_code=$?
    echo ""
    echo "==========================="
    echo "Exit code: $exit_code"
    
    return $exit_code
}

main() {
    local cmd="${1:-run}"  # Either "build", "run", or "run-only"
    
    case "$cmd" in
        build)
            build_components_test
            ;;
        run)
            ensure_cdftt_binary
            build_components_test
            run_components_test
            ;;
        run-only)
            run_components_test
            ;;
        *)
            usage
            exit 1
            ;;
    esac
}

main "$@"