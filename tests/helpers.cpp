#include "helpers.hpp"

#include <cmath>
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>
#include <sys/wait.h>
#include <thread>
#include <unistd.h>
#include <regex>

namespace cdftt_tests {

// Resolve the cdftt binary path.
// Priority to check for it: 
//        1)  CDFTT_BINARY env var (explicit path) 
//   else 2) > LIBCDFTTSRC env var
//   else 3) PATH fallback.
static std::string find_cdftt_binary() {
    const char* bin = std::getenv("CDFTT_BINARY");
    if (bin && bin[0] != '\0')
        return bin;
    const char* src = std::getenv("LIBCDFTTSRC");
    if (src && src[0] != '\0')
        return std::string(src) + "/cdftt/cdftt";
    return "cdftt";
}

// Converts a live output stream (from a live process: "file descriptor") into a string until EOF.
static std::string drain_fd(int fd) {
    std::string out;
    char buf[4096];  // Arbitrary buffer size for reading; the loop will continue until EOF.
    ssize_t n;
    while ((n = ::read(fd, buf, sizeof(buf))) > 0)
        out.append(buf, static_cast<std::size_t>(n));
    return out;
}

// 1) Spawn cdftt with input_file_path
// 2) wait for termination
// 3) return
//        exit code 
//      + captured stdout/stderr
RunResult run_cdftt(const std::string& input_file_path) {

    // Path to the cdftt binary to execute
    const std::string binary = find_cdftt_binary();

    // Create pipes for capturing stdout and stderr from the child process (the cdftt process).
    int out_pipe[2], err_pipe[2];
    if (::pipe(out_pipe) != 0 || ::pipe(err_pipe) != 0)
        throw std::runtime_error(std::string("pipe() failed: ") + std::strerror(errno));

    // Fork the process to run cdftt in the child.
    // This allows us to capture stdout/stderr via the pipes and get the exit code on completion.
    // Using popen instead would be easier to read, but it would only capture stdout
    // and not stderr, which is important for diagnosing test failures.
    const pid_t pid = ::fork();
    if (pid < 0)
        throw std::runtime_error(std::string("fork() failed: ") + std::strerror(errno));

    if (pid == 0) {
        // ! We are now inside the child process.
        // The child's only job: set up the pipes, then become cdftt.

        // Each pipe has two ends: [0]=read, [1]=write.
        // The child will only write, so close the "read" ends as they are unused here.
        ::close(out_pipe[0]);
        ::close(err_pipe[0]);

        // Redirect the child's stdout and stderr into the "write" ends of the pipes.
        // After this, anything cdftt prints goes into the pipes instead of the terminal.
        ::dup2(out_pipe[1], STDOUT_FILENO);
        ::dup2(err_pipe[1], STDERR_FILENO);

        // The original write-end descriptors are now redundant (dup2 made copies).
        ::close(out_pipe[1]);
        ::close(err_pipe[1]);

        // The expected argv is: { "/path/to/cdftt", "/path/to/input.txt", NULL }
        // Where the "NULL" signals the end of the argument list.
        const char* argv[] = { binary.c_str(), input_file_path.c_str(), nullptr };

        // ! Execute the cdftt binary. If successful, this does not return; the child process becomes cdftt.
        ::execvp(binary.c_str(), const_cast<char* const*>(argv));
        // if everything went well, we never reach this point because execvp replaces the process image with cdftt.

        // execvp only returns if something went wrong and it failed (binary not found, not executable, etc.).
        // Report the failure to stderr so tests see a concrete cause, then exit immediately.
        const int exec_errno = errno;
        std::string message = "execvp() failed for \"" + binary + "\": " + std::strerror(exec_errno) + "\n";
        (void) ::write(STDERR_FILENO, message.c_str(), message.size());
        ::_exit(exec_errno == ENOENT ? 127 : 126);
    }

    // ! Only the parent will arrive here

    // Close write ends as the parent only needs to read
    // drain both pipes concurrently to avoid blocking (can apparently occur when one pipe buffer fills
    // while the other is unread)
    ::close(out_pipe[1]);
    ::close(err_pipe[1]);

    // Store the results in a custom struct to return after the child process finishes.
    RunResult result;
    // Thread to read stderr concurrently 
    std::thread err_thread([&]() {
        result.stderr_text = drain_fd(err_pipe[0]);  // stops when the child process closes the write end (on exit)
        ::close(err_pipe[0]);
    });
    // Read stdout in the main thread
    result.stdout_text = drain_fd(out_pipe[0]);  // stops when the child process closes the write end (on exit)
    ::close(out_pipe[0]);

    // Merges the stderr thread back into the main thread to ensure all output is captured before we return.
    err_thread.join();

    // Wait for the child process to finish and capture its exit code
    int status = 0;
    ::waitpid(pid, &status, 0);
    result.exit_code = WIFEXITED(status)     // Did it exit normally (not killed by a signal)?
                     ? WEXITSTATUS(status)   // If so, get the exit code (0 for success, non-zero for failure).
                     : -1;                   // Else (e.g., killed by signal), we return -1.

    // ! NB: Non-zero exits are not thrown here.
    // The caller will assert on the exit code as needed. We just return whatever we got.
    return result;
}

// Write a cdftt input file to disk at the given path.
// Each test calls this to create its own input before calling run_cdftt
// so tests don't depend on pre-existing files and can control every parameter.
// If the file already exists, it is overwritten from scratch.
void write_input_file(const std::string& path, const std::string& content) {
    std::ofstream f(path, std::ios::out | std::ios::trunc);
    if (!f.is_open())
        throw std::runtime_error("Cannot open file for writing: " + path);
    f << content;
    if (f.fail())
        throw std::runtime_error("Write failed for file: " + path);
}

// Fail the test with a clear message when path does not denote an existing
// regular file. Used after grid-generating jobs to confirm output was written.
void assert_file_exists(const std::string& path) {
    struct stat st{};
    if (::stat(path.c_str(), &st) != 0 || !S_ISREG(st.st_mode))
        throw std::runtime_error("Expected file does not exist: " + path);
}

static void assert_contains(const std::string& haystack,
                            const std::string& fragment,
                            const std::string& stream_name) {
    if (haystack.find(fragment) == std::string::npos)
        throw std::runtime_error(stream_name + ": expected to find \"" + fragment
                                 + "\" but it was absent.\nFull text:\n" + haystack);
}

void assert_zero_exit(const RunResult& result) {
    if (result.exit_code != 0)
        // throw std::runtime_error("Expected exit code 0 but got " + std::to_string(result.exit_code)
        //                          + "\nstdout:\n" + result.stdout_text
        //                          + "\nstderr:\n" + result.stderr_text);
        throw std::runtime_error("Expected exit code 0 but got " + std::to_string(result.exit_code)
                                 + "\nstderr:\n" + result.stderr_text);
}

void assert_nonzero_exit(const RunResult& result) {
    if (result.exit_code == 0)
        // throw std::runtime_error("Expected non-zero exit code but got 0.\nstdout:\n"
        //                          + result.stdout_text + "\nstderr:\n" + result.stderr_text);
        throw std::runtime_error("Expected non-zero exit code but got 0.\nstderr:\n"
                                 + result.stderr_text);
}

void assert_stdout_contains(const RunResult& result, const std::string& fragment) {
    assert_contains(result.stdout_text, fragment, "stdout");
    if (result.stdout_text.find(fragment) == std::string::npos)
        throw std::runtime_error("stdout: expected to find \"" + fragment
                                 + "\" but it was absent.");
}

void assert_stdout_does_not_contain(const RunResult& result, const std::string& fragment) {
    if (result.stdout_text.find(fragment) != std::string::npos)
        // throw std::runtime_error("stdout: expected not to find \"" + fragment
        //                          + "\" but it was present.\nFull text:\n" + result.stdout_text);
        throw std::runtime_error("stdout: expected not to find \"" + fragment
                                 + "\" but it was present.");
}

void assert_stderr_contains(const RunResult& result, const std::string& fragment) {
    assert_contains(result.stderr_text, fragment, "stderr");
}

struct CubeHeader {
    int atom_lines = 0;
    int nx = 0;
    int ny = 0;
    int nz = 0;
};

static CubeHeader parse_cube_header(std::ifstream& f, const std::string& path) {
    std::string line;

    // Two title/comment lines.
    if (!std::getline(f, line) || !std::getline(f, line))
        throw std::runtime_error("Cube parse error: missing title/comment lines in " + path);

    // Line 3: number of atoms + origin.
    if (!std::getline(f, line))
        throw std::runtime_error("Cube parse error: missing atom/origin line in " + path);
    std::istringstream atom_ss(line);
    int atom_count = 0;
    double ox = 0.0, oy = 0.0, oz = 0.0;
    if (!(atom_ss >> atom_count >> ox >> oy >> oz))
        throw std::runtime_error("Cube parse error: invalid atom/origin line in " + path);

    CubeHeader h;
    h.atom_lines = std::abs(atom_count);

    // Lines 4-6: voxel counts + axis vectors.
    auto parse_axis_line = [&](int& n_voxels, const char* axis_name) {
        if (!std::getline(f, line))
            throw std::runtime_error(std::string("Cube parse error: missing ") + axis_name + " axis line in " + path);
        std::istringstream axis_ss(line);
        double ax = 0.0, ay = 0.0, az = 0.0;
        if (!(axis_ss >> n_voxels >> ax >> ay >> az))
            throw std::runtime_error(std::string("Cube parse error: invalid ") + axis_name + " axis line in " + path);
        if (n_voxels <= 0)
            throw std::runtime_error(std::string("Cube parse error: non-positive ") + axis_name + " voxel count in " + path);
    };

    parse_axis_line(h.nx, "X");
    parse_axis_line(h.ny, "Y");
    parse_axis_line(h.nz, "Z");

    // Skip atom specification lines.
    for (int i = 0; i < h.atom_lines; ++i) {
        if (!std::getline(f, line))
            throw std::runtime_error("Cube parse error: truncated atom specification block in " + path);
    }

    return h;
}

void assert_cube_header_parseable(const std::string& path) {
    assert_file_exists(path);

    std::ifstream f(path);
    if (!f.is_open())
        throw std::runtime_error("Cube parse error: cannot open file " + path);

    (void) parse_cube_header(f, path);
}

void assert_cube_has_finite_values(const std::string& path) {
    assert_file_exists(path);

    std::ifstream f(path);
    if (!f.is_open())
        throw std::runtime_error("Cube parse error: cannot open file " + path);

    const CubeHeader h = parse_cube_header(f, path);
    const long long expected_values = static_cast<long long>(h.nx) * h.ny * h.nz;

    std::string token;
    long long numeric_values = 0;
    while (f >> token) {
        char* end = nullptr;
        const double v = std::strtod(token.c_str(), &end);
        if (end == token.c_str() || *end != '\0')
            throw std::runtime_error("Cube parse error: non-numeric voxel token \"" + token + "\" in " + path);
        if (!std::isfinite(v))
            throw std::runtime_error("Cube value error: non-finite voxel value in " + path);
        ++numeric_values;
    }

    if (numeric_values < expected_values) {
        throw std::runtime_error("Cube parse error: insufficient voxel values in " + path
                                 + " (expected at least " + std::to_string(expected_values)
                                 + ", got " + std::to_string(numeric_values) + ")");
    }
}


std::string repo_root() {
    const char* from_env = std::getenv("CDFTT_REPO_ROOT");
    if (from_env && from_env[0] != '\0')
        return from_env;

    return ".";
}


void assert_file_size_less_than(const std::string& path, std::size_t max_bytes) {
    std::ifstream f(path, std::ios::binary | std::ios::ate);
    if (!f.is_open())
        throw std::runtime_error("File size error: cannot open file " + path);
    const std::size_t actual_size = static_cast<std::size_t>(f.tellg());
    if (actual_size >= max_bytes)
        throw std::runtime_error("File size error: file " + path + " is not less than " + std::to_string(max_bytes) + " bytes");
}

int count_lines_matching_regex(const std::string& text, const std::string& regex_pattern) {
    std::regex re(regex_pattern);
    std::istringstream ss(text);
    std::string line;
    int matches = 0;
    while (std::getline(ss, line)) {
        if (std::regex_search(line, re))
            ++matches;
    }
    return matches;
}

double extract_double_after_label(const std::string& text, const std::string& label) {
    const std::size_t pos = text.find(label);
    if (pos == std::string::npos)
        throw std::runtime_error("Label not found: " + label);

    // Find the substring after the label and attempt to parse the first double.
    const std::string after = text.substr(pos + label.size());
    std::istringstream ss(after);
    double v;
    if (!(ss >> v))
        throw std::runtime_error("Failed to parse double after label: " + label);
    return v;
}

std::vector<double> read_first_n_voxels(const std::string& path, std::size_t n) {
    std::ifstream f(path);
    if (!f.is_open())
        throw std::runtime_error("Cube parse error: cannot open file " + path);

    std::string line;
    // Skip title/comment two lines
    if (!std::getline(f, line) || !std::getline(f, line))
        throw std::runtime_error("Cube parse error: missing title/comment lines in " + path);

    // atom/origin line
    if (!std::getline(f, line))
        throw std::runtime_error("Cube parse error: missing atom/origin line in " + path);
    std::istringstream atom_ss(line);
    int atom_count = 0;
    double ox, oy, oz;
    if (!(atom_ss >> atom_count >> ox >> oy >> oz))
        throw std::runtime_error("Cube parse error: invalid atom/origin line in " + path);
    const int atom_lines = std::abs(atom_count);

    // skip three axis lines
    for (int i = 0; i < 3; ++i) {
        if (!std::getline(f, line))
            throw std::runtime_error("Cube parse error: missing axis line in " + path);
    }

    // skip atom spec lines
    for (int i = 0; i < atom_lines; ++i) {
        if (!std::getline(f, line))
            throw std::runtime_error("Cube parse error: truncated atom specification block in " + path);
    }

    std::vector<double> vals;
    vals.reserve(n);
    std::string token;
    while (vals.size() < n && (f >> token)) {
        char* end = nullptr;
        const double v = std::strtod(token.c_str(), &end);
        if (end == token.c_str() || *end != '\0')
            throw std::runtime_error("Cube parse error: non-numeric voxel token \"" + token + "\" in " + path);
        vals.push_back(v);
    }

    return vals;
}

} // namespace cdftt_tests
