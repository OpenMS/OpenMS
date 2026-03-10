// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/ExternalProcess.h>

#include <OpenMS/DATASTRUCTURES/String.h>

#include <algorithm>
#include <sstream>
#include <thread>
#include <utility>
#include <vector>
#include <string>

#ifdef _WIN32
  // Force Boost.Process to use std::filesystem instead of boost::filesystem to avoid linking Boost.Filesystem on Windows
  // This prevents unresolved externals like boost::filesystem::path_traits::convert(...)
  // when only std::filesystem is desired and available with MSVC.
  #ifndef BOOST_PROCESS_USE_STD_FS
  #define BOOST_PROCESS_USE_STD_FS 1
  #endif
  #include <boost/process.hpp>
  namespace bp = boost::process;
#else
  #include <unistd.h>
  #include <sys/wait.h>
  #include <sys/types.h>
  #include <sys/select.h>
  #include <fcntl.h>
  #include <signal.h>
#endif

namespace OpenMS
{

  /// default Ctor; callbacks for stdout/stderr are empty
  ExternalProcess::ExternalProcess()
    : ExternalProcess([&](const String& /*out*/) {}, [&](const String& /*out*/) {}) // call other Ctor to set callbacks!
  {
  }

  ExternalProcess::ExternalProcess(std::function<void(const String&)> callbackStdOut, std::function<void(const String&)> callbackStdErr)
    : callbackStdOut_(std::move(callbackStdOut)),
      callbackStdErr_(std::move(callbackStdErr))
  {
  }

  ExternalProcess::~ExternalProcess() = default;

  /// re-wire the callbacks used using run()
  void ExternalProcess::setCallbacks(std::function<void(const String&)> callbackStdOut, std::function<void(const String&)> callbackStdErr)
  {
    callbackStdOut_ = std::move(callbackStdOut);
    callbackStdErr_ = std::move(callbackStdErr);
  }

  static inline std::string join_cmd_(const std::string& exe, const std::vector<std::string>& args)
  {
    std::ostringstream oss;
    oss << exe;
    for (const auto& a : args)
    {
      oss << " " << a;
    }
    return oss.str();
  }

  ExternalProcess::RETURNSTATE ExternalProcess::run(const std::string& exe,
                                                    const std::vector<std::string>& args,
                                                    const std::string& working_dir,
                                                    const bool verbose,
                                                    IO_MODE io_mode,
                                                    const std::map<std::string, std::string>& env)
  {
    String error_msg;
    return run(exe, args, working_dir, verbose, error_msg, io_mode, env);
  }

  ExternalProcess::RETURNSTATE ExternalProcess::run(const std::string& exe,
                                                    const std::vector<std::string>& args,
                                                    const std::string& working_dir,
                                                    const bool verbose,
                                                    String& error_msg,
                                                    IO_MODE io_mode,
                                                    const std::map<std::string, std::string>& env)
  {
    error_msg.clear();

    if (verbose)
    {
      callbackStdOut_(String("Running: ") + String(join_cmd_(exe, args)) + String("\n"));
    }

#ifdef _WIN32
    // Windows implementation using boost::process
    bp::environment bp_env = boost::this_process::environment();
    // Add custom environment variables
    for (const auto& kv : env)
    {
      bp_env[kv.first] = kv.second;
    }
    const bool read_io = (io_mode == IO_MODE::READ_ONLY || io_mode == IO_MODE::READ_WRITE);

    bp::ipstream child_out;
    bp::ipstream child_err;
    std::unique_ptr<bp::child> child_ptr;

    try
    {
      if (working_dir.empty())
      {
        if (read_io)
        {
          child_ptr = std::make_unique<bp::child>(
            exe,
            bp::args(args),
            bp_env,
            bp::std_out > child_out,
            bp::std_err > child_err
          );
        }
        else
        {
          child_ptr = std::make_unique<bp::child>(
            exe,
            bp::args(args),
            bp_env,
            bp::std_out > bp::null,
            bp::std_err > bp::null
          );
        }
      }
      else
      {
        if (read_io)
        {
          child_ptr = std::make_unique<bp::child>(
            exe,
            bp::args(args),
            bp_env,
            bp::start_dir = working_dir,
            bp::std_out > child_out,
            bp::std_err > child_err
          );
        }
        else
        {
          child_ptr = std::make_unique<bp::child>(
            exe,
            bp::args(args),
            bp_env,
            bp::start_dir = working_dir,
            bp::std_out > bp::null,
            bp::std_err > bp::null
          );
        }
      }
    }
    catch (const bp::process_error& /*e*/)
    {
      // Failed to start process at all (equivalent to QProcess::FailedToStart)
      error_msg = "Process '" + String(exe) + "' failed to start. Does it exist? Is it executable?";
      if (verbose)
      {
        callbackStdErr_(error_msg + '\n');
      }
      return RETURNSTATE::FAILED_TO_START;
    }

    bp::child& c = *child_ptr;

    std::thread t_out;
    std::thread t_err;

    if (read_io)
    {
      t_out = std::thread([this, &child_out]()
      {
        std::string line;
        while (std::getline(child_out, line))
        {
          callbackStdOut_(String(line) + '\n');
        }
        if (child_out.good())
        {
          std::ostringstream oss;
          oss << child_out.rdbuf();
          auto rest = oss.str();
          if (!rest.empty()) callbackStdOut_(String(rest));
        }
      });

      t_err = std::thread([this, &child_err]()
      {
        std::string line;
        while (std::getline(child_err, line))
        {
          callbackStdErr_(String(line) + '\n');
        }
        if (child_err.good())
        {
          std::ostringstream oss;
          oss << child_err.rdbuf();
          auto rest = oss.str();
          if (!rest.empty()) callbackStdErr_(String(rest));
        }
      });
    }

    c.wait();

    if (t_out.joinable()) t_out.join();
    if (t_err.joinable()) t_err.join();

    int exit_code = c.exit_code();

    if (exit_code != 0)
    {
      // Map common 127 to FAILED_TO_START; otherwise NONZERO_EXIT
      if (exit_code == 127)
      {
        error_msg = "Process '" + String(exe) + "' failed to start (exit code 127). Does it exist? Is it executable?";
        if (verbose) { callbackStdErr_(error_msg + '\n'); }
        return RETURNSTATE::FAILED_TO_START;
      }
      error_msg = "Process '" + String(exe) + "' did not finish successfully (exit code: " + String(exit_code) + "). Please check the log.";
      if (verbose) { callbackStdErr_(error_msg + '\n'); }
      return RETURNSTATE::NONZERO_EXIT;
    }

    if (verbose) { callbackStdOut_("Executed '" + String(exe) + "' successfully!\n"); }
    return RETURNSTATE::SUCCESS;

#else
    // POSIX implementation using fork/exec and pipes (no boost link requirements)
    const bool want_read = (io_mode == IO_MODE::READ_ONLY || io_mode == IO_MODE::READ_WRITE);

    int pipe_out[2] = {-1, -1};
    int pipe_err[2] = {-1, -1};

    if (want_read)
    {
      if (pipe(pipe_out) == -1 || pipe(pipe_err) == -1)
      {
        error_msg = "Failed to create pipes for process communication";
        if (verbose) { callbackStdErr_(error_msg + '\n'); }
        // cleanup
        if (pipe_out[0] != -1) close(pipe_out[0]);
        if (pipe_out[1] != -1) close(pipe_out[1]);
        if (pipe_err[0] != -1) close(pipe_err[0]);
        if (pipe_err[1] != -1) close(pipe_err[1]);
        return RETURNSTATE::FAILED_TO_START;
      }
      // set non-blocking for better select/read behavior
      fcntl(pipe_out[0], F_SETFL, O_NONBLOCK);
      fcntl(pipe_err[0], F_SETFL, O_NONBLOCK);
    }

    pid_t pid = fork();
    if (pid == -1)
    {
      error_msg = "Failed to fork process";
      if (verbose) { callbackStdErr_(error_msg + '\n'); }
      if (want_read)
      {
        if (pipe_out[0] != -1) close(pipe_out[0]);
        if (pipe_out[1] != -1) close(pipe_out[1]);
        if (pipe_err[0] != -1) close(pipe_err[0]);
        if (pipe_err[1] != -1) close(pipe_err[1]);
      }
      return RETURNSTATE::FAILED_TO_START;
    }
    else if (pid == 0)
    {
      // Child
      if (want_read)
      {
        dup2(pipe_out[1], STDOUT_FILENO);
        dup2(pipe_err[1], STDERR_FILENO);
        close(pipe_out[0]); close(pipe_out[1]);
        close(pipe_err[0]); close(pipe_err[1]);
      }

      if (!working_dir.empty())
      {
        chdir(working_dir.c_str());
      }

      // Set custom environment variables in the child process
      for (const auto& kv : env)
      {
        setenv(kv.first.c_str(), kv.second.c_str(), 1);
      }

      std::vector<char*> argv_vec;
      argv_vec.push_back(const_cast<char*>(exe.c_str()));
      for (const auto& a : args) argv_vec.push_back(const_cast<char*>(a.c_str()));
      argv_vec.push_back(nullptr);

      execvp(exe.c_str(), argv_vec.data());
      _exit(127); // If we reach here, exec failed
    }

    // Parent
    if (want_read)
    {
      // close write ends
      close(pipe_out[1]);
      close(pipe_err[1]);

      bool child_finished = false;
      int status_from_loop = 0;

      fd_set fds;
      char buffer[4096];

      while (true)
      {
        FD_ZERO(&fds);
        FD_SET(pipe_out[0], &fds);
        FD_SET(pipe_err[0], &fds);

        struct timeval timeout;
        timeout.tv_sec = 0;
        timeout.tv_usec = 50000; // 50ms

        int max_fd = std::max(pipe_out[0], pipe_err[0]) + 1;
        int result = select(max_fd, &fds, nullptr, nullptr, &timeout);

        if (result > 0)
        {
          if (FD_ISSET(pipe_out[0], &fds))
          {
            ssize_t bytes_read = read(pipe_out[0], buffer, sizeof(buffer));
            if (bytes_read > 0)
            {
              callbackStdOut_(String(std::string(buffer, buffer + bytes_read)));
            }
          }
          if (FD_ISSET(pipe_err[0], &fds))
          {
            ssize_t bytes_read = read(pipe_err[0], buffer, sizeof(buffer));
            if (bytes_read > 0)
            {
              callbackStdErr_(String(std::string(buffer, buffer + bytes_read)));
            }
          }
        }

        int status_tmp = 0;
        pid_t pr = waitpid(pid, &status_tmp, WNOHANG);
        if (pr == pid)
        {
          child_finished = true;
          status_from_loop = status_tmp;
          break;
        }
      }

      // drain any remaining
      ssize_t bytes_read;
      while ((bytes_read = read(pipe_out[0], buffer, sizeof(buffer))) > 0)
      {
        callbackStdOut_(String(std::string(buffer, buffer + bytes_read)));
      }
      while ((bytes_read = read(pipe_err[0], buffer, sizeof(buffer))) > 0)
      {
        callbackStdErr_(String(std::string(buffer, buffer + bytes_read)));
      }

      close(pipe_out[0]);
      close(pipe_err[0]);

      int status = 0;
      if (!child_finished)
      {
        waitpid(pid, &status, 0);
      }
      else
      {
        status = status_from_loop;
      }

      if (WIFEXITED(status))
      {
        int exit_code = WEXITSTATUS(status);
        if (exit_code == 0)
        {
          if (verbose) { callbackStdOut_("Executed '" + String(exe) + "' successfully!\n"); }
          return RETURNSTATE::SUCCESS;
        }
        if (exit_code == 127)
        {
          error_msg = "Process '" + String(exe) + "' failed to start. Does it exist? Is it executable?";
          if (verbose) { callbackStdErr_(error_msg + '\n'); }
          return RETURNSTATE::FAILED_TO_START;
        }
        error_msg = "Process '" + String(exe) + "' did not finish successfully (exit code: " + String(exit_code) + "). Please check the log.";
        if (verbose) { callbackStdErr_(error_msg + '\n'); }
        return RETURNSTATE::NONZERO_EXIT;
      }
      else if (WIFSIGNALED(status))
      {
        error_msg = "Process '" + String(exe) + "' crashed (signal: " + String(WTERMSIG(status)) + ").";
        if (verbose) { callbackStdErr_(error_msg + '\n'); }
        return RETURNSTATE::CRASH;
      }

      if (verbose) { callbackStdOut_("Executed '" + String(exe) + "' successfully!\n"); }
      return RETURNSTATE::SUCCESS;
    }
    else
    {
      // No IO requested: just wait for child
      int status = 0;
      waitpid(pid, &status, 0);
      if (WIFEXITED(status))
      {
        int exit_code = WEXITSTATUS(status);
        if (exit_code == 0)
        {
          if (verbose) { callbackStdOut_("Executed '" + String(exe) + "' successfully!\n"); }
          return RETURNSTATE::SUCCESS;
        }
        if (exit_code == 127)
        {
          error_msg = "Process '" + String(exe) + "' failed to start. Does it exist? Is it executable?";
          if (verbose) { callbackStdErr_(error_msg + '\n'); }
          return RETURNSTATE::FAILED_TO_START;
        }
        error_msg = "Process '" + String(exe) + "' did not finish successfully (exit code: " + String(exit_code) + "). Please check the log.";
        if (verbose) { callbackStdErr_(error_msg + '\n'); }
        return RETURNSTATE::NONZERO_EXIT;
      }
      else if (WIFSIGNALED(status))
      {
        error_msg = "Process '" + String(exe) + "' crashed (signal: " + String(WTERMSIG(status)) + ").";
        if (verbose) { callbackStdErr_(error_msg + '\n'); }
        return RETURNSTATE::CRASH;
      }
      if (verbose) { callbackStdOut_("Executed '" + String(exe) + "' successfully!\n"); }
      return RETURNSTATE::SUCCESS;
    }
#endif
  }

} // namespace OpenMS
