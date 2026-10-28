// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/SYSTEM/ExternalProcess.h>

#ifdef _WIN32
  #ifndef NOMINMAX
    #define NOMINMAX
  #endif
  #ifndef WIN32_LEAN_AND_MEAN
    #define WIN32_LEAN_AND_MEAN
  #endif
  // Asio requires Winsock2 before windows.h can include the legacy Winsock API.
  #include <winsock2.h>
  #include <windows.h>
#endif

#include <boost/version.hpp>
#include <boost/asio/io_context.hpp>

// Boost.Process v1 compatibility shims removed in Boost 1.88; use v1/ prefix for 1.88+
#if BOOST_VERSION >= 108800
  #include <boost/process/v1/args.hpp>
  #include <boost/process/v1/async_pipe.hpp>
  #include <boost/process/v1/child.hpp>
  #include <boost/process/v1/env.hpp>
  #include <boost/process/v1/group.hpp>
  #include <boost/process/v1/io.hpp>
  #include <boost/process/v1/search_path.hpp>
  #include <boost/process/v1/start_dir.hpp>
  #include <boost/process/v1/extend.hpp>
  #ifdef _WIN32
    #include <boost/process/v1/error.hpp>
  #endif
#else
  #include <boost/process/args.hpp>
  #include <boost/process/async_pipe.hpp>
  #include <boost/process/child.hpp>
  #include <boost/process/env.hpp>
  #include <boost/process/group.hpp>
  #include <boost/process/io.hpp>
  #include <boost/process/search_path.hpp>
  #include <boost/process/start_dir.hpp>
  #include <boost/process/extend.hpp>
  #ifdef _WIN32
    #include <boost/process/error.hpp>
  #endif
#endif

#include <array>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <memory>
#include <optional>
#include <system_error>
#include <thread>
#include <utility>

#ifndef _WIN32
  #include <cerrno>
  #include <fcntl.h>
  #include <sys/resource.h>
  #include <sys/wait.h> // for WIFSIGNALED
  #include <unistd.h>
  #ifdef __linux__
    #include <sys/syscall.h>
    // Sysroots with old kernel headers (e.g. conda-forge, manylinux) lack the
    // number, but it is 436 on every Linux architecture except alpha.
    #if ! defined(SYS_close_range) && ! defined(__alpha__)
      #define SYS_close_range 436
    #endif
  #endif
#endif

#if BOOST_VERSION >= 108800
namespace bp = boost::process::v1;
#else
namespace bp = boost::process;
#endif

namespace OpenMS
{
#ifdef _WIN32
namespace
{
  // Windows passes a command-line string, not argv. Boost.Process v1 does not
  // escape trailing backslashes, empty arguments or tabs according to the CRT
  // rules, so encode it ourselves before CreateProcessW is called.
  std::wstring windowsCommandLine(const std::wstring& executable, const std::vector<std::wstring>& arguments)
  {
    std::wstring command = L"\"" + executable + L"\"";
    for (const auto& argument : arguments)
    {
      if (! argument.empty() && argument.find_first_of(L" \t\n\v\"") == std::wstring::npos)
      {
        command += L" " + argument;
        continue;
      }
      command += L" \"";
      std::size_t backslashes = 0;
      for (const wchar_t character : argument)
      {
        if (character == L'\\')
        {
          ++backslashes;
          continue;
        }
        // Before a quote, 2N+1 backslashes encode N backslashes and a literal
        // quote. Elsewhere backslashes have no special meaning.
        command.append(character == L'\"' ? 2 * backslashes + 1 : backslashes, L'\\');
        command += character;
        backslashes = 0;
      }
      // Double trailing backslashes so they cannot escape the closing quote.
      command.append(2 * backslashes, L'\\');
      command += L'\"';
    }
    return command;
  }

  // OpenMS strings are UTF-8. Match PathUtils::to_path's fallback for callers
  // that still pass an argument encoded in the native Windows code page.
  std::wstring wideString(const std::string& text)
  {
    if (text.empty()) return {};
    if (text.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
      throw bp::process_error(std::make_error_code(std::errc::value_too_large), "Process argument is too long");
    }
    const int length = static_cast<int>(text.size());
    UINT code_page = CP_UTF8;
    DWORD flags = MB_ERR_INVALID_CHARS;
    int size = MultiByteToWideChar(code_page, flags, text.data(), length, nullptr, 0);
    if (size == 0 && GetLastError() == ERROR_NO_UNICODE_TRANSLATION)
    {
      code_page = CP_ACP;
      flags = 0;
      size = MultiByteToWideChar(code_page, flags, text.data(), length, nullptr, 0);
    }
    if (size == 0)
    {
      throw bp::process_error(std::error_code(static_cast<int>(GetLastError()), std::system_category()), "Cannot convert process argument to UTF-16");
    }
    std::wstring result(size, L'\0');
    if (MultiByteToWideChar(code_page, flags, text.data(), length, result.data(), size) == 0)
    {
      throw bp::process_error(std::error_code(static_cast<int>(GetLastError()), std::system_category()), "Cannot convert process argument to UTF-16");
    }
    return result;
  }
} // namespace
#endif

/// default Ctor; callbacks for stdout/stderr are empty
ExternalProcess::ExternalProcess(): ExternalProcess([](const std::string& /*out*/) {}, [](const std::string& /*out*/) {})
{
}

  ExternalProcess::ExternalProcess(std::function<void(const std::string&)> callbackStdOut, std::function<void(const std::string&)> callbackStdErr)
    : callbackStdOut_(std::move(callbackStdOut)),
      callbackStdErr_(std::move(callbackStdErr))
  {
  }

  ExternalProcess::~ExternalProcess() = default;

  /// re-wire the callbacks used using run()
  void ExternalProcess::setCallbacks(std::function<void(const std::string&)> callbackStdOut, std::function<void(const std::string&)> callbackStdErr)
  {
    callbackStdOut_ = std::move(callbackStdOut);
    callbackStdErr_ = std::move(callbackStdErr);
  }

  ExternalProcess::RETURNSTATE ExternalProcess::run(const std::string& exe, const std::vector<std::string>& args, const std::string& working_dir, const bool verbose,
                                                     IO_MODE io_mode, const std::map<std::string, std::string>& env, std::function<void()> idle_callback)
  {
    std::string error_msg;
    return run(exe, args, working_dir, verbose, error_msg, io_mode, env, std::move(idle_callback));
  }

  ExternalProcess::RETURNSTATE ExternalProcess::run(const std::string& exe, const std::vector<std::string>& args, const std::string& working_dir, const bool verbose,
                                                     std::string& error_msg, IO_MODE io_mode, const std::map<std::string, std::string>& env, std::function<void()> idle_callback)
  {
    Result result = runWithResult(exe, args, working_dir, verbose, io_mode, env, std::move(idle_callback));
    error_msg = std::move(result.error_message);
    return result.state;
  }

  ExternalProcess::Result ExternalProcess::runWithResult(const std::string& exe,
                                                         const std::vector<std::string>& args,
                                                         const std::string& working_dir,
                                                         bool verbose,
                                                         IO_MODE io_mode,
                                                         const std::map<std::string, std::string>& env,
                                                         std::function<void()> idle_callback,
                                                         const std::atomic_bool* cancel)
  {
    Result result;
    auto cancellation_requested = [&]() { return cancel != nullptr && cancel->load(std::memory_order_relaxed); };
    if (cancellation_requested())
    {
      result.state = RETURNSTATE::CANCELLED;
      result.error_message = "Process '" + exe + "' was cancelled before launch.";
      return result;
    }

    if (verbose)
    {
      std::string command = "Running: " + exe;
      for (const auto& argument : args)
      {
        command += " " + argument;
      }
      callbackStdOut_(command + "\n");
    }

    // Only cancellable calls own a process group (POSIX) or job object (Windows).
    // Ordinary calls keep the child in the caller's group: terminal signals such
    // as Ctrl+C still reach it, and an enclosing cancellable call (e.g. a workflow
    // running an adapter which runs a search engine) can still stop the whole tree.
    const bool own_group = cancel != nullptr;
    std::optional<bp::group> group; // created in the try block: on Windows, creating a job can throw
    bp::child child;
    auto stop_and_wait = [&]() {
      std::error_code ignored;
      bool group_terminated = false;
      if (group && group->valid())
      {
        std::error_code termination_error;
        group->terminate(termination_error);
        group_terminated = ! termination_error;
      }
      if (child.valid())
      {
        // Without a group, or if terminating it failed, stop the direct child.
        if (! group_terminated && child.running(ignored)) { child.terminate(ignored); }
        child.wait(ignored);
      }
      if (group) { group->detach(); }
    };

    try
    {
      if (own_group) { group.emplace(); }
      // Boost.Process does not search PATH when passed a plain executable name.
#ifdef _WIN32
      std::wstring executable = wideString(exe);
      const auto resolved = bp::search_path(executable);
      if (! resolved.empty()) executable = resolved.wstring();
      std::vector<std::wstring> process_arguments;
      process_arguments.reserve(args.size());
      for (const auto& argument : args)
        process_arguments.push_back(wideString(argument));
      const std::wstring start_dir = working_dir.empty() ? L"." : wideString(working_dir);
      bp::wenvironment process_environment = boost::this_process::wenvironment();
      for (const auto& [key, value] : env)
      {
        const auto name = wideString(key);
        // Boost v1's environment container compares keys case-sensitively,
        // while Windows does not (notably inherited Path versus supplied PATH).
        std::vector<std::wstring> existing_names;
        for (const auto& entry : process_environment)
        {
          const auto existing = entry.get_name();
          if (CompareStringOrdinal(existing.c_str(), -1, name.c_str(), -1, TRUE) == CSTR_EQUAL) { existing_names.push_back(existing); }
        }
        for (const auto& existing : existing_names)
          process_environment.erase(existing);
        process_environment[name] = wideString(value);
      }
      auto command_line = windowsCommandLine(executable, process_arguments);
      auto set_command_line = bp::extend::on_setup([&](auto& executor) { executor.cmd_line = command_line.data(); });

      // Boost v1 assigns a child to its job in group_ref::on_success. Launch
      // suspended so it cannot spawn descendants before that assignment. Built
      // group initializers run before explicit extend handlers (execute_impl.hpp),
      // so the resume handler below runs only after the group has been assigned.
      auto suspend = bp::extend::on_setup([](auto& executor) { executor.creation_flags |= CREATE_SUSPENDED; });
      auto resume = bp::extend::on_success([](auto& executor) {
        if (ResumeThread(executor.proc_info.hThread) == static_cast<DWORD>(-1))
        {
          executor.set_error(std::error_code(static_cast<int>(GetLastError()), std::system_category()), "ResumeThread failed");
        }
      });
      std::error_code launch_error;
      auto cleanup_failed_launch = bp::extend::on_error([](auto& executor, const std::error_code&) noexcept {
        // Assignment/resume can fail before bp::child is returned to us. The
        // executor still owns these handles: terminate and wait, never close them.
        const auto process = executor.proc_info.hProcess;
        if (process != nullptr && process != INVALID_HANDLE_VALUE)
        {
          TerminateProcess(process, EXIT_FAILURE);
          WaitForSingleObject(process, INFINITE);
        }
      });
#else
      std::string executable = exe;
      const auto resolved = bp::search_path(executable);
      if (! resolved.empty()) { executable = resolved.string(); }
      // Use a separate environment: native_environment writes through to the
      // parent process, which would also race with other workflow workers.
      bp::environment process_environment = boost::this_process::environment();
      for (const auto& [key, value] : env)
      {
        process_environment[key] = value;
      }
      const auto& process_arguments = args;
      const std::string start_dir = working_dir.empty() ? "." : working_dir;
      struct rlimit descriptor_limits {};
      if (::getrlimit(RLIMIT_NOFILE, &descriptor_limits) != 0)
      {
        throw bp::process_error(std::error_code(errno, std::system_category()), "Cannot determine descriptor limit");
      }
      const auto max_descriptor = static_cast<rlim_t>(std::numeric_limits<int>::max());
      const auto descriptor_limit = static_cast<int>(descriptor_limits.rlim_cur < max_descriptor ? descriptor_limits.rlim_cur : max_descriptor);
      auto restrict_descriptors = bp::extend::on_exec_setup([descriptor_limit](auto& executor) {
        // Mark after fork: parallel launches may create more pipes at any time.
        // CLOEXEC retains Boost's launch-error pipe until a successful exec, and
        // leaves the child's redirected standard streams and parent untouched.
        auto mark = [&executor](int descriptor) {
          int status;
          do
          {
            status = ::fcntl(descriptor, F_SETFD, FD_CLOEXEC);
          } while (status == -1 && errno == EINTR);
          if (status == -1 && errno != EBADF)
          {
            executor.set_error(std::error_code(errno, std::system_category()), "Cannot restrict inherited descriptors");
            ::_exit(EXIT_FAILURE);
          }
        };
#if defined(__linux__) && defined(SYS_close_range)
        constexpr unsigned close_range_cloexec = 1U << 2; // Linux CLOSE_RANGE_CLOEXEC ABI (5.11+)
        if (::syscall(SYS_close_range, 3U, std::numeric_limits<unsigned>::max(), close_range_cloexec) == 0) return;
#endif
#ifdef __linux__
        // Older kernels: mark only the open descriptors. Containers can have a
        // descriptor limit of 2^30, which would make the scan below take minutes.
        // open, getdents64, fcntl and close are async-signal-safe and allocate nothing.
        const int listing = ::open("/proc/self/fd", O_RDONLY | O_DIRECTORY | O_CLOEXEC);
        if (listing >= 0)
        {
          alignas(8) char buffer[4096];
          bool complete = true;
          for (;;)
          {
            const long count = ::syscall(SYS_getdents64, listing, buffer, sizeof(buffer));
            if (count == 0) break;
            if (count < 0)
            {
              if (errno == EINTR) continue;
              complete = false;
              break;
            }
            for (long offset = 0; offset < count;)
            {
              // struct linux_dirent64: d_ino (8), d_off (8), d_reclen (2), d_type (1), d_name
              unsigned short length = 0;
              std::memcpy(&length, buffer + offset + 16, sizeof(length));
              const char* name = buffer + offset + 19;
              int descriptor = 0;
              bool numeric = *name != '\0';
              for (const char* c = name; *c != '\0'; ++c)
              {
                if (*c < '0' || *c > '9' || descriptor > (std::numeric_limits<int>::max() - 9) / 10)
                {
                  numeric = false;
                  break;
                }
                descriptor = descriptor * 10 + (*c - '0');
              }
              if (numeric && descriptor >= 3 && descriptor != listing) { mark(descriptor); }
              offset += length;
            }
          }
          ::close(listing);
          if (complete) return;
        }
#endif
        // Portable fallback: scan the configured fd limit.
        for (int descriptor = 3; descriptor < descriptor_limit; ++descriptor)
        {
          mark(descriptor);
        }
      });
#endif

      const bool can_read = io_mode == IO_MODE::READ_ONLY || io_mode == IO_MODE::READ_WRITE;
      boost::asio::io_context io_context;
      std::unique_ptr<bp::async_pipe> stdout_pipe, stderr_pipe;
      auto launch = [&](auto&&... properties) {
        return bp::child(executable, bp::args(process_arguments), bp::start_dir(start_dir), process_environment,
                         std::forward<decltype(properties)>(properties)...,
#ifdef _WIN32
                         launch_error, set_command_line, suspend, resume, cleanup_failed_launch
#else
                         restrict_descriptors
#endif
        );
      };
      // A child in its own POSIX process group is a background job of the
      // caller's terminal: reading the terminal would stop it (SIGTTIN) and
      // block this call. Give it an empty stdin, as QProcess-based runners did.
      if (can_read)
      {
        stdout_pipe = std::make_unique<bp::async_pipe>(io_context);
        stderr_pipe = std::make_unique<bp::async_pipe>(io_context);
        child = own_group ? launch(bp::std_out > *stdout_pipe, bp::std_err > *stderr_pipe, bp::std_in < bp::null, *group)
                          : launch(bp::std_out > *stdout_pipe, bp::std_err > *stderr_pipe);
      }
      else
      {
        child = own_group ? launch(bp::std_out > bp::null, bp::std_err > bp::null, bp::std_in < bp::null, *group)
                          : launch(bp::std_out > bp::null, bp::std_err > bp::null);
      }

#ifdef _WIN32
      if (launch_error) { throw bp::process_error(launch_error, "Cannot launch or initialize process group"); }
#endif

      if (! child.valid())
      {
        result.error_message = "Process '" + exe + "' failed to start. Does it exist? Is it executable?";
        stop_and_wait();
      }
      else
      {
        std::array<char, 4096> stdout_buffer {}, stderr_buffer {};
        bool stdout_done = ! can_read, stderr_done = ! can_read;
        std::function<void(const boost::system::error_code&, std::size_t)> read_stdout, read_stderr;
        read_stdout = [&](const boost::system::error_code& error, std::size_t count) {
          if (count != 0) { callbackStdOut_(std::string(stdout_buffer.data(), count)); }
          stdout_done = bool(error);
          if (! error) { stdout_pipe->async_read_some(boost::asio::buffer(stdout_buffer), read_stdout); }
        };
        read_stderr = [&](const boost::system::error_code& error, std::size_t count) {
          if (count != 0) { callbackStdErr_(std::string(stderr_buffer.data(), count)); }
          stderr_done = bool(error);
          if (! error) { stderr_pipe->async_read_some(boost::asio::buffer(stderr_buffer), read_stderr); }
        };
        if (can_read)
        {
          stdout_pipe->async_read_some(boost::asio::buffer(stdout_buffer), read_stdout);
          stderr_pipe->async_read_some(boost::asio::buffer(stderr_buffer), read_stderr);
        }

        bool cancelled = false;
        auto poll = [&]() {
          if (stdout_done && stderr_done) { std::this_thread::sleep_for(std::chrono::milliseconds(25)); }
          else
          {
            io_context.restart();
            io_context.run_for(std::chrono::milliseconds(25));
          }
          if (idle_callback) { idle_callback(); }
        };
        while (child.running())
        {
          if (cancellation_requested())
          {
            cancelled = true;
            stop_and_wait();
            break;
          }
          poll();
        }
        child.wait();
        result.exit_code = child.exit_code();
#ifdef _WIN32
        const bool crashed = result.exit_code < 0 || static_cast<unsigned int>(result.exit_code) > 0x80000000u;
#else
        // Darwin's wait-status macros require an lvalue.
        int native_status = child.native_exit_code();
        const bool crashed = WIFSIGNALED(native_status);
#endif
        // A failed wrapper can leave its workers alive. In a cancellable call, stop
        // them before draining their pipes or returning failure to the workflow,
        // which may remove their temporary directory or start a new run using the
        // same output files. Without an owned group, this only reaps the child.
        if (! cancelled && (crashed || result.exit_code != 0)) { stop_and_wait(); }

        // A descendant may inherit a pipe and outlive the direct child. Never wait
        // indefinitely for that descendant to close the pipe after the child exits.
        const auto drain_deadline = std::chrono::steady_clock::now() + std::chrono::seconds(1);
        while (! (stdout_done && stderr_done) && std::chrono::steady_clock::now() < drain_deadline)
        {
          if (! cancelled && cancellation_requested())
          {
            cancelled = true;
            stop_and_wait();
          }
          poll();
        }
        if (can_read)
        {
          boost::system::error_code ignored;
          stdout_pipe->close(ignored);
          stderr_pipe->close(ignored);
          io_context.restart();
          io_context.poll();
        }
        if (group) { group->detach(); }

        if (cancelled)
        {
          result.state = RETURNSTATE::CANCELLED;
          result.error_message = "Process '" + exe + "' was cancelled.";
        }
        else if (crashed)
        {
          result.state = RETURNSTATE::CRASH;
          result.error_message = "Process '" + exe + "' crashed hard (segfault-like). Please check the log.";
        }
        else if (result.exit_code != 0)
        {
          result.state = RETURNSTATE::NONZERO_EXIT;
          result.error_message
            = "Process '" + exe + "' did not finish successfully (exit code: " + StringUtils::toStr(result.exit_code) + "). Please check the log.";
        }
        else
        {
          result.state = RETURNSTATE::SUCCESS;
        }
      }
    }
    catch (const bp::process_error& error)
    {
      const bool started = child.valid();
      stop_and_wait();
      result.state = started ? RETURNSTATE::CRASH : RETURNSTATE::FAILED_TO_START;
      result.error_message = "Process '" + exe + (started ? "' could not be monitored: " : "' failed to start: ") + error.what();
    }
    catch (...)
    {
      // Exceptions from caller callbacks must not leave a child or its descendants running.
      stop_and_wait();
      throw;
    }

    if (verbose)
    {
      if (result.state == RETURNSTATE::SUCCESS) { callbackStdOut_("Executed '" + exe + "' successfully!\n"); }
      else
      {
        callbackStdErr_(result.error_message + "\n");
      }
    }
    return result;
  }

} // namespace OpenMS
