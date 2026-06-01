// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/ExternalProcess.h>

#include <OpenMS/DATASTRUCTURES/String.h>

#include <boost/version.hpp>
#include <boost/asio/io_context.hpp>

// Boost.Process v1 compatibility shims removed in Boost 1.88; use v1/ prefix for 1.88+
#if BOOST_VERSION >= 108800
#include <boost/process/v1/child.hpp>
#include <boost/process/v1/args.hpp>
#include <boost/process/v1/io.hpp>
#include <boost/process/v1/async_pipe.hpp>
#include <boost/process/v1/start_dir.hpp>
#include <boost/process/v1/env.hpp>
#include <boost/process/v1/search_path.hpp>
#else
#include <boost/process/child.hpp>
#include <boost/process/args.hpp>
#include <boost/process/io.hpp>
#include <boost/process/async_pipe.hpp>
#include <boost/process/start_dir.hpp>
#include <boost/process/env.hpp>
#include <boost/process/search_path.hpp>
#endif

#include <array>
#include <chrono>
#include <utility>

#ifndef _WIN32
#include <sys/wait.h> // for WIFSIGNALED
#endif

#if BOOST_VERSION >= 108800
namespace bp = boost::process::v1;
#else
namespace bp = boost::process;
#endif

namespace OpenMS
{

  /// default Ctor; callbacks for stdout/stderr are empty
  ExternalProcess::ExternalProcess()
    : ExternalProcess([](const String& /*out*/) {}, [](const String& /*out*/) {})
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

  ExternalProcess::RETURNSTATE ExternalProcess::run(const String& exe, const std::vector<String>& args, const String& working_dir, const bool verbose,
                                                     IO_MODE io_mode, const std::map<String, String>& env, std::function<void()> idle_callback)
  {
    String error_msg;
    return run(exe, args, working_dir, verbose, error_msg, io_mode, env, std::move(idle_callback));
  }

  ExternalProcess::RETURNSTATE ExternalProcess::run(const String& exe, const std::vector<String>& args, const String& working_dir, const bool verbose,
                                                     String& error_msg, IO_MODE io_mode, const std::map<String, String>& env, std::function<void()> idle_callback)
  {
    error_msg.clear();

    if (verbose)
    {
      String cmd = "Running: " + exe;
      for (const auto& a : args)
      {
        cmd += " " + a;
      }
      cmd += "\n";
      callbackStdOut_(cmd);
    }

    // Convert args from OpenMS::String to std::string
    std::vector<std::string> std_args;
    std_args.reserve(args.size());
    for (const auto& a : args)
    {
      std_args.push_back(static_cast<std::string>(a));
    }

    // Resolve executable via PATH search (bp::child with a plain string does not search PATH)
    std::string exe_str = static_cast<std::string>(exe);
    auto resolved = bp::search_path(exe_str);
    if (!resolved.empty())
    {
      exe_str = resolved.string();
    }

    // Build environment: inherit current + add custom vars
    auto proc_env = boost::this_process::environment();
    for (const auto& kv : env)
    {
      proc_env[static_cast<std::string>(kv.first)] = static_cast<std::string>(kv.second);
    }

    try
    {
      // Determine whether we can read stdout/stderr
      bool can_read = (io_mode == IO_MODE::READ_ONLY || io_mode == IO_MODE::READ_WRITE);

      if (can_read)
      {
        // Async mode: use io_context + async_pipe for streaming output
        boost::asio::io_context io_ctx;
        bp::async_pipe ap_out(io_ctx);
        bp::async_pipe ap_err(io_ctx);

        std::string start_dir_str = working_dir.empty() ? "." : static_cast<std::string>(working_dir);

        bp::child child(
          exe_str,
          bp::args(std_args),
          bp::std_out > ap_out,
          bp::std_err > ap_err,
          bp::start_dir(start_dir_str),
          proc_env
        );

        if (!child.valid())
        {
          error_msg = "Process '" + exe + "' failed to start. Does it exist? Is it executable?";
          if (verbose)
          {
            callbackStdErr_(error_msg + "\n");
          }
          return RETURNSTATE::FAILED_TO_START;
        }

        // Async read buffers
        std::array<char, 4096> buf_out{};
        std::array<char, 4096> buf_err{};

        // Declare read handlers with std::function to allow recursive calls
        std::function<void(const boost::system::error_code&, std::size_t)> read_out;
        std::function<void(const boost::system::error_code&, std::size_t)> read_err;

        read_out = [&](const boost::system::error_code& ec, std::size_t n) {
          if (!ec && n > 0)
          {
            callbackStdOut_(String(std::string(buf_out.data(), n)));
          }
          if (!ec)
          {
            ap_out.async_read_some(boost::asio::buffer(buf_out), read_out);
          }
        };

        read_err = [&](const boost::system::error_code& ec, std::size_t n) {
          if (!ec && n > 0)
          {
            callbackStdErr_(String(std::string(buf_err.data(), n)));
          }
          if (!ec)
          {
            ap_err.async_read_some(boost::asio::buffer(buf_err), read_err);
          }
        };

        // Start async reads
        ap_out.async_read_some(boost::asio::buffer(buf_out), read_out);
        ap_err.async_read_some(boost::asio::buffer(buf_err), read_err);

        // Poll loop: run pending I/O handlers, call idle_callback for GUI event processing
        while (child.running())
        {
          io_ctx.run_for(std::chrono::milliseconds(50));
          if (idle_callback)
          {
            idle_callback();
          }
        }
        // Drain remaining I/O after child exits
        io_ctx.run();

        child.wait(); // collect exit status

        int exit_code = child.exit_code();

#ifdef _WIN32
        // On Windows, abnormal termination often results in specific exit codes
        // like 0xC0000005 (access violation), etc.
        bool crashed = (exit_code < 0 || static_cast<unsigned int>(exit_code) > 0x80000000u);
#else
        // On POSIX, use waitpid status macros on the native exit code
        auto native_code = child.native_exit_code();
        bool crashed = WIFSIGNALED(native_code);
#endif
        if (crashed)
        {
          error_msg = "Process '" + exe + "' crashed hard (segfault-like). Please check the log.";
          if (verbose)
          {
            callbackStdErr_(error_msg + "\n");
          }
          return RETURNSTATE::CRASH;
        }
        else if (exit_code != 0)
        {
          error_msg = "Process '" + exe + "' did not finish successfully (exit code: " + String(exit_code) + "). Please check the log.";
          if (verbose)
          {
            callbackStdErr_(error_msg + "\n");
          }
          return RETURNSTATE::NONZERO_EXIT;
        }
      }
      else
      {
        // No I/O mode or write-only: just run and wait
        std::string start_dir_str = working_dir.empty() ? "." : static_cast<std::string>(working_dir);

        bp::child child(
          exe_str,
          bp::args(std_args),
          bp::std_out > bp::null,
          bp::std_err > bp::null,
          bp::start_dir(start_dir_str),
          proc_env
        );

        if (!child.valid())
        {
          error_msg = "Process '" + exe + "' failed to start. Does it exist? Is it executable?";
          if (verbose)
          {
            callbackStdErr_(error_msg + "\n");
          }
          return RETURNSTATE::FAILED_TO_START;
        }

        // Poll loop to allow idle_callback
        while (child.running())
        {
          std::this_thread::sleep_for(std::chrono::milliseconds(50));
          if (idle_callback)
          {
            idle_callback();
          }
        }
        child.wait();

        int exit_code = child.exit_code();
        if (exit_code != 0)
        {
          error_msg = "Process '" + exe + "' did not finish successfully (exit code: " + String(exit_code) + "). Please check the log.";
          if (verbose)
          {
            callbackStdErr_(error_msg + "\n");
          }
          return RETURNSTATE::NONZERO_EXIT;
        }
      }
    }
    catch (const bp::process_error& /*e*/)
    {
      error_msg = "Process '" + exe + "' failed to start. Does it exist? Is it executable?";
      if (verbose)
      {
        callbackStdErr_(error_msg + "\n");
      }
      return RETURNSTATE::FAILED_TO_START;
    }

    if (verbose)
    {
      callbackStdOut_("Executed '" + exe + "' successfully!\n");
    }
    return RETURNSTATE::SUCCESS;
  }

} // namespace OpenMS
