// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/RWrapper.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/SYSTEM/File.h>

#include <boost/version.hpp>

// Boost.Process v1 compatibility shims removed in Boost 1.88; use v1/ prefix for 1.88+
#if BOOST_VERSION >= 108800
#include <boost/process/v1/child.hpp>
#include <boost/process/v1/args.hpp>
#include <boost/process/v1/io.hpp>
#include <boost/process/v1/search_path.hpp>
#else
#include <boost/process/child.hpp>
#include <boost/process/args.hpp>
#include <boost/process/io.hpp>
#include <boost/process/search_path.hpp>
#endif

#include <string>
#include <thread>

#if BOOST_VERSION >= 108800
namespace bp = boost::process::v1;
#else
namespace bp = boost::process;
#endif

namespace OpenMS
{

  bool RWrapper::runScript(const String& script_file, const std::vector<String>& cmd_args, const String& executable /*= "Rscript"*/, bool find_R /*= false */, bool verbose /*= true */)
  {
    if (find_R && !findR(executable, verbose))
    {
      return false;
    }

    String fullscript;
    try
    {
      fullscript = findScript(script_file, verbose);
    }
    catch (...)
    {
      return false;
    }

    if (verbose)
    {
      OPENMS_LOG_INFO << "Running R script '" << fullscript << "' ...";
    }

    // Build argument list: --vanilla --quiet <script> <cmd_args...>
    std::vector<std::string> args;
    args.push_back("--vanilla");
    args.push_back("--quiet");
    args.push_back(static_cast<std::string>(fullscript));
    for (const auto& a : cmd_args)
    {
      args.push_back(static_cast<std::string>(a));
    }

    try
    {
      bp::ipstream pipe_out;
      bp::ipstream pipe_err;
      bp::child child(
        bp::search_path(static_cast<std::string>(executable)),
        bp::args(args),
        bp::std_out > pipe_out,
        bp::std_err > pipe_err
      );

      // Drain both pipes concurrently to prevent deadlock
      // (child may block if pipe buffer fills before we read)
      std::string stdout_content, stderr_content;
      auto drain_out = [&pipe_out, &stdout_content]() {
        std::string line;
        while (std::getline(pipe_out, line)) { stdout_content += line + "\n"; }
      };
      auto drain_err = [&pipe_err, &stderr_content]() {
        std::string line;
        while (std::getline(pipe_err, line)) { stderr_content += line + "\n"; }
      };
      std::thread reader_err(drain_err); // start stderr reader first
      drain_out();                        // read stdout on this thread
      reader_err.join();

      child.wait();

      int exit_code = child.exit_code();
      if (exit_code != 0)
      {
        if (verbose)
        {
          OPENMS_LOG_INFO << " failed" << std::endl;
          OPENMS_LOG_ERROR << "\n--- ERROR MESSAGES ---\n" << stderr_content;
          OPENMS_LOG_ERROR << "\n--- OTHER MESSAGES ---\n" << stdout_content;
          OPENMS_LOG_ERROR << "\n\nScript failed. See above for an error description. " << std::endl;
        }
        return false;
      }
    }
    catch (const bp::process_error&)
    {
      if (verbose)
      {
        OPENMS_LOG_INFO << " failed" << std::endl;
        OPENMS_LOG_ERROR << "Error: Could not run '" << executable << "'. Is it installed and in PATH?" << std::endl;
      }
      return false;
    }

    if (verbose)
    {
      OPENMS_LOG_INFO << " success" << std::endl;
    }
    return true;
  }

  bool RWrapper::findR(const String& executable /*= "Rscript"*/, bool verbose /*= true*/)
  {
    if (verbose) OPENMS_LOG_INFO << "Finding R interpreter 'Rscript' ...";

    std::vector<std::string> args = {"--vanilla", "-e", "sessionInfo()"};

    try
    {
      // Use both pipes to match Qt MergedChannels behavior (stderr merged into output)
      bp::ipstream pipe_out;
      bp::ipstream pipe_err;
      bp::child child(
        bp::search_path(static_cast<std::string>(executable)),
        bp::args(args),
        bp::std_out > pipe_out,
        bp::std_err > pipe_err
      );

      // Drain both pipes before wait to prevent deadlock
      std::string stdout_content, stderr_content;
      std::thread reader_err([&pipe_err, &stderr_content]() {
        std::string line;
        while (std::getline(pipe_err, line)) { stderr_content += line + "\n"; }
      });
      {
        std::string line;
        while (std::getline(pipe_out, line)) { stdout_content += line + "\n"; }
      }
      reader_err.join();
      // Merge like Qt MergedChannels
      std::string merged_content = stdout_content + stderr_content;

      child.wait();

      if (child.exit_code() != 0)
      {
        if (verbose)
        {
          OPENMS_LOG_INFO << " failed" << std::endl;
          // Construct args string for error message
          std::string args_str;
          for (const auto& a : args)
          {
            if (!args_str.empty()) args_str += " ";
            args_str += a;
          }
          OPENMS_LOG_ERROR << "Error: 'Rscript' executable returned with error (command: 'Rscript " << args_str << "')\n";
          OPENMS_LOG_ERROR << "Output was:\n------>\n" << merged_content
                    << "\n<------\n"
                    << "Make sure 'Rscript' is installed properly." << std::endl;
        }
        return false;
      }
    }
    catch (const bp::process_error&)
    {
      if (verbose)
      {
        OPENMS_LOG_INFO << " failed" << std::endl;
        OPENMS_LOG_ERROR << "Error: Could not find or run '" << executable << "' executable (FailedToStart).\n";
        OPENMS_LOG_ERROR << "Please install 'Rscript', make sure it's in PATH and is flagged as executable." << std::endl;
      }
      return false;
    }

    if (verbose)
    {
      OPENMS_LOG_INFO << " success" << std::endl;
    }
    if (verbose)
    {
      OPENMS_LOG_INFO << "Trying to invoke 'Rscript' ... success" << std::endl;
    }
    return true;
  }

  OpenMS::String RWrapper::findScript(const String& script_file, bool verbose /*= true*/)
  {
    String s;
    try
    {
      s = File::find(script_file, StringList(1, File::getOpenMSDataPath().ensureLastChar('/') + "SCRIPTS"));
    }
    catch (...)
    {
      if (verbose)
      {
        OPENMS_LOG_ERROR << "\n\nCould not find R script '" << script_file << "'!\n" << std::endl;
      }
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, script_file);
    }
    return s;
  }

} // namespace OpenMS
