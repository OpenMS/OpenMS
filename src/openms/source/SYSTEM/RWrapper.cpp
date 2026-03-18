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

#include <boost/process/child.hpp>
#include <boost/process/args.hpp>
#include <boost/process/io.hpp>
#include <boost/process/search_path.hpp>

#include <string>

namespace bp = boost::process;

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
        static_cast<std::string>(executable),
        bp::args(args),
        bp::std_out > pipe_out,
        bp::std_err > pipe_err
      );

      child.wait();

      int exit_code = child.exit_code();
      if (exit_code != 0)
      {
        if (verbose)
        {
          OPENMS_LOG_INFO << " failed" << std::endl;
          OPENMS_LOG_ERROR << "\n--- ERROR MESSAGES ---\n";
          std::string line;
          while (std::getline(pipe_err, line))
          {
            OPENMS_LOG_ERROR << line << "\n";
          }
          OPENMS_LOG_ERROR << "\n--- OTHER MESSAGES ---\n";
          while (std::getline(pipe_out, line))
          {
            OPENMS_LOG_ERROR << line << "\n";
          }
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
      bp::ipstream pipe_out;
      bp::child child(
        static_cast<std::string>(executable),
        bp::args(args),
        bp::std_out > pipe_out,
        bp::std_err > bp::null // merge-like: just discard stderr
      );

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
          std::string line;
          OPENMS_LOG_ERROR << "Output was:\n------>\n";
          while (std::getline(pipe_out, line))
          {
            OPENMS_LOG_ERROR << line << "\n";
          }
          OPENMS_LOG_ERROR << "\n<------\n"
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
