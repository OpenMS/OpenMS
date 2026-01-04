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
#include <OpenMS/SYSTEM/ExternalProcess.h>

#include <sstream>

namespace OpenMS
{

  bool RWrapper::runScript( const String& script_file, const std::vector<String>& cmd_args, const String& executable /*= "Rscript"*/, bool find_R /*= false */, bool verbose /*= true */)
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

    std::vector<String> args;
    args.push_back("--vanilla");
    args.push_back("--quiet");
    args.push_back(fullscript);
    args.insert(args.end(), cmd_args.begin(), cmd_args.end());

    String proc_stdout, proc_stderr;
    auto lam_out = [&](const String& out) { proc_stdout += out; };
    auto lam_err = [&](const String& out) { proc_stderr += out; };
    ExternalProcess ep(lam_out, lam_err);

    auto rt = ep.run(executable, args, "", false); // verbose=false for ExternalProcess

    if (rt != ExternalProcess::RETURNSTATE::SUCCESS)
    {
      if (verbose)
      {
        OPENMS_LOG_INFO << " failed" << std::endl;
        OPENMS_LOG_ERROR << "\n--- ERROR MESSAGES ---\n";
        OPENMS_LOG_ERROR << proc_stderr;
        OPENMS_LOG_ERROR << "\n--- OTHER MESSAGES ---\n";
        OPENMS_LOG_ERROR << proc_stdout;
        OPENMS_LOG_ERROR << "\n\nScript failed. See above for an error description. " << std::endl;
      }
      return false;
    }
    if (verbose)
    {
      OPENMS_LOG_INFO << " success" << std::endl;
    }
    return true;
  }

  bool RWrapper::findR( const String& executable /*= "Rscript"*/, bool verbose /*= true*/ )
  {
    if (verbose) OPENMS_LOG_INFO << "Finding R interpreter 'Rscript' ...";

    std::vector<String> args = {"--vanilla", "-e", "sessionInfo()"};

    String proc_stdout, proc_stderr;
    auto lam_out = [&](const String& out) { proc_stdout += out; };
    auto lam_err = [&](const String& out) { proc_stderr += out; };
    ExternalProcess ep(lam_out, lam_err);

    auto rt = ep.run(executable, args, "", false); // verbose=false for ExternalProcess

    if (rt == ExternalProcess::RETURNSTATE::FAILED_TO_START)
    {
      if (verbose)
      {
        OPENMS_LOG_INFO << " failed" << std::endl;
        OPENMS_LOG_ERROR << "Error: Could not find or run '" << executable << "' executable (FailedToStart).\n";
        if (!proc_stdout.empty())
        {
          OPENMS_LOG_ERROR << "Output was:\n------>\n"
                    << proc_stdout
                    << "\n<------\n";
        }
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
      OPENMS_LOG_INFO << "Trying to invoke 'Rscript' ...";
    }
    if (rt != ExternalProcess::RETURNSTATE::SUCCESS)
    {
      if (verbose)
      {
        OPENMS_LOG_INFO << " failed" << std::endl;

        // Build command string for error message
        std::ostringstream cmd_stream;
        for (size_t i = 0; i < args.size(); ++i) {
          if (i > 0) cmd_stream << " ";
          cmd_stream << args[i];
        }

        OPENMS_LOG_ERROR << "Error: 'Rscript' executable returned with error (command: 'Rscript " << cmd_stream.str() << "')\n"
                  << "Output was:\n------>\n"
                  << proc_stdout
                  << "\n<------\n"
                  << "Make sure 'Rscript' is installed properly." << std::endl;
      }
      return false;
    }
    if (verbose)
    {
      OPENMS_LOG_INFO << " success" << std::endl;
    }
    return true;
  }

  OpenMS::String RWrapper::findScript( const String& script_file, bool verbose /*= true*/ )
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
