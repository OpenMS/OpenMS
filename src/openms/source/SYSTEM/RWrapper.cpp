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

namespace OpenMS
{

  bool RWrapper::runScript( const String& script_file, const std::vector<std::string>& cmd_args, const std::string& executable /*= "Rscript"*/, bool find_R /*= false */, bool verbose /*= true */)
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
    
    std::vector<std::string> args = {"--vanilla", "--quiet", static_cast<std::string>(fullscript)};
    // Add command arguments
    args.insert(args.end(), cmd_args.begin(), cmd_args.end());

    ExternalProcess process;
    String error_msg;
    ExternalProcess::RETURNSTATE result = process.run(executable, args, "", verbose, error_msg);

    if (result != ExternalProcess::RETURNSTATE::SUCCESS)
    {
      if (verbose)
      {
        OPENMS_LOG_INFO << " failed" << std::endl;
        OPENMS_LOG_ERROR << "\nScript failed. Error: " << error_msg << std::endl;
      }
      return false;
    }
    if (verbose)
    {
      OPENMS_LOG_INFO << " success" << std::endl;
    }
    return true;
  }

  bool RWrapper::findR( const std::string& executable /*= "Rscript"*/, bool verbose /*= true*/ )
  {
    if (verbose) OPENMS_LOG_INFO << "Finding R interpreter 'Rscript' ...";

    std::vector<std::string> args = {"--vanilla", "-e", "sessionInfo()"};
    
    ExternalProcess process;
    String error_msg;
    ExternalProcess::RETURNSTATE result = process.run(executable, args, "", false, error_msg);

    if (result == ExternalProcess::RETURNSTATE::FAILED_TO_START)
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
      OPENMS_LOG_INFO << "Trying to invoke 'Rscript' ...";
    }
    if (result != ExternalProcess::RETURNSTATE::SUCCESS)
    {
      if (verbose)
      {
        OPENMS_LOG_INFO << " failed" << std::endl;
        OPENMS_LOG_ERROR << "Error: 'Rscript' executable returned with error\n"
                  << "Error: " << error_msg << "\n"
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
