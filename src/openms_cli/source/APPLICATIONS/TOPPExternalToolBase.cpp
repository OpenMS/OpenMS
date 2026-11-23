// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Clemens Groepl, Johannes Junker, Stephan Aiche, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPExternalToolBase.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/SYSTEM/ExternalProcess.h>
#include <OpenMS/SYSTEM/File.h>

using namespace std;

namespace OpenMS
{

  // constructors are inherited from TOPPBase (see header); the out-of-line dtor anchors the vtable here
  TOPPExternalToolBase::~TOPPExternalToolBase() = default;

  StringList TOPPExternalToolBase::getThirdPartyToolLocations()
  {
    return getThirdPartyToolLocations(File::getOpenMSDataPath());
  }

  StringList TOPPExternalToolBase::getThirdPartyToolLocations(const std::string& data_path)
  {
    StringList locations = File::listDirectories(data_path + "/THIRDPARTY");
    for (std::string& location : locations)
    {
      StringUtils::ensureLastChar(location, '/');
    }
    return locations;
  }

  namespace
  {
    // A name with a directory part is meant as it is, and an empty name is not a name.
    bool isPlainFileName(const std::string& exe_filename)
    {
      return !exe_filename.empty() && exe_filename.find_first_of("/\\") == std::string::npos;
    }
  }

  bool TOPPExternalToolBase::findThirdPartyExecutable(std::string& exe_filename)
  {
    // checked here as well, so that such a name does not need the shared-data directory
    if (!isPlainFileName(exe_filename))
    {
      return false;
    }
    return findThirdPartyExecutable(exe_filename, getThirdPartyToolLocations());
  }

  bool TOPPExternalToolBase::findThirdPartyExecutable(std::string& exe_filename, const StringList& locations)
  {
    if (!isPlainFileName(exe_filename))
    {
      return false;
    }
    for (const std::string& location : locations)
    {
      const std::string candidate = location + exe_filename;
      if (File::exists(candidate) && !File::isDirectory(candidate))
      {
        exe_filename = candidate;
        return true;
      }
    }
    return false;
  }

  bool TOPPExternalToolBase::findExecutable_(std::string& executable) const
  {
    // The PATH first, then the third-party tools that ship with OpenMS: the Linux and macOS
    // packages install them under share/OpenMS/THIRDPARTY without putting them on the PATH.
    return File::findExecutable(executable) || findThirdPartyExecutable(executable);
  }

  TOPPBase::ExitCodes TOPPExternalToolBase::runExternalProcess_(const std::string& executable, const std::vector<std::string>& arguments, const std::string& workdir, const std::map<std::string, std::string>& env) const
  {
    std::string proc_stdout, proc_stderr; // collect all output (might be useful if program crashes, see below)
    return runExternalProcess_(executable, arguments, proc_stdout, proc_stderr, workdir, env);
  }

  TOPPBase::ExitCodes TOPPExternalToolBase::runExternalProcess_(const std::string& executable, const std::vector<std::string>& arguments, std::string& proc_stdout, std::string& proc_stderr, const std::string& workdir, const std::map<std::string, std::string>& env) const
  {
    proc_stdout.clear();
    proc_stderr.clear();

    // callbacks: invoked whenever output is available.
    auto lam_out = [&](const std::string& out) { proc_stdout += out; if (debug_level_ >= 4) OPENMS_LOG_INFO << out; };
    auto lam_err = [&](const std::string& out) { proc_stderr += out; if (debug_level_ >= 4) OPENMS_LOG_INFO << out; };
    ExternalProcess ep(lam_out, lam_err);

    const auto& rt = ep.run(executable, arguments, workdir, true, ExternalProcess::IO_MODE::READ_WRITE, env); // does automatic escaping etc... start
    if (debug_level_ < 4 && rt != ExternalProcess::RETURNSTATE::SUCCESS)
    { // error occurred: if not written already in callback, do it now
      writeLogError_("Standard output: " + proc_stdout);
      writeLogError_("Standard error: " + proc_stderr);
    }
    switch (rt)
    {
      case ExternalProcess::RETURNSTATE::SUCCESS:
        return EXECUTION_OK;
      case ExternalProcess::RETURNSTATE::NONZERO_EXIT:
      case ExternalProcess::RETURNSTATE::CRASH:
        return EXTERNAL_PROGRAM_ERROR;
      case ExternalProcess::RETURNSTATE::FAILED_TO_START:
        return EXTERNAL_PROGRAM_NOTFOUND;
      default:
        throw Exception::InternalToolError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown return state of external process.");
    }
  }

} // namespace OpenMS
