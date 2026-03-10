// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/PythonInfo.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/ExternalProcess.h>

#include <sstream>
#include <filesystem>

using namespace std;

namespace OpenMS
{
  bool PythonInfo::canRun(String& python_executable, String& error_msg)
  {
    stringstream ss;
    String py_original = python_executable;
    if (!File::findExecutable(python_executable))
    {
      ss << "  Python not found at '" << python_executable << "'!\n"
         << "  Make sure Python is installed and this location is correct.\n";
      if (std::filesystem::path(static_cast<std::string>(python_executable)).is_relative())
      {
        static String path;
        if (path.empty())
        {
          path = getenv("PATH");
        }
        ss << "  You might need to add the Python binary to your PATH variable\n"
           << "  or use an absolute path+filename pointing to Python.\n"
           << "  The current SYSTEM PATH is: '" << path << "'.\n\n";
#ifdef __APPLE__
        ss << "  On MacOSX, application bundles change the system PATH; Open your executable (e.g. KNIME/TOPPAS/TOPPView) from within the bundle (e.g. ./TOPPAS.app/Contents/MacOS/TOPPAS) to preserve the system PATH or use an absolute path to Python!\n";
#endif
      }
      error_msg = ss.str();
      return false;
    }

    if (python_executable != py_original)
    {
      ss << "Python executable ('" << py_original << "') resolved to '" << python_executable << "'\n";
    }

    ExternalProcess process;
    String proc_error_msg;
    std::vector<std::string> args = {"--version"};
    ExternalProcess::RETURNSTATE result = process.run(static_cast<std::string>(python_executable), args, "", false, proc_error_msg);
    
    bool success = (result == ExternalProcess::RETURNSTATE::SUCCESS);
    if (!success)
    {
      if (result == ExternalProcess::RETURNSTATE::FAILED_TO_START)
      {
        ss << "  Python found at '" << python_executable << "' but failed to run!\n"
           << "  Make sure you have the rights to execute this binary file.\n";
      }
      else
      {
        ss << "  Error executing '" << python_executable << "'!\n"
           << "  Error description: '" << proc_error_msg << "'.\n";
      }
    }

    error_msg = ss.str();
    return success;
  }

  bool PythonInfo::isPackageInstalled(const String& python_executable, const String& package_name)
  {
    ExternalProcess process;
    String error_msg;
    std::vector<std::string> args = {"-c", static_cast<std::string>(String("import ") + package_name)};
    ExternalProcess::RETURNSTATE result = process.run(static_cast<std::string>(python_executable), args, "", false, error_msg);
    
    return (result == ExternalProcess::RETURNSTATE::SUCCESS);
  }

  String PythonInfo::getVersion(const String& python_executable)
  {
    String v;
    ExternalProcess process;
    String error_msg;
    String output;
    
    // Set up callbacks to capture output
    process.setCallbacks(
      [&output](const String& stdout_line) { output += stdout_line; },
      [&output](const String& stderr_line) { output += stderr_line; }
    );
    
    std::vector<std::string> args = {"--version"};
    ExternalProcess::RETURNSTATE result = process.run(static_cast<std::string>(python_executable), args, "", false, error_msg);
    
    if (result == ExternalProcess::RETURNSTATE::SUCCESS)
    {
      v = output;
      v.trim(); // remove '\n'
    }
    return v;
  }

} // namespace OpenMS
