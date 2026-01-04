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

#include <cstdio>
#include <array>
#include <sstream>

#ifdef OPENMS_WINDOWSPLATFORM
#include <io.h>
#define popen _popen
#define pclose _pclose
#else
#include <sys/wait.h>
#endif

using namespace std;

namespace OpenMS
{
  namespace
  {
    // Helper function to check if a path is relative
    bool isRelativePath(const String& path)
    {
      if (path.empty()) return true;

#ifdef OPENMS_WINDOWSPLATFORM
      // Windows: absolute paths start with drive letter (C:) or UNC path (\\)
      if (path.size() >= 2 && path[1] == ':') return false; // C:\path
      if (path.size() >= 2 && path[0] == '\\' && path[1] == '\\') return false; // \\server\path
      return true;
#else
      // Unix/Mac: absolute paths start with /
      return path[0] != '/';
#endif
    }

    // Helper function to execute a command and capture output
    bool executeCommand(const String& command, String& output, int& exit_code)
    {
      std::array<char, 128> buffer;
      output.clear();

      FILE* pipe = popen(command.c_str(), "r");
      if (!pipe)
      {
        return false; // Failed to start
      }

      // Read output
      while (fgets(buffer.data(), buffer.size(), pipe) != nullptr)
      {
        output += buffer.data();
      }

      int status = pclose(pipe);
      if (status == -1)
      {
        exit_code = -1; // pclose failed
        return false;
      }

#ifdef OPENMS_WINDOWSPLATFORM
      // On Windows, pclose returns the exit code directly
      exit_code = status;
#else
      // On Unix, decode waitpid-style status
      if (WIFEXITED(status))
      {
        exit_code = WEXITSTATUS(status);
      }
      else
      {
        exit_code = -1; // Process did not exit normally
      }
#endif
      return true;
    }
  }

  bool PythonInfo::canRun(String& python_executable, String& error_msg)
  {
    stringstream ss;
    String py_original = python_executable;
    if (!File::findExecutable(python_executable))
    {
      ss << "  Python not found at '" << python_executable << "'!\n"
         << "  Make sure Python is installed and this location is correct.\n";
      if (isRelativePath(python_executable))
      {
        static String path;
        if (path.empty())
        {
          path = getenv("PATH") ? getenv("PATH") : "";
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

    String command = python_executable + " --version 2>&1";
    String output;
    int exit_code = 0;
    bool started = executeCommand(command, output, exit_code);
    bool success = started && exit_code == 0;

    if (!success)
    {
      if (!started)
      {
        ss << "  Python found at '" << python_executable << "' but failed to run!\n"
           << "  Make sure you have the rights to execute this binary file.\n";
      }
      else if (exit_code != 0)
      {
        ss << "  Error executing '" << python_executable << "'!\n"
           << "  Exit code: " << exit_code << "\n"
           << "  Output: " << output << "\n";
      }
    }

    error_msg = ss.str();
    return success;
  }

  bool PythonInfo::isPackageInstalled(const String& python_executable, const String& package_name)
  {
    String command = python_executable + " -c \"import " + package_name + "\" 2>&1";
    String output;
    int exit_code = 0;
    bool started = executeCommand(command, output, exit_code);
    return (started && exit_code == 0);
  }

  String PythonInfo::getVersion(const String& python_executable)
  {
    String command = python_executable + " --version 2>&1";
    String output;
    int exit_code = 0;
    bool started = executeCommand(command, output, exit_code);

    if (started && exit_code == 0)
    {
      output.trim(); // remove '\n'
      return output;
    }
    return "";
  }

} // namespace OpenMS
