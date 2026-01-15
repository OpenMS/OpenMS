// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/JavaInfo.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/SYSTEM/File.h>

#include <cstdio>
#include <array>
#include <memory>

#ifdef OPENMS_WINDOWSPLATFORM
#include <io.h>
#define popen _popen
#define pclose _pclose
#else
#include <sys/wait.h>
#endif

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

  bool JavaInfo::canRun(const String& java_executable, bool verbose_on_error)
  {
    // Build command: java_executable -version
    String command = java_executable + " -version 2>&1"; // Redirect stderr to stdout

    String output;
    int exit_code = 0;
    bool started = executeCommand(command, output, exit_code);

    bool success = started && exit_code == 0;

    if (!success && verbose_on_error)
    {
        OPENMS_LOG_ERROR << "Java-Check:\n";
        if (!started)
        {
          OPENMS_LOG_ERROR
            << "  Java not found at '" << java_executable << "'!\n"
            << "  Make sure Java is installed and this location is correct.\n";
          if (isRelativePath(java_executable))
          {
            static String path;
            if (path.empty())
            {
              path = getenv("PATH") ? getenv("PATH") : "";
            }
            OPENMS_LOG_ERROR << "  You might need to add the Java binary to your PATH variable\n"
              << "  or use an absolute path+filename pointing to Java.\n"
              << "  The current SYSTEM PATH is: '" << path << "'.\n\n"
  #ifdef __APPLE__
              << "  On MacOSX, application bundles change the system PATH; Open your executable (e.g. KNIME/TOPPAS/TOPPView) from within the bundle (e.g. ./TOPPAS.app/Contents/MacOS/TOPPAS) to preserve the system PATH or use an absolute path to Java!\n"
  #endif
              << std::endl;
          }
          else
          {
            OPENMS_LOG_ERROR << "  You gave an absolute path to Java. Please check if it's correct.\n"
              << "  You can also try 'java' if your system path is correctly configured.\n"
              << std::endl;
          }
        }
        else if (exit_code != 0)
        {
          OPENMS_LOG_ERROR << "  Error executing '" << java_executable << "'!\n"
                    << "  Exit code: " << exit_code << "\n"
                    << "  Output: " << output << std::endl;
        }
    }
    return success;
  }

} // namespace OpenMS
