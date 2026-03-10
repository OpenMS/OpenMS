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
#include <OpenMS/SYSTEM/ExternalProcess.h>

#include <filesystem>

namespace OpenMS
{

  bool JavaInfo::canRun(const String& java_executable, bool verbose_on_error)
  {
    ExternalProcess process;
    String error_msg;
    std::vector<std::string> args = {"-version"};
    ExternalProcess::RETURNSTATE result = process.run(static_cast<std::string>(java_executable), args, "", false, error_msg);
    
    bool success = (result == ExternalProcess::RETURNSTATE::SUCCESS);
    if (!success && verbose_on_error)
    {
        OPENMS_LOG_ERROR << "Java-Check:\n";
        if (result == ExternalProcess::RETURNSTATE::FAILED_TO_START)
        {
          OPENMS_LOG_ERROR
            << "  Java not found at '" << java_executable << "'!\n"
            << "  Make sure Java is installed and this location is correct.\n";
          if (std::filesystem::path(static_cast<std::string>(java_executable)).is_relative())
          {
            static String path;
            if (path.empty())
            {
              path = getenv("PATH");
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
        else
        {
          OPENMS_LOG_ERROR << "  Error executing '" << java_executable << "'!\n"
                    << "  Error description: '" << error_msg << "'.\n";
        }
    }
    return success;
  }

} // namespace OpenMS
