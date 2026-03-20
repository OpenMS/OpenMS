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
#include <OpenMS/SYSTEM/PathUtils.h>

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

#include <chrono>
#include <filesystem>

#if BOOST_VERSION >= 108800
namespace bp = boost::process::v1;
#else
namespace bp = boost::process;
#endif

namespace OpenMS
{

  bool JavaInfo::canRun(const String& java_executable, bool verbose_on_error)
  {
    try
    {
      bp::ipstream pipe_out;
      bp::ipstream pipe_err;
      bp::child child(
        bp::search_path(static_cast<std::string>(java_executable)),
        bp::args({"-version"}),
        bp::std_out > pipe_out,
        bp::std_err > pipe_err
      );

      bool finished = child.wait_for(std::chrono::seconds(30));
      if (!finished)
      {
        child.terminate();
        child.wait();
        if (verbose_on_error)
        {
          OPENMS_LOG_ERROR << "Java-Check:\n";
          OPENMS_LOG_ERROR
            << "  Java was found at '" << java_executable << "' but the process timed out (can happen on very busy systems).\n"
            << "  Please free some resources or if you want to run the TOPP tool nevertheless set the TOPP tools 'force' flag in order to avoid this check." << std::endl;
        }
        return false;
      }

      if (child.exit_code() != 0)
      {
        if (verbose_on_error)
        {
          OPENMS_LOG_ERROR << "Java-Check:\n";
          OPENMS_LOG_ERROR << "  Error executing '" << java_executable << "'!\n"
                           << "  Java returned a non-zero exit code (" << child.exit_code() << ").\n";
        }
        return false;
      }

      return true;
    }
    catch (const bp::process_error& /*e*/)
    {
      if (verbose_on_error)
      {
        OPENMS_LOG_ERROR << "Java-Check:\n";
        OPENMS_LOG_ERROR
          << "  Java not found at '" << java_executable << "'!\n"
          << "  Make sure Java is installed and this location is correct.\n";
        if (to_path(java_executable).is_relative())
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
      return false;
    }
  }

} // namespace OpenMS
