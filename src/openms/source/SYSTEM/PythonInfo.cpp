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
#include <sstream>

#if BOOST_VERSION >= 108800
namespace bp = boost::process::v1;
#else
namespace bp = boost::process;
#endif

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
      if (to_path(python_executable).is_relative())
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

    try
    {
      bp::ipstream pipe_out;
      bp::ipstream pipe_err;
      bp::child child(
        bp::search_path(static_cast<std::string>(python_executable)),
        bp::args({"--version"}),
        bp::std_out > pipe_out,
        bp::std_err > pipe_err
      );

      bool finished = child.wait_for(std::chrono::seconds(30));
      if (!finished)
      {
        child.terminate();
        child.wait(); // collect after terminate
        ss << "  Python was found at '" << python_executable << "' but the process timed out (can happen on very busy systems).\n"
           << "  Please free some resources or if you want to run the TOPP tool nevertheless set the TOPP tools 'force' flag in order to avoid this check.\n";
        error_msg = ss.str();
        return false;
      }

      error_msg = ss.str();
      return (child.exit_code() == 0);
    }
    catch (const bp::process_error& /*e*/)
    {
      ss << "  Python found at '" << python_executable << "' but failed to run!\n"
         << "  Make sure you have the rights to execute this binary file.\n";
      error_msg = ss.str();
      return false;
    }
  }

  bool PythonInfo::isPackageInstalled(const String& python_executable, const String& package_name)
  {
    try
    {
      bp::child child(
        bp::search_path(static_cast<std::string>(python_executable)),
        bp::args({"-c", "import " + static_cast<std::string>(package_name)}),
        bp::std_out > bp::null,
        bp::std_err > bp::null
      );

      bool finished = child.wait_for(std::chrono::seconds(30));
      if (!finished)
      {
        child.terminate();
        child.wait();
        return false;
      }
      return (child.exit_code() == 0);
    }
    catch (const bp::process_error&)
    {
      return false;
    }
  }

  String PythonInfo::getVersion(const String& python_executable)
  {
    String v;
    try
    {
      bp::ipstream pipe_out;
      bp::ipstream pipe_err;
      bp::child child(
        bp::search_path(static_cast<std::string>(python_executable)),
        bp::args({"--version"}),
        bp::std_out > pipe_out,
        bp::std_err > pipe_err
      );

      bool finished = child.wait_for(std::chrono::seconds(30));
      if (finished && child.exit_code() == 0)
      {
        std::string line;
        while (std::getline(pipe_out, line))
        {
          v += line;
        }
        while (std::getline(pipe_err, line))
        {
          v += line; // some pythons report version on stderr
        }
        v.trim(); // remove '\n'
      }
      else if (!finished)
      {
        child.terminate();
        child.wait();
      }
    }
    catch (const bp::process_error&)
    {
      // return empty string
    }
    return v;
  }

} // namespace OpenMS
