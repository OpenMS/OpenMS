// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

namespace OpenMS
{
  class String;
  /**
    @brief Detect Python and retrieve information.

    Similar classes exist for other external tools, e.g. JavaInfo .

    @ingroup System
  */
  class OPENMS_DLLAPI PythonInfo
  {
  public:
    /**
      @brief Determine if Python is installed and executable

      The call fails if either Python is not installed or if a relative location is given and Python is not on the search PATH.
      If Python is found, the executable name will be modified to the absolute path.
      If Python is not found, an error message will be put into @p error_msg

      @param python_executable Path to Python executable. Can be absolute, relative or just a filename
      @param error_msg On error, contains detailed error description (e.g. 
      @return Returns false if Python executable can not be called; true if Python executable can be executed
    **/
    static bool canRun(String& python_executable, String& error_msg);


    /**
     @brief Determine if the Python given in @p python_executable has the package @p package_name already installed

     If Python cannot be found, the function will just return false.
     Thus, make sure that PythonInfo::canRun() succeeds before calling this function.

     @param python_executable As determined by canRun()...
     @param package_name The package you want to test (mind lower/upper case!)
     @return true if package is installed
    */
    static bool isPackageInstalled(const String& python_executable, const String& package_name);

    /**
     @brief Determine the version of Python given in @p python_executable by calling '--version'

     If Python cannot be found, the function will return the empty string.
     Thus, make sure that PythonInfo::canRun() succeeds before calling this function.

     @param python_executable As determined by canRun()...
     @return the output of 'python --version'
    */
    static String getVersion(const String& python_executable);
  };

}

