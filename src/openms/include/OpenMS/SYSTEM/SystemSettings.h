// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow, Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <string>

namespace OpenMS
{
  class Param;

  /**
    @brief Per-user OpenMS configuration: the OpenMS.ini system parameters and the
    home, temporary and database directories derived from them.

    Every lookup honours the environment overrides (OPENMS_HOME_PATH, OPENMS_TMPDIR,
    XDG_CONFIG_HOME) before consulting OpenMS.ini and the platform defaults.

    This class depends on File and Param; File does not depend on it. These
    lookups used to be static members of File.

    @ingroup System
  */
  class OPENMS_DLLAPI SystemSettings
  {
public:
    /// Returns the OpenMS home path (environment variable overwrites the default home path)
    static std::string getOpenMSHomePath();

    /// @brief Returns the per-user OpenMS configuration directory (the directory that holds OpenMS.ini)
    ///
    /// Follows the XDG base directory specification on unix-like systems
    /// (&lt;XDG_CONFIG_HOME&gt;/OpenMS or &lt;home&gt;/.config/OpenMS) and uses &lt;home&gt;/.OpenMS otherwise.
    /// The returned path has no trailing separator and the directory is not guaranteed to exist.
    /// @return String containing the per-user OpenMS configuration directory path
    static std::string getOpenMSConfigDir();

    /// The current OpenMS temporary data path (for temporary files).
    /// Looks up the following locations, taking the first one which is set:
    ///   - the OPENMS_TMPDIR environment variable
    ///   - a non-empty 'temp_dir' entry in the OpenMS.ini read by getSystemParameters()
    ///   - the system temp directory (usually defined by environment 'TMP' or 'TEMP')
    /// The path is returned as configured; it is not checked for existence.
    static std::string getTempDirectory();

    /// The current OpenMS user data path (for result files), with a trailing '/'.
    /// Resolved in the following order, taking the first one which is set:
    ///   1. the OPENMS_HOME_PATH environment variable
    ///   2. a non-empty "home_dir" entry in the OpenMS.ini read by getSystemParameters()
    ///   3. the user's home directory (HOME, or USERPROFILE on Windows)
    static std::string getUserDirectory();

    /// Returns the OpenMS.ini system parameters, read from getOpenMSConfigDir() + "/OpenMS.ini".
    /// If that file does not exist, the built-in defaults are returned. If it exists but its
    /// 'version' is missing or outdated, missing entries are filled in from the defaults (the
    /// file itself is not rewritten).
    static Param getSystemParameters();

    /// uses File::find() to search for a file names @p db_name
    /// in the 'id_db_dir' param of the OpenMS system parameters
    /// @exception FileNotFound is thrown, if the file is not found
    static std::string findDatabase(const std::string& db_name);

private:
    /// get defaults for the system's Temp-path, user home directory etc.
    static Param getSystemParameterDefaults_();
  };
}
