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

    This class depends on File and Param; File does not depend on it. The former
    File entry points remain as deprecated forwarders.

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
    /// Looks up the following locations, taking the first one which is non-null:
    ///   - environment variable OPENMS_TMPDIR
    ///   - 'temp_dir' in the ~/OpenMS.ini file
    ///   - System temp directory (usually defined by environment 'TMP' or 'TEMP'
    static std::string getTempDirectory();

    /// The current OpenMS user data path (for result files)
    /// Tries to set the user directory in following order:
    ///   1. OPENMS_HOME_DIR if environmental variable set
    ///   2. "home_dir" entry in OpenMS.ini
    ///   3. user home directory
    static std::string getUserDirectory();

    /// get the system's default OpenMS.ini file in the users home directory (&lt;home&gt;/OpenMS/OpenMS.ini)
    /// or create/repair it if required
    /// order:
    ///   1. &lt;OPENMS_HOME_DIR&gt;/OpenMS/OpenMS.ini if environmental variable set
    ///   2. user home directory &lt;home&gt;/OpenMS/OpenMS.ini
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
