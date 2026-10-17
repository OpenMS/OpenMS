// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Clemens Groepl, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/APPLICATIONS/OpenMS_CLIConfig.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <map>
#include <string>
#include <vector>

namespace OpenMS
{
  /**
    @brief Base class for TOPP tools that invoke an external executable.

    Adds @ref runExternalProcess_ - a thin adapter over @ref ExternalProcess that maps its
    RETURNSTATE to @ref TOPPBase::ExitCodes and logs stdout/stderr on failure (or when the
    debug level is >= 4). TOPP tools that shell out to an external program should derive from
    this class instead of @ref TOPPBase directly.

    Their executable parameters (tagged @c is_executable) are resolved on the @c PATH first and
    then among the third-party tools that ship with OpenMS (@ref findThirdPartyExecutable). The
    Linux and macOS packages install those under @c share/OpenMS/THIRDPARTY without putting
    them on the @c PATH, so an adapter finds e.g. the bundled Comet without @c -comet_executable.
  */
  class OPENMS_CLI_DLLAPI TOPPExternalToolBase : public TOPPBase
  {
  public:
    /// No default constructor
    TOPPExternalToolBase() = delete;

    /// No default copy constructor.
    TOPPExternalToolBase(const TOPPExternalToolBase&) = delete;

    /// Inherit @ref TOPPBase's constructor (this class adds behaviour, not state).
    using TOPPBase::TOPPBase;

    /// Destructor
    ~TOPPExternalToolBase() override;

    /**
      @brief The folders of the third-party tools that ship with OpenMS

      Every folder of @c THIRDPARTY in the OpenMS shared-data directory
      (File::getOpenMSDataPath()), e.g. @c share/OpenMS/THIRDPARTY/Comet/, sorted by name.
      The installers put the search engines and other tools that OpenMS bundles there; the
      Windows installer also puts these folders on the @c PATH, the Linux and macOS packages
      do not.

      @return The folders, with '/' as separator and ending in '/'. Empty if there are none, as
              in a build tree: its @c THIRDPARTY folder holds only a ReadMe, and the packaging
              fills it.
    */
    static StringList getThirdPartyToolLocations();

    /**
      @brief The folders of the third-party tools in @c THIRDPARTY of the shared-data directory @p data_path

      Same as getThirdPartyToolLocations(), for the given directory (e.g. for testing).
    */
    static StringList getThirdPartyToolLocations(const std::string& data_path);

    /**
      @brief Looks for an executable among the third-party tools that ship with OpenMS

      Searches the folders of getThirdPartyToolLocations(), in that order. Only a plain file
      name such as @c comet.exe or @c sage is looked up; a name with a directory part, such as
      @c /opt/x/percolator or @c sub/tool, is left as it is, and so is an empty name.

      @param[in,out] exe_filename File name of the executable; replaced by its full path if found
      @return true if the executable was found
    */
    static bool findThirdPartyExecutable(std::string& exe_filename);

    /**
      @brief Looks for an executable in the given folders of third-party tools

      Same as findThirdPartyExecutable(std::string&), searching @p locations (folders ending in
      '/', as getThirdPartyToolLocations() returns them) instead.
    */
    static bool findThirdPartyExecutable(std::string& exe_filename, const StringList& locations);

  protected:
    /**
      @brief Searches the PATH (File::findExecutable()), then the third-party tools that ship with OpenMS (findThirdPartyExecutable())

      An executable on the PATH takes precedence, so a user's own version of a bundled tool wins.
    */
    bool findExecutable_(std::string& executable) const override;

    /// Runs an external process via ExternalProcess and prints its stderr output on failure or if debug_level > 4
    ExitCodes runExternalProcess_(const std::string& executable, const std::vector<std::string>& arguments, const std::string& workdir = "", const std::map<std::string, std::string>& env = {}) const;

    /// Runs an external process via ExternalProcess and prints its stderr output on failure or if debug_level > 4
    /// Additionally returns the process' stdout and stderr
    ExitCodes runExternalProcess_(const std::string& executable, const std::vector<std::string>& arguments, std::string& proc_stdout, std::string& proc_stderr, const std::string& workdir = "", const std::map<std::string, std::string>& env = {}) const;
  };

} // namespace OpenMS
