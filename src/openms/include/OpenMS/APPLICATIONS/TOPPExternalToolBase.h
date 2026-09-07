// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Clemens Groepl, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

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
  */
  class OPENMS_DLLAPI TOPPExternalToolBase : public TOPPBase
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

  protected:
    /// Runs an external process via ExternalProcess and prints its stderr output on failure or if debug_level > 4
    ExitCodes runExternalProcess_(const std::string& executable, const std::vector<std::string>& arguments, const std::string& workdir = "", const std::map<std::string, std::string>& env = {}) const;

    /// Runs an external process via ExternalProcess and prints its stderr output on failure or if debug_level > 4
    /// Additionally returns the process' stdout and stderr
    ExitCodes runExternalProcess_(const std::string& executable, const std::vector<std::string>& arguments, std::string& proc_stdout, std::string& proc_stderr, const std::string& workdir = "", const std::map<std::string, std::string>& env = {}) const;
  };

} // namespace OpenMS
