// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Clemens Groepl, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS
{
  class String;

  /**
      @brief Version information class.

      The OpenMS release version and revision data can be retrieved as a string
      or as integers.

      Note that the term <i>"version"</i> refers to releases (such as 1.0,
      1.1.1-alpha, 1.2, ...), which follows the https://semver.org/ definition
      of version numbers using major, minor, patch and pre-release identifiers
      (separated by a dash).
      The term <i>"revision"</i> refers to a revision control system
      such as git and is mainly of interest for developers. The term <i>"branch"</i>
      refers to the git branch that this build of OpenMS is based on.

      The VersionInfo class contains only static methods.

      @ingroup Concept
  */
  class OPENMS_DLLAPI VersionInfo
  {
public:

    struct OPENMS_DLLAPI VersionDetails
    {
      Int version_major = 0;
      Int version_minor = 0;
      Int version_patch = 0;
      String pre_release_identifier;

      VersionDetails() = default;

      /// Copy constructor
      VersionDetails(const VersionDetails & other) = default;

      /// Copy assignment
      VersionDetails& operator=(const VersionDetails& other) = default;

      /**
        @brief parse String and return as proper struct

        @returns VersionInfo::empty on failure
      */
      static VersionDetails create(const String & version);

      bool operator<(const VersionDetails & rhs) const;
      bool operator==(const VersionDetails & rhs) const;
      bool operator!=(const VersionDetails & rhs) const;
      bool operator>(const VersionDetails & rhs) const;

      static const VersionDetails EMPTY; // 0.0.0 version for comparison
    };

    /// Return the build time of OpenMS
    static String getTime();

    /// Return the version number of OpenMS
    static String getVersion();

    /// Return the version number of OpenMS
    static VersionDetails getVersionStruct();

    /**
      @brief Return the revision number from revision control system, e.g. git.

      On released versions of OpenMS, the result is "exported".
      The result can be possibly be "" on some platforms, which means that
      revision info is unavailable.  You should check for both cases in your
      code.

      @internal The current git version is queried by the build system regularly and
      the result is written as a header file which is included by VersionInfo.cpp.
    */
    static String getRevision();

    /**
      @brief Return the branch name from revision control system, e.g. git.

      On released versions of OpenMS the result is "exported".
      The result can be possibly be "" on some platforms, which means that
      revision info is unavailable.  You should check for both cases in your
      code.

      @internal The current git branch is queried by the build system regularly and
      the result is written as a header file which is included by VersionInfo.cpp.
    */
    static String getBranch();

  };

}

