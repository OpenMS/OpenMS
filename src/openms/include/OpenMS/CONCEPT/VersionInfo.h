// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
      Int version_major;
      Int version_minor;
      Int version_patch;
      String pre_release_identifier;

      VersionDetails() :
        version_major(0), version_minor(0), version_patch(0), pre_release_identifier("")
      {
      }


      VersionDetails(const VersionDetails & other):
        version_major(other.version_major),
        version_minor(other.version_minor),
        version_patch(other.version_patch),
        pre_release_identifier(other.pre_release_identifier)
      {
      }

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

