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

    /**
      @brief Parsed semver-style version: @c major.minor.patch with an optional pre-release identifier.

      Result type of @ref VersionInfo::getVersionStruct (and of @ref create when called with
      a free-form version string). Default-constructed instances compare equal to
      @ref EMPTY and represent both the "0.0.0 with no pre-release" version and the
      parse-failure sentinel returned by @ref create.
    */
    struct OPENMS_DLLAPI VersionDetails
    {
      Int version_major = 0;             ///< Major number ahead of the first @c '.'.
      Int version_minor = 0;             ///< Minor number between the first and (optional) second @c '.'.
      Int version_patch = 0;             ///< Patch number after the second @c '.'; left at @c 0 when the version string had only two components.
      String pre_release_identifier;     ///< Pre-release suffix after a trailing @c '-' (everything to the end of the string); empty when no @c '-' is present.

      /// Default-construct to @c 0.0.0 with an empty pre-release identifier (equivalent to @ref EMPTY).
      VersionDetails() = default;

      /// Copy constructor.
      VersionDetails(const VersionDetails & other) = default;

      /// Copy assignment.
      VersionDetails& operator=(const VersionDetails& other) = default;

      /**
        @brief Parse a semver-style @c major[.minor[.patch[-pre_release]]] string into the struct.

        Splits on @c '.' and the optional trailing @c '-':
          - At least one @c '.' is required — strings without a dot return @ref EMPTY.
          - The major number must parse as an integer.
          - The minor number must parse as an integer; if no second @c '.' follows, @c patch
            is left at @c 0 and parsing succeeds.
          - The patch number must parse as an integer; if no trailing @c '-' is present,
            @c pre_release_identifier is left empty and parsing succeeds.
          - Everything after the trailing @c '-' (verbatim, no further parsing) becomes
            @c pre_release_identifier.

        Any integer conversion failure (caught @c OpenMS::Exception::ConversionError) yields
        @ref EMPTY; the function does not propagate the exception.

        @param[in] version Version string in @c "X.Y[.Z[-PRE]]" form.
        @return Parsed struct on success, @ref EMPTY on any of the failure paths above.
      */
      static VersionDetails create(const String & version);

      /**
        @brief Compare two versions in @em descending pre-release order.

        Compares lexicographically on @c (major, minor, patch). Tie-break on the pre-release
        identifier follows the convention that "pre-release < release of the same triple",
        but only checks presence — i.e. @c "1.0.0-alpha @c < @c 1.0.0" holds, while
        @c "1.0.0-alpha @c < @c 1.0.0-beta" is @c false (both have a pre-release, so the
        triples are compared as equal).

        @param[in] rhs Version to compare against.
        @return @c true if @c *this is strictly less than @p rhs by the rule above.
      */
      bool operator<(const VersionDetails & rhs) const;
      /**
        @brief Field-wise equality on all four members, including string equality on the pre-release identifier.

        @param[in] rhs Version to compare against.
        @return @c true if every field of @c *this equals the corresponding field of @p rhs.
      */
      bool operator==(const VersionDetails & rhs) const;
      /**
        @brief Field-wise inequality (the logical negation of @ref operator==).

        @param[in] rhs Version to compare against.
        @return @c true if at least one field differs.
      */
      bool operator!=(const VersionDetails & rhs) const;
      /**
        @brief Equivalent to @c !(*this @c < @c rhs @c || @c *this @c == @c rhs); inherits the @c operator< caveat about pre-release ties.

        @param[in] rhs Version to compare against.
        @return @c true if @c *this is strictly greater than @p rhs.
      */
      bool operator>(const VersionDetails & rhs) const;

      /// Sentinel @c 0.0.0 version used both as the default-constructed value and as the parse-failure return of @ref create.
      static const VersionDetails EMPTY;
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

