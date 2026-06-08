// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/ToolDescription.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>

#include <map>


namespace OpenMS
{
  /**
    @brief Registry of TOPP tools and their (TOPPAS / KNIME) categories.

    Tool descriptions come from two layers:

      -# A hard-coded list of @em official TOPP tools maintained in @ref ToolHandler::getTOPPToolList
         (one @ref Internal::ToolDescription per tool, including its KNIME-style category string).
         This is what @ref ToolHandler::getTOPPToolList, @ref ToolHandler::getTypes
         and @ref ToolHandler::getCategory consult.
      -# Optional @em internal-tool @c .ttd config files loaded lazily by
         @c getInternalToolConfigFiles_ from @c [OpenMS]/share/TOOLS/INTERNAL plus an
         OS-specific subdirectory (@c .../LINUX on Mac and Linux, @c .../WINDOWS on Windows).
         The search path can be augmented via the @c OPENMS_TTD_INTERNAL_PATH environment
         variable.

    @note @ref ToolHandler::getExternalToolsPath returns the @c [OpenMS]/share/TOOLS/EXTERNAL location, but
          this directory is not currently scanned for @c .ttd entries by this class — only the
          @c TOOLS/INTERNAL tree is loaded.

    Used by TOPPAS for the visual workflow editor's tool palette and by the TOPP runtime to
    look up tools and their categories.

    @note The hard-coded list (rather than driving it from files at load time) is a deliberate
          trade-off: adding a tool requires a recompile anyway, and parsing all share files on
          library load would slow start-up.

    @ingroup System
  */

  /// Map: TOPP tool name -> its @ref Internal::ToolDescription (category + per-type configuration).
  typedef std::map<String, Internal::ToolDescription> ToolListType;

  class OPENMS_DLLAPI ToolHandler
  {
public:

    /**
      @brief List every official TOPP tool shipped with this OpenMS release, keyed by tool name.

      Includes only the hard-coded "internal" tool registry — external @c .ttd entries are not
      merged here. Each value carries the tool's KNIME-style category string (e.g.
      @c "Quantitation", @c "File Converter") used by TOPPAS for grouping.

      @return Map @c toolname -> @ref Internal::ToolDescription for every tool.
    */
    static ToolListType getTOPPToolList();

    /**
      @brief Return the alternative "types" / sub-commands a tool supports, or an empty list if it has none.

      Most tools have a single behaviour; a small number expose multiple sub-modes via @c -type
      (e.g. @c FeatureFinderCentroided vs. @c FeatureFinderIsotopeWavelet sharing infrastructure).

      @param[in] toolname Name of the TOPP tool to query.
      @return Type names (may be empty); empty also when the tool is unknown.
    */
    static StringList getTypes(const String& toolname);

    /**
      @brief Return the KNIME-style category string of a tool.

      @param[in] toolname Name of the TOPP tool to query.
      @return Category string (e.g. @c "Quantitation") or an empty string if @p toolname is unknown.
    */
    static String getCategory(const String& toolname);

    /**
      @brief Resolved file-system path of the external-tool config directory.
      @return @c File::getOpenMSDataPath() + @c "/TOOLS/EXTERNAL".
    */
    static String getExternalToolsPath();

    /**
      @brief Resolved file-system path of the internal-tool config directory (root of the @c .ttd search).
      @return @c File::getOpenMSDataPath() + @c "/TOOLS/INTERNAL".
    */
    static String getInternalToolsPath();

private:

    /// Lazily load the internal tool registry from @c .ttd files under @ref getInternalToolsPath
    static std::vector<Internal::ToolDescription> getInternalTools_();

    /// Enumerate @c .ttd config file paths under @ref getInternalToolsPath (and its OS-specific subdir)
    static StringList getInternalToolConfigFiles_();

    /// Parse the @c .ttd config files (idempotent — sets @c tools_internal_loaded_ on success)
    static void loadInternalToolConfig_();

    /// Cached internal tool registry; populated lazily by @ref loadInternalToolConfig_
    static std::vector<Internal::ToolDescription> tools_internal_;

    /// Whether @c tools_internal_ has been initialised this run
    static bool tools_internal_loaded_;
  };

} // namespace OpenMS
