// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/APPLICATIONS/OpenMS_CLIConfig.h>

#include <OpenMS/DATASTRUCTURES/ToolDescription.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>

#include <map>


namespace OpenMS
{
  /**
    @brief Registry of TOPP tools and their TOPPAS categories.

    The registry is @em data, not code: it is assembled from the tab-separated @c .tsv files
    under @ref ToolHandler::getToolRegistryPath (@c [OpenMS share]/TOOLS plus an OS-specific
    subdirectory, @c .../LINUX on Mac and Linux and @c .../WINDOWS on Windows). The search
    path can be augmented through the @c OPENMS_TOOL_REGISTRY_PATH environment variable
    (@c OPENMS_TTD_INTERNAL_PATH, its name in earlier releases, is still read as well). A
    directory that the search reaches twice contributes its files once.

    An entry carries only what the registry is for -- the tool's name and the category
    TOPPAS groups it under, plus the tool's @c -type sub-modes where it has any:

    @code
    # <tool name>	<TOPPAS category>
    FileConverter	File Converter
    @endcode

    Everything else about a tool (its description, parameters, valid formats, citations) is
    self-described by its binary through @c -write_ctd, so a tool does not have to be built
    from the same source tree, or even the same repository, as the library: installing its
    binary and a @c .tsv naming it is enough to register it.

    OpenMS' own tools are registered by @c share/OpenMS/TOOLS/OpenMS.tsv, which the build
    generates from the @c openms_topp_tool() declarations in @c src/topp/executables.cmake
    (see @c cmake/topp_tool_macros.cmake). That file is therefore never edited by hand and
    lists exactly the tools the build produced, so a tool that a build option did not build
    is not registered either. A build tree that was never installed has no
    @c share/OpenMS/TOOLS of its own; there the generated file is read from the build tree.

    Used by TOPPAS for the visual workflow editor's tool palette and by the TOPP runtime to
    look up tools and their categories.

    @note The assembled registry is parsed once per process and cached, so the repeated
          lookups @ref OpenMS::TOPPBase makes during start-up cost one map lookup each.

    @ingroup System
  */

  /// Map: TOPP tool name -> its @ref Internal::ToolDescription (category + per-type configuration).
  typedef std::map<std::string, Internal::ToolDescription> ToolListType;

  class OPENMS_CLI_DLLAPI ToolHandler
  {
public:

    /**
      @brief List every TOPP tool enabled in this build, keyed by tool name.

      Assembled from the @c .tsv files under @ref getToolRegistryPath (a name collision
      between two of them throws Exception::InvalidValue, naming both files). Each value
      carries the tool's category string (e.g. @c "Quantitation", @c "File Converter") used
      by TOPPAS for grouping. Every @c .tsv of the directory is read, so a tool of another
      repository is listed as soon as its file is installed there.

      Prefer @ref getTOPPToolListRef when a copy is not needed.

      @return Map @c toolname -> @ref Internal::ToolDescription for every tool enabled in this build.
    */
    static ToolListType getTOPPToolList();

    /**
      @brief The registry of @ref getTOPPToolList without copying it.

      Same content as @ref getTOPPToolList, served from the process-wide cache. Use this for
      membership tests and lookups; the reference stays valid for the lifetime of the process.

      @return Reference to the cached map @c toolname -> @ref Internal::ToolDescription.
    */
    static const ToolListType& getTOPPToolListRef();

    /**
      @brief Return the alternative "types" / sub-commands a tool supports, or an empty list if it has none.

      Most tools have a single behaviour; a small number expose multiple sub-modes via @c -type
      (e.g. @c FeatureFinderCentroided vs. @c FeatureFinderIsotopeWavelet sharing infrastructure).

      @param[in] toolname Name of the TOPP tool to query.
      @return Type names (may be empty); empty also when the tool is unknown, so that a tool
              which is not in the registry can still write its CTD/CWL description.
    */
    static StringList getTypes(const std::string& toolname);

    /**
      @brief Return the category string of a tool.

      @param[in] toolname Name of the TOPP tool to query.
      @return Category string (e.g. @c "Quantitation") or an empty string if @p toolname is unknown.
    */
    static std::string getCategory(const std::string& toolname);

    /**
      @brief Resolved file-system path of the tool registry directory (root of the @c .tsv search).
      @return @c File::getOpenMSDataPath() + @c "/TOOLS".
    */
    static std::string getToolRegistryPath();

private:

    /// Parse every file of @ref getToolRegistryFiles_ into the registry map.
    /// Throws Exception::InvalidValue when two files register the same tool name.
    static ToolListType loadRegistry_();

    /// Enumerate the @c .tsv files to read: everything under @ref getToolRegistryPath (and its
    /// OS-specific subdirectory), or under the build tree's generated registry when this build
    /// was never installed, plus the directories the environment variables add. A directory the
    /// search reaches more than once contributes its files once.
    static StringList getToolRegistryFiles_();
  };

} // namespace OpenMS
