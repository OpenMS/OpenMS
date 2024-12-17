// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/ToolDescription.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <map>

#include <QtCore/qcontainerfwd.h> // for QStringList


namespace OpenMS
{
  /**
     @brief Handles lists of TOPP tools and their categories (for TOPPAS)

     Path's were *.ttd files are searched for:

     Default:
       The OpenMS share directory ([OpenMS]/share/TOOLS/EXTERNAL)
       OS specific directories
        - [OpenMS]/share/TOOLS/EXTERNAL/LINUX   (for Mac and Linux)
        - [OpenMS]/share/TOOLS/EXTERNAL/WINDOWS (for Windows)
     Environment:
       OPENMS_TTD_PATH  (use only one path here!)
  */

  /*
    internal details & discussion:

    We could create the list of TOPP tools from a set of files in /share
    instead of hard coding them here.
    Advantage: - no recompile if new tool is added (but the new tool will necessitate that anyway)
               - quickly changing a tool's category (e.g. from PreProcessing to Quantitation) and thus its place in TOPPAS
                   even users could rearrange tools themselves...
    Disadvantage:
               - when to library loads, we'd need to parse all the files. Making our start-up time even longer...
               - when files are broken/missing, we will have a hard time initializing the lib
  */

  /// map each tool to its ToolDescription
  typedef std::map<String, Internal::ToolDescription> ToolListType;

  class OPENMS_DLLAPI ToolHandler
  {
public:

    /// Returns the list of official TOPP tools contained in the OpenMS/TOPP release.
    static ToolListType getTOPPToolList(const bool includeGenericWrapper = false);

    /// get all types of a tool (empty if none)
    static StringList getTypes(const String& toolname);

    /// Returns the category string from TOPP tools
    /// @return empty string if tool was not found
    static String getCategory(const String& toolname);

    /// get getOpenMSDataPath() + "/TOOLS/EXTERNAL"
    static String getExternalToolsPath();

    /// get File::getOpenMSDataPath() + "/TOOLS/INTERNAL"
    static String getInternalToolsPath();

private:

    static Internal::ToolDescription getExternalTools_();
    static QStringList getExternalToolConfigFiles_();
    static void loadExternalToolConfig_();
    static Internal::ToolDescription tools_external_;
    static bool tools_external_loaded_;

    static std::vector<Internal::ToolDescription> getInternalTools_();
    static QStringList getInternalToolConfigFiles_();
    static void loadInternalToolConfig_();
    static std::vector<Internal::ToolDescription> tools_internal_;
    static bool tools_internal_loaded_;
  };

} // namespace OpenMS

