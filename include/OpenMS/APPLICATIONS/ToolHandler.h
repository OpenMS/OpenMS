// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_APPLICATIONS_TOOLHANDLER_H
#define OPENMS_APPLICATIONS_TOOLHANDLER_H

#include <OpenMS/DATASTRUCTURES/ToolDescription.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS
{
  /**
     @brief Handles lists of TOPP and UTILS tools and their categories (for TOPPAS)

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

    We could create the list of TOPP tools and UTILS from a set of files in /share
    instead of hard coding them here.
    Advantage: - no recompile if new tool is added (but the new tool will necessitate that anyway)
               - quickly changing a tool's category (e.g. from PreProcessing to Quantitation) and thus its place in TOPPAS
                   even users could rearrange tools themselves...
    Disadvantage:
               - when to library loads, we'd need to parse all the files. Making our start-up time even longer...
               - when files are broken/missing, we will have a hard time initializing the lib
  */

  /// map each TOPP/UTIL to its ToolDescription
  typedef Map<String, Internal::ToolDescription> ToolListType;

  class OPENMS_DLLAPI ToolHandler
  {
public:

    /// Returns the list of official TOPP tools contained in the OpenMS/TOPP release.
    static ToolListType getTOPPToolList(const bool includeGenericWrapper = false);

    /// Returns the list of official UTIL tools contained in the OpenMS/TOPP release.
    static ToolListType getUtilList();

    /// get all types of a tool (empty if none)
    static StringList getTypes(const String & toolname);

    /// Returns the category string from TOPP or UTIL tools
    /// @return empty string if tool was not found
    static String getCategory(const String & toolname);

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

#endif //OPENMS_APPLICATIONS_TOOLHANDLER_H
