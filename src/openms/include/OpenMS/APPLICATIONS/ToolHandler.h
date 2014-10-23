// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_APPLICATIONS_TOOLHANDLER_H
#define OPENMS_APPLICATIONS_TOOLHANDLER_H

#include <OpenMS/DATASTRUCTURES/ToolDescription.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/DATASTRUCTURES/String.h>

class QStringList;

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
