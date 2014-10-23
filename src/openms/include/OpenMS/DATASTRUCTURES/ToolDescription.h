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
// $Maintainer: $
// $Authors: Chris Bielow, Mathias Walzer $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_TOOLDESCRIPTION_H
#define OPENMS_DATASTRUCTURES_TOOLDESCRIPTION_H

#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <map>

#include <OpenMS/OpenMSConfig.h>

namespace OpenMS
{

  namespace Internal
  {
    /**
      @brief Maps input/output files to filenames for the external program

    */
    struct FileMapping
    {
      String location; // a regex/macro mix; to be expanded by tool;
      String target;   // TOPP parameter that determines the desired name
      // thus: move location -> target

      FileMapping & operator=(const FileMapping & rhs);
    };

    /**
      @brief Filename mappings for all input/output files

    */
    struct MappingParam
    {
      std::map<Int, String> mapping;
      std::vector<FileMapping> pre_moves;
      std::vector<FileMapping> post_moves;

      OPENMS_DLLAPI MappingParam & operator=(const MappingParam & rhs);

    };

    /**
        @brief ToolDescription Class.

        This class represents a ToolDescription.

        @ingroup Datastructures
    */
    struct OPENMS_DLLAPI ToolDescriptionInternal
    {
      bool is_internal;
      String name;
      String category;
      StringList types; // -types of the tool (if any, e.g. ['centroided','wavelet'])

      // default C'Tor
      ToolDescriptionInternal();

      // C'Tor with arguments
      ToolDescriptionInternal(const bool p_is_internal, const String & p_name, const String & p_category, const StringList & p_types);

      // short C'Tor
      ToolDescriptionInternal(const String & p_name, const StringList & p_types);

      ToolDescriptionInternal & operator=(const ToolDescriptionInternal & rhs);

      bool operator==(const ToolDescriptionInternal & rhs) const;

      bool operator<(const ToolDescriptionInternal & rhs) const;
    };

    struct OPENMS_DLLAPI ToolExternalDetails
    {
      String text_startup;
      String text_fail;
      String text_finish;
      String category;
      String commandline;
      String path; //< filename to external tool
      String working_directory; //< folder where the command will be executed from
      MappingParam tr_table;
      Param param;
    };

    /**
      Used for internal and external tools
    */
    struct OPENMS_DLLAPI ToolDescription :
      ToolDescriptionInternal
    {
      // additional details for external tools (one entry for each 'type')
      std::vector<ToolExternalDetails> external_details;

      // default CTor
      ToolDescription();

      // C'Tor for internal TOPP tools
      ToolDescription(const String & p_name, const String & p_category, const StringList & p_types = StringList());

      void addExternalType(const String & type, const ToolExternalDetails & details);

      void append(const ToolDescription & other);

      ToolDescription & operator=(const ToolDescription & rhs);
    };

  } // namespace Internal

} // namespace OPENMS

#endif // OPENMS_DATASTRUCTURES_TOOLDESCRIPTION_H
