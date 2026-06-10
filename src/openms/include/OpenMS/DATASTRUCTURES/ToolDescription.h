// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow, Mathias Walzer $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <OpenMS/OpenMSConfig.h>

#include <map>

namespace OpenMS
{
  namespace Internal
  {
    /**
      @brief Maps input/output files to filenames for the external program

    */
    struct FileMapping
    {
      std::string location; ///< a regex/macro mix; to be expanded by tool;
      std::string target; ///< TOPP parameter that determines the desired name
      // thus: move location -> target

      /// Default constructor
      FileMapping() = default;

      /// Copy constructor
      FileMapping(const FileMapping& other) = default;

      /// Copy assignment
      FileMapping& operator=(const FileMapping& rhs) = default;
    };

    /**
      @brief Filename mappings for all input/output files

    */
    struct MappingParam
    {
      std::map<Int, std::string> mapping;
      std::vector<FileMapping> pre_moves;
      std::vector<FileMapping> post_moves;

      /// Default constructor
      MappingParam() = default;

      /// Copy constructor
      MappingParam(const MappingParam& other) = default;

      /// Copy assignment
      MappingParam& operator=(const MappingParam& other) = default;  
    };

    /**
        @brief ToolDescription Class.

        This class represents a ToolDescription.

        @ingroup Datastructures
    */
    struct OPENMS_DLLAPI ToolDescriptionInternal
    {
      bool is_internal = false;
      std::string name;
      std::string category;
      StringList types; ///< -types of the tool

      /// default C'Tor
      ToolDescriptionInternal() = default;

      /// C'Tor with arguments
      ToolDescriptionInternal(const bool p_is_internal, const std::string& p_name, const std::string& p_category, const StringList& p_types);

      /// short C'Tor
      ToolDescriptionInternal(const std::string& p_name, const StringList& p_types);

      /// Copy assignment
      ToolDescriptionInternal& operator=(const ToolDescriptionInternal& rhs) = default;

      bool operator==(const ToolDescriptionInternal& rhs) const;

      bool operator<(const ToolDescriptionInternal& rhs) const;
    };

    struct OPENMS_DLLAPI ToolExternalDetails
    {
      std::string text_startup;
      std::string text_fail;
      std::string text_finish;
      std::string category;
      std::string commandline;
      std::string path; ///< filename to external tool
      std::string working_directory; ///< folder where the command will be executed from
      MappingParam tr_table;
      Param param;
    };

    /**
      Used for internal and external tools
    */
    struct OPENMS_DLLAPI ToolDescription :
      ToolDescriptionInternal
    {
      /// additional details for external tools (one entry for each 'type')
      std::vector<ToolExternalDetails> external_details;

      /// default CTor
      ToolDescription() = default;

      /// Copy C'Tor
      ToolDescription(const ToolDescription& other) = default;

      /// C'Tor for internal TOPP tools
      ToolDescription(const std::string& p_name, const std::string& p_category, const StringList& p_types = StringList());

      /// Copy assignment
      ToolDescription& operator=(const ToolDescription& rhs) = default;

    };
  } // namespace Internal
} // namespace OPENMS
