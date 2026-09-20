// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow, Mathias Walzer $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/DATASTRUCTURES/TypeAliases.h>

#include <OpenMS/OpenMSConfig.h>

namespace OpenMS
{
  namespace Internal
  {
    /**
        @brief What the TOPP tool registry knows about one tool.

        The name ToolHandler looks it up by, the category TOPPAS groups it under, and the
        @c -type sub-modes it offers. Everything else about a tool comes from the tool's own
        binary (@c -write_ctd), not from here.

        @ingroup Datastructures
    */
    struct OPENMS_DLLAPI ToolDescriptionInternal
    {
      std::string name;
      std::string category;
      StringList types; ///< -types of the tool

      /// default C'Tor
      ToolDescriptionInternal() = default;

      /// C'Tor with arguments
      ToolDescriptionInternal(const std::string& p_name, const std::string& p_category, const StringList& p_types);

      /// short C'Tor
      ToolDescriptionInternal(const std::string& p_name, const StringList& p_types);

      /// Copy assignment
      ToolDescriptionInternal& operator=(const ToolDescriptionInternal& rhs) = default;

      bool operator==(const ToolDescriptionInternal& rhs) const;

      bool operator<(const ToolDescriptionInternal& rhs) const;
    };

    /**
      A tool as the registry describes it.
    */
    struct OPENMS_DLLAPI ToolDescription :
      ToolDescriptionInternal
    {
      /// default CTor
      ToolDescription() = default;

      /// Copy C'Tor
      ToolDescription(const ToolDescription& other) = default;

      /// C'Tor from a registry entry
      ToolDescription(const std::string& p_name, const std::string& p_category, const StringList& p_types = StringList());

      /// Copy assignment
      ToolDescription& operator=(const ToolDescription& rhs) = default;

    };
  } // namespace Internal
} // namespace OPENMS
