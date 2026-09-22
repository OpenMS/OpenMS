// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow, Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/ToolDescription.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

using namespace std;

namespace OpenMS
{

  namespace Internal
  {
    // C'Tor with arguments
    ToolDescriptionInternal::ToolDescriptionInternal(const std::string& p_name, const std::string& p_category, const StringList& p_types) :
      name(p_name),
      category(p_category),
      types(p_types)
    {
    }

    ToolDescriptionInternal::ToolDescriptionInternal(const std::string& p_name, const StringList& p_types) :
      
      name(p_name),
      category(),
      types(p_types)
    {
    }

    bool ToolDescriptionInternal::operator==(const ToolDescriptionInternal& rhs) const
    {
      if (this == &rhs)
        return true;

      return name == rhs.name
             && category == rhs.category
             && types == rhs.types;
    }

    bool ToolDescriptionInternal::operator<(const ToolDescriptionInternal& rhs) const
    {
      if (this == &rhs)
        return false;

      return name + "." + ListUtils::concatenate(types, ",") < rhs.name + "." + ListUtils::concatenate(rhs.types, ",");
    }
    
    // C'Tor for internal TOPP tools
    ToolDescription::ToolDescription(const std::string& p_name, const std::string& p_category, const StringList& p_types) :
      ToolDescriptionInternal(p_name, p_category, p_types)
    {
    }

  }

} // namespace OpenMS
