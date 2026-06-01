// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow, Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/ToolDescription.h>

using namespace std;

namespace OpenMS
{

  namespace Internal
  {
    // C'Tor with arguments
    ToolDescriptionInternal::ToolDescriptionInternal(const bool p_is_internal, const String& p_name, const String& p_category, const StringList& p_types) :
      is_internal(p_is_internal),
      name(p_name),
      category(p_category),
      types(p_types)
    {
    }

    ToolDescriptionInternal::ToolDescriptionInternal(const String& p_name, const StringList& p_types) :
      
      name(p_name),
      category(),
      types(p_types)
    {
    }

    bool ToolDescriptionInternal::operator==(const ToolDescriptionInternal& rhs) const
    {
      if (this == &rhs)
        return true;

      return is_internal == rhs.is_internal
             && name == rhs.name
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
    ToolDescription::ToolDescription(const String& p_name, const String& p_category, const StringList& p_types) :
      ToolDescriptionInternal(true, p_name, p_category, p_types)
    {
    }

  }

} // namespace OpenMS
