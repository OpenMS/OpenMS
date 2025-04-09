// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow, Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/LogStream.h>
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

    void ToolDescription::addExternalType(const String& type, const ToolExternalDetails& details)
    {
      types.push_back(type);
      external_details.push_back(details);
    }

    void ToolDescription::append(const ToolDescription& other)
    {
      // sanity check
      if (is_internal != other.is_internal
         || name != other.name
          //|| category != other.category
         || (is_internal && !external_details.empty())
         || (other.is_internal && !other.external_details.empty())
         || (!is_internal && external_details.size() != types.size())
         || (!other.is_internal && other.external_details.size() != other.types.size())
          )
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Extending (external) ToolDescription failed!", "");
      }

      // append types and external information
      types.insert(types.end(), other.types.begin(), other.types.end());
      external_details.insert(external_details.end(), other.external_details.begin(), other.external_details.end());

      // check that types are unique
      std::set<String> unique_check;
      unique_check.insert(types.begin(), types.end());
      if (unique_check.size() != types.size())
      {
        OPENMS_LOG_ERROR << "A type appears at least twice for the TOPP tool '" << name << "'. Types given are '" << ListUtils::concatenate(types, ", ") << "'\n";
        if (name == "GenericWrapper")
        {
          OPENMS_LOG_ERROR << "Check the .ttd files in your share/ folder and remove duplicate types!\n";
        }
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "see above!", "");
      }
    }
  }

  Internal::ToolDescription bla;

} // namespace OpenMS
