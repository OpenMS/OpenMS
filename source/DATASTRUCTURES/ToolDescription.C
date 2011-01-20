// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: $
// $Authors: Chris Bielow, Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/ToolDescription.h>

using namespace std;

namespace OpenMS
{

  namespace Internal
  {

    // default C'Tor
    ToolDescriptionInternal::ToolDescriptionInternal()
      : is_internal(false),
        name(),
        category(),
        types()
    {
    }

    // C'Tor with arguments
    ToolDescriptionInternal::ToolDescriptionInternal(const bool p_is_internal, const String& p_name, const String& p_category, const StringList& p_types)
      : is_internal(p_is_internal),
        name(p_name),
        category(p_category),
        types(p_types)
    {
    }

    ToolDescriptionInternal& ToolDescriptionInternal::operator=(const ToolDescriptionInternal& rhs)
    {
      if (this==&rhs) return *this;
      
      is_internal = rhs.is_internal;
      name = rhs.name;
      category = rhs.category;
      types = rhs.types;
      return *this;
    }




    // default CTor
    ToolDescription::ToolDescription()
      : external_details()
    {
    }

    // C'Tor for internal TOPP tools
    ToolDescription::ToolDescription(const String& p_name, const String& p_category, const StringList& p_types)
      : ToolDescriptionInternal(true, p_name, p_category, p_types)
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
      if ( is_internal != other.is_internal
        || name != other.name
        //|| category != other.category
        || (is_internal && external_details.size()>0)
        || (other.is_internal && other.external_details.size()>0)
        || (!is_internal && external_details.size() != types.size())
        || (!other.is_internal && other.external_details.size() != other.types.size())
        )
      {
        throw Exception::InvalidValue(__FILE__,__LINE__, __PRETTY_FUNCTION__, "Extending (external) ToolDescription failed!", "");
      }
      
      // append types and external information
      types.insert(types.end(), other.types.begin(), other.types.end());
      external_details.insert(external_details.end(), other.external_details.begin(), other.external_details.end());

      // check that types are unique
      std::set<String> unique_check;
      unique_check.insert(types.begin(),types.end());
      if (unique_check.size() != types.size())
      {
        LOG_ERROR << "A type appears at least twice for the TOPP/UTIL '" << name << "'. Types given are '" << types.concatenate(", ") << "'\n";
        if (name =="GenericWrapper") LOG_ERROR << "Check the .ttd files in your share/ folder and remove duplicate types!\n";
        Exception::InvalidValue(__FILE__,__LINE__, __PRETTY_FUNCTION__, "see above!", "");
      }
    }

    ToolDescription& ToolDescription::operator=(const ToolDescription& rhs)
    {
      if (this==&rhs) return *this;
      
      ToolDescriptionInternal::operator=(rhs);
      external_details = rhs.external_details;
      return *this;
    }
  }

  Internal::ToolDescription bla;
	
} // namespace OpenMS


