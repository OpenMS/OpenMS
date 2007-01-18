// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <sstream>

#include <OpenMS/METADATA/MetaInfoRegistry.h>
#include <OpenMS/DATASTRUCTURES/String.h>

using namespace std;

namespace OpenMS
{

	MetaInfoRegistry::MetaInfoRegistry(): next_index_(1024), name_to_index_(), index_to_name_(), index_to_description_(), index_to_unit_()
	{
		name_to_index_["isotopic_range"] = 1;
		index_to_name_[1] = "isotopic_range";
		index_to_description_[1] = "consecutive numbering of the peaks in an isotope pattern. 0 is the monoisotopic peak";
		index_to_unit_[1] = "";
		
		name_to_index_["cluster_id"] = 2;
		index_to_name_[2] = "cluster_id";
		index_to_description_[2] = "consecutive numbering of isotope clusters in a spectrum";
		index_to_unit_[2] = "";

		name_to_index_["label"] = 3;
		index_to_name_[3] = "label";
		index_to_description_[3] = "label e.g. shown in visialization";
		index_to_unit_[3] = "";

		name_to_index_["icon"] = 4;
		index_to_name_[4] = "icon";
		index_to_description_[4] = "icon shown in visialization";
		index_to_unit_[4] = "";
		
		name_to_index_["color"] = 5;
		index_to_name_[5] = "color";
		index_to_description_[5] = "color used for visialization e.g. #FF00FF for purple";
		index_to_unit_[5] = "";
	}
	

	MetaInfoRegistry::MetaInfoRegistry(const MetaInfoRegistry& rhs)
	{
		*this = rhs;
	}
	
	MetaInfoRegistry::~MetaInfoRegistry()
	{
		
	}
	
	MetaInfoRegistry& MetaInfoRegistry::operator = (const MetaInfoRegistry& rhs)
	{
		if (this==&rhs) return *this;
		
		next_index_ = rhs.next_index_;
		name_to_index_ = rhs.name_to_index_;
		index_to_name_ = rhs.index_to_name_;
		index_to_description_ = rhs.index_to_description_;
		index_to_unit_ = rhs.index_to_unit_;
		
		return *this;
	}
		

	UnsignedInt MetaInfoRegistry::registerName(const string& name, const string& description, const string& unit) 
	{
		map<string,UnsignedInt>::iterator it = name_to_index_.find(name);
		if (it == name_to_index_.end())
		{
			name_to_index_[name] = next_index_;
			index_to_name_[next_index_] = name;
			index_to_description_[next_index_] = description;
			index_to_unit_[next_index_] = unit;
			return next_index_++;
		}
		else
		{
			return it->second;
		}
	}
	
	UnsignedInt MetaInfoRegistry::getIndex(const string& name) const throw(Exception::InvalidValue)
	{
		map<string,UnsignedInt>::const_iterator it = name_to_index_.find(name);
		if (it != name_to_index_.end())
		{
			return it->second;
		}
		throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Unregistered Name!",name);
	}


	string MetaInfoRegistry::getDescription(UnsignedInt index) const throw(Exception::InvalidValue)
	{
		map<UnsignedInt,string>::const_iterator it = index_to_description_.find(index);
		if (it == index_to_description_.end())
		{
			throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Unregistered index!",String(index));
		}
		return it->second;
	}
 
	string MetaInfoRegistry::getDescription(const std::string& name) const throw(Exception::InvalidValue)
	{
		UnsignedInt index = getIndex(name);
		if  (index==0)
		{
			throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Unregistered Name!",name);
		}
		return (index_to_description_.find(index))->second;
	}

	string MetaInfoRegistry::getUnit(UnsignedInt index) const throw(Exception::InvalidValue)
	{
		map<UnsignedInt,string>::const_iterator it = index_to_unit_.find(index);
		if (it == index_to_description_.end())
		{
			throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Unregistered index!",String(index));
		}
		return it->second;
	}
	
	string MetaInfoRegistry::getUnit(const string& name) const throw(Exception::InvalidValue)
	{
		UnsignedInt index = getIndex(name);
		if  (index==0)
		{
			throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Unregistered Name!",name);
		}
		return (index_to_unit_.find(index))->second;
	}

  string  MetaInfoRegistry::getName(UnsignedInt index) const throw(Exception::InvalidValue)
  {
		map<UnsignedInt,string>::const_iterator it = index_to_name_.find(index);
		if (it != index_to_name_.end())
		{
			return it->second;
		}
		throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Unregistered index!",String(index));
  }
	
} //namespace
