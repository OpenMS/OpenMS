// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/MetaInfo.h>

using namespace std;

namespace OpenMS
{
	
	MetaInfoRegistry MetaInfo::registry_ = MetaInfoRegistry();
	
	MetaInfo::MetaInfo()
	{
		
	}
	
	MetaInfo::MetaInfo(const MetaInfo& rhs)
	{
		*this = rhs;
	}
	
	MetaInfo::~MetaInfo()
	{
		
	}
	
	MetaInfo& MetaInfo::operator = (const MetaInfo& rhs)
	{
		if (this==&rhs) return *this;
		
		index_to_value_ = rhs.index_to_value_;
		
		return *this;
	}

  bool MetaInfo::operator== (const MetaInfo& rhs) const
  {
  	return index_to_value_ == rhs.index_to_value_;
  }
  
  bool MetaInfo::operator!= (const MetaInfo& rhs) const
  {
  	return !(operator==(rhs));
 	}
	
	const DataValue& MetaInfo::getValue(const String& name) const
	{
		map<UInt,DataValue>::const_iterator it = index_to_value_.find(registry_.getIndex(name));
		if (it != index_to_value_.end())
		{
			return it->second;
		}
		return DataValue::EMPTY;
	}
	
	const DataValue& MetaInfo::getValue(UInt index) const
	{
		map<UInt,DataValue>::const_iterator it = index_to_value_.find(index);
		if (it != index_to_value_.end())
		{
			return it->second;
		}
		return DataValue::EMPTY;
	}

	void MetaInfo::setValue(const String& name, const DataValue& value)
	{
		index_to_value_[registry_.getIndex(name)] = value;
	}
	
	void MetaInfo::setValue(UInt index, const DataValue& value)
	{
		index_to_value_[index] = value;
	}	
	
	MetaInfoRegistry& MetaInfo::registry()
	{
		return registry_;
	}

	bool MetaInfo::exists(const String& name) const
	{
		try
		{
			if (index_to_value_.find(registry_.getIndex(name))==index_to_value_.end())
			{
				return false;
			}
		}
		catch (Exception::InvalidValue)
		{
			return false;
		}
		return true;		
	}
	
	bool MetaInfo::exists(UInt index) const
	{
		if (index_to_value_.find(index)==index_to_value_.end())
		{
			return false;
		}
		return true;
	}

	void MetaInfo::removeValue(const String& name)
	{
		map<UInt,DataValue>::iterator it = index_to_value_.find(registry_.getIndex(name));
		if (it != index_to_value_.end())
		{
			index_to_value_.erase(it);
		}
	}

	void MetaInfo::removeValue(UInt index)
	{
		map<UInt,DataValue>::iterator it = index_to_value_.find(index);
		if (it != index_to_value_.end())
		{
			index_to_value_.erase(it);
		}
	}

  void MetaInfo::getKeys(vector<String>& keys) const
  {
    keys.resize(index_to_value_.size());
    UInt i =0;
		for (map<UInt,DataValue>::const_iterator it = index_to_value_.begin(); it != index_to_value_.end(); ++it)
		{
			keys[i++]=registry_.getName(it->first);
		}
  }

	void MetaInfo::getKeys(vector<UInt>& keys) const
  {
    keys.resize(index_to_value_.size());
    UInt i =0;
		for (map<UInt,DataValue>::const_iterator it = index_to_value_.begin(); it != index_to_value_.end(); ++it)
		{
			keys[i++]=it->first;
		}
  }

  bool MetaInfo::empty() const
	{
		return index_to_value_.empty();
	}
	
  void MetaInfo::clear()
	{
		index_to_value_.clear();
	}

} //namespace
