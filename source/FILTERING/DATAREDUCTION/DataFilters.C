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

#include <OpenMS/FILTERING/DATAREDUCTION/DataFilters.h>
#include <OpenMS/KERNEL/Feature.h>

#include <iostream>

using namespace std;

namespace OpenMS
{
	
	String DataFilters::DataFilter::toString() const
	{
		String out;
		//field
		if (field==INTENSITY) out = "Intensity ";
		else if (field==QUALITY) out = "Quality ";
		else if (field==CHARGE) out = "Charge ";
		else if (field==META_DATA) out = "Meta::" + meta_name + " ";
		//operation
		if (op==GREATER_EQUAL) out += ">= ";
		else if (op==EQUAL) out += "= ";
		else if (op==LESS_EQUAL) out += "<= ";
		else if (op==EXISTS) out += "exists";
		//value
		if (field == META_DATA)
		{
			if (op != EXISTS)
			{
				if (value_is_numerical) out = out + value;
				else out = out + value_string;
			}
			return out;
		}
		out = out + value;
		return out;
	}
	
	void DataFilters::DataFilter::fromString(const String& filter) throw (Exception::InvalidValue)
	{
		bool meta = false;
		String tmp = filter;
		tmp.trim();
		vector<String> parts;
		tmp.split(' ',parts);
		int size = parts.size();
		if (size!=2 && size!=3) throw Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__,tmp,"Invalid filter format");
		//field
		tmp=parts[0];
		tmp.toLower();
		if (tmp=="intensity") field = INTENSITY;
		else if (tmp=="charge") field = CHARGE;
		else if (tmp=="quality") field = QUALITY;
		else if (tmp.hasPrefix(String("meta::")))
		{
			meta = true;
			field = META_DATA;
			meta_name = tmp.suffix(tmp.size() - 6);
		}
		else Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__,tmp,"Invalid field name");
		//operation
		tmp=parts[1];
		if (tmp==">=") op = GREATER_EQUAL;
		else if (tmp=="=") op = EQUAL;
		else if (tmp=="<=") op = LESS_EQUAL;
		else if (tmp=="exists" && meta)
		{
			op = EXISTS;
			return;
		}
		else Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__,tmp,"Invalid operator");
		//value
		tmp = parts[2];
		try
		{
			value = tmp.toDouble();
			value_is_numerical = true;
		}
		catch (Exception::ConversionError)
		{
			value_is_numerical = false;
			if (!meta) // non meta values must be numerical
			{
				Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__,tmp,"Invalid value");
			}
			else
			{
				value_string = tmp;
			}
		}
	}
	
	void DataFilters::add(const DataFilter& filter)
	{
		filters_.push_back(filter);
	}
	
	void DataFilters::remove(UInt index) throw (Exception::IndexOverflow)
	{
		if (index>=filters_.size()) throw Exception::IndexOverflow(__FILE__,__LINE__,__PRETTY_FUNCTION__,index,filters_.size());
		filters_.erase(filters_.begin()+index);
	}

	void DataFilters::replace(UInt index, const DataFilter& filter) throw (Exception::IndexOverflow)
	{
		filters_[index] = filter;
	}


	void DataFilters::clear()
	{
		filters_.clear();
	}


	UInt DataFilters::size() const
	{
		return filters_.size();
	}
	
	const DataFilters::DataFilter& DataFilters::operator[](UInt index) const
	{
		return filters_[index];
	}
	
	bool DataFilters::passes(const Feature& feature) const
	{
		for (vector<DataFilter>::const_iterator it=filters_.begin(); it!=filters_.end(); ++it)
		{
			if (it->field==INTENSITY)
			{
				if (it->op==GREATER_EQUAL && feature.getIntensity()<it->value) return false;
				else if (it->op==LESS_EQUAL && feature.getIntensity()>it->value) return false;
				else if (it->op==EQUAL && feature.getIntensity()!=it->value) return false;
			}
			else if (it->field==QUALITY)
			{
				if (it->op==GREATER_EQUAL && feature.getOverallQuality()<it->value) return false;
				else if (it->op==LESS_EQUAL && feature.getOverallQuality()>it->value) return false;
				else if (it->op==EQUAL && feature.getOverallQuality()!=it->value) return false;
			}
			else if (it->field==CHARGE)
			{
				if (it->op==EQUAL && feature.getCharge()!=it->value) return false;
				else if (it->op==GREATER_EQUAL && feature.getCharge()<it->value) return false;
				else if (it->op==LESS_EQUAL && feature.getCharge()>it->value) return false;
			}
			else if (it->field==META_DATA)
			{
				const MetaInfoInterface& mii = static_cast<MetaInfoInterface>(feature);
				if(!metaPasses(mii,it)) return false;
			}
		}
		return true;
	}

}//Namespace

