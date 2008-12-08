// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/DATASTRUCTURES/StringList.h>

using namespace std;

namespace OpenMS
{

	StringList::StringList()
	{
	}
	
	StringList::StringList(const StringList& rhs)
		: vector<String>(rhs)
	{
		
	}
	
	StringList::StringList(const vector<String>& rhs)
		: vector<String>(rhs)
	{
		
	}
	
	StringList::StringList(const vector<string>& rhs)
		: vector<String>(rhs.begin(),rhs.end())
	{
		
	}

	StringList& StringList::operator=(const StringList& rhs)
	{
		vector<String>::operator=(rhs);
		return *this;
	}
	
	StringList& StringList::operator=(const vector<String>& rhs)
	{
		vector<String>::operator=(rhs);
		return *this;
	}
	
	StringList& StringList::operator=(const vector<string>& rhs)
	{
		this->resize(rhs.size());
		for (UInt i=0; i<rhs.size(); ++i)
		{
			this->operator[](i)= rhs[i];
		}
		return *this;
	}

	StringList StringList::create(const String& list)
	{
		StringList out;
		if (!list.split(',',out) && list!="")
		{
			out.push_back(list);
		}
		return out;
	}
	
	bool StringList::contains(const String& s) const
	{
		for (UInt i=0; i<this->size(); ++i)
		{
			if (this->operator[](i)==s) return true;
		}
		return false;
	}
	
	void StringList::toUpper()
	{
		for (UInt i=0; i<this->size(); ++i)
		{
			this->operator[](i).toUpper();
		}
	}
	
	void StringList::toLower()
	{
		for (UInt i=0; i<this->size(); ++i)
		{
			this->operator[](i).toLower();
		}
	}

		String StringList::concatenate(const String& glue) const
		{
			
			if (size()==0) return "";
			
			String output = *(begin());
			for (const_iterator it=begin()+1; it!=end(); ++it)
			{
				output += (glue + *it) ;
			}
			
			return output;
		}

	// ----------------- Output operator ----------------------

	ostream& operator<<(std::ostream& os, const StringList& p)
	{
		os << "[";
		if (p.size()>0)
		{
			os << p[0];
		}
		
		for (UInt i=1; i<p.size(); ++i)
		{
			os << ", " << p[i];
		}
		os << "]";
		return os;
	}


} // namespace OpenMS


