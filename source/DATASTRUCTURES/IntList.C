// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: David Wojnar $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/IntList.h>

using namespace std;

namespace OpenMS
{

	IntList::IntList()
	{
	}
	
	IntList::IntList(const IntList& rhs)
		: vector<Int>(rhs)
	{
		
	}
	
	IntList::IntList(const vector<Int>& rhs)
		: vector<Int>(rhs)
	{
		
	}

	IntList::IntList(const vector<UInt>& rhs)
	{
		this->resize(rhs.size());
		for (Size i=0;i<rhs.size();++i)
		{
			(*this)[i]=(Int)rhs[i];
		}
	}
	
	IntList& IntList::operator=(const IntList& rhs)
	{
		vector<Int>::operator=(rhs);
		return *this;
	}
	
	IntList& IntList::operator=(const vector<UInt>& rhs)
	{
		this->resize(rhs.size());
		for (Size i=0;i<rhs.size();++i)
		{
			(*this)[i]=(Int)rhs[i];
		}
		return *this;
	}

	IntList& IntList::operator=(const vector<Int>& rhs)
	{
		vector<Int>::operator=(rhs);
		return *this;
	}

	IntList IntList::create(const String& list)
	{
		IntList ret;
		vector<String> out;
		list.split(',', out);
		ret.resize(out.size());
		for (Size i = 0; i < out.size(); ++i)
		{
			ret[i] = out[i].toInt();
		}
		return ret;
	}
	
	IntList IntList::create(const StringList& list)
	{
		IntList ret;
		for(UInt i = 0 ; i < list.size(); ++i)
		{
			ret.push_back(list[i].toInt());
		}
		return ret;
	}
	
	bool IntList::contains(Int s) const
	{
		for (Size i=0; i<this->size(); ++i)
		{
			if (this->operator[](i)==s) return true;
		}
		return false;
	}

	// ----------------- Output operator ----------------------

	ostream& operator<<(std::ostream& os, const IntList& p)
	{
		os << "[";
		if (p.size()>0)
		{
			os << p[0];
		}
		
		for (Size i=1; i<p.size(); ++i)
		{
			os << ", " << p[i];
		}
		os << "]";
		return os;
	}


} // namespace OpenMS



