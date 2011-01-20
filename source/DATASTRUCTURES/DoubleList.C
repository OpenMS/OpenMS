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
// $Maintainer: David Wojnar $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/DoubleList.h>
using namespace std;

namespace OpenMS
{

	DoubleList::DoubleList()
	{
	}
	
	DoubleList::DoubleList(const DoubleList& rhs)
		: vector<DoubleReal>(rhs)
	{
		
	}
	
	DoubleList::DoubleList(const vector<DoubleReal>& rhs)
		: vector<DoubleReal>(rhs)
	{
		
	}

	DoubleList::DoubleList(const vector<Real>& rhs)
	{
		this->resize(rhs.size());
		for (Size i=0;i<rhs.size();++i)
		{
			(*this)[i]=(DoubleReal)rhs[i];
		}
	}
	
	DoubleList& DoubleList::operator=(const DoubleList& rhs)
	{
		vector<DoubleReal>::operator=(rhs);
		return *this;
	}
	
	DoubleList& DoubleList::operator=(const vector<Real>& rhs)
	{
		this->resize(rhs.size());
		for (Size i=0;i<rhs.size();++i)
		{
			(*this)[i]=(DoubleReal)rhs[i];
		}
		return *this;
	}

	DoubleList& DoubleList::operator=(const vector<DoubleReal>& rhs)
	{
		vector<DoubleReal>::operator=(rhs);
		return *this;
	}

	DoubleList DoubleList::create(const String& list)
	{
		DoubleList ret;
		vector<String> out;
		list.split(',', out);
		ret.resize(out.size());
		for (Size i = 0; i < out.size(); ++i)
		{
			ret[i]=out[i].toDouble();
		}
		return ret;
	}
	
	DoubleList DoubleList::create(const StringList& list)
	{
		DoubleList ret;
		for(UInt i = 0 ; i < list.size(); ++i)
		{
			ret.push_back(list[i].toDouble());
		}
		return ret;
	}
	
	bool DoubleList::contains(DoubleReal s, DoubleReal tolerance) const
	{
		for (Size i=0; i<this->size(); ++i)
		{
			if (std::fabs(this->operator[](i)-s)<tolerance) return true;
		}
		return false;
	}

	// ----------------- Output operator ----------------------

	ostream& operator<<(std::ostream& os, const DoubleList& p)
	{
		os << "[";
		if (p.size()>0)
		{
			os << precisionWrapper(p[0]);
		}
		
		for (Size i=1; i<p.size(); ++i)
		{
			os << ", " << precisionWrapper(p[i]);
		}
		os << "]";
		return os;
	}


} // namespace OpenMS



