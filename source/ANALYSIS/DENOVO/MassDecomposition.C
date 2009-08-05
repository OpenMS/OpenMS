// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/DENOVO/MassDecomposition.h>

using namespace std;

namespace OpenMS
{

	MassDecomposition::MassDecomposition(const String& deco)
		: number_of(0),
			number_of_max_aa(0),
			min_number(1000),
			max_number(0)
	{
		String tmp = deco;
		vector<String> split;
		
		// some more info per line
		if (deco.has('('))
		{
			Size pos = tmp.find('(', 0);
			tmp = tmp.substr(0, pos);
			tmp.trim();
		}

		tmp.split(' ', split);
		UInt sum = 0;
		number_of_max_aa = 0;
		// only one aa type?
		if (split.size() == 0 && tmp.size() != 0)
		{
			split.push_back(tmp);
		}
		if (split.size() != 0)
		{
			for (UInt i = 0; i != split.size(); ++i)
			{
				char aa = split[i][0];
				String s = split[i];
				s.erase(0, 1);
				UInt n = (UInt)s.toInt();
				if (number_of_max_aa < n)
				{
					number_of_max_aa = n;
				}
				sum += n;
				decomp[aa] = n;
			}
			number_of = sum;
			if (number_of < min_number)
			{
				min_number = number_of;
			}
			if (number_of > max_number)
			{
				max_number = number_of;
			}
		}
	}

	MassDecomposition::MassDecomposition(const MassDecomposition& rhs)
		: decomp(rhs.decomp),
			number_of(rhs.number_of),
			number_of_max_aa(rhs.number_of_max_aa),
			min_number(rhs.min_number),
			max_number(rhs.max_number)
	{
	}
	
	MassDecomposition& MassDecomposition::operator = (const MassDecomposition& rhs)
	{
		if (&rhs != this)
		{
			decomp = rhs.decomp;
			number_of = rhs.number_of;
			number_of_max_aa = rhs.number_of_max_aa;
			min_number = rhs.min_number;
			max_number = rhs.max_number;
		}
		return *this;
	}

	MassDecomposition& MassDecomposition::operator += (const MassDecomposition& d)
	{
		for (Map<char, UInt>::const_iterator it = d.decomp.begin(); it != d.decomp.end(); ++it)
		{
			if (decomp.find(it->first) == decomp.end())
			{
				decomp.insert(*it);
			}
			else
			{
				decomp[it->first] += it->second;
			}
		}

		return *this;
	}

	bool MassDecomposition::operator < (const MassDecomposition& rhs) const
	{
		return decomp < rhs.decomp;
	}	

	bool MassDecomposition::operator == (const String& deco) const
	{
		Map<char, UInt> tmp;
		for (String::ConstIterator it = deco.begin(); it != deco.end(); ++it)
		{
			char aa = *it;
			if (decomp.find(aa) == decomp.end())
			{
				return false;
			}

			if (tmp.find(aa) != tmp.end())
			{
				tmp[aa]++;
			}
			else
			{
				tmp[aa] = 1;
			}
		}

		return tmp == decomp;
	}

	String MassDecomposition::toString() const
	{
		String s;
		for (Map<char, UInt>::const_iterator it = decomp.begin(); it != decomp.end(); ++it)
		{
			s += it->first + String(it->second) + String(" ");
		}
		return s;
	}

	String MassDecomposition::toExpandedString() const
	{
		String s;
		for (Map<char, UInt>::const_iterator it = decomp.begin(); it != decomp.end(); ++it)
		{
			s += String(it->second, it->first);
		}
		return s;
	}

	bool MassDecomposition::containsTag(const String& tag) const
	{
		Map<char, UInt> tmp;
		for (String::ConstIterator it = tag.begin(); it != tag.end(); ++it)
		{
			char aa = *it;
			if (decomp.find(aa) == decomp.end())
			{
				return false;
			}
			if (tmp.find(aa) != tmp.end())
			{
				tmp[aa]++;
			}
			else
			{
				tmp[aa] = 1;
			}
		}

		// check if tag decomp is compatible with decomp
		for (Map<char, UInt>::const_iterator it = tmp.begin(); it != tmp.end(); ++it)
		{
			if (decomp.find(it->first)->second < it->second)
			{
				return false;
			}
		}

		return true;
	}

	bool MassDecomposition::compatible(const MassDecomposition& deco) const
	{
		for (Map<char, UInt>::const_iterator it = deco.decomp.begin(); it != deco.decomp.end(); ++it)
		{
			if (decomp.find(it->first) == decomp.end() || decomp.find(it->first)->second != it->second)
			{
				return false;
			}
		}
		return true;
	}

	MassDecomposition MassDecomposition::operator + (const MassDecomposition& rhs) const
	{
		MassDecomposition d(*this);
		for (Map<char, UInt>::const_iterator it = rhs.decomp.begin(); it != rhs.decomp.end(); ++it)
  	{
  		if (!d.decomp.has(it->first))
    	{
      	d.decomp.insert(*it);
    	}
    	else
    	{
      	d.decomp[it->first] += it->second;
   	 	}
  	}

  	return d;
	}

	UInt MassDecomposition::getNumberOfMaxAA() const
	{
		return number_of_max_aa;
	}

} // namespace OpenMS

