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
// $Maintainer: Marc Sturm $
// $Authors: $
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
		for (Size i=0; i<rhs.size(); ++i)
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
	
	StringList StringList::create( const char * const * list, UInt size )
	{
		StringList out;
		for ( UInt i = 0; i < size; ++i )
		{
			out.push_back(list[i]);
		}
		return out;
	}
	
	bool StringList::contains(const String& s) const
	{
		for (Size i=0; i<this->size(); ++i)
		{
			if (this->operator[](i)==s) return true;
		}
		return false;
	}
	
	void StringList::toUpper()
	{
		for (Size i=0; i<this->size(); ++i)
		{
			this->operator[](i).toUpper();
		}
	}
	
	void StringList::toLower()
	{
		for (Size i=0; i<this->size(); ++i)
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
		
		for (Size i=1; i<p.size(); ++i)
		{
			os << ", " << p[i];
		}
		os << "]";
		return os;
	}

	StringList::Iterator StringList::search(const String& text, bool trim)
	{
		return search(begin(),text,trim);
	}
	
	StringList::ConstIterator StringList::search(const String& text, bool trim) const
  {
    return search(begin(),text,trim);
  }
	
	StringList::Iterator StringList::search(const Iterator& start, const String& text, bool trim)
	{
		String pattern = text;
		if (trim)
		{
			pattern.trim();
		}		
		
		String tmp;
		
		Iterator it = start;

		while(it!=end())
		{
			tmp = *it;
			if (trim)
			{
				if (tmp.trim().hasPrefix(pattern)) return it;	
			}
			else
			{
				if (tmp.hasPrefix(pattern)) return it;
			}
			++it;
		}
		
		//nothing found
		return end();
	}

  StringList::ConstIterator StringList::search(const ConstIterator& start, const String& text, bool trim) const
  {
    String pattern = text;
    if (trim)
    {
      pattern.trim();
    }

    String tmp;

    ConstIterator it = start;

    while(it!=end())
    {
      tmp = *it;
      if (trim)
      {
        if (tmp.trim().hasPrefix(pattern)) return it;
      }
      else
      {
        if (tmp.hasPrefix(pattern)) return it;
      }
      ++it;
    }

    //nothing found
    return end();
  }

	
	StringList::Iterator StringList::searchSuffix(const Iterator& start, const String& text, bool trim)
	{
		String pattern = text;
		if (trim)
		{
			pattern.trim();
		}		
		
		String tmp;
		
		Iterator it = start;
		
		while(it!=end())
		{
			tmp = *it;
			if (trim)
			{
				if (tmp.trim().hasSuffix(pattern)) return it;	
			}
			else
			{
				if (tmp.hasSuffix(pattern)) return it;
			}
			++it;
		}
		
		//nothing found
		return end();		
	}

  StringList::ConstIterator StringList::searchSuffix(const ConstIterator& start, const String& text, bool trim) const
  {
    String pattern = text;
    if (trim)
    {
      pattern.trim();
    }

    String tmp;

    ConstIterator it = start;

    while(it!=end())
    {
      tmp = *it;
      if (trim)
      {
        if (tmp.trim().hasSuffix(pattern)) return it;
      }
      else
      {
        if (tmp.hasSuffix(pattern)) return it;
      }
      ++it;
    }

    //nothing found
    return end();
  }


	StringList::Iterator StringList::searchSuffix(const String& text, bool trim)
	{
		return searchSuffix(begin(),text,trim);		
	}

  StringList::ConstIterator StringList::searchSuffix(const String& text, bool trim) const
  {
    return searchSuffix(begin(),text,trim);
  }

} // namespace OpenMS


