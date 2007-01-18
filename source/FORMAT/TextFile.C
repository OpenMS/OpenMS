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

#include <OpenMS/FORMAT/TextFile.h>

#include <fstream>
#include <iostream>

using namespace std;

namespace OpenMS 
{
	
	TextFile::TextFile()
		: vector<String>()
	{
		
	}

	TextFile::~TextFile()
	{
	}
	
	TextFile::TextFile(const String& filename, bool trim_lines) 
		throw (Exception::FileNotFound)
		: vector<String>()
	{
		load(filename, trim_lines);
	}
  
  
	void TextFile::load(const String& filename, bool trim_lines) 
		throw (Exception::FileNotFound)
	{
		ifstream is(filename.c_str());
    if (!is)
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
    }

		clear();
    
    String str;
    while(getline(is,str,'\n'))
    {
    	if (trim_lines)
    	{
    		push_back(str.trim());
    	}
    	else
    	{
    		push_back(str);
    	}
    }		
	}


	void TextFile::save(const String& filename) 
		throw (Exception::UnableToCreateFile)
	{
		ofstream os;
		os.open (filename.c_str(), ofstream::out);
		
		if(!os)
		{
			 throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
		}
		
		for (Iterator it = begin(); it!=end(); ++it)
		{
			if (it->hasSuffix("\n"))
			{
				if (it->hasSuffix("\r\n"))
				{
					os << it->substr(0,-2)<< endl;
				}
				else
				{
					os << it->substr(0,-1) << endl;
				}
			}
			else
			{
				os << *it << endl;
			}
		}
		os.close();
	}
	
	
	TextFile::Iterator TextFile::search(const String& text, bool trim)
	{
		return search(begin(),text,trim);
	}
	
	TextFile::ConstIterator TextFile::search(const String& text, bool trim) const
  {
    return search(begin(),text,trim);
  }
	
	TextFile::Iterator TextFile::search(const Iterator& start, const String& text, bool trim)
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

  TextFile::ConstIterator TextFile::search(const ConstIterator& start, const String& text, bool trim) const
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

	
	TextFile::Iterator TextFile::searchSuffix(const Iterator& start, const String& text, bool trim)
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

  TextFile::ConstIterator TextFile::searchSuffix(const ConstIterator& start, const String& text, bool trim) const
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


	TextFile::Iterator TextFile::searchSuffix(const String& text, bool trim)
	{
		return searchSuffix(begin(),text,trim);		
	}

  TextFile::ConstIterator TextFile::searchSuffix(const String& text, bool trim) const
  {
    return searchSuffix(begin(),text,trim);
  }


	String TextFile::asString() const
	{
		String tmp;
		tmp.implode(this->begin(),this->end());
		return tmp;
	}

} // namespace OpenMS

