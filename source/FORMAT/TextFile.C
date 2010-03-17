// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <fstream>
#include <iostream>

using namespace std;

namespace OpenMS 
{
	
	TextFile::TextFile()
		: StringList()
	{
		
	}

	TextFile::~TextFile()
	{
	}
	
	TextFile::TextFile(const String& filename, bool trim_lines, Int first_n) 
		: StringList()
	{
		load(filename, trim_lines, first_n);
	}
  
  
	void TextFile::load(const String& filename, bool trim_lines, Int first_n) 
	{
    ifstream is(filename.c_str(),ios_base::in | ios_base::binary);
    if (!is)
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
    }

		clear();
		
    String str;
    bool had_enough=false;
    while(getline(is,str,'\n') && !had_enough)
    {
      // platform specific line endings: 
      // Windows LE: \r\n
      //    we now have a line with \r at the end: get rid of it
      if (str.size()>=1 && *str.rbegin()=='\r') str = str.substr(0,str.size()-1);

      // Mac (OS<=9): \r
      //    we just read the whole file into a string: split it
      StringList lines;
      if (str.hasSubstring("\r")) lines = StringList::create(str,'\r');
      else lines.push_back(str);

      // Linux&MacOSX: \n
      //    nothing to do

      for (Size i=0;i<lines.size();++i)
      {
        str = lines[i];
    	  if (trim_lines)
    	  {
    		  push_back(str.trim());
    	  }
    	  else
    	  {
    		  push_back(str);
    	  }
      	
    	  if (first_n>-1 && (Int)(size())==first_n)
    	  {
          had_enough=true;
    		  break;
    	  }
      }
    }		
	}


	void TextFile::store(const String& filename) 
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
					os << it->substr(0,-2)<< "\n";
				}
				else
				{
					os << it->substr(0,-1) << "\n";
				}
			}
			else
			{
				os << *it << "\n";
			}
		}
		os.close();
	}
	
} // namespace OpenMS

