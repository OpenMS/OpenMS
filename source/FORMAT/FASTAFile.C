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
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/TextFile.h>

#include <fstream>

using namespace std;

namespace OpenMS 
{
  FASTAFile::FASTAFile()
  {
  	
  }
  
  FASTAFile::~FASTAFile()
  {
    
  }

  void FASTAFile::load(const String& filename, FASTAType& data) throw (Exception::FileNotFound,Exception::ParseError)
  {
  	data.clear();
  	
   	TextFile file(filename, true);
    TextFile::iterator running_iterator;
    TextFile::iterator end_iterator;
    String actual_tag = "";
    String actual_sequence = "";
    vector<String> parts;

		if (file.size() == 0)
		{
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, 
				"File is empty!", filename);
		}
	
  	running_iterator = file.search(">");
  	while(running_iterator != file.end())
  	{
  		end_iterator = running_iterator;
  		end_iterator++;
  		end_iterator = file.search(end_iterator, ">");
			running_iterator->suffix('>').trim().split(' ', parts);
			if (parts.size() > 0)
			{
				actual_tag = parts[0];
			}
			else
			{
				actual_tag = running_iterator->suffix('>');
			}
			running_iterator++;
  		if (end_iterator != file.end())
  		{
  			while(running_iterator != end_iterator 
  						&& running_iterator != file.end())
  			{
  				actual_sequence = actual_sequence + (*running_iterator);
  				running_iterator++;
  			}
  		}
  		else
  		{
  			while(running_iterator != file.end())
  			{
  				actual_sequence = actual_sequence + (*running_iterator);
  				running_iterator++;
  			}  			  			
  		}
			data.push_back(make_pair(actual_tag, actual_sequence));
	  	actual_sequence = "";  			  			
  	}		    
  }

	void FASTAFile::store(const String& filename, const FASTAType& data) const throw (Exception::UnableToCreateFile)
	{
		ofstream outfile;
		outfile.open(filename.c_str(), ofstream::out);
		
		if (!outfile.good())
		{
			throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
		}
		
		for (FASTAType::const_iterator it = data.begin(); it!=data.end(); ++it)
		{
			outfile << ">" <<it->first <<endl<< it->second<< endl <<endl;
		}
		outfile.close();
	}


} // namespace OpenMS
