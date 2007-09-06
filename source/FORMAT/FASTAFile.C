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
#include <OpenMS/SYSTEM/File.h>

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
 
		if (!File::exists(filename))
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
		}
		
		if (!File::readable(filename))
		{
			throw Exception::FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
		}
	
		String tag, seq;
	
		ifstream in(filename.c_str());
		String line;
		while (getline(in, line, '\n'))
		{
			if (line.size() > 0)
			{
				if (line[0] == '>')
				{
					if (tag != "" && seq != "")
					{
						FASTAEntry entry;
						entry.first = tag;
						entry.second = seq;
						data.push_back(entry);
						tag = "";
						seq = "";
					}
								
					line.erase(line.begin());
					tag = line.trim();
				}
				else
				{
					seq += line.trim();
				}
			}
		}

		if (tag != "" && seq != "")
		{
			FASTAEntry entry;
      entry.first = tag;
      entry.second = seq;
      data.push_back(entry);
		}
		in.close();
		
		return;
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
			
			outfile << ">" << it->first << endl;
			
			String tmp(it->second);
			while (tmp.size() > 80)
			{
				outfile << tmp.prefix(80) << endl;
				tmp.erase(0, 80);
			}

			if (tmp.size() > 0)
			{
				outfile << tmp << endl;
			}
		}
		outfile.close();
	}


} // namespace OpenMS
