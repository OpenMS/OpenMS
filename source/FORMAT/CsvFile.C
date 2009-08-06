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
// $Maintainer: David Wojnar $
// $Authors: $
// --------------------------------------------------------------------------
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>

using namespace std;

namespace OpenMS 
{
	
	CsvFile::CsvFile()
		: TextFile(),itemseperator(','),itemenclosed(false)
	{
		
	}

	CsvFile::~CsvFile()
	{
	}
	
	CsvFile::CsvFile(const String& filename, char is,bool ie, Int first_n)
		: TextFile(),itemseperator(is),itemenclosed(ie)
	{
		load(filename, false, first_n);
	}
  
  
	void CsvFile::fload(const String& filename, char is,bool ie, Int first_n)
	{
		itemseperator = is;
		itemenclosed = ie;
		load(filename,true,first_n);
	}

	bool CsvFile::getRow(UInt row,StringList &list)
	{
		if(row > this->size())
		{
			throw Exception::InvalidIterator(__FILE__, __LINE__,__PRETTY_FUNCTION__);
		}
		bool splitted = this->operator[](row).split(itemseperator, list);
		if(!splitted)
		{
			return splitted;
		}
		for(UInt i = 0; i< list.size(); i++)
		{
			if(itemenclosed)
			{
				list[i] = list[i].substr(1, list[i].size()-2);
			}
		}
			return true;
	}

} // namespace OpenMS
