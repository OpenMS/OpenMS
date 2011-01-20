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
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>
#include <OpenMS/CONCEPT/Exception.h>

using namespace std;

namespace OpenMS 
{
	
	CsvFile::CsvFile()
		: TextFile(),itemseperator_(','),itemenclosed_(false)
	{
		
	}

	CsvFile::~CsvFile()
	{
	}
	
	CsvFile::CsvFile(const String& filename, char is,bool ie, Int first_n)
		: TextFile(),itemseperator_(is),itemenclosed_(ie)
	{
		load(filename, false, first_n);
	}
  
  
	void CsvFile::fload(const String& filename, char is,bool ie, Int first_n)
	{
		itemseperator_ = is;
		itemenclosed_ = ie;
		load(filename,true,first_n);
	}

	bool CsvFile::getRow(Size row,StringList &list)
	{
		if(row > this->size())
		{
			throw Exception::InvalidIterator(__FILE__, __LINE__,__PRETTY_FUNCTION__);
		}
		bool splitted = this->operator[](row).split(itemseperator_, list);
		if(!splitted)
		{
			return splitted;
		}
		for(Size i = 0; i< list.size(); i++)
		{
			if(itemenclosed_)
			{
				list[i] = list[i].substr(1, list[i].size()-2);
			}
		}
			return true;
	}

} // namespace OpenMS
