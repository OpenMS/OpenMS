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
// $Maintainer: Cornelia Friedle $
// --------------------------------------------------------------------------
//

#include <OpenMS/FILTERING/DATAREDUCTION/DataReducer.h>
//derived classes
#include<OpenMS/FILTERING/DATAREDUCTION/MaxReducer.h>
#include<OpenMS/FILTERING/DATAREDUCTION/SumReducer.h>
#include<OpenMS/CONCEPT/Factory.h>

using namespace std;

namespace OpenMS
{
	DataReducer::DataReducer()
		: FactoryProduct("DataRaducer")
	{

	} 
	
	DataReducer::~DataReducer()
	{
	} 
	
	DataReducer::DataReducer(const DataReducer& source)
    : FactoryProduct(source)
  {
  }
	
	DataReducer& DataReducer::operator=(const DataReducer& source)
  {
		if (&source == this) return *this;
	   
	  FactoryProduct::operator=(source);
	   
		return *this;
  }
	
	void DataReducer::registerChildren()
	{
		Factory<DataReducer>::registerProduct(MaxReducer::getProductName(), &MaxReducer::create);
		Factory<DataReducer>::registerProduct(SumReducer::getProductName(), &SumReducer::create);
	}

}// namespace openms

