// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
	DataReducer::DataReducer(): FactoryProduct()
	{
	} 
	
	DataReducer::~DataReducer()
	{
	} 
	
	DataReducer::DataReducer(const DataReducer& source)
    : FactoryProduct(source)
  {
		name_ = source.getName();
  }
	
	DataReducer& DataReducer::operator=(const DataReducer& source)
  {
   FactoryProduct::operator=(source);
	 return *this;
  }
	
	void DataReducer::setParameter(const Param& p)
	{
		if(checkParameter(p))
		{
			param_ = p;
		}
		else
		{
			throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,"The Parmater is empty");
		}
		}
	
	bool DataReducer::checkParameter(const Param& p)
	{
		if(p.empty())
		{
			throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,"The Parmater is empty");
 			return false;
		}
		else
		{
			return true;
		}
	}
	
	void DataReducer::registerChildren()
	{
		Factory<DataReducer>::registerProduct(MaxReducer::getName(), &MaxReducer::create);
		Factory<DataReducer>::registerProduct(SumReducer::getName(), &SumReducer::create);
	}

}// namespace openms

