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
// $Maintainer: Clemens Groepl, Andreas Bertsch $
// $Authors: Chris Bauer $
// --------------------------------------------------------------------------


#include <OpenMS/CHEMISTRY/PepIterator.h>
#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/FORMAT/FastaIterator.h>
#include <OpenMS/CHEMISTRY/EdwardsLippertIterator.h>
#include <OpenMS/FORMAT/FastaIteratorIntern.h>
#include <OpenMS/CHEMISTRY/EdwardsLippertIteratorTryptic.h>
#include <OpenMS/CHEMISTRY/TrypticIterator.h>

namespace OpenMS
{
	
	typedef std::pair <String, String> FASTAEntry;
	
	PepIterator::PepIterator()
	{
	}

	PepIterator::PepIterator(const PepIterator& /*source*/) 
	{
	}
	
	PepIterator::~PepIterator()
	{
	}

	void PepIterator::registerChildren()
	{
		//register new products here
		Factory<PepIterator>::registerProduct(EdwardsLippertIterator::getProductName(),& EdwardsLippertIterator::create);
		Factory<PepIterator>::registerProduct(EdwardsLippertIteratorTryptic::getProductName(),& EdwardsLippertIteratorTryptic::create);	
		Factory<PepIterator>::registerProduct(FastaIterator::getProductName(),& FastaIterator::create);	
		Factory<PepIterator>::registerProduct(FastaIteratorIntern::getProductName(),& FastaIteratorIntern::create);
		Factory<PepIterator>::registerProduct(TrypticIterator::getProductName(),& TrypticIterator::create);
	}
}
