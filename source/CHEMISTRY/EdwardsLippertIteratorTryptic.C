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
// $Maintainer: Clemens Groepl,Andreas Bertsch$
// $Authors: $
// --------------------------------------------------------------------------
 

#include <OpenMS/CHEMISTRY/EdwardsLippertIteratorTryptic.h>
namespace OpenMS 
{
	EdwardsLippertIteratorTryptic::EdwardsLippertIteratorTryptic()
		: EdwardsLippertIterator()
	{
	}

	EdwardsLippertIteratorTryptic::EdwardsLippertIteratorTryptic(const EdwardsLippertIteratorTryptic& rhs)
		: EdwardsLippertIterator(rhs)
	{
	}
	
	EdwardsLippertIteratorTryptic::~EdwardsLippertIteratorTryptic()
	{
	}
	
	EdwardsLippertIteratorTryptic& EdwardsLippertIteratorTryptic::operator = (const EdwardsLippertIteratorTryptic& rhs)
	{
		if (this != &rhs)
		{
			this->EdwardsLippertIterator::operator = (rhs);
		}
		return *this;
	}
	
	bool EdwardsLippertIteratorTryptic::isDigestingEnd(char aa1, char aa2)
	{
		return (aa1 == 'K' || aa1 == 'R') && aa2 != 'P';
	}

}
