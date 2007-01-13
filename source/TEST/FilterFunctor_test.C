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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/Identification.h>
///////////////////////////

#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>

///////////////////////////

#include <vector>
#include <iostream>

///////////////////////////
START_TEST(FilterFunctor, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

CHECK(FilterFunctor())
	// nothing to test
RESULT

CHECK(FilterFunctor(const FilterFunctor& source))
	// nothing to test
RESULT

CHECK(FilterFunctor& operator = (const FilterFunctor& source))
	// nothing to test
RESULT

CHECK(static void registerChildren())
	// TODO
RESULT

CHECK(template<typename SpectrumType> double apply(SpectrumType& /* spectrum */))
	// nothing to test
RESULT

CHECK(~FilterFunctor())
	// nothing to test
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
