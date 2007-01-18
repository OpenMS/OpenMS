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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiModule.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>

#include <OpenMS/CONCEPT/Exception.h>


///////////////////////////

START_TEST(FeaFiModule, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

// default ctor
FeaFiModule* ptr = 0;
CHECK(FeaFiModule())
	ptr = new FeaFiModule();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

// destructor
CHECK(~FeaFiModule())
	delete ptr;
RESULT

// assignment operator
CHECK(FeaFiModule& operator = (const FeaFiModule& source))
	FeaFiModule fm1;
	FeaFiTraits fft;
  fm1.setTraits(&fft);
  FeaFiModule fm2;
  fm2 = fm1;

  FeaFiModule fm3;
	fm3.setTraits(&fft);
	
  fm1 = FeaFiModule();
	TEST_EQUAL(fm3,fm2)
RESULT

// copy constructor
CHECK(TestModel(const FeaFiModule& source))
	FeaFiModule fm1;	
  FeaFiTraits fft;
  fm1.setTraits(&fft);

  FeaFiModule fm2(fm1);

  FeaFiModule fm3;
  fm3.setTraits(&fft);

  fm1 = FeaFiModule();
	TEST_EQUAL(fm2, fm3)
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
