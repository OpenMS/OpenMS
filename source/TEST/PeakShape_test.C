// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Eva Lange  $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakShape.h>

///////////////////////////

START_TEST(PeakShape, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

PeakShape* peakshape_ptr=0;
CHECK((PeakShape()))
  peakshape_ptr = new PeakShape;
  TEST_NOT_EQUAL(peakshape_ptr, 0)
RESULT

CHECK((virtual ~PeakShape()))
		delete peakshape_ptr;
RESULT
	
CHECK((PeakShape& operator = (const PeakShape& rhs)))
		PeakShape peakshape;
    peakshape.height = 10003.232;
    peakshape.mz_position = 0.323;
    peakshape.left_width = 2.998;
    peakshape.right_width = 2.776;
    peakshape.area = 8329832.141;
    peakshape.type = PeakShape::LORENTZ_PEAK;
    
    PeakShape peakshape_copy;
    peakshape_copy = peakshape;

    TEST_REAL_EQUAL(peakshape_copy.height, 10003.232) 
    TEST_REAL_EQUAL(peakshape_copy.mz_position, 0.323)
    TEST_REAL_EQUAL(peakshape_copy.left_width, 2.998)
    TEST_REAL_EQUAL(peakshape_copy.right_width, 2.776)
    TEST_REAL_EQUAL(peakshape_copy.area, 8329832.141)
    TEST_EQUAL(peakshape_copy.type, PeakShape::LORENTZ_PEAK)
RESULT

CHECK((PeakShape(const PeakShape& rhs)))
    PeakShape peakshape;
    peakshape.height = 10003.232;
    peakshape.mz_position = 0.323;
    peakshape.left_width = 2.998;
    peakshape.right_width = 2.776;
    peakshape.area = 8329832.141;
    peakshape.type = PeakShape::LORENTZ_PEAK;
    
    PeakShape peakshape_copy(peakshape);
   
    TEST_REAL_EQUAL(peakshape.height,10003.232) 
    TEST_REAL_EQUAL(peakshape.mz_position, 0.323)
    TEST_REAL_EQUAL(peakshape.left_width, 2.998)
    TEST_REAL_EQUAL(peakshape.right_width,2.776)
    TEST_REAL_EQUAL(peakshape.area, 8329832.141)
    TEST_EQUAL(peakshape.type, PeakShape::LORENTZ_PEAK)
RESULT

CHECK((PeakShape(DoubleReal height_, DoubleReal mz_position_, DoubleReal left_width_, DoubleReal right_width_, DoubleReal area_, PeakIterator left_, PeakIterator right_, Type type_)))
    DoubleReal height = 100.0;
    DoubleReal mz_position = 0.0;
    DoubleReal left_width = 3.0;
    DoubleReal right_width = 3.0;
    DoubleReal area = 309.23292;
    PeakShape::Type type = PeakShape::LORENTZ_PEAK;

    PeakShape peakshape(height,
    mz_position,
    left_width,
    right_width,
		area,
		type);

    TEST_REAL_EQUAL(peakshape.height,height) 
	  TEST_REAL_EQUAL(peakshape.mz_position, mz_position)
		TEST_REAL_EQUAL(peakshape.left_width, left_width)
    TEST_REAL_EQUAL(peakshape.right_width, right_width)
    TEST_REAL_EQUAL(peakshape.area, area)
    TEST_REAL_EQUAL(peakshape.r_value, 0.0)
		TEST_EQUAL(peakshape.type, PeakShape::LORENTZ_PEAK)
RESULT

CHECK((DoubleReal getSymmetricMeasure() const))
		DoubleReal height = 100.0;
    DoubleReal mz_position = 0.0;
    DoubleReal left_width = 3.0;
    DoubleReal right_width = 9.0;
    DoubleReal area = 309.23292;
    PeakShape::Type type = PeakShape::SECH_PEAK;

    PeakShape peakshape(height,
												mz_position,
												left_width,
												right_width,
												area,
												type);

    DoubleReal sym_value = peakshape.getSymmetricMeasure();
    TEST_REAL_EQUAL(sym_value,3.0/9.0)
RESULT

CHECK((DoubleReal operator() (DoubleReal x) const))
    DoubleReal height = 100.0;
    DoubleReal mz_position = 0.0;
    DoubleReal left_width = 4.0;
    DoubleReal right_width = 4.0;
    DoubleReal area = 100;
    PeakShape::Type type = PeakShape::LORENTZ_PEAK;
    
    PeakShape peakshape(height,
												mz_position,
												left_width,
												right_width,
												area,
												type);
   
    TEST_REAL_EQUAL(peakshape.getFWHM(),.5)
RESULT

CHECK((DoubleReal getFWHM() const))
  DoubleReal height = 100.0;
  DoubleReal mz_position = 0.0;
  DoubleReal left_width = 4.0;
  DoubleReal right_width = 4.0;
  DoubleReal area = 100;
  PeakShape::Type type = PeakShape::LORENTZ_PEAK;
    
  PeakShape p(height,
							mz_position,
							left_width,
							right_width,
							area,
							type);


  TEST_REAL_EQUAL(p.getFWHM(),1/right_width + 1/left_width)
RESULT

CHECK(bool operator==(const PeakShape &rhs) const)
	PeakShape p1,p2;
	TEST_EQUAL(p1==p2,true)
	
	p1.mz_position = 14.4;
	TEST_EQUAL(p1==p2,false)
	
	p2.mz_position = 14.4;
	TEST_EQUAL(p1==p2,true)
RESULT

CHECK(bool operator!=(const PeakShape &rhs) const)
	PeakShape p1,p2;
	TEST_EQUAL(p1!=p2,false)
	
	p1.mz_position = 14.4;
	TEST_EQUAL(p1!=p2,true)
	
	p2.mz_position = 14.4;
	TEST_EQUAL(p1!=p2,false)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

