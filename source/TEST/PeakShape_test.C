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
// $Maintainer: Alexandra Zerck $
// $Authors: Eva Lange $
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
START_SECTION((PeakShape()))
  peakshape_ptr = new PeakShape;
  TEST_NOT_EQUAL(peakshape_ptr, 0)
END_SECTION

START_SECTION((virtual ~PeakShape()))
		delete peakshape_ptr;
END_SECTION
	
START_SECTION((PeakShape& operator = (const PeakShape& rhs)))
		PeakShape peakshape;
    peakshape.height = 10003.232;
    peakshape.mz_position = 0.323;
    peakshape.left_width = 2.998;
    peakshape.right_width = 2.776;
    peakshape.area = 8329832.141;
    peakshape.type = PeakShape::LORENTZ_PEAK;
    
    PeakShape peakshape_copy;
    peakshape_copy = peakshape;

    TEST_REAL_SIMILAR(peakshape_copy.height, 10003.232) 
    TEST_REAL_SIMILAR(peakshape_copy.mz_position, 0.323)
    TEST_REAL_SIMILAR(peakshape_copy.left_width, 2.998)
    TEST_REAL_SIMILAR(peakshape_copy.right_width, 2.776)
    TEST_REAL_SIMILAR(peakshape_copy.area, 8329832.141)
    TEST_EQUAL(peakshape_copy.type, PeakShape::LORENTZ_PEAK)
END_SECTION

START_SECTION((PeakShape(const PeakShape& rhs)))
    PeakShape peakshape;
    peakshape.height = 10003.232;
    peakshape.mz_position = 0.323;
    peakshape.left_width = 2.998;
    peakshape.right_width = 2.776;
    peakshape.area = 8329832.141;
    peakshape.type = PeakShape::LORENTZ_PEAK;
    
    PeakShape peakshape_copy(peakshape);
   
    TEST_REAL_SIMILAR(peakshape.height,10003.232) 
    TEST_REAL_SIMILAR(peakshape.mz_position, 0.323)
    TEST_REAL_SIMILAR(peakshape.left_width, 2.998)
    TEST_REAL_SIMILAR(peakshape.right_width,2.776)
    TEST_REAL_SIMILAR(peakshape.area, 8329832.141)
    TEST_EQUAL(peakshape.type, PeakShape::LORENTZ_PEAK)
END_SECTION
MSSpectrum<> spec;
spec.resize(100);
for(Int i = 0; i<100;++i)
{
  spec[i].setMZ(i*0.1);
  spec[i].setIntensity(100.); 
}
START_SECTION((PeakShape(DoubleReal height_, DoubleReal mz_position_, DoubleReal left_width_, DoubleReal right_width_, DoubleReal area_, PeakIterator left_, PeakIterator right_, Type type_)))
    DoubleReal height = 100.0;
    DoubleReal mz_position = 0.0;
    DoubleReal left_width = 3.0;
    DoubleReal right_width = 3.0;
    DoubleReal area = 309.23292;
    PeakShape::Type type = PeakShape::LORENTZ_PEAK;

    
    PeakShape::PeakIterator it1 = spec.begin()+2;
    PeakShape::PeakIterator it2 = spec.begin()+30;
    PeakShape peakshape(height,
    mz_position,
    left_width,
    right_width,
		area,it1,it2,
		type);

    TEST_EQUAL(peakshape.iteratorsSet(), true)
    TEST_REAL_SIMILAR(peakshape.height,height) 
	  TEST_REAL_SIMILAR(peakshape.mz_position, mz_position)
		TEST_REAL_SIMILAR(peakshape.left_width, left_width)
    TEST_REAL_SIMILAR(peakshape.right_width, right_width)
    TEST_REAL_SIMILAR(peakshape.area, area)
    TEST_REAL_SIMILAR(peakshape.r_value, 0.0)
		TEST_EQUAL(peakshape.type, PeakShape::LORENTZ_PEAK)
END_SECTION

START_SECTION((PeakShape(DoubleReal height_, DoubleReal mz_position_, DoubleReal left_width_, DoubleReal right_width_, DoubleReal area_, Type type_)))
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

    TEST_EQUAL(peakshape.iteratorsSet(), false)
    TEST_REAL_SIMILAR(peakshape.height,height) 
	  TEST_REAL_SIMILAR(peakshape.mz_position, mz_position)
		TEST_REAL_SIMILAR(peakshape.left_width, left_width)
    TEST_REAL_SIMILAR(peakshape.right_width, right_width)
    TEST_REAL_SIMILAR(peakshape.area, area)
    TEST_REAL_SIMILAR(peakshape.r_value, 0.0)
		TEST_EQUAL(peakshape.type, PeakShape::LORENTZ_PEAK)
END_SECTION

START_SECTION((bool iteratorsSet() const))
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

    
    PeakShape::PeakIterator it1 = spec.begin()+2;
    PeakShape::PeakIterator it2 = spec.begin()+30;
    PeakShape peakshape2(height,
    mz_position,
    left_width,
    right_width,
		area,it1,it2,
		type);

    TEST_EQUAL(peakshape2.iteratorsSet(), true)
    TEST_EQUAL(peakshape.iteratorsSet(), false)
END_SECTION

START_SECTION((PeakIterator getRightEndpoint() const))
    DoubleReal height = 100.0;
    DoubleReal mz_position = 4.0;
    DoubleReal left_width = 3.0;
    DoubleReal right_width = 3.0;
    DoubleReal area = 309.23292;
    PeakShape::Type type = PeakShape::LORENTZ_PEAK;

    PeakShape::PeakIterator it1 = spec.begin()+2;
    PeakShape::PeakIterator it2 = spec.begin()+30;
    PeakShape peakshape(height,
    mz_position,
    left_width,
    right_width,
		area,it1,it2,
		type);

    TEST_REAL_SIMILAR(peakshape.getRightEndpoint()->getMZ(), (spec.begin()+30)->getMZ())
    TEST_REAL_SIMILAR(peakshape.getRightEndpoint()->getIntensity(), (spec.begin()+30)->getIntensity())
END_SECTION


START_SECTION((void setRightEndpoint(PeakIterator right_endpoint)))
    DoubleReal height = 100.0;
    DoubleReal mz_position = 4.0;
    DoubleReal left_width = 3.0;
    DoubleReal right_width = 3.0;
    DoubleReal area = 309.23292;
    PeakShape::Type type = PeakShape::LORENTZ_PEAK;

    PeakShape::PeakIterator it1 = spec.begin()+2;
    PeakShape::PeakIterator it2 = spec.begin()+30;
    PeakShape peakshape(height,
    mz_position,
    left_width,
    right_width,
		area,
		type);

    peakshape.setLeftEndpoint(it1);
    peakshape.setRightEndpoint(it2); 
    TEST_EQUAL(peakshape.iteratorsSet(), true)
    TEST_REAL_SIMILAR(peakshape.getRightEndpoint()->getMZ(), (spec.begin()+30)->getMZ())
    TEST_REAL_SIMILAR(peakshape.getRightEndpoint()->getIntensity(), (spec.begin()+30)->getIntensity())
END_SECTION

START_SECTION((PeakIterator getLeftEndpoint() const))
    DoubleReal height = 100.0;
    DoubleReal mz_position = 4.0;
    DoubleReal left_width = 3.0;
    DoubleReal right_width = 3.0;
    DoubleReal area = 309.23292;
    PeakShape::Type type = PeakShape::LORENTZ_PEAK;

    PeakShape::PeakIterator it1 = spec.begin()+2;
    PeakShape::PeakIterator it2 = spec.begin()+30;
    PeakShape peakshape(height,
    mz_position,
    left_width,
    right_width,
		area,it1,it2,
		type);

    TEST_REAL_SIMILAR(peakshape.getLeftEndpoint()->getMZ(), (spec.begin()+2)->getMZ())
    TEST_REAL_SIMILAR(peakshape.getLeftEndpoint()->getIntensity(), (spec.begin()+2)->getIntensity())
END_SECTION


START_SECTION((void setLeftEndpoint(PeakIterator left_endpoint)))
    DoubleReal height = 100.0;
    DoubleReal mz_position = 4.0;
    DoubleReal left_width = 3.0;
    DoubleReal right_width = 3.0;
    DoubleReal area = 309.23292;
    PeakShape::Type type = PeakShape::LORENTZ_PEAK;

    PeakShape::PeakIterator it1 = spec.begin()+2;
    PeakShape::PeakIterator it2 = spec.begin()+30;
    PeakShape peakshape(height,
    mz_position,
    left_width,
    right_width,
		area,
		type);
    peakshape.setLeftEndpoint(it1);

    TEST_EQUAL(peakshape.iteratorsSet(), false)
    TEST_REAL_SIMILAR(peakshape.getLeftEndpoint()->getMZ(), (spec.begin()+2)->getMZ())
    TEST_REAL_SIMILAR(peakshape.getLeftEndpoint()->getIntensity(), (spec.begin()+2)->getIntensity())
END_SECTION


START_SECTION((DoubleReal getSymmetricMeasure() const))
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
    TEST_REAL_SIMILAR(sym_value,3.0/9.0)
END_SECTION

START_SECTION((DoubleReal operator() (DoubleReal x) const))
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
   
    TEST_REAL_SIMILAR(peakshape.getFWHM(),.5)
END_SECTION

START_SECTION((DoubleReal getFWHM() const))
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


  TEST_REAL_SIMILAR(p.getFWHM(),1/right_width + 1/left_width)
END_SECTION

START_SECTION(bool operator==(const PeakShape &rhs) const)
	PeakShape p1,p2;
	TEST_EQUAL(p1==p2,true)
	
	p1.mz_position = 14.4;
	TEST_EQUAL(p1==p2,false)
	
	p2.mz_position = 14.4;
	TEST_EQUAL(p1==p2,true)
END_SECTION

START_SECTION(bool operator!=(const PeakShape &rhs) const)
	PeakShape p1,p2;
	TEST_EQUAL(p1!=p2,false)
	
	p1.mz_position = 14.4;
	TEST_EQUAL(p1!=p2,true)
	
	p2.mz_position = 14.4;
	TEST_EQUAL(p1!=p2,false)
END_SECTION

START_SECTION(([PeakShape::PositionLess] bool operator()(const PeakShape &a, const PeakShape &b)))
  PeakShape p1(0.,123.,0.,0.,0.,PeakShape::LORENTZ_PEAK);
  PeakShape p2(0.,124.,0.,0.,0.,PeakShape::LORENTZ_PEAK);
  PeakShape::PositionLess comp;
  TEST_EQUAL(comp(p1,p2),true);
  TEST_EQUAL(comp(p2,p1),false);
END_SECTION

  
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

