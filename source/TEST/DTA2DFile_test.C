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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FORMAT/DTA2DFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/Peak2D.h>

using namespace OpenMS;
using namespace std;

DRange<1> makeRange(float a, float b)
{
	DPosition<1> pa(a), pb(b);
	return DRange<1>(pa, pb);
}

///////////////////////////

START_TEST(DTAFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DTA2DFile* ptr = 0;
CHECK((DTA2DFile()))
	ptr = new DTA2DFile;
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~DTA2DFile()))
	delete ptr;
RESULT

CHECK(const PeakFileOptions& getOptions() const)
	DTA2DFile file;
	TEST_EQUAL(file.getOptions().hasMSLevels(),false)
RESULT

CHECK(PeakFileOptions& getOptions())
	DTA2DFile file;
	file.getOptions().addMSLevel(1);
	TEST_EQUAL(file.getOptions().hasMSLevels(),true);
RESULT

CHECK((template<typename MapType> void load(const String& filename, MapType& map) ))
	PRECISION(0.01)

	MSExperiment<> e;
	DTA2DFile dta;

	//test exception
	TEST_EXCEPTION( Exception::FileNotFound , dta.load("dummy/dummy.dta2d",e) )

	// real test
	dta.load("data/DTA2DFile_test_1.dta2d",e);
	TEST_EQUAL(e.size(), 9);
	ABORT_IF(e.size() != 9)

	MSExperiment<>::const_iterator it(e.begin());
	TEST_REAL_EQUAL(it->getContainer()[0].getPosition()[0], 230.02)
	TEST_REAL_EQUAL(it->getRT(), 4711.1)
	TEST_REAL_EQUAL(it->getContainer()[0].getIntensity(), 47218.89)
	++it;

	TEST_REAL_EQUAL(it->getContainer()[0].getPosition()[0], 231.51)
	TEST_REAL_EQUAL(it->getRT(), 4711.2)
	TEST_REAL_EQUAL(it->getContainer()[0].getIntensity(), 89935.22)
	++it;
		
	TEST_REAL_EQUAL(it->getContainer()[0].getPosition()[0], 139.42)
	TEST_REAL_EQUAL(it->getRT(), 4711.3)
	TEST_REAL_EQUAL(it->getContainer()[0].getIntensity(), 318.52)
	++it;

	TEST_REAL_EQUAL(it->getContainer()[0].getPosition()[0], 149.93)
	TEST_REAL_EQUAL(it->getRT(), 4711.4)
	TEST_REAL_EQUAL(it->getContainer()[0].getIntensity(), 61870.99)
	++it;
		
	TEST_REAL_EQUAL(it->getContainer()[0].getPosition()[0], 169.65)
	TEST_REAL_EQUAL(it->getRT(), 4711.5)
	TEST_REAL_EQUAL(it->getContainer()[0].getIntensity(), 62074.22)
	++it;
		
	TEST_REAL_EQUAL(it->getContainer()[0].getPosition()[0], 189.30)
	TEST_REAL_EQUAL(it->getRT(), 4711.6)
	TEST_REAL_EQUAL(it->getContainer()[0].getIntensity(), 53737.85)
	++it;

	TEST_REAL_EQUAL(it->getContainer()[0].getPosition()[0], 202.28)
	TEST_REAL_EQUAL(it->getRT(), 4711.7)
	TEST_REAL_EQUAL(it->getContainer()[0].getIntensity(), 49410.25)
	++it;
		
	TEST_REAL_EQUAL(it->getContainer()[0].getPosition()[0], 207.82)
	TEST_REAL_EQUAL(it->getRT(), 4711.8)
	TEST_REAL_EQUAL(it->getContainer()[0].getIntensity(), 17038.71)
	++it;
		
	TEST_REAL_EQUAL(it->getContainer()[0].getPosition()[0], 219.72)
	TEST_REAL_EQUAL(it->getRT(), 4711.9)
	TEST_REAL_EQUAL(it->getContainer()[0].getIntensity(), 73629.98)


	dta.load("data/DTA2DFile_test_2.dta2d",e);
	DPeakArray<Peak2D> array;
	e.get2DData(array);
	TEST_EQUAL(array.size(), 11);
	ABORT_IF(array.size() != 11)

	DPeakArray<Peak2D>::ConstIterator it2 = array.begin();

	TEST_REAL_EQUAL(it2->getMZ(), 230.02)
	TEST_REAL_EQUAL(it2->getRT(), 4711.1)
	TEST_REAL_EQUAL(it2->getIntensity(), 47218.89)
	++it2;

	TEST_REAL_EQUAL(it2->getMZ(), 430.02)
	TEST_REAL_EQUAL(it2->getRT(), 4711.1)
	TEST_REAL_EQUAL(it2->getIntensity(), 47219.89)
	++it2;

	TEST_REAL_EQUAL(it2->getMZ(), 630.02)
	TEST_REAL_EQUAL(it2->getRT(), 4711.1)
	TEST_REAL_EQUAL(it2->getIntensity(), 47210.89)
	++it2;

	TEST_REAL_EQUAL(it2->getMZ(), 231.51)
	TEST_REAL_EQUAL(it2->getRT(), 4711.2)
	TEST_REAL_EQUAL(it2->getIntensity(), 89935.22)
	++it2;

	TEST_REAL_EQUAL(it2->getMZ(), 139.42)
	TEST_REAL_EQUAL(it2->getRT(), 4711.3)
	TEST_REAL_EQUAL(it2->getIntensity(), 318.52)
	++it2;

	TEST_REAL_EQUAL(it2->getMZ(), 149.93)
	TEST_REAL_EQUAL(it2->getRT(), 4711.4)
	TEST_REAL_EQUAL(it2->getIntensity(), 61870.99)
	++it2;

	TEST_REAL_EQUAL(it2->getMZ(), 169.65)
	TEST_REAL_EQUAL(it2->getRT(), 4711.5)
	TEST_REAL_EQUAL(it2->getIntensity(), 62074.22)
	++it2;

	TEST_REAL_EQUAL(it2->getMZ(), 189.30)
	TEST_REAL_EQUAL(it2->getRT(), 4711.6)
	TEST_REAL_EQUAL(it2->getIntensity(), 53737.85)
	++it2;

	TEST_REAL_EQUAL(it2->getMZ(), 202.28)
	TEST_REAL_EQUAL(it2->getRT(), 4711.7)
	TEST_REAL_EQUAL(it2->getIntensity(), 49410.25)
	++it2;

	TEST_REAL_EQUAL(it2->getMZ(), 207.82)
	TEST_REAL_EQUAL(it2->getRT(), 4711.8)
	TEST_REAL_EQUAL(it2->getIntensity(), 17038.71)
	++it2;

	TEST_REAL_EQUAL(it2->getMZ(), 219.72)
	TEST_REAL_EQUAL(it2->getRT(), 4711.9)
	TEST_REAL_EQUAL(it2->getIntensity(), 73629.98)

	// Test with DPeak 

	MSExperiment<Peak1D > e3;
	dta.load("data/DTA2DFile_test_1.dta2d",e3);
	TEST_EQUAL(e3.size(), 9);
	ABORT_IF(e3.size() != 9)

	MSExperiment<Peak1D >::const_iterator it3(e3.begin());
	TEST_EQUAL(it3->size(), 3);
	ABORT_IF(it3->size() != 3)
	TEST_REAL_EQUAL(it3->getRT(), 4711.1)
	TEST_REAL_EQUAL(it3->getContainer()[0].getPosition()[0], 230.02)
	TEST_REAL_EQUAL(it3->getContainer()[0].getIntensity(), 47218.89)
	TEST_REAL_EQUAL(it3->getContainer()[1].getPosition()[0], 430.02)
	TEST_REAL_EQUAL(it3->getContainer()[1].getIntensity(), 47219.89)
	TEST_REAL_EQUAL(it3->getContainer()[2].getPosition()[0], 630.02)
	TEST_REAL_EQUAL(it3->getContainer()[2].getIntensity(), 47210.89)
	++it3;

	TEST_REAL_EQUAL(it3->getContainer()[0].getPosition()[0], 231.51)
	TEST_REAL_EQUAL(it3->getRT(), 4711.2)
	TEST_REAL_EQUAL(it3->getContainer()[0].getIntensity(), 89935.22)
	++it3;
		
	TEST_REAL_EQUAL(it3->getContainer()[0].getPosition()[0], 139.42)
	TEST_REAL_EQUAL(it3->getRT(), 4711.3)
	TEST_REAL_EQUAL(it3->getContainer()[0].getIntensity(), 318.52)
	++it3;

	TEST_REAL_EQUAL(it3->getContainer()[0].getPosition()[0], 149.93)
	TEST_REAL_EQUAL(it3->getRT(), 4711.4)
	TEST_REAL_EQUAL(it3->getContainer()[0].getIntensity(), 61870.99)
	++it3;
		
	TEST_REAL_EQUAL(it3->getContainer()[0].getPosition()[0], 169.65)
	TEST_REAL_EQUAL(it3->getRT(), 4711.5)
	TEST_REAL_EQUAL(it3->getContainer()[0].getIntensity(), 62074.22)
	++it3;
		
	TEST_REAL_EQUAL(it3->getContainer()[0].getPosition()[0], 189.30)
	TEST_REAL_EQUAL(it3->getRT(), 4711.6)
	TEST_REAL_EQUAL(it3->getContainer()[0].getIntensity(), 53737.85)
	++it3;

	TEST_REAL_EQUAL(it3->getContainer()[0].getPosition()[0], 202.28)
	TEST_REAL_EQUAL(it3->getRT(), 4711.7)
	TEST_REAL_EQUAL(it3->getContainer()[0].getIntensity(), 49410.25)
	++it3;
		
	TEST_REAL_EQUAL(it3->getContainer()[0].getPosition()[0], 207.82)
	TEST_REAL_EQUAL(it3->getRT(), 4711.8)
	TEST_REAL_EQUAL(it3->getContainer()[0].getIntensity(), 17038.71)
	++it3;
		
	TEST_REAL_EQUAL(it3->getContainer()[0].getPosition()[0], 219.72)
	TEST_REAL_EQUAL(it3->getRT(), 4711.9)
	TEST_REAL_EQUAL(it3->getContainer()[0].getIntensity(), 73629.98)


  MSExperiment<> e4;

  dta.load("data/DTA2DFile_test_3.dta2d",e4);
  TEST_EQUAL(e4.size(),9)
	TEST_REAL_EQUAL(e4[0].getRT(), 282666)
	TEST_REAL_EQUAL(e4[1].getRT(), 282672)
	TEST_REAL_EQUAL(e4[2].getRT(), 282678)
	TEST_REAL_EQUAL(e4[3].getRT(), 282684)
	TEST_REAL_EQUAL(e4[4].getRT(), 282690)
	
RESULT

CHECK((template<typename MapType> void store(const String& filename, const MapType& map) const ))
	PRECISION(0.1)
	std::string tmp_filename;
  MSExperiment<> e;
  DTA2DFile f;

  NEW_TMP_FILE(tmp_filename);
  f.load("data/DTA2DFile_test_1.dta2d",e);
	f.store(tmp_filename,e);

	MSExperiment<> e2;
	f.load(tmp_filename,e2);
	DPeakArray<Peak2D> array;
	e2.get2DData(array);
	TEST_EQUAL(array.size(), 11);
	ABORT_IF(array.size() != 11)

	DPeakArray<Peak2D>::ConstIterator it2 = array.begin();

	TEST_REAL_EQUAL(it2->getMZ(), 230.02)
	TEST_REAL_EQUAL(it2->getRT(), 4711.1)
	TEST_REAL_EQUAL(it2->getIntensity(), 47218.89)
	++it2;

	TEST_REAL_EQUAL(it2->getMZ(), 430.02)
	TEST_REAL_EQUAL(it2->getRT(), 4711.1)
	TEST_REAL_EQUAL(it2->getIntensity(), 47219.89)
	++it2;

	TEST_REAL_EQUAL(it2->getMZ(), 630.02)
	TEST_REAL_EQUAL(it2->getRT(), 4711.1)
	TEST_REAL_EQUAL(it2->getIntensity(), 47210.89)
	++it2;

	TEST_REAL_EQUAL(it2->getMZ(), 231.51)
	TEST_REAL_EQUAL(it2->getRT(), 4711.2)
	TEST_REAL_EQUAL(it2->getIntensity(), 89935.22)
	++it2;

	TEST_REAL_EQUAL(it2->getMZ(), 139.42)
	TEST_REAL_EQUAL(it2->getRT(), 4711.3)
	TEST_REAL_EQUAL(it2->getIntensity(), 318.52)
	++it2;

	TEST_REAL_EQUAL(it2->getMZ(), 149.93)
	TEST_REAL_EQUAL(it2->getRT(), 4711.4)
	TEST_REAL_EQUAL(it2->getIntensity(), 61870.99)
	++it2;

	TEST_REAL_EQUAL(it2->getMZ(), 169.65)
	TEST_REAL_EQUAL(it2->getRT(), 4711.5)
	TEST_REAL_EQUAL(it2->getIntensity(), 62074.22)
	++it2;

	TEST_REAL_EQUAL(it2->getMZ(), 189.30)
	TEST_REAL_EQUAL(it2->getRT(), 4711.6)
	TEST_REAL_EQUAL(it2->getIntensity(), 53737.85)
	++it2;

	TEST_REAL_EQUAL(it2->getMZ(), 202.28)
	TEST_REAL_EQUAL(it2->getRT(), 4711.7)
	TEST_REAL_EQUAL(it2->getIntensity(), 49410.25)
	++it2;

	TEST_REAL_EQUAL(it2->getMZ(), 207.82)
	TEST_REAL_EQUAL(it2->getRT(), 4711.8)
	TEST_REAL_EQUAL(it2->getIntensity(), 17038.71)
	++it2;

	TEST_REAL_EQUAL(it2->getMZ(), 219.72)
	TEST_REAL_EQUAL(it2->getRT(), 4711.9)
	TEST_REAL_EQUAL(it2->getIntensity(), 73629.98)

	MSExperiment< Peak1D > e3;
	f.load(tmp_filename,e3);
	DPeakArray<Peak2D > array2;
	e2.get2DData(array2);
	TEST_EQUAL(array2.size(), 11);
	ABORT_IF(array2.size() != 11)

	DPeakArray<Peak2D >::ConstIterator it3 = array2.begin();

	TEST_REAL_EQUAL(it3->getMZ(), 230.02)
	TEST_REAL_EQUAL(it3->getRT(), 4711.1)
	TEST_REAL_EQUAL(it3->getIntensity(), 47218.89)
	++it3;

	TEST_REAL_EQUAL(it3->getMZ(), 430.02)
	TEST_REAL_EQUAL(it3->getRT(), 4711.1)
	TEST_REAL_EQUAL(it3->getIntensity(), 47219.89)
	++it3;

	TEST_REAL_EQUAL(it3->getMZ(), 630.02)
	TEST_REAL_EQUAL(it3->getRT(), 4711.1)
	TEST_REAL_EQUAL(it3->getIntensity(), 47210.89)
	++it3;

	TEST_REAL_EQUAL(it3->getMZ(), 231.51)
	TEST_REAL_EQUAL(it3->getRT(), 4711.2)
	TEST_REAL_EQUAL(it3->getIntensity(), 89935.22)
	++it3;

	TEST_REAL_EQUAL(it3->getMZ(), 139.42)
	TEST_REAL_EQUAL(it3->getRT(), 4711.3)
	TEST_REAL_EQUAL(it3->getIntensity(), 318.52)
	++it3;

	TEST_REAL_EQUAL(it3->getMZ(), 149.93)
	TEST_REAL_EQUAL(it3->getRT(), 4711.4)
	TEST_REAL_EQUAL(it3->getIntensity(), 61870.99)
	++it3;

	TEST_REAL_EQUAL(it3->getMZ(), 169.65)
	TEST_REAL_EQUAL(it3->getRT(), 4711.5)
	TEST_REAL_EQUAL(it3->getIntensity(), 62074.22)
	++it3;

	TEST_REAL_EQUAL(it3->getMZ(), 189.30)
	TEST_REAL_EQUAL(it3->getRT(), 4711.6)
	TEST_REAL_EQUAL(it3->getIntensity(), 53737.85)
	++it3;

	TEST_REAL_EQUAL(it3->getMZ(), 202.28)
	TEST_REAL_EQUAL(it3->getRT(), 4711.7)
	TEST_REAL_EQUAL(it3->getIntensity(), 49410.25)
	++it3;

	TEST_REAL_EQUAL(it3->getMZ(), 207.82)
	TEST_REAL_EQUAL(it3->getRT(), 4711.8)
	TEST_REAL_EQUAL(it3->getIntensity(), 17038.71)
	++it3;

	TEST_REAL_EQUAL(it3->getMZ(), 219.72)
	TEST_REAL_EQUAL(it3->getRT(), 4711.9)
	TEST_REAL_EQUAL(it3->getIntensity(), 73629.98)


RESULT

CHECK(([EXTRA] load with RT range))
	PRECISION(0.01)

	MSExperiment<> e;
	DTA2DFile dta;

	dta.getOptions().setRTRange(makeRange(282670, 282685));
	dta.load("data/DTA2DFile_test_3.dta2d",e);
	TEST_REAL_EQUAL(e[0].getRT(), 282672)
	TEST_REAL_EQUAL(e[1].getRT(), 282678)
	TEST_REAL_EQUAL(e[2].getRT(), 282684)
RESULT

CHECK(([EXTRA] load with MZ range))
	PRECISION(0.01)

	MSExperiment<> e;
	DTA2DFile dta;

	dta.getOptions().setMZRange(makeRange(150, 220));
	dta.load("data/DTA2DFile_test_1.dta2d",e);

	MSExperiment<>::const_iterator it(e.begin());
	TEST_REAL_EQUAL(it->getContainer()[0].getPosition()[0], 169.65)
	++it;

	TEST_REAL_EQUAL(it->getContainer()[0].getPosition()[0], 189.30)
	++it;

	TEST_REAL_EQUAL(it->getContainer()[0].getPosition()[0], 202.28)
	++it;

	TEST_REAL_EQUAL(it->getContainer()[0].getPosition()[0], 207.82)
	++it;

	TEST_REAL_EQUAL(it->getContainer()[0].getPosition()[0], 219.72)
RESULT

CHECK(([EXTRA] load with intensity range))
	PRECISION(0.01)

	MSExperiment<> e;
	DTA2DFile dta;

	dta.getOptions().setIntensityRange(makeRange(30000, 70000));
	dta.load("data/DTA2DFile_test_1.dta2d",e);

	MSExperiment<>::const_iterator it(e.begin());
	TEST_REAL_EQUAL(it->getContainer()[0].getIntensity(), 47218.89)
	++it;

	TEST_REAL_EQUAL(it->getContainer()[0].getIntensity(), 61870.99)
	++it;

	TEST_REAL_EQUAL(it->getContainer()[0].getIntensity(), 62074.22)
	++it;

	TEST_REAL_EQUAL(it->getContainer()[0].getIntensity(), 53737.85)
	++it;

	TEST_REAL_EQUAL(it->getContainer()[0].getIntensity(), 49410.25)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
