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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/KERNEL/DSpectrum.h>
#include <OpenMS/KERNEL/DPeakArray.h>

///////////////////////////

START_TEST(DTAFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

DTAFile* ptr = 0;
CHECK(DTAFile())
	ptr = new DTAFile;
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~DTAFile())
	delete ptr;
RESULT

CHECK(template<typename SpectrumType> void load(const String& filename, SpectrumType& spectrum) throw(Exception::FileNotFound, Exception::ParseError))
	PRECISION(0.01)
	DSpectrum<1> s;
	DTAFile f1;
	
	TEST_EXCEPTION(Exception::FileNotFound, f1.load("data_Idontexist",s);)

	f1.load("data/DTAFile_test.dta",s);
	
	TEST_EQUAL(s.size(), 25);
	TEST_REAL_EQUAL(s.getPrecursorPeak().getPosition()[0], 582.40666)
	TEST_REAL_EQUAL(s.getPrecursorPeak().getCharge(), 3)

	ABORT_IF(s.size() != 25)
	DSpectrum<1>::ConstIterator it(s.begin());
	
	TEST_REAL_EQUAL(it->getPosition()[0], 139.42)
	TEST_REAL_EQUAL(it->getIntensity(), 318.52)
	++it;

	TEST_REAL_EQUAL(it->getPosition()[0], 149.93)
	TEST_REAL_EQUAL(it->getIntensity(), 61870.99)
	++it;
		
	TEST_REAL_EQUAL(it->getPosition()[0], 169.65)
	TEST_REAL_EQUAL(it->getIntensity(), 62074.22)
	++it;
		
	TEST_REAL_EQUAL(it->getPosition()[0], 189.30)
	TEST_REAL_EQUAL(it->getIntensity(), 53737.85)
	++it;
		
	TEST_REAL_EQUAL(it->getPosition()[0], 202.28)
	TEST_REAL_EQUAL(it->getIntensity(), 49410.25)
	++it;
		
	TEST_REAL_EQUAL(it->getPosition()[0], 207.82)
	TEST_REAL_EQUAL(it->getIntensity(), 17038.71)
	++it;
		
	TEST_REAL_EQUAL(it->getPosition()[0], 219.72)
	TEST_REAL_EQUAL(it->getIntensity(), 73629.98)
	++it;
		
	TEST_REAL_EQUAL(it->getPosition()[0], 230.02)
	TEST_REAL_EQUAL(it->getIntensity(), 47218.89)
	++it;
		
	TEST_REAL_EQUAL(it->getPosition()[0], 231.51)
	TEST_REAL_EQUAL(it->getIntensity(), 89935.22)
	++it;
		
	TEST_REAL_EQUAL(it->getPosition()[0], 263.88)
	TEST_REAL_EQUAL(it->getIntensity(), 81685.67)
	++it;
		
	TEST_REAL_EQUAL(it->getPosition()[0], 318.38)
	TEST_REAL_EQUAL(it->getIntensity(), 11876.49)
	++it;
		
	TEST_REAL_EQUAL(it->getPosition()[0], 362.22)
	TEST_REAL_EQUAL(it->getIntensity(), 91984.30)
	++it;
		
	TEST_REAL_EQUAL(it->getPosition()[0], 389.84)
	TEST_REAL_EQUAL(it->getIntensity(), 35049.17)
	++it;
		
	TEST_REAL_EQUAL(it->getPosition()[0], 489.86)
	TEST_REAL_EQUAL(it->getIntensity(), 55259.42)
	++it;
		
	TEST_REAL_EQUAL(it->getPosition()[0], 508.47)
	TEST_REAL_EQUAL(it->getIntensity(), 63411.68)
	++it;
		
	TEST_REAL_EQUAL(it->getPosition()[0], 521.44)
	TEST_REAL_EQUAL(it->getIntensity(), 42096.49)
	++it;
		
	TEST_REAL_EQUAL(it->getPosition()[0], 546.98)
	TEST_REAL_EQUAL(it->getIntensity(), 7133.92)
	++it;
		
	TEST_REAL_EQUAL(it->getPosition()[0], 562.72)
	TEST_REAL_EQUAL(it->getIntensity(), 47540.11)
	++it;
		
	TEST_REAL_EQUAL(it->getPosition()[0], 579.61)
	TEST_REAL_EQUAL(it->getIntensity(), 63350.29)
	++it;
		
	TEST_REAL_EQUAL(it->getPosition()[0], 600.60)
	TEST_REAL_EQUAL(it->getIntensity(), 29075.85)
	++it;
		
	TEST_REAL_EQUAL(it->getPosition()[0], 628.64)
	TEST_REAL_EQUAL(it->getIntensity(), 41030.62)
	++it;
		
	TEST_REAL_EQUAL(it->getPosition()[0], 629.66)
	TEST_REAL_EQUAL(it->getIntensity(), 51153.79)
	++it;
		
	TEST_REAL_EQUAL(it->getPosition()[0], 630.37)
	TEST_REAL_EQUAL(it->getIntensity(), 31411.63)
	++it;
		
	TEST_REAL_EQUAL(it->getPosition()[0], 652.64)
	TEST_REAL_EQUAL(it->getIntensity(), 26967.56)
	++it;
		
	TEST_REAL_EQUAL(it->getPosition()[0], 712.18)
	TEST_REAL_EQUAL(it->getIntensity(), 29235.72)



	//TEST WITH DPeakArray and DRawDataPoint

	DSpectrum<1, DPeakArray<1,DRawDataPoint<1> > > s2;
	f1.load("data/DTAFile_test.dta",s2);
	
	TEST_EQUAL(s2.size(), 25);
	TEST_REAL_EQUAL(s2.getPrecursorPeak().getPosition()[0], 582.4066)
	TEST_REAL_EQUAL(s2.getPrecursorPeak().getCharge(), 3)

	ABORT_IF(s2.size() != 25)
	DSpectrum<1, DPeakArray<1,DRawDataPoint<1> > >::ConstIterator it2(s2.begin());
	
	TEST_REAL_EQUAL(it2->getPosition()[0], 139.42)
	TEST_REAL_EQUAL(it2->getIntensity(), 318.52)
	++it2;

	TEST_REAL_EQUAL(it2->getPosition()[0], 149.93)
	TEST_REAL_EQUAL(it2->getIntensity(), 61870.99)
	++it2;
		
	TEST_REAL_EQUAL(it2->getPosition()[0], 169.65)
	TEST_REAL_EQUAL(it2->getIntensity(), 62074.22)
	++it2;
		
	TEST_REAL_EQUAL(it2->getPosition()[0], 189.30)
	TEST_REAL_EQUAL(it2->getIntensity(), 53737.85)
	++it2;
		
	TEST_REAL_EQUAL(it2->getPosition()[0], 202.28)
	TEST_REAL_EQUAL(it2->getIntensity(), 49410.25)
	++it2;
		
	TEST_REAL_EQUAL(it2->getPosition()[0], 207.82)
	TEST_REAL_EQUAL(it2->getIntensity(), 17038.71)
	++it2;
		
	TEST_REAL_EQUAL(it2->getPosition()[0], 219.72)
	TEST_REAL_EQUAL(it2->getIntensity(), 73629.98)
	++it2;
		
	TEST_REAL_EQUAL(it2->getPosition()[0], 230.02)
	TEST_REAL_EQUAL(it2->getIntensity(), 47218.89)
	++it2;
		
	TEST_REAL_EQUAL(it2->getPosition()[0], 231.51)
	TEST_REAL_EQUAL(it2->getIntensity(), 89935.22)
	++it2;
		
	TEST_REAL_EQUAL(it2->getPosition()[0], 263.88)
	TEST_REAL_EQUAL(it2->getIntensity(), 81685.67)
	++it2;
		
	TEST_REAL_EQUAL(it2->getPosition()[0], 318.38)
	TEST_REAL_EQUAL(it2->getIntensity(), 11876.49)
	++it2;
		
	TEST_REAL_EQUAL(it2->getPosition()[0], 362.22)
	TEST_REAL_EQUAL(it2->getIntensity(), 91984.30)
	++it2;
		
	TEST_REAL_EQUAL(it2->getPosition()[0], 389.84)
	TEST_REAL_EQUAL(it2->getIntensity(), 35049.17)
	++it2;
		
	TEST_REAL_EQUAL(it2->getPosition()[0], 489.86)
	TEST_REAL_EQUAL(it2->getIntensity(), 55259.42)
	++it2;
		
	TEST_REAL_EQUAL(it2->getPosition()[0], 508.47)
	TEST_REAL_EQUAL(it2->getIntensity(), 63411.68)
	++it2;
		
	TEST_REAL_EQUAL(it2->getPosition()[0], 521.44)
	TEST_REAL_EQUAL(it2->getIntensity(), 42096.49)
	++it2;
		
	TEST_REAL_EQUAL(it2->getPosition()[0], 546.98)
	TEST_REAL_EQUAL(it2->getIntensity(), 7133.92)
	++it2;
		
	TEST_REAL_EQUAL(it2->getPosition()[0], 562.72)
	TEST_REAL_EQUAL(it2->getIntensity(), 47540.11)
	++it2;
		
	TEST_REAL_EQUAL(it2->getPosition()[0], 579.61)
	TEST_REAL_EQUAL(it2->getIntensity(), 63350.29)
	++it2;
		
	TEST_REAL_EQUAL(it2->getPosition()[0], 600.60)
	TEST_REAL_EQUAL(it2->getIntensity(), 29075.85)
	++it2;
		
	TEST_REAL_EQUAL(it2->getPosition()[0], 628.64)
	TEST_REAL_EQUAL(it2->getIntensity(), 41030.62)
	++it2;
		
	TEST_REAL_EQUAL(it2->getPosition()[0], 629.66)
	TEST_REAL_EQUAL(it2->getIntensity(), 51153.79)
	++it2;
		
	TEST_REAL_EQUAL(it2->getPosition()[0], 630.37)
	TEST_REAL_EQUAL(it2->getIntensity(), 31411.63)
	++it2;
		
	TEST_REAL_EQUAL(it2->getPosition()[0], 652.64)
	TEST_REAL_EQUAL(it2->getIntensity(), 26967.56)
	++it2;
		
	TEST_REAL_EQUAL(it2->getPosition()[0], 712.18)
	TEST_REAL_EQUAL(it2->getIntensity(), 29235.72)

RESULT

CHECK(template<typename SpectrumType> void store(const String& filename, const SpectrumType& spectrum) const throw(Exception::UnableToCreateFile))
	String filename;
	NEW_TMP_FILE(filename);
	
	DTAFile dta;
	DSpectrum<1> spec, spec2;
	DSpectrum<1>::PeakType peak;
	
	spec.getPrecursorPeak().getPosition()[0] = 582.40666;
	spec.getPrecursorPeak().setCharge(3);
	
	peak.getPosition()[0] = 11.4;
	peak.setIntensity(11.5);
	spec.getContainer().push_back(peak);
	
	peak.getPosition()[0] = 12.4;
	peak.setIntensity(12.5);
	spec.getContainer().push_back(peak);
	
	peak.getPosition()[0] = 13.4;
	peak.setIntensity(13.5);
	spec.getContainer().push_back(peak);
	
	//store file
	dta.store(filename,spec);
	//load file
	dta.load(filename,spec2);

	TEST_REAL_EQUAL(spec.getPrecursorPeak().getPosition()[0],582.40666)
	TEST_REAL_EQUAL(spec.getPrecursorPeak().getCharge(),3)
	
	ABORT_IF(spec2.getContainer().size() != 3)
	
	DSpectrum<1>::ConstIterator it = spec2.begin();

	TEST_REAL_EQUAL(it->getPosition()[0], 11.4)
	TEST_REAL_EQUAL(it->getIntensity(), 11.5)
	++it;
	
	TEST_REAL_EQUAL(it->getPosition()[0], 12.4)
	TEST_REAL_EQUAL(it->getIntensity(), 12.5)
	++it;

	TEST_REAL_EQUAL(it->getPosition()[0], 13.4)
	TEST_REAL_EQUAL(it->getIntensity(), 13.5)


	//TEST WITH DPeakArray and DRawDataPoint
	
	DSpectrum<1, DPeakArray<1,DRawDataPoint<1> > > raw_spec, raw_spec2;
	DSpectrum<1, DPeakArray<1,DRawDataPoint<1> > >::PeakType raw_peak;
	
	raw_peak.getPosition()[0] = 11.4;
	raw_peak.setIntensity(11.5);
	raw_spec.getContainer().push_back(raw_peak);
	
	raw_peak.getPosition()[0] = 12.4;
	raw_peak.setIntensity(12.5);
	raw_spec.getContainer().push_back(raw_peak);
	
	raw_peak.getPosition()[0] = 13.4;
	raw_peak.setIntensity(13.5);
	raw_spec.getContainer().push_back(raw_peak);
	
	//store file
	dta.store(filename,raw_spec);
	//load file
	dta.load(filename,raw_spec2);
	
	ABORT_IF(raw_spec2.getContainer().size() != 3)
	
	DSpectrum<1, DPeakArray<1,DRawDataPoint<1> > >::ConstIterator it2 = raw_spec2.begin();

	TEST_REAL_EQUAL(it2->getPosition()[0], 11.4)
	TEST_REAL_EQUAL(it2->getIntensity(), 11.5)
	++it2;
	
	TEST_REAL_EQUAL(it2->getPosition()[0], 12.4)
	TEST_REAL_EQUAL(it2->getIntensity(), 12.5)
	++it2;

	TEST_REAL_EQUAL(it2->getPosition()[0], 13.4)
	TEST_REAL_EQUAL(it2->getIntensity(), 13.5)

RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
