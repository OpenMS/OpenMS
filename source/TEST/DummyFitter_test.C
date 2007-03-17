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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/DummyFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>
#include <OpenMS/CONCEPT/Exception.h>

///////////////////////////

START_TEST(DummyFitter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using std::stringstream;

enum
{
	RT = RawDataPoint2D::RT,
	MZ = RawDataPoint2D::MZ
};

// default ctor
DummyFitter* ptr = 0;
CHECK((DummyFitter()))
	ptr = new DummyFitter();
  TEST_EQUAL(ptr->getName(), "DummyFitter")
	TEST_NOT_EQUAL(ptr, 0)
RESULT

// destructor
CHECK((virtual ~DummyFitter()))
	delete ptr;
RESULT

CHECK((static const String getProductName()))
	TEST_EQUAL(DummyFitter::getProductName(),"DummyFitter")
	TEST_EQUAL(DummyFitter().getName(),"DummyFitter")
RESULT

CHECK((static BaseModelFitter* create()))
	TEST_NOT_EQUAL(DummyFitter::create(),0)
RESULT

CHECK((DummyFitter& operator=(const DummyFitter &rhs)))
	DummyFitter ms1;
	DummyFitter ms2;
	
	ms1 = ms2;
	
	TEST_EQUAL(ms1 == ms2, true)
RESULT

CHECK((DummyFitter(const DummyFitter &rhs)))
	DummyFitter ms1;
	DummyFitter ms2(ms1);
		
	TEST_EQUAL(ms1 == ms2, true)
RESULT

CHECK(([EXTRA]void DummyFitter::setParameters(const Param& param)))
	DummyFitter* fitter = new DummyFitter();
	
	// check default params
	Param p1 = fitter->getParameters();
	TEST_EQUAL(p1.getValue("min_num_peaks:final"),DataValue(5))
	TEST_EQUAL(p1.getValue("min_num_peaks:extended"),DataValue(10))
	TEST_EQUAL(p1.getValue("use_fwhm_intensity"),DataValue(0))
	
	// change default settings
	Param p2;
	p2.setValue("min_num_peaks:final",20);
	p2.setValue("min_num_peaks:extended",34);
	fitter->setParameters(p2);
	
	// check changes
	Param p3 = fitter->getParameters();
	TEST_EQUAL(p3.getValue("min_num_peaks:final"),DataValue(20))
	TEST_EQUAL(p3.getValue("min_num_peaks:extended"),DataValue(34))

RESULT

CHECK((Feature fit(const IndexSet &range)))
	
	// Test feature construction
	PRECISION(0.1)

	FeaFiTraits traits;
	double mzs[] = {675, 675.5, 676, 676.5, 677, 677.5, 678};
	const UInt mz_num = 7;
	double rts[] = { 1260, 1260.5, 1261, 1261.5, 1262, 1262.5,
											1263, 1263.5, 1264, 1264.5, 1265};
	const UInt rt_num = 11;

	double intens[] = { 1.65879841, 6.652431187, 19.59411554, 42.38668296, 67.34288093,
		78.58007608, 67.34288093, 42.38668296, 19.59411554, 6.652431187, 1.65879841,
		20.20830161, 81.04320276, 238.7051942, 516.3755092, 820.4042402, 957.3013023,
		820.4042402, 516.3755092, 238.7051942, 81.04320276, 20.20830161, 90.56732447,
		363.210436, 1069.80246, 2314.234476, 3676.796717, 4290.326784, 3676.796717,
		2314.234476, 1069.80246, 363.210436, 90.56732447, 149.3202743, 598.8327716,
		1763.806071, 3815.527605, 6062.012955, 7073.553026, 6062.012955, 3815.527605,
		1763.806071, 598.8327716, 149.3202743, 90.56732447, 363.210436, 1069.80246,
		2314.234476, 3676.796717, 4290.326784, 3676.796717, 2314.234476, 1069.80246,
		363.210436, 90.56732447, 20.20830161, 81.04320276, 238.7051942, 516.3755092,
		820.4042402, 957.3013023, 820.4042402, 516.3755092, 238.7051942, 81.04320276,
		20.20830161, 1.65879841, 6.652431187, 19.59411554, 42.38668296, 67.34288093,
		78.58007608, 67.34288093, 42.38668296, 19.59411554, 6.652431187, 1.65879841};

	Peak2D p;
	DPeakArray<2, Peak2D> peak_array;
	for (UInt mz=0; mz<mz_num; mz++) 
	{
		for (UInt rt=0; rt<rt_num; rt++)
		{
			p.setMZ(mzs[mz]);
			p.setRT(rts[rt]);
			p.setIntensity(intens[mz*rt_num+rt]);
			peak_array.push_back(p);
		}
	}
	peak_array.sortByPosition();
		
	MSExperimentExtern<Peak1D > exp;
	exp.set2DData(peak_array);
		
	traits.setData(exp.begin(), exp.end(),100);
	
	DummyFitter fitter;
	fitter.setTraits(&traits);
	
	FeaFiModule::IndexSet  set;
	for (UInt mz=0; mz<mz_num; mz++) 
	{
		for (UInt rt=0; rt<rt_num; rt++)
		{
			set.insert(std::make_pair(rt,mz));
		}
	}
	Feature feature = fitter.fit(set);

	// rt and m/z should be determined by point with max. ion count
	TEST_REAL_EQUAL(feature.getMZ(), 676.5);	
	TEST_REAL_EQUAL(feature.getRT(), 1262.5);
	
	// intensity should be sum of ion counts by default
	TEST_REAL_EQUAL(feature.getIntensity(), 79820.9);
	// charge estimate is always zero
	TEST_EQUAL(feature.getCharge(), 0);
	
	// test quality
	PRECISION(0.01)
	TEST_REAL_EQUAL(feature.getOverallQuality(), 1.0);
	TEST_REAL_EQUAL(feature.getQuality(0), 0.0);
	TEST_REAL_EQUAL(feature.getQuality(1), 0.0);
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
