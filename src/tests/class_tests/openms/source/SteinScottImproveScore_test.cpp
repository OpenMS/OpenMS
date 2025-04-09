// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Vipul Patel $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <iostream>

#include <OpenMS/COMPARISON/SteinScottImproveScore.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/PROCESSING/SCALING/Normalizer.h>

///////////////////////////

START_TEST(SteinScottImproveScore, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

SteinScottImproveScore* ptr = nullptr;
SteinScottImproveScore* nullPointer = nullptr;

START_SECTION(SteinScottImproveScore())
	ptr = new SteinScottImproveScore();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(virtual ~SteinScottImproveScore())
	delete ptr;
END_SECTION

ptr = new SteinScottImproveScore();

START_SECTION(SteinScottImproveScore(const SteinScottImproveScore& source))
SteinScottImproveScore copy(*ptr);
	TEST_EQUAL(copy.getName(), ptr->getName());
	TEST_EQUAL(copy.getParameters(), ptr->getParameters());
END_SECTION

START_SECTION(SteinScottImproveScore& operator = (const SteinScottImproveScore& source))
SteinScottImproveScore copy;
	copy = *ptr;
	TEST_EQUAL(copy.getName(), ptr->getName());
	TEST_EQUAL(copy.getParameters(), ptr->getParameters());
END_SECTION

START_SECTION(double operator () (const PeakSpectrum& spec) const)
	
	MSSpectrum spectrum;
	spectrum.setRT(1);
	
		spectrum.setMSLevel(1);
		
		for (float mz=500.0; mz<=900; mz+=100.0)
		    { 
		      Peak1D peak;
		      peak.setMZ(mz);
		      peak.setIntensity(mz);
		      spectrum.push_back(peak);
		      
		    }
  double score = (*ptr)(spectrum);
	if(score >0.99) score =1;  
	TEST_REAL_SIMILAR(score, 1);
END_SECTION

START_SECTION(double operator () (const PeakSpectrum& spec1, const PeakSpectrum& spec2) const)
	MSSpectrum spectrum1,spectrum2;
	spectrum1.setRT(1);
	spectrum2.setRT(1);
	spectrum1.setMSLevel(1);
	spectrum2.setMSLevel(1);
		
	for (float mz=500.0; mz<=900; mz+=100.0)
	    { 
	      Peak1D peak;
	      peak.setMZ(mz);
	      peak.setIntensity(mz);
	      spectrum1.push_back(peak);
	      spectrum2.push_back(peak);
	    }
	
  double score = (*ptr)(spectrum1, spectrum2);
	if(score >0.99) score =1;
  TEST_REAL_SIMILAR(score, 1.0)
END_SECTION

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
