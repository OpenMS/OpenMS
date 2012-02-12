// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Erhan Kenar $
// $Authors: Vipul Patel $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <iostream>

#include <OpenMS/COMPARISON/SPECTRA/SteinScottImproveScore.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>

///////////////////////////

START_TEST(SteinScottImproveScore, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

SteinScottImproveScore* ptr = 0;
SteinScottImproveScore* nullPointer = 0;

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
	
	MSSpectrum<> spectrum;
	spectrum.setRT(1);
	
		spectrum.setMSLevel(1);
		
		for (Real mz=500.0; mz<=900; mz+=100.0)
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
	MSSpectrum<> spectrum1,spectrum2;
	spectrum1.setRT(1);
	spectrum2.setRT(1);
	spectrum1.setMSLevel(1);
	spectrum2.setMSLevel(1);
		
	for (Real mz=500.0; mz<=900; mz+=100.0)
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

START_SECTION(static PeakSpectrumCompareFunctor* create())
	PeakSpectrumCompareFunctor* psf = SteinScottImproveScore::create();
	SteinScottImproveScore stein;
	TEST_EQUAL(psf->getParameters(), stein.getParameters())
	TEST_EQUAL(psf->getName(), stein.getName())
END_SECTION

START_SECTION(static const String getProductName())
	TEST_EQUAL(ptr->getProductName(), "SteinScottImproveScore")
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
