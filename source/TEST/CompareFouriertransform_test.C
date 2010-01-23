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
// $Maintainer: Vipul Patel $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <iostream>

#include <OpenMS/COMPARISON/SPECTRA/CompareFouriertransform.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>

///////////////////////////

START_TEST(CompareFouriertransform, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

CompareFouriertransform* ptr = 0;

START_SECTION(CompareFouriertransform())
	ptr = new CompareFouriertransform();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION(virtual ~CompareFouriertransform())
	delete ptr;
END_SECTION

ptr = new CompareFouriertransform();

START_SECTION(CompareFouriertransform(const CompareFouriertransform& source))
CompareFouriertransform copy(*ptr);
	TEST_EQUAL(copy.getName(), ptr->getName());
	TEST_EQUAL(copy.getParameters(), ptr->getParameters());
END_SECTION

START_SECTION(CompareFouriertransform& operator = (const CompareFouriertransform& source))
CompareFouriertransform copy;
	copy = *ptr;
	TEST_EQUAL(copy.getName(), ptr->getName());
	TEST_EQUAL(copy.getParameters(), ptr->getParameters());
END_SECTION

START_SECTION(double operator () (const PeakSpectrum& ) const)
	
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
	  
	TEST_REAL_SIMILAR(score, 0);
END_SECTION

START_SECTION(void transform(PeakSpectrum & spec) )
	
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
	ptr->transform(spectrum);
	MSSpectrum<>::FloatDataArrays& temp = spectrum.getFloatDataArrays();
	TEST_STRING_SIMILAR("Fouriertransformation", temp[temp.size()-1].getName())  
	
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
	ptr->transform(spectrum1);
	ptr->transform(spectrum2);

   double score = ptr->operator()(spectrum1, spectrum2);
   TEST_REAL_SIMILAR(score, 1.0)
END_SECTION

START_SECTION(static PeakSpectrumCompareFunctor* create())
	PeakSpectrumCompareFunctor* psf = CompareFouriertransform::create();
  CompareFouriertransform cft;
	TEST_EQUAL(psf->getParameters(), cft.getParameters())
	TEST_EQUAL(psf->getName(), cft.getName())
END_SECTION

START_SECTION(static const String getProductName())
	TEST_EQUAL(ptr->getProductName(), "CompareFouriertransform")
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
