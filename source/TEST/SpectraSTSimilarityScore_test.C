// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: David Wojnar $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <iostream>

#include <OpenMS/COMPARISON/SPECTRA/SpectraSTSimilarityScore.h>
#include <OpenMS/FORMAT/MSPFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/DATASTRUCTURES/String.h>

///////////////////////////

START_TEST(SpectraSTSimilarityScore, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

SpectraSTSimilarityScore* ptr = 0;

START_SECTION(SpectraSTSimilarityScore())
	ptr = new SpectraSTSimilarityScore();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION(~SpectraSTSimilarityScore())
	delete ptr;
END_SECTION
TOLERANCE_ABSOLUTE(0.01)
ptr = new SpectraSTSimilarityScore();

START_SECTION(SpectraSTSimilarityScore(const SpectraSTSimilarityScore& source))
	SpectraSTSimilarityScore copy(*ptr);
	TEST_EQUAL(copy.getName(), ptr->getName());
	TEST_EQUAL(copy.getParameters(), ptr->getParameters());
END_SECTION

START_SECTION(SpectraSTSimilarityScore& operator = (const SpectraSTSimilarityScore& source))
	SpectraSTSimilarityScore copy;
	copy = *ptr;
	TEST_EQUAL(copy.getName(), ptr->getName());
	TEST_EQUAL(copy.getParameters(), ptr->getParameters());
END_SECTION

START_SECTION(DoubleReal operator () (const PeakSpectrum& spec) const)
	RichPeakMap exp;
	PeakSpectrum s1;
	MSPFile msp;
	std::vector< PeptideIdentification > ids;
  const String filename(OPENMS_GET_TEST_DATA_PATH("SpectraSTSimilarityScore_1.MSP"));
  msp.load(filename, ids, exp);
	for(Size k = 0; k < exp[0].size(); ++k)
	{
			Peak1D peak;
			peak.setIntensity(exp[0][k].getIntensity());
			peak.setMZ(exp[0][k].getMZ());
			peak.setPosition(exp[0][k].getPosition());
			s1.push_back(peak);
	}	
  DoubleReal score = (*ptr)(s1);
  TEST_REAL_SIMILAR(score, 1);
END_SECTION

START_SECTION(DoubleReal operator () (const PeakSpectrum& spec1, const PeakSpectrum& spec2) const)
  PeakSpectrum s1, s2, s3;
	RichPeakMap exp;
	MSPFile msp;
	std::vector< PeptideIdentification > ids;
  const String filename(OPENMS_GET_TEST_DATA_PATH("SpectraSTSimilarityScore_1.MSP"));
  msp.load(filename, ids, exp);
	for(Size k = 0; k < exp[0].size(); ++k)
	{
			Peak1D peak;
			peak.setIntensity(exp[0][k].getIntensity());
			peak.setMZ(exp[0][k].getMZ());
			peak.setPosition(exp[0][k].getPosition());
			s1.push_back(peak);
	}	
	for(Size k = 0; k < exp[1].size(); ++k)
	{
			Peak1D peak;
			peak.setIntensity(exp[1][k].getIntensity());
			peak.setMZ(exp[1][k].getMZ());
			peak.setPosition(exp[1][k].getPosition());
			s2.push_back(peak);
	}	
  TOLERANCE_ABSOLUTE(0.01)

  DoubleReal score = (*ptr)(s1, s2);
  TEST_REAL_SIMILAR(score, 1)
  
  for(Size k = 0; k < exp[2].size(); ++k)
	{
			Peak1D peak;
			peak.setIntensity(exp[2][k].getIntensity());
			peak.setMZ(exp[2][k].getMZ());
			peak.setPosition(exp[2][k].getPosition());
			s3.push_back(peak);
	}	
	score = (*ptr)(s1, s3);
  TEST_REAL_SIMILAR(score, 0)
END_SECTION

START_SECTION((DoubleReal operator()(const BinnedSpectrum &bin1, const BinnedSpectrum &bin2) const))
  PeakSpectrum s1, s2, s3;
	RichPeakMap exp;
	MSPFile msp;
	std::vector< PeptideIdentification > ids;
  const String filename(OPENMS_GET_TEST_DATA_PATH("SpectraSTSimilarityScore_1.MSP"));
  msp.load(filename, ids, exp);
	for(Size k = 0; k < exp[0].size(); ++k)
	{
			Peak1D peak;
			peak.setIntensity(exp[0][k].getIntensity());
			peak.setMZ(exp[0][k].getMZ());
			peak.setPosition(exp[0][k].getPosition());
			s1.push_back(peak);
	}	
	for(Size k = 0; k < exp[1].size(); ++k)
	{
			Peak1D peak;
			peak.setIntensity(exp[1][k].getIntensity());
			peak.setMZ(exp[0][k].getMZ());
			peak.setPosition(exp[1][k].getPosition());
			s2.push_back(peak);
	}	
  TOLERANCE_ABSOLUTE(0.01)

  DoubleReal score = (*ptr)(ptr->transform(s1), ptr->transform(s2));
  TEST_REAL_SIMILAR(score, 1)
  
  	for(Size k = 0; k < exp[2].size(); ++k)
	{
			Peak1D peak;
			peak.setIntensity(exp[2][k].getIntensity());
			peak.setMZ(exp[2][k].getMZ());
			peak.setPosition(exp[2][k].getPosition());
			s3.push_back(peak);
	}	
	   score = (*ptr)(ptr->transform(s1), ptr->transform(s3));
  TEST_REAL_SIMILAR(score, 0)
END_SECTION

START_SECTION(bool preprocess(PeakSpectrum &spec, Real remove_peak_intensity_threshold=2.01, UInt cut_peaks_below=1000, Size min_peak_number=5, Size max_peak_number=150))
	PeakSpectrum s1, s2, s3;
	RichPeakMap exp;
	MSPFile msp;
	std::vector< PeptideIdentification > ids;
  const String filename(OPENMS_GET_TEST_DATA_PATH("SpectraSTSimilarityScore_1.MSP"));
  msp.load(filename, ids, exp);
	for(Size k = 0; k < exp[0].size(); ++k)
	{
			Peak1D peak;
			peak.setIntensity(exp[0][k].getIntensity());
			peak.setMZ(exp[0][k].getMZ());
			peak.setPosition(exp[0][k].getPosition());
			s1.push_back(peak);
	}	
	for(Size k = 0; k < exp[1].size(); ++k)
	{
			Peak1D peak;
			peak.setIntensity(exp[1][k].getIntensity());
			peak.setMZ(exp[0][k].getMZ());
			peak.setPosition(exp[1][k].getPosition());
			s2.push_back(peak);
	}	
	  	for(Size k = 0; k < exp[2].size(); ++k)
	{
			Peak1D peak;
			peak.setIntensity(exp[2][k].getIntensity());
			peak.setMZ(exp[2][k].getMZ());
			peak.setPosition(exp[2][k].getPosition());
			s3.push_back(peak);
	}	
  TOLERANCE_ABSOLUTE(0.01)
	ptr->preprocess(s1,2,10000);
	TEST_EQUAL(s1.size(),6)
	
	//min_peaks
	TEST_EQUAL(ptr->preprocess(s2,2,1000,12),false)
	//max_peaks
	ptr->preprocess(s3,1,10000,5,8);
	TEST_EQUAL(s3.size(),8)
END_SECTION



START_SECTION(static PeakSpectrumCompareFunctor* create())
	PeakSpectrumCompareFunctor* psf = SpectraSTSimilarityScore::create();
	SpectraSTSimilarityScore spectrast;
	TEST_EQUAL(psf->getParameters(), spectrast.getParameters())
	TEST_EQUAL(psf->getName(), spectrast.getName())
END_SECTION

START_SECTION(static const String getProductName())
	TEST_EQUAL(ptr->getProductName(), "SpectraSTSimilarityScore")
END_SECTION

START_SECTION(DoubleReal delta_D(DoubleReal top_hit, DoubleReal runner_up))
SpectraSTSimilarityScore spectrast;
TEST_EXCEPTION( Exception::DivisionByZero, spectrast.delta_D(0,5))
TEST_REAL_SIMILAR(spectrast.delta_D(5,4),0.2)
TEST_REAL_SIMILAR(spectrast.delta_D(25,1),0.96)
END_SECTION

START_SECTION((DoubleReal compute_F(DoubleReal dot_product, DoubleReal delta_D, DoubleReal dot_bias)))
//pretty straightforward function
NOT_TESTABLE
END_SECTION

START_SECTION(DoubleReal dot_bias(const BinnedSpectrum &bin1, const BinnedSpectrum &bin2, DoubleReal dot_product=-1) const)
	
END_SECTION
START_SECTION(BinnedSpectrum transform(const PeakSpectrum& spec))


END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
