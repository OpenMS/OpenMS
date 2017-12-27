// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

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

SpectraSTSimilarityScore* ptr = nullptr;
SpectraSTSimilarityScore* nullPointer = nullptr;

START_SECTION(SpectraSTSimilarityScore())
	ptr = new SpectraSTSimilarityScore();
	TEST_NOT_EQUAL(ptr, nullPointer)
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

START_SECTION(double operator () (const PeakSpectrum& spec) const)
	PeakMap exp;
	PeakSpectrum s1;
	MSPFile msp;
	std::vector< PeptideIdentification > ids;
  const String filename(OPENMS_GET_TEST_DATA_PATH("SpectraSTSimilarityScore_1.msp"));
  msp.load(filename, ids, exp);
	for(Size k = 0; k < exp[0].size(); ++k)
	{
			Peak1D peak;
			peak.setIntensity(exp[0][k].getIntensity());
			peak.setMZ(exp[0][k].getMZ());
			peak.setPosition(exp[0][k].getPosition());
			s1.push_back(peak);
	}	
  double score = (*ptr)(s1);
  TEST_REAL_SIMILAR(score, 1);
END_SECTION

START_SECTION(double operator () (const PeakSpectrum& spec1, const PeakSpectrum& spec2) const)
  PeakSpectrum s1, s2, s3;
	PeakMap exp;
	MSPFile msp;
	std::vector< PeptideIdentification > ids;
  const String filename(OPENMS_GET_TEST_DATA_PATH("SpectraSTSimilarityScore_1.msp"));
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

  double score = (*ptr)(s1, s2);
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

START_SECTION((double operator()(const BinnedSpectrum &bin1, const BinnedSpectrum &bin2) const))
  PeakSpectrum s1, s2, s3;
	PeakMap exp;
	MSPFile msp;
	std::vector< PeptideIdentification > ids;
  const String filename(OPENMS_GET_TEST_DATA_PATH("SpectraSTSimilarityScore_1.msp"));
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

  double score = (*ptr)(ptr->transform(s1), ptr->transform(s2));
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

START_SECTION(bool preprocess(PeakSpectrum &spec, float remove_peak_intensity_threshold=2.01, UInt cut_peaks_below=1000, Size min_peak_number=5, Size max_peak_number=150))
	PeakSpectrum s1, s2, s3;
	PeakMap exp;
	MSPFile msp;
	std::vector< PeptideIdentification > ids;
  const String filename(OPENMS_GET_TEST_DATA_PATH("SpectraSTSimilarityScore_1.msp"));
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

START_SECTION(double delta_D(double top_hit, double runner_up))
SpectraSTSimilarityScore spectrast;
TEST_EXCEPTION( Exception::DivisionByZero, spectrast.delta_D(0,5))
TEST_REAL_SIMILAR(spectrast.delta_D(5,4),0.2)
TEST_REAL_SIMILAR(spectrast.delta_D(25,1),0.96)
END_SECTION

START_SECTION((double compute_F(double dot_product, double delta_D, double dot_bias)))
//pretty straightforward function
NOT_TESTABLE
END_SECTION

START_SECTION(double dot_bias(const BinnedSpectrum &bin1, const BinnedSpectrum &bin2, double dot_product=-1) const)
	PeakSpectrum s1,s2;
	Peak1D peak;
	peak.setIntensity(1);
	peak.setMZ(1);
	s1.push_back(peak);
		peak.setIntensity(0);
	peak.setMZ(2);
	s1.push_back(peak);
		peak.setIntensity(2);
	peak.setMZ(3);
	s1.push_back(peak);
		peak.setIntensity(3);
	peak.setMZ(4);
	s1.push_back(peak);
	
		peak.setIntensity(0);
	peak.setMZ(1);
	s2.push_back(peak);
		peak.setIntensity(4);
	peak.setMZ(2);
	s2.push_back(peak);
		peak.setIntensity(5);
	peak.setMZ(3);
	s2.push_back(peak);
		peak.setIntensity(6);
	peak.setMZ(4);
	s2.push_back(peak);
		peak.setIntensity(0);
	peak.setMZ(5);
	s2.push_back(peak);
	BinnedSpectrum bin(1,1,s1);

	BinnedSpectrum bin2(1,1,s2);

		
	TEST_REAL_SIMILAR(ptr->dot_bias(bin,bin2,1), 98.585 );
	TEST_REAL_SIMILAR(ptr->dot_bias(bin2,bin,1), 98.585);
END_SECTION
START_SECTION(BinnedSpectrum transform(const PeakSpectrum& spec))
	PeakSpectrum s1;
	Peak1D peak;
	peak.setIntensity(1);
	peak.setMZ(1);
	s1.push_back(peak);
		peak.setIntensity(0);
	peak.setMZ(2);
	s1.push_back(peak);
		peak.setIntensity(2);
	peak.setMZ(3);
	s1.push_back(peak);
		peak.setIntensity(3);
	peak.setMZ(4);
	s1.push_back(peak);
	BinnedSpectrum bin = ptr->transform(s1);
	
	SparseVector<float>::SparseVectorIterator iter = bin.getBins().begin();
	TEST_REAL_SIMILAR((double)*iter,0.1205);
		iter++;
	TEST_REAL_SIMILAR((double)*iter,0.3614);
		iter++;
			TEST_REAL_SIMILAR((double)*iter,0.602);
		iter++;
			TEST_REAL_SIMILAR((double)*iter,0.602);
	delete ptr;
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
