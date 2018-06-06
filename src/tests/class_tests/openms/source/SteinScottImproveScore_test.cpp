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
// $Authors: Vipul Patel $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

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
