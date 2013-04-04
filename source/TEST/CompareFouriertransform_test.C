// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Erhan Kenar $
// $Authors: Vipul Patel $
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
CompareFouriertransform* nullPointer = 0;

START_SECTION(CompareFouriertransform())
	ptr = new CompareFouriertransform();
	TEST_NOT_EQUAL(ptr, nullPointer)
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
