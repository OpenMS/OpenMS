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
// $Authors: Volker Mosthaf, Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/COMPARISON/SPECTRA/SpectrumCheapDPCorr.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(SpectrumCheapDPCorr, "$Id$")

/////////////////////////////////////////////////////////////

SpectrumCheapDPCorr* e_ptr = nullptr;
SpectrumCheapDPCorr* e_nullPointer = nullptr;
START_SECTION(SpectrumCheapDPCorr())
	e_ptr = new SpectrumCheapDPCorr;
	TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION(~SpectrumCheapDPCorr())
	delete e_ptr;
END_SECTION

e_ptr = new SpectrumCheapDPCorr();

START_SECTION(SpectrumCheapDPCorr(const SpectrumCheapDPCorr& source))
	SpectrumCheapDPCorr copy(*e_ptr);
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION(SpectrumCheapDPCorr& operator = (const SpectrumCheapDPCorr& source))
	SpectrumCheapDPCorr copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION(double operator () (const PeakSpectrum& a, const PeakSpectrum& b) const)
	DTAFile dta_file;
	PeakSpectrum spec1;
	dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec1);

	DTAFile dta_file2;
	PeakSpectrum spec2;
	dta_file2.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests_2.dta"), spec2);

	double score = (*e_ptr)(spec1, spec2);

	TOLERANCE_ABSOLUTE(0.1)
	TEST_REAL_SIMILAR(score, 10145.4)

	score = (*e_ptr)(spec1, spec1);

	TEST_REAL_SIMILAR(score, 12295.5)
	
	SpectrumCheapDPCorr corr;
	score = corr(spec1, spec2);
	TEST_REAL_SIMILAR(score, 10145.4)

	score = corr(spec1, spec1);

	TEST_REAL_SIMILAR(score, 12295.5)
	
END_SECTION

START_SECTION(const PeakSpectrum& lastconsensus() const)
	TEST_EQUAL(e_ptr->lastconsensus().size(), 121)
END_SECTION

START_SECTION((Map<UInt, UInt> getPeakMap() const))
	TEST_EQUAL(e_ptr->getPeakMap().size(), 121)
END_SECTION

START_SECTION(double operator () (const PeakSpectrum& a) const)
  DTAFile dta_file;
  PeakSpectrum spec1;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec1);

  double score = (*e_ptr)(spec1);

  TEST_REAL_SIMILAR(score, 12295.5)

END_SECTION

START_SECTION(static PeakSpectrumCompareFunctor* create())
	PeakSpectrumCompareFunctor* cf = SpectrumCheapDPCorr::create();
	SpectrumCheapDPCorr corr;
	TEST_EQUAL(cf->getParameters(), corr.getParameters())
	TEST_EQUAL(cf->getName(), corr.getName())
END_SECTION

START_SECTION(static const String getProductName())
	TEST_EQUAL(SpectrumCheapDPCorr::getProductName(), "SpectrumCheapDPCorr")
END_SECTION

START_SECTION(void setFactor(double f))
	e_ptr->setFactor(0.3);

	TEST_EXCEPTION(Exception::OutOfRange, e_ptr->setFactor(1.1))
END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
