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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <iostream>

#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignmentScore.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>
#include <OpenMS/FORMAT/DTAFile.h>

///////////////////////////

START_TEST(SpectrumAlignmentScore, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

SpectrumAlignmentScore* ptr = nullptr;
SpectrumAlignmentScore* nullPointer = nullptr;

START_SECTION(SpectrumAlignmentScore())
	ptr = new SpectrumAlignmentScore();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(double operator () (const PeakSpectrum& spec1, const PeakSpectrum& spec2) const)
	PeakSpectrum s1, s2;
	DTAFile().load(OPENMS_GET_TEST_DATA_PATH("PILISSequenceDB_DFPIANGER_1.dta"), s1);
	DTAFile().load(OPENMS_GET_TEST_DATA_PATH("PILISSequenceDB_DFPIANGER_1.dta"), s2);

	Normalizer normalizer;
	Param p(normalizer.getParameters());
	p.setValue("method", "to_one");
	normalizer.setParameters(p);
	normalizer.filterSpectrum(s1);
	normalizer.filterSpectrum(s2);
	
	TOLERANCE_ABSOLUTE(0.01)

	double score = (*ptr)(s1, s2);
	TEST_REAL_SIMILAR(score, 1.48268)

	s2.resize(100);

	score = (*ptr)(s1, s2);

	normalizer.filterSpectrum(s2);
	TEST_REAL_SIMILAR(score, 3.82472)
END_SECTION

START_SECTION(virtual ~SpectrumAlignmentScore())
	delete ptr;
END_SECTION

ptr = new SpectrumAlignmentScore();

START_SECTION(SpectrumAlignmentScore(const SpectrumAlignmentScore &source))
	SpectrumAlignmentScore sas1;
	Param p(sas1.getParameters());
	p.setValue("tolerance", 0.2);
	sas1.setParameters(p);

	SpectrumAlignmentScore sas2(sas1);

	TEST_EQUAL(sas1.getName(), sas2.getName())
	TEST_EQUAL(sas1.getParameters(), sas2.getParameters())

END_SECTION


START_SECTION(SpectrumAlignmentScore& operator=(const SpectrumAlignmentScore &source))
  SpectrumAlignmentScore sas1;
  Param p(sas1.getParameters());
  p.setValue("tolerance", 0.2);
  sas1.setParameters(p);

  SpectrumAlignmentScore sas2;

	sas2 = sas1;

  TEST_EQUAL(sas1.getName(), sas2.getName())
  TEST_EQUAL(sas1.getParameters(), sas2.getParameters())
END_SECTION

START_SECTION(double operator()(const PeakSpectrum &spec) const)
	PeakSpectrum s1;
	DTAFile().load(OPENMS_GET_TEST_DATA_PATH("PILISSequenceDB_DFPIANGER_1.dta"), s1);

  Normalizer normalizer;
  Param p(normalizer.getParameters());
  p.setValue("method", "to_one");
  normalizer.setParameters(p);
  normalizer.filterSpectrum(s1);

	double score = (*ptr)(s1);
	TEST_REAL_SIMILAR(score, 1.48268);
	
END_SECTION


START_SECTION(static PeakSpectrumCompareFunctor* create())
	PeakSpectrumCompareFunctor* pscf = SpectrumAlignmentScore::create();
	SpectrumAlignmentScore sas;
	TEST_EQUAL(pscf->getParameters(), sas.getParameters())
	TEST_EQUAL(pscf->getName(), sas.getName())
END_SECTION

START_SECTION(static const String getProductName())
	TEST_STRING_EQUAL(SpectrumAlignmentScore::getProductName(), "SpectrumAlignmentScore")
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
