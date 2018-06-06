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

#include <OpenMS/COMPARISON/SPECTRA/ZhangSimilarityScore.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>

///////////////////////////

START_TEST(ZhangSimilarityScore, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

ZhangSimilarityScore* ptr = nullptr;
ZhangSimilarityScore* nullPointer = nullptr;

START_SECTION(ZhangSimilarityScore())
	ptr = new ZhangSimilarityScore();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~ZhangSimilarityScore())
	delete ptr;
END_SECTION

ptr = new ZhangSimilarityScore();

START_SECTION(ZhangSimilarityScore(const ZhangSimilarityScore& source))
	ZhangSimilarityScore copy(*ptr);
	TEST_EQUAL(copy.getName(), ptr->getName());
	TEST_EQUAL(copy.getParameters(), ptr->getParameters());
END_SECTION

START_SECTION(ZhangSimilarityScore& operator = (const ZhangSimilarityScore& source))
	ZhangSimilarityScore copy;
	copy = *ptr;
	TEST_EQUAL(copy.getName(), ptr->getName());
	TEST_EQUAL(copy.getParameters(), ptr->getParameters());
END_SECTION

START_SECTION(double operator () (const PeakSpectrum& spec) const)
	PeakSpectrum s1;
  DTAFile().load(OPENMS_GET_TEST_DATA_PATH("PILISSequenceDB_DFPIANGER_1.dta"), s1);

  Normalizer normalizer;
  Param p(normalizer.getParameters());
  p.setValue("method", "to_one");
  normalizer.setParameters(p);
  normalizer.filterSpectrum(s1);

  double score = (*ptr)(s1);
  TEST_REAL_SIMILAR(score, 1.82682);
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
  TEST_REAL_SIMILAR(score, 1.82682)

  s2.resize(100);

  score = (*ptr)(s1, s2);

  normalizer.filterSpectrum(s2);
  TEST_REAL_SIMILAR(score, 0.328749)
END_SECTION

START_SECTION(static PeakSpectrumCompareFunctor* create())
	PeakSpectrumCompareFunctor* psf = ZhangSimilarityScore::create();
	ZhangSimilarityScore zhang;
	TEST_EQUAL(psf->getParameters(), zhang.getParameters())
	TEST_EQUAL(psf->getName(), zhang.getName())
END_SECTION

START_SECTION(static const String getProductName())
	TEST_EQUAL(ptr->getProductName(), "ZhangSimilarityScore")
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
