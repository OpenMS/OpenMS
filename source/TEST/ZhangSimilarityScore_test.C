// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

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

ZhangSimilarityScore* ptr = 0;

CHECK(ZhangSimilarityScore())
	ptr = new ZhangSimilarityScore();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~ZhangSimilarityScore())
	delete ptr;
RESULT

ptr = new ZhangSimilarityScore();

CHECK(ZhangSimilarityScore(const ZhangSimilarityScore& source))
	ZhangSimilarityScore copy(*ptr);
	TEST_EQUAL(copy.getName(), ptr->getName());
	TEST_EQUAL(copy.getParameters(), ptr->getParameters());
RESULT

CHECK(ZhangSimilarityScore& operator = (const ZhangSimilarityScore& source))
	ZhangSimilarityScore copy;
	copy = *ptr;
	TEST_EQUAL(copy.getName(), ptr->getName());
	TEST_EQUAL(copy.getParameters(), ptr->getParameters());
RESULT

CHECK(double operator () (const PeakSpectrum& spec) const)
	PeakSpectrum s1;
  DTAFile().load("data/PILISSequenceDB_DFPIANGER_1.dta", s1);

  Normalizer normalizer;
  Param p(normalizer.getParameters());
  p.setValue("method", "to_one");
  normalizer.setParameters(p);
  normalizer.filterSpectrum(s1);

  double score = (*ptr)(s1);
  TEST_REAL_EQUAL(score, 1.82682);
RESULT

CHECK(double operator () (const PeakSpectrum& spec1, const PeakSpectrum& spec2) const)
  PeakSpectrum s1, s2;
  DTAFile().load("data/PILISSequenceDB_DFPIANGER_1.dta", s1);
  DTAFile().load("data/PILISSequenceDB_DFPIANGER_1.dta", s2);

  Normalizer normalizer;
  Param p(normalizer.getParameters());
  p.setValue("method", "to_one");
  normalizer.setParameters(p);
  normalizer.filterSpectrum(s1);
  normalizer.filterSpectrum(s2);

  PRECISION(0.01)

  double score = (*ptr)(s1, s2);
  TEST_REAL_EQUAL(score, 1.82682)

  s2.getContainer().resize(100);

  score = (*ptr)(s1, s2);

  normalizer.filterSpectrum(s2);
  TEST_REAL_EQUAL(score, 0.328749)
RESULT

CHECK(static PeakSpectrumCompareFunctor* create())
	PeakSpectrumCompareFunctor* psf = ZhangSimilarityScore::create();
	ZhangSimilarityScore zhang;
	TEST_EQUAL(psf->getParameters(), zhang.getParameters())
	TEST_EQUAL(psf->getName(), zhang.getName())
RESULT

CHECK(static const String getProductName())
	TEST_EQUAL(ptr->getProductName(), "ZhangSimilarityScore")
RESULT
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
