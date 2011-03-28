// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

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

SpectrumAlignmentScore* ptr = 0;
SpectrumAlignmentScore* nullPointer = 0;

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
