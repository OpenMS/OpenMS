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
// $Maintainer: Andreas Bertsch $
// $Authors: Volker Mosthaf, Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/COMPARISON/SPECTRA/SpectrumCheapDPCorr.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(SpectrumCheapDPCorr, "$Id$")

/////////////////////////////////////////////////////////////

SpectrumCheapDPCorr* e_ptr = 0;
START_SECTION(SpectrumCheapDPCorr())
	e_ptr = new SpectrumCheapDPCorr;
	TEST_NOT_EQUAL(e_ptr, 0)
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
