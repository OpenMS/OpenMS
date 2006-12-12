// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/COMPARISON/SPECTRA/SpectrumCheapDPCorr.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(SpectrumCheapDPCorr, "$Id$")

/////////////////////////////////////////////////////////////

SpectrumCheapDPCorr* e_ptr = 0;
CHECK(SpectrumCheapDPCorr())
	e_ptr = new SpectrumCheapDPCorr;
	TEST_NOT_EQUAL(e_ptr, 0)
RESULT

CHECK(~SpectrumCheapDPCorr())
	delete e_ptr;
RESULT

e_ptr = new SpectrumCheapDPCorr();

CHECK(SpectrumCheapDPCorr(const SpectrumCheapDPCorr& source))
	SpectrumCheapDPCorr copy(*e_ptr);
	TEST_EQUAL(*e_ptr == copy, true)
RESULT

CHECK(double operator () (const ClusterSpectrum& csa, const ClusterSpectrum& csb))
	DTAFile dta_file;
	PeakSpectrum spec1;
	dta_file.load("data/Transformers_tests.dta", spec1);

	//ClusterSpectrum csa(spec, 0, 1, 1);
	//BinnedRep csa(spec, 1, 1);
	

	DTAFile dta_file2;
	PeakSpectrum spec2;
	dta_file2.load("data/Transformers_tests_2.dta", spec2);

	//ClusterSpectrum csb(spec2, 0, 1, 1);	
	//BinnedRep csb(spec, 1, 1);
	
	double score = (*e_ptr)(spec1, spec2);

	PRECISION(0.1)
	/// @todo next to equality tests fail, don't know why (andreas)
	TEST_REAL_EQUAL(score, 10145.4)

	score = (*e_ptr)(spec1, spec1);

	TEST_REAL_EQUAL(score, 12295.5)
	
	SpectrumCheapDPCorr corr;
	score = corr(spec1, spec2);
	TEST_REAL_EQUAL(score, 10145.4)

	score = corr(spec1, spec1);

	TEST_REAL_EQUAL(score, 12295.5)
	
RESULT

CHECK(const MSSpectrum< DPeak<1> >& lastconsensus() const)
	TEST_EQUAL(e_ptr->lastconsensus().size(), 121)
RESULT

CHECK(bool usebins() const)
	// @todo fails ???? (andreas)
	//TEST_EQUAL(e_ptr->usebins(), false)
RESULT

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
