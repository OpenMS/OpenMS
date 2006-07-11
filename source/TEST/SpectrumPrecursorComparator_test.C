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
// $Id: $
// $Author: marc_sturm $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/COMPARISON/SPECTRA/SpectrumPrecursorComparator.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(SpectrumPrecursorComparator, "$Id: $")

/////////////////////////////////////////////////////////////

SpectrumPrecursorComparator* e_ptr = 0;
CHECK(SpectrumPrecursorComparator())
	e_ptr = new SpectrumPrecursorComparator;
	TEST_NOT_EQUAL(e_ptr, 0)
RESULT

CHECK(~SpectrumPrecursorComparator())
	delete e_ptr;
RESULT

e_ptr = new SpectrumPrecursorComparator();

CHECK(SpectrumPrecursorComparator(const SpectrumPrecursorComparator& source))
	SpectrumPrecursorComparator copy(*e_ptr);
	TEST_EQUAL(*e_ptr == copy, true)
RESULT

CHECK(double operator () (const ClusterSpectrum& csa, const ClusterSpectrum& csb) const)
	DTAFile dta_file;
	PeakSpectrum spec;
	dta_file.load("data/spectrum.dta", spec);

	ClusterSpectrum csa(spec, 0, 1, 1);

	DTAFile dta_file2;
	PeakSpectrum spec2;
	dta_file2.load("data/spectrum2.dta", spec2);

	ClusterSpectrum csb(spec2, 0, 1, 1);

	double score = (*e_ptr)(csa, csb);

	TEST_REAL_EQUAL(score, 1.7685)

	score = (*e_ptr)(csa, csa);

	TEST_REAL_EQUAL(score, 2)
RESULT

CHECK(bool usebins() const)
	TEST_REAL_EQUAL(e_ptr->usebins(), false)
RESULT

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
