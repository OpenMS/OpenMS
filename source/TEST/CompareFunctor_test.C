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
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/COMPARISON/SPECTRA/CompareFunctor.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterSpectrum.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(CompareFunctor, "$Id: $")

/////////////////////////////////////////////////////////////

CompareFunctor* e_ptr = 0;
CHECK(CompareFunctor())
	e_ptr = new CompareFunctor;
	TEST_NOT_EQUAL(e_ptr, 0)
RESULT

CHECK(~CompareFunctor())
	delete e_ptr;
RESULT

e_ptr = new CompareFunctor();

CHECK(CompareFunctor(const CompareFunctor& source))
	CompareFunctor copy(*e_ptr);
	TEST_EQUAL(*e_ptr == copy, true)
RESULT

CHECK(double operator () (const ClusterSpectrum& csa, const ClusterSpectrum& csb))
	DTAFile dta_file;
	PeakSpectrum spec;
	dta_file.load("data/spectrum.dta", spec);

	DTAFile dta_file2;
	PeakSpectrum spec2;
	dta_file2.load("data/spectrum2.dta", spec2);

	double filter = e_ptr->filter(spec, spec);
	
	TEST_REAL_EQUAL(filter, 1.0)

	filter = e_ptr->filter(spec, spec2);

	TEST_REAL_EQUAL(filter, 1.0)
RESULT

CHECK(bool usebins() const)
	TEST_EQUAL(e_ptr->usebins(), false)
RESULT

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
