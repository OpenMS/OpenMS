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

#include <OpenMS/FILTERING/TRANSFORMERS/IntensityDistBins.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterSpectrum.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(IntensityDistBins, "$Id$")

/////////////////////////////////////////////////////////////

IntensityDistBins* e_ptr = 0;
CHECK(IntensityDistBins())
	e_ptr = new IntensityDistBins;
	TEST_NOT_EQUAL(e_ptr, 0)
RESULT

CHECK(~IntensityDistBins())
	delete e_ptr;
RESULT

e_ptr = new IntensityDistBins();

CHECK(IntensityDistBins(const IntensityDistBins& source))
	IntensityDistBins copy(*e_ptr);
	TEST_EQUAL(*e_ptr == copy, true)
RESULT

CHECK(std::vector<double> operator () (const ClusterSpectrum& spec))
	DTAFile dta_file;
	PeakSpectrum spec;
	dta_file.load("data/spectrum.dta", spec);
	
	vector<double> filter = (*e_ptr)(spec);

	TEST_EQUAL(filter.size(), 10)

	TEST_REAL_EQUAL(filter[0], 0.595041)

RESULT

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
