// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FORMAT/PeakTypeEstimator.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <iostream>
#include <vector>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(String, "$Id: PeakTypeEstimator_test.C 6139 2009-10-20 07:37:40Z andreas_bertsch $")

/////////////////////////////////////////////////////////////

PeakTypeEstimator* ptr = 0;

START_SECTION(([EXTRA]PeakTypeEstimator()))
	ptr = new PeakTypeEstimator();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION(([EXTRA] ~PeakTypeEstimator()))
	delete ptr;
END_SECTION

START_SECTION((template<typename PeakConstIterator> SpectrumSettings::SpectrumType estimateType(const PeakConstIterator& begin, const PeakConstIterator& end) const))
	DTAFile file;
	MSExperiment<> exp;
	exp.resize(4);
	PeakTypeEstimator pte;
	// raw data (with zeros)
	file.load(OPENMS_GET_TEST_DATA_PATH("PeakTypeEstimator_raw.dta"),exp[0]);
	// TOF raw data (without zeros)
	file.load(OPENMS_GET_TEST_DATA_PATH("PeakTypeEstimator_rawTOF.dta"),exp[1]);
	//peak data
	file.load(OPENMS_GET_TEST_DATA_PATH("PeakTypeEstimator_peak.dta"),exp[2]);
	//too few data points
	exp[3].resize(4);
	
	TEST_EQUAL(pte.estimateType(exp[0].begin(),exp[0].end()), SpectrumSettings::RAWDATA);
	TEST_EQUAL(pte.estimateType(exp[1].begin(),exp[1].end()), SpectrumSettings::RAWDATA);
	TEST_EQUAL(pte.estimateType(exp[2].begin(),exp[2].end()), SpectrumSettings::PEAKS);
	TEST_EQUAL(pte.estimateType(exp[3].begin(),exp[3].end()), SpectrumSettings::UNKNOWN);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
