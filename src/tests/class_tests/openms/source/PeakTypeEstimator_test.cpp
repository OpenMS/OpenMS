// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/PeakTypeEstimator.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <iostream>
#include <vector>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(String, "$Id$")

/////////////////////////////////////////////////////////////

PeakTypeEstimator* ptr = nullptr;
PeakTypeEstimator* nullPointer = nullptr;

START_SECTION(([EXTRA]PeakTypeEstimator()))
	ptr = new PeakTypeEstimator();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(([EXTRA] ~PeakTypeEstimator()))
	delete ptr;
END_SECTION

START_SECTION((template<typename PeakConstIterator> SpectrumSettings::SpectrumType estimateType(const PeakConstIterator& begin, const PeakConstIterator& end) const))
	DTAFile file;
	MSSpectrum spec;
	// raw data (with zeros)
	file.load(OPENMS_GET_TEST_DATA_PATH("PeakTypeEstimator_raw.dta"), spec);
  TEST_EQUAL(PeakTypeEstimator::estimateType(spec.begin(), spec.end()), SpectrumSettings::PROFILE);
  // TOF raw data (without zeros)
	file.load(OPENMS_GET_TEST_DATA_PATH("PeakTypeEstimator_rawTOF.dta"), spec);
  TEST_EQUAL(PeakTypeEstimator::estimateType(spec.begin(), spec.end()), SpectrumSettings::PROFILE);
  // peak data
	file.load(OPENMS_GET_TEST_DATA_PATH("PeakTypeEstimator_peak.dta"), spec);
  TEST_EQUAL(PeakTypeEstimator::estimateType(spec.begin(), spec.end()), SpectrumSettings::CENTROID);
  // too few data points
  spec.resize(4);
	TEST_EQUAL(PeakTypeEstimator::estimateType(spec.begin(), spec.end()), SpectrumSettings::UNKNOWN);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
