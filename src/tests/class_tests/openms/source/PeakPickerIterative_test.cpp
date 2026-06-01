// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerIterative.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PeakPickerIterative, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeakPickerIterative* ptr = nullptr;
PeakPickerIterative* null_ptr = nullptr;
START_SECTION(PeakPickerIterative())
{
	ptr = new PeakPickerIterative();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~PeakPickerIterative())
{
	delete ptr;
}
END_SECTION

START_SECTION((void updateMembers_()))
{
  // TODO
}
END_SECTION

START_SECTION((~PeakPickerIterative()))
{
  // TODO
}
END_SECTION

START_SECTION((template < typename PeakType > void pick(const MSSpectrum &input, MSSpectrum &output)))
{
  // TODO
}
END_SECTION

START_SECTION((template < typename PeakType > void pickExperiment(const MSExperiment< PeakType > &input, MSExperiment< PeakType > &output)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



