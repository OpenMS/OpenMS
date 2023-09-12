// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerSH.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PeakPickerSH, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeakPickerSH* ptr = nullptr;
PeakPickerSH* null_ptr = nullptr;
START_SECTION(PeakPickerSH())
{
	ptr = new PeakPickerSH();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~PeakPickerSH())
{
	delete ptr;
}
END_SECTION

START_SECTION((virtual ~PeakPickerSH()))
{
  // TODO
}
END_SECTION

START_SECTION((template < typename PeakType > void pick(const MSSpectrum &input, MSSpectrum &output, float fWindowWidth)))
{
  // TODO
}
END_SECTION

START_SECTION((void pickExperiment(const PeakMap &input, PeakMap &output)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



