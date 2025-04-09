// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/FORMAT/DATAACCESS/NoopMSDataConsumer.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(NoopMSDataConsumer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

NoopMSDataConsumer* ptr = nullptr;
NoopMSDataConsumer* null_ptr = nullptr;
START_SECTION(NoopMSDataConsumer())
{
	ptr = new NoopMSDataConsumer();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~NoopMSDataConsumer())
{
	delete ptr;
}
END_SECTION

START_SECTION((void setExperimentalSettings(const ExperimentalSettings &)))
{
  // TODO
}
END_SECTION

START_SECTION((void setExpectedSize(Size, Size)))
{
  // TODO
}
END_SECTION

START_SECTION((void consumeSpectrum(SpectrumType &)))
{
  // TODO
}
END_SECTION

START_SECTION((void consumeChromatogram(ChromatogramType &)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



