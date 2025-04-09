// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MSDataWritingConsumer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MSDataWritingConsumer* ptr = 0;
MSDataWritingConsumer* null_ptr = 0;
START_SECTION(MSDataWritingConsumer())
{
	ptr = new MSDataWritingConsumer();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~MSDataWritingConsumer())
{
	delete ptr;
}
END_SECTION

START_SECTION((virtual void setExperimentalSettings(const ExperimentalSettings &exp)))
{
  // TODO
}
END_SECTION

START_SECTION((virtual void setExpectedSize(Size expectedSpectra, Size expectedChromatograms)))
{
  // TODO
}
END_SECTION

START_SECTION((virtual void consumeSpectrum(SpectrumType &s)))
{
  // TODO
}
END_SECTION

START_SECTION((virtual void consumeChromatogram(ChromatogramType &c)))
{
  // TODO
}
END_SECTION

START_SECTION((MSDataWritingConsumer(String filename)))
{
  // TODO
}
END_SECTION

START_SECTION((virtual ~MSDataWritingConsumer()))
{
  // TODO
}
END_SECTION

START_SECTION((virtual void addDataProcessing(DataProcessing d)))
{
  // TODO
}
END_SECTION

START_SECTION((virtual Size getNrSpectraWritten()))
{
  // TODO
}
END_SECTION

START_SECTION((virtual Size getNrChromatogramsWritten()))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



