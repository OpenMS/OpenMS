// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

///////////////////////////

#include <OpenMS/CONCEPT/Types.h>

#include "H5Cpp.h"
#include "blosc_filter.h"

using namespace std;

#include <OpenMS/FORMAT/Connection_mzMLb.h>

//using MzMLb = pwiz::msdata::mzmlb::Connection_mzMLb;

#include <string>
#include <iostream>

START_TEST(MzMLb, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

using MzMLb = OpenMS::Connection_mzMLb; // TODO: maybe rename? MzMLbStream?

START_SECTION((MzMLb()))
{
  // open mzMLb file
  const std::string filename( OPENMS_GET_TEST_DATA_PATH("small.mzMLb") );
  auto mzMLb = MzMLb(filename);

  // Read from the stream and print to std::cout
  char* buffer;
  mzMLb.read(buffer, (std::streamsize)1e3);
  std::cout << buffer << std::endl;
}
END_SECTION

END_TEST