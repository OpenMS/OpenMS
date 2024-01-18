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
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>

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
  // load blosc plugin TODO: maybe move to some other part
  char *version, *date;
  auto return_code = register_blosc(&version, &date);
  TEST_EQUAL(return_code >= 0, true);
  std::cout << "Blosc version info: " << version << " " << date << std::endl;

  // open mzMLb file
  const std::string filename( OPENMS_GET_TEST_DATA_PATH("small.mzMLb") ); // file converted with pwiz
  auto mzMLb = MzMLb(filename);
  std::streamsize xml_size = mzMLb.size("mzML");
  std::cout << xml_size << std::endl; // size of XML part?
  
  // Allocate the buffer (plus one for the null terminator)
  std::string xml_buffer(xml_size, '\0');

  // Read the XML blob
  mzMLb.read(&xml_buffer[0], xml_size);
  std::cout << xml_buffer << std::endl;
  
  MSExperiment run;
  MzMLFile().loadBuffer(xml_buffer, run);
}
END_SECTION

END_TEST