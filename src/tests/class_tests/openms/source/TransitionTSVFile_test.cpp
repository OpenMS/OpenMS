// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/TraMLFile.h>

#include <boost/assign/std/vector.hpp>

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(TransitionTSVFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TransitionTSVFile* ptr = nullptr;
TransitionTSVFile* nullPointer = nullptr;

START_SECTION(TransitionTSVFile())
{
  ptr = new TransitionTSVFile();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~TransitionTSVFile())
{
  delete ptr;
}
END_SECTION

START_SECTION( void convertTargetedExperimentToTSV(const char * filename, OpenMS::TargetedExperiment & targeted_exp))
{
  // see TOPP tool test
  NOT_TESTABLE
}
END_SECTION

START_SECTION( void convertTSVToTargetedExperiment(const char * filename, OpenMS::TargetedExperiment & targeted_exp))
{
  // see TOPP tool test
  NOT_TESTABLE
}
END_SECTION

START_SECTION( void validateTargetedExperiment(OpenMS::TargetedExperiment & targeted_exp))
{
  NOT_TESTABLE
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



