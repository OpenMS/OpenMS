// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger, Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/TraMLFile.h>

#include <boost/assign/std/vector.hpp>

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(TransitionPQPFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TransitionPQPFile* ptr = nullptr;
TransitionPQPFile* nullPointer = nullptr;

START_SECTION(TransitionPQPFile())
{
  ptr = new TransitionPQPFile();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~TransitionPQPFile())
{
  delete ptr;
}
END_SECTION

START_SECTION( void convertTargetedExperimentToPQP(const char * filename, OpenMS::TargetedExperiment & targeted_exp))
{
  // see TOPP tool test
  NOT_TESTABLE
}
END_SECTION

START_SECTION( void convertPQPToTargetedExperiment(const char * filename, OpenMS::TargetedExperiment & targeted_exp))
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



