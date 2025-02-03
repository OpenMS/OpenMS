// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Witold Wolski $
// --------------------------------------------------------------------------

#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionHelper.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/Transitions.h>

#include <OpenMS/CONCEPT/ClassTest.h>
using namespace OpenMS;

using namespace std;

///////////////////////////

START_TEST(ITrans2Trans, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(initializeXCorrMatrix)
{
  OpenSwath::LightTargetedExperiment  lte;
  OpenSwath::TargetedExperiment  te;
  //OpenSwath::convert(lte, te);
  //TODO (wolski): write tests
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
