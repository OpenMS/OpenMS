// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include "OpenSwathTestHelper.h"

#include <OpenMS/FEATUREFINDER/EmgScoring.h>
#include <OpenMS/KERNEL/MRMTransitionGroup.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>

///////////////////////////

using namespace OpenMS;

START_TEST(EmgScoring, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

EmgScoring* ptr = nullptr;
EmgScoring* nullPointer = nullptr;

START_SECTION(EmgScoring())
{
  ptr = new EmgScoring();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~EmgScoring())
{
  delete ptr;
}
END_SECTION

START_SECTION(Param getDefaults())
{
  EmgScoring emgscore;
  Param p = emgscore.getDefaults();
  TEST_NOT_EQUAL(&p, nullPointer)
}
END_SECTION

START_SECTION(void setFitterParam(Param param))
{
  EmgScoring emgscore;
  Param p = emgscore.getDefaults();
  TEST_NOT_EQUAL(&p, nullPointer)
  emgscore.setFitterParam(p);
}
END_SECTION

START_SECTION(( template < typename SpectrumType, class TransitionT > double calcElutionFitScore(MRMFeature &mrmfeature, MRMTransitionGroup< SpectrumType, TransitionT > &transition_group)))
{
  // test a set of feature (belonging to the same peptide)
  double elution_model_fit_score;
  EmgScoring emgscore;

  MRMFeature feature = OpenSWATH_Test::createMockFeature();
  OpenSWATH_Test::MRMTransitionGroupType transition_group = OpenSWATH_Test::createMockTransitionGroup();

  elution_model_fit_score = emgscore.calcElutionFitScore(feature, transition_group);
  TEST_REAL_SIMILAR(elution_model_fit_score, 0.924365639)

}
END_SECTION

START_SECTION( double elutionModelFit(ConvexHull2D::PointArrayType current_section, bool smooth_data) )
{
  // test a single feature
  double elution_model_fit_score;
  EmgScoring emgscore;

  MRMFeature feature = OpenSWATH_Test::createMockFeature();
  Feature f = feature.getFeature("tr1");

  elution_model_fit_score = emgscore.elutionModelFit(f.getConvexHulls()[0].getHullPoints() , false);
  TEST_REAL_SIMILAR(elution_model_fit_score, 0.981013417243958) 
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

