// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include "OpenSwathTestHelper.h"

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgScoring.h>
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

