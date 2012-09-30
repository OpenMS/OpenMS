// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
#include <OpenMS/KERNEL/MRMTransitionGroup.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h>

// for setup
#include "OpenSwathTestHelper.h"

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgScoring.h>
///////////////////////////

using namespace OpenMS;
using namespace std;
using namespace std;

typedef OpenSWATH_Test::RichPeakChromatogram RichPeakChromatogram;
typedef OpenSWATH_Test::MRMTransitionGroupType MRMTransitionGroupType;

START_TEST(EmgScoring, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

EmgScoring* ptr = 0;
EmgScoring* nullPointer = 0;

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

START_SECTION(double calcElutionFitScore(MRMFeature & mrmfeature, MRMTransitionGroup<SpectrumType, PeakType, TransitionT> & transition_group) )
{
  MRMTransitionGroupType transition_group;
  std::vector< RichPeakChromatogram > picked_chroms;

  OpenSWATH_Test::setup_MRMFeatureFinderScoring(transition_group, picked_chroms);

  // create the corresponding mrm feature
  int chr_idx = 0, peak_idx = 0;
  MRMFeature mrmfeature = MRMTransitionGroupPicker().createMRMFeature(transition_group, picked_chroms, chr_idx, peak_idx);

  EmgScoring emgscore;
  double elution_model_fit_score = emgscore.calcElutionFitScore(mrmfeature, transition_group);
  TEST_REAL_SIMILAR(elution_model_fit_score, 0.924365639)

}
END_SECTION

START_SECTION (void setFitterParam(Param& param))
{
  // just handing the parameters to the internal fitter
  NOT_TESTABLE
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

