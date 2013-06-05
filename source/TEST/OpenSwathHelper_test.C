// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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

///////////////////////////

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>
#include <boost/assign/std/vector.hpp>
#include <boost/assign/list_of.hpp>

///////////////////////////

START_TEST(OpenSwathHelper, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace std;
using namespace OpenMS;
using namespace OpenSwath;

OpenSwathHelper* ptr = 0;
OpenSwathHelper* nullPointer = 0;

START_SECTION(OpenSwathHelper())
{
  ptr = new OpenSwathHelper();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~OpenSwathHelper())
{
  delete ptr;
}
END_SECTION

START_SECTION(static void selectSwathTransitions(const OpenMS::TargetedExperiment &targeted_exp, OpenMS::TargetedExperiment &transition_exp_used, double min_upper_edge_dist, double lower, double upper))
{
  TargetedExperiment exp1;
  TargetedExperiment exp2;

  ReactionMonitoringTransition tr1;
  ReactionMonitoringTransition tr2;
  ReactionMonitoringTransition tr3;

  tr1.setPrecursorMZ(100.0);
  tr2.setPrecursorMZ(200.0);
  tr3.setPrecursorMZ(300.0);

  std::vector<ReactionMonitoringTransition> transitions;
  transitions.push_back(tr1);
  transitions.push_back(tr2);
  transitions.push_back(tr3);

  exp1.setTransitions(transitions);

  // select all transitions between 200 and 500
  OpenSwathHelper::selectSwathTransitions(exp1, exp2, 1.0, 199.9, 500);
  TEST_EQUAL(exp2.getTransitions().size(), 2)
}
END_SECTION

START_SECTION(static void selectSwathTransitions(const OpenSwath::LightTargetedExperiment &targeted_exp, OpenSwath::LightTargetedExperiment &transition_exp_used, double min_upper_edge_dist, double lower, double upper))
{
  LightTargetedExperiment exp1;
  LightTargetedExperiment exp2;

  LightTransition tr1;
  LightTransition tr2;
  LightTransition tr3;

  tr1.precursor_mz = 100.0;
  tr2.precursor_mz = 200.0;
  tr3.precursor_mz = 300.0;

  std::vector<LightTransition> transitions;
  transitions.push_back(tr1);
  transitions.push_back(tr2);
  transitions.push_back(tr3);

  exp1.transitions = transitions;

  // select all transitions between 200 and 500
  OpenSwathHelper::selectSwathTransitions(exp1, exp2, 1.0, 199.9, 500);
  TEST_EQUAL(exp2.getTransitions().size(), 2)
}
END_SECTION

START_SECTION( (template < class TargetedExperimentT > static bool checkSwathMapAndSelectTransitions(const OpenMS::MSExperiment< Peak1D > &exp, const TargetedExperimentT &targeted_exp, TargetedExperimentT &transition_exp_used, double min_upper_edge_dist)))
{
  // tested above already
  NOT_TESTABLE
}
END_SECTION

START_SECTION(static void checkSwathMap(const OpenMS::MSExperiment< Peak1D > &swath_map, double &lower, double &upper))
{
  OpenMS::MSExperiment<Peak1D> swath_map;
  OpenMS::MSSpectrum<Peak1D> spectrum;
  OpenMS::Precursor prec;
  std::vector<Precursor> precursors;
  prec.setIsolationWindowLowerOffset(200);
  prec.setIsolationWindowUpperOffset(300);
  precursors.push_back(prec);
  spectrum.setPrecursors(precursors);
  swath_map.addSpectrum(spectrum);

  double lower, upper;
  OpenSwathHelper::checkSwathMap(swath_map, lower, upper);

  TEST_REAL_SIMILAR(lower, 200);
  TEST_REAL_SIMILAR(upper, 300);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


