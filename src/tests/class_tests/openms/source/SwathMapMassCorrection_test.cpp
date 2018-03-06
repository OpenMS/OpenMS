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

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/SwathMapMassCorrection.h>
///////////////////////////

#include <OpenMS/KERNEL/MRMTransitionGroup.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/SwathMap.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>

using namespace OpenMS;

typedef OpenSwath::LightTransition TransitionType;

OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType getData()
{
  OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType map;
  return map;
}

void addTransitions( OpenMS::MRMFeatureFinderScoring::MRMTransitionGroupType & transition_group)
{
  {
    String native_id = "tr2";
    TransitionType tr;
    tr.product_mz = 500.00;
    tr.precursor_mz = 412;
    tr.transition_name = native_id;
    transition_group.addTransition(tr, native_id );
  }

  {
    String native_id = "tr2";
    TransitionType tr;
    tr.product_mz = 600.00;
    tr.precursor_mz = 412;
    tr.transition_name = native_id;
    transition_group.addTransition(tr, native_id );
  }

  {
    String native_id = "tr3";
    TransitionType tr;
    tr.product_mz = 700.00;
    tr.precursor_mz = 412;
    tr.transition_name = native_id;
    transition_group.addTransition(tr, native_id );
  }

  {
    String native_id = "tr4";
    TransitionType tr;
    tr.product_mz = 800.00;
    tr.precursor_mz = 412;
    tr.transition_name = native_id;
    transition_group.addTransition(tr, native_id );
  }

}

START_TEST(SwathMapMassCorrection, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SwathMapMassCorrection* ptr = nullptr;
SwathMapMassCorrection* nullPointer = nullptr;

START_SECTION(SwathMapMassCorrection())
  ptr = new SwathMapMassCorrection;
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(virtual ~SwathMapMassCorrection())
    delete ptr;
END_SECTION

START_SECTION( static void correctMZ(OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType & transition_group_map, std::vector< OpenSwath::SwathMap > & swath_maps, std::string corr_type))
{

  // targets for correction are : 500.00, 600.00, 700.00, 800.00
  // "measured data" as input   : 500.02, 600.00, 699.97, 800.02

  MRMFeature feature;
  feature.setRT(3120);
  OpenMS::MRMFeatureFinderScoring::MRMTransitionGroupType transition_group;
  transition_group.addFeature(feature);
  addTransitions(transition_group);

  // Add one group to the map
  OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType transition_group_map;
  transition_group_map["group1"] = transition_group;

  // Create a mock spectrum fitting to the transition group
  boost::shared_ptr<PeakMap > exp(new PeakMap);
  {
    MSSpectrum spec;
    Peak1D p;

    p.setMZ(500.02);
    p.setIntensity(150);
    spec.push_back(p);
    p.setMZ(600.00);
    p.setIntensity(150);
    spec.push_back(p);
    p.setMZ(699.97);
    p.setIntensity(22500.01); // double the weight of all other data
    spec.push_back(p);
    p.setMZ(800.02);
    p.setIntensity(150);
    spec.push_back(p);
    spec.setRT(3121); // 3120 is the feature
    exp->addSpectrum(spec);
  }
  OpenSwath::SpectrumAccessPtr sptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp);

  OpenSwath::SwathMap map;
  map.sptr = sptr; 
  map.lower = 400;
  map.upper = 425;
  map.center = 412.5;
  map.ms1 = false;

  // should work with empty maps
  std::vector< OpenSwath::SwathMap > empty_swath_maps;
  SwathMapMassCorrection::correctMZ(transition_group_map, empty_swath_maps, "none");
  SwathMapMassCorrection::correctMZ(transition_group_map, empty_swath_maps, "unweighted_regression");

  std::vector<double> data;

  {
      std::vector< OpenSwath::SwathMap > swath_maps;
      swath_maps.push_back(map);
      SwathMapMassCorrection::correctMZ(transition_group_map, swath_maps, "none", 1.0);
      data = swath_maps[0].sptr->getSpectrumById(0)->getMZArray()->data;
      TEST_REAL_SIMILAR(data[0], 500.02)
      TEST_REAL_SIMILAR(data[1], 600.00)
      TEST_REAL_SIMILAR(data[2], 699.97)
      TEST_REAL_SIMILAR(data[3], 800.02)
  }

  {
      std::vector< OpenSwath::SwathMap > swath_maps;
      swath_maps.push_back(map);
      SwathMapMassCorrection::correctMZ(transition_group_map, swath_maps, "unweighted_regression", 0.05);
      data = swath_maps[0].sptr->getSpectrumById(0)->getMZArray()->data;
      TEST_REAL_SIMILAR(data[0], -0.00428216 + 0.999986 * 500.02) // 500.00857204075
      TEST_REAL_SIMILAR(data[1], -0.00428216 + 0.999986 * 600.00) // 599.987143224553
      TEST_REAL_SIMILAR(data[2], -0.00428216 + 0.999986 * 699.97) // 699.955714551266
      TEST_REAL_SIMILAR(data[3], -0.00428216 + 0.999986 * 800.02) // 800.004284734697
  }

  {
      std::vector< OpenSwath::SwathMap > swath_maps;
      swath_maps.push_back(map);
      SwathMapMassCorrection::correctMZ(transition_group_map, swath_maps, "unweighted_regression", 1.0);
      data = swath_maps[0].sptr->getSpectrumById(0)->getMZArray()->data;
      TEST_REAL_SIMILAR(data[0], -0.0219795 + 1.00003 * 500.02) // 500.01300527988
      TEST_REAL_SIMILAR(data[1], -0.0219795 + 1.00003 * 600.00) // 599.99600151022
      TEST_REAL_SIMILAR(data[2], -0.0219795 + 1.00003 * 699.97) // 699.96899744088
      TEST_REAL_SIMILAR(data[3], -0.0219795 + 1.00003 * 800.02) // 800.02199576900
  }

  {
      std::vector< OpenSwath::SwathMap > swath_maps;
      swath_maps.push_back(map);
      SwathMapMassCorrection::correctMZ(transition_group_map, swath_maps, "weighted_regression", 1.0);
      data = swath_maps[0].sptr->getSpectrumById(0)->getMZArray()->data;
      TEST_REAL_SIMILAR(data[0], -0.0315101 + 1.00005 * 500.02) // 500.01539273402
      TEST_REAL_SIMILAR(data[1], -0.0315101 + 1.00005 * 600.00) // 600.00077200650
      TEST_REAL_SIMILAR(data[2], -0.0315101 + 1.00005 * 699.97) // 699.97615074094
      TEST_REAL_SIMILAR(data[3], -0.0315101 + 1.00005 * 800.02) // 800.03153377967
  }

  {
    std::vector< OpenSwath::SwathMap > swath_maps;
    swath_maps.push_back(map);
    SwathMapMassCorrection::correctMZ(transition_group_map, swath_maps, "quadratic_regression", 1.0);
    data = swath_maps[0].sptr->getSpectrumById(0)->getMZArray()->data;
    TEST_REAL_SIMILAR(data[0], -0.7395987927448004 + 1.002305255194642 * 500.02 -1.750157412772069e-06 * 500.02 * 500.02) // 499.995500552639
    TEST_REAL_SIMILAR(data[1], -0.7395987927448004 + 1.002305255194642 * 600.00 -1.750157412772069e-06 * 600.00 * 600.00) // 600.013497655443
    TEST_REAL_SIMILAR(data[2], -0.7395987927448004 + 1.002305255194642 * 699.97 -1.750157412772069e-06 * 699.97 * 699.97) // 699.986507058627
    TEST_REAL_SIMILAR(data[3], -0.7395987927448004 + 1.002305255194642 * 800.02 -1.750157412772069e-06 * 800.02 * 800.02) // 800.004494718161
  }

  {
    std::vector< OpenSwath::SwathMap > swath_maps;
    swath_maps.push_back(map);
    SwathMapMassCorrection::correctMZ(transition_group_map, swath_maps, "weighted_quadratic_regression", 1.0);
    data = swath_maps[0].sptr->getSpectrumById(0)->getMZArray()->data;
    TEST_REAL_SIMILAR(data[0], -0.8323316718451679 + 1.002596944948891 * 500.02 -1.967834556637627e-06 * 500.02 * 500.02) // 499.994194744862
    TEST_REAL_SIMILAR(data[1], -0.8323316718451679 + 1.002596944948891 * 600.00 -1.967834556637627e-06 * 600.00 * 600.00) // 600.0174148571
    TEST_REAL_SIMILAR(data[2], -0.8323316718451679 + 1.002596944948891 * 699.97 -1.967834556637627e-06 * 699.97 * 699.97) // 699.991295598558
    TEST_REAL_SIMILAR(data[3], -0.8323316718451679 + 1.002596944948891 * 800.02 -1.967834556637627e-06 * 800.02 * 800.02) // 800.005799138426
  }

  {
    std::vector< OpenSwath::SwathMap > swath_maps;
    swath_maps.push_back(map);
    SwathMapMassCorrection::correctMZ(transition_group_map, swath_maps, "quadratic_regression_delta_ppm", 1.0);
    data = swath_maps[0].sptr->getSpectrumById(0)->getMZArray()->data;
    TEST_REAL_SIMILAR(data[0], 499.997160932778)
    TEST_REAL_SIMILAR(data[1], 600.010219722383)
    TEST_REAL_SIMILAR(data[2], 699.988081672119)
    TEST_REAL_SIMILAR(data[3], 800.004537672719)
  }

  {
    std::vector< OpenSwath::SwathMap > swath_maps;
    swath_maps.push_back(map);
    SwathMapMassCorrection::correctMZ(transition_group_map, swath_maps, "weighted_quadratic_regression_delta_ppm", 1.0);
    data = swath_maps[0].sptr->getSpectrumById(0)->getMZArray()->data;
    TEST_REAL_SIMILAR(data[0], 499.996336995751)
    TEST_REAL_SIMILAR(data[1], 600.013185628794)
    TEST_REAL_SIMILAR(data[2], 699.992311403648)
    TEST_REAL_SIMILAR(data[3], 800.005854568825)
  }

}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
