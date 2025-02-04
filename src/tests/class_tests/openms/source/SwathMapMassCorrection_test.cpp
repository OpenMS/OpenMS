// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>
#include <OpenMS/IONMOBILITY/IMDataConverter.h>

using namespace OpenMS;

typedef OpenSwath::LightTransition TransitionType;

typedef std::map<String, OpenMS::MRMFeatureFinderScoring::MRMTransitionGroupType *> TransitionGroupMapPtrType;

OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType getData()
{
  OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType map;
  return map;
}

OpenSwath::LightTargetedExperiment addTransitions( OpenMS::MRMFeatureFinderScoring::MRMTransitionGroupType & transition_group)
{
  OpenSwath::LightTargetedExperiment exp;
  {
    String native_id = "tr1";
    TransitionType tr;
    tr.product_mz = 500.00;
    tr.precursor_mz = 412;
    tr.peptide_ref = "pep1";
    tr.precursor_im = 11;
    tr.transition_name = native_id;
    transition_group.addTransition(tr, native_id );
    exp.transitions.push_back(tr);
  }

  {
    String native_id = "tr2";
    TransitionType tr;
    tr.product_mz = 600.00;
    tr.precursor_mz = 412;
    tr.peptide_ref = "pep1";
    tr.precursor_im = 11;
    tr.transition_name = native_id;
    transition_group.addTransition(tr, native_id );
    exp.transitions.push_back(tr);
  }

  {
    String native_id = "tr3";
    TransitionType tr;
    tr.product_mz = 700.00;
    tr.precursor_mz = 412;
    tr.peptide_ref = "pep1";
    tr.precursor_im = 11;
    tr.transition_name = native_id;
    transition_group.addTransition(tr, native_id );
    exp.transitions.push_back(tr);
  }

  {
    String native_id = "tr4";
    TransitionType tr;
    tr.product_mz = 800.00;
    tr.precursor_mz = 412;
    tr.peptide_ref = "pep1";
    tr.precursor_im = 11;
    tr.transition_name = native_id;
    transition_group.addTransition(tr, native_id );
    exp.transitions.push_back(tr);
  }

  OpenSwath::LightCompound cmp;
  cmp.id = "pep1";
  cmp.drift_time = 11;
  exp.compounds.push_back(cmp);
  return exp;
}

void addTransitionsPep2( OpenMS::MRMFeatureFinderScoring::MRMTransitionGroupType & transition_group, OpenSwath::LightTargetedExperiment& exp)
{
  String native_id = "tr5";
  TransitionType tr;
  tr.product_mz = 900.00;
  tr.precursor_mz = 500.0;
  tr.precursor_im = 15;
  tr.peptide_ref = "pep2";
  tr.transition_name = native_id;
  transition_group.addTransition(tr, native_id );
  exp.transitions.push_back(tr);

  OpenSwath::LightCompound cmp;
  cmp.id = "pep2";
  cmp.drift_time = 15;
  exp.compounds.push_back(cmp);
}
void addTransitionsPep3( OpenMS::MRMFeatureFinderScoring::MRMTransitionGroupType & transition_group, OpenSwath::LightTargetedExperiment& exp)
{
  String native_id = "tr6";
  TransitionType tr;
  tr.product_mz = 950.00;
  tr.precursor_mz = 600.0;
  tr.precursor_im = 20;
  tr.peptide_ref = "pep3";
  tr.transition_name = native_id;
  transition_group.addTransition(tr, native_id );
  exp.transitions.push_back(tr);

  OpenSwath::LightCompound cmp;
  cmp.id = "pep3";
  cmp.drift_time = 20;
  exp.compounds.push_back(cmp);
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

START_SECTION( void correctMZ(OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType & transition_group_map, std::vector< OpenSwath::SwathMap > & swath_maps, const std::string& corr_type, const bool pasef))
{
  // None of the tests here test pasef flag
  bool pasef = false;

  // targets for correction are : 500.00, 600.00, 700.00, 800.00
  // "measured data" as input   : 500.02, 600.00, 699.97, 800.02
  // IM of feature is 11

  MRMFeature feature;
  feature.setRT(3120);
  OpenMS::MRMFeatureFinderScoring::MRMTransitionGroupType transition_group;
  transition_group.addFeature(feature);
  OpenSwath::LightTargetedExperiment targ_exp = addTransitions(transition_group);

  // Add one group to the map
  std::map<String, OpenMS::MRMFeatureFinderScoring::MRMTransitionGroupType *> transition_group_map;
  transition_group_map["group1"] = &transition_group;
  transition_group_map["group2"] = &transition_group;
  transition_group_map["group3"] = &transition_group;

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

  // Create secondary mock spectrum for testing PASEF flag, this spectrum should never be used

  boost::shared_ptr<PeakMap > exp2(new PeakMap);
  {
    MSSpectrum spec;
    Peak1D p;

    p.setMZ(499.97);
    p.setIntensity(1500000); // This has the highest weight in this spectrum
    spec.push_back(p);
    p.setMZ(600.00);
    p.setIntensity(150);
    spec.push_back(p);
    p.setMZ(699.97);
    p.setIntensity(150);
    spec.push_back(p);
    p.setMZ(800.02);
    p.setIntensity(150);
    spec.push_back(p);
    spec.setRT(3122); // 3120 is the feature
    exp2->addSpectrum(spec);
  }

  OpenSwath::SpectrumAccessPtr sptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp);
  OpenSwath::SpectrumAccessPtr sptr2 = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp2);

  OpenSwath::SwathMap map;
  map.sptr = sptr;
  map.lower = 400;
  map.upper = 425;
  map.center = 412.5;
  map.imLower = 10; // only used for pasef results
  map.imUpper = 20; // only used for pasef results
  map.ms1 = false;

  OpenSwath::SwathMap mapPASEF;
  mapPASEF.sptr = sptr2;
  mapPASEF.lower = 400;
  mapPASEF.center = 412.5;
  mapPASEF.ms1 = false;
  mapPASEF.imLower=20;
  mapPASEF.imUpper=30;

  SwathMapMassCorrection mc;

  // should work with empty maps
  std::vector< OpenSwath::SwathMap > empty_swath_maps;
  mc.correctMZ(transition_group_map, targ_exp, empty_swath_maps, pasef);

  auto p = mc.getDefaults();
  p.setValue("mz_correction_function", "unweighted_regression");
  mc.setParameters(p);
  mc.correctMZ(transition_group_map, targ_exp, empty_swath_maps, pasef);

  std::vector<double> data;

  {
      auto p = mc.getDefaults();
      p.setValue("mz_correction_function", "none");
      p.setValue("mz_extraction_window", 1.0);
      mc.setParameters(p);

      std::vector< OpenSwath::SwathMap > swath_maps;
      swath_maps.push_back(map);
      mc.correctMZ(transition_group_map, targ_exp, swath_maps, pasef);
      data = swath_maps[0].sptr->getSpectrumById(0)->getMZArray()->data;
      TEST_REAL_SIMILAR(data[0], 500.02)
      TEST_REAL_SIMILAR(data[1], 600.00)
      TEST_REAL_SIMILAR(data[2], 699.97)
      TEST_REAL_SIMILAR(data[3], 800.02)
  }

  {
      auto p = mc.getDefaults();
      p.setValue("mz_correction_function", "unweighted_regression");
      p.setValue("mz_extraction_window", 0.05);
      mc.setParameters(p);

      std::vector< OpenSwath::SwathMap > swath_maps;
      swath_maps.push_back(map);
      mc.correctMZ(transition_group_map, targ_exp, swath_maps, pasef);
      data = swath_maps[0].sptr->getSpectrumById(0)->getMZArray()->data;
      TEST_REAL_SIMILAR(data[0], -0.00428216 + 0.999986 * 500.02) // 500.00857204075
      TEST_REAL_SIMILAR(data[1], -0.00428216 + 0.999986 * 600.00) // 599.987143224553
      TEST_REAL_SIMILAR(data[2], -0.00428216 + 0.999986 * 699.97) // 699.955714551266
      TEST_REAL_SIMILAR(data[3], -0.00428216 + 0.999986 * 800.02) // 800.004284734697
  }

  {
      auto p = mc.getDefaults();
      p.setValue("mz_correction_function", "unweighted_regression");
      p.setValue("mz_extraction_window", 1.0);
      mc.setParameters(p);

      std::vector< OpenSwath::SwathMap > swath_maps;
      swath_maps.push_back(map);
      mc.correctMZ(transition_group_map, targ_exp, swath_maps, pasef);
      data = swath_maps[0].sptr->getSpectrumById(0)->getMZArray()->data;
      TEST_REAL_SIMILAR(data[0], -0.0219795 + 1.00003 * 500.02) // 500.01300527988
      TEST_REAL_SIMILAR(data[1], -0.0219795 + 1.00003 * 600.00) // 599.99600151022
      TEST_REAL_SIMILAR(data[2], -0.0219795 + 1.00003 * 699.97) // 699.96899744088
      TEST_REAL_SIMILAR(data[3], -0.0219795 + 1.00003 * 800.02) // 800.02199576900
  }

  {
    auto p = mc.getDefaults();
    p.setValue("mz_correction_function", "unweighted_regression");
    p.setValue("mz_extraction_window", 1.0);
    mc.setParameters(p);

    std::vector< OpenSwath::SwathMap > swath_maps;
    swath_maps.push_back(map);
    mc.correctMZ(transition_group_map, targ_exp, swath_maps, pasef);
    data = swath_maps[0].sptr->getSpectrumById(0)->getMZArray()->data;
    TEST_REAL_SIMILAR(data[0], -0.0315101 + 1.00005 * 500.02) // 500.01539273402
    TEST_REAL_SIMILAR(data[1], -0.0315101 + 1.00005 * 600.00) // 600.00077200650
    TEST_REAL_SIMILAR(data[2], -0.0315101 + 1.00005 * 699.97) // 699.97615074094
    TEST_REAL_SIMILAR(data[3], -0.0315101 + 1.00005 * 800.02) // 800.03153377967
  }

  {
    auto p = mc.getDefaults();
    p.setValue("mz_correction_function", "quadratic_regression");
    p.setValue("mz_extraction_window", 1.0);
    mc.setParameters(p);

    std::vector< OpenSwath::SwathMap > swath_maps;
    swath_maps.push_back(map);
    mc.correctMZ(transition_group_map, targ_exp, swath_maps, pasef);
    data = swath_maps[0].sptr->getSpectrumById(0)->getMZArray()->data;
    TEST_REAL_SIMILAR(data[0], -0.7395987927448004 + 1.002305255194642 * 500.02 -1.750157412772069e-06 * 500.02 * 500.02) // 499.995500552639
    TEST_REAL_SIMILAR(data[1], -0.7395987927448004 + 1.002305255194642 * 600.00 -1.750157412772069e-06 * 600.00 * 600.00) // 600.013497655443
    TEST_REAL_SIMILAR(data[2], -0.7395987927448004 + 1.002305255194642 * 699.97 -1.750157412772069e-06 * 699.97 * 699.97) // 699.986507058627
    TEST_REAL_SIMILAR(data[3], -0.7395987927448004 + 1.002305255194642 * 800.02 -1.750157412772069e-06 * 800.02 * 800.02) // 800.004494718161
  }

  {
    auto p = mc.getDefaults();
    p.setValue("mz_correction_function", "weighted_quadratic_regression");
    p.setValue("mz_extraction_window", 1.0);
    mc.setParameters(p);

    std::vector< OpenSwath::SwathMap > swath_maps;
    swath_maps.push_back(map);
    mc.correctMZ(transition_group_map, targ_exp, swath_maps, pasef);
    data = swath_maps[0].sptr->getSpectrumById(0)->getMZArray()->data;
    TEST_REAL_SIMILAR(data[0], -0.8323316718451679 + 1.002596944948891 * 500.02 -1.967834556637627e-06 * 500.02 * 500.02) // 499.994194744862
    TEST_REAL_SIMILAR(data[1], -0.8323316718451679 + 1.002596944948891 * 600.00 -1.967834556637627e-06 * 600.00 * 600.00) // 600.0174148571
    TEST_REAL_SIMILAR(data[2], -0.8323316718451679 + 1.002596944948891 * 699.97 -1.967834556637627e-06 * 699.97 * 699.97) // 699.991295598558
    TEST_REAL_SIMILAR(data[3], -0.8323316718451679 + 1.002596944948891 * 800.02 -1.967834556637627e-06 * 800.02 * 800.02) // 800.005799138426
  }

  {
    auto p = mc.getDefaults();
    p.setValue("mz_correction_function", "quadratic_regression_delta_ppm");
    p.setValue("mz_extraction_window", 1.0);
    mc.setParameters(p);

    std::vector< OpenSwath::SwathMap > swath_maps;
    swath_maps.push_back(map);
    mc.correctMZ(transition_group_map, targ_exp, swath_maps, pasef);
    data = swath_maps[0].sptr->getSpectrumById(0)->getMZArray()->data;
    TEST_REAL_SIMILAR(data[0], 499.997160932778)
    TEST_REAL_SIMILAR(data[1], 600.010219722383)
    TEST_REAL_SIMILAR(data[2], 699.988081672119)
    TEST_REAL_SIMILAR(data[3], 800.004537672719)
  }

  {
    auto p = mc.getDefaults();
    p.setValue("mz_correction_function", "weighted_quadratic_regression_delta_ppm");
    p.setValue("mz_extraction_window", 1.0);
    mc.setParameters(p);

    std::vector< OpenSwath::SwathMap > swath_maps;
    swath_maps.push_back(map);
    mc.correctMZ(transition_group_map, targ_exp, swath_maps, pasef);
    data = swath_maps[0].sptr->getSpectrumById(0)->getMZArray()->data;
    TEST_REAL_SIMILAR(data[0], 499.996336995751)
    TEST_REAL_SIMILAR(data[1], 600.013185628794)
    TEST_REAL_SIMILAR(data[2], 699.992311403648)
    TEST_REAL_SIMILAR(data[3], 800.005854568825)
  }

  {
      auto p = mc.getDefaults();
      p.setValue("mz_correction_function", "unweighted_regression");
      p.setValue("mz_extraction_window", 0.05);
      mc.setParameters(p);

      std::vector< OpenSwath::SwathMap > swath_maps;
      swath_maps.push_back(map);
      mc.correctMZ(transition_group_map, targ_exp, swath_maps, pasef);
      data = swath_maps[0].sptr->getSpectrumById(0)->getMZArray()->data;
      TEST_REAL_SIMILAR(data[0], -0.00428216 + 0.999986 * 500.02) // 500.00857204075
      TEST_REAL_SIMILAR(data[1], -0.00428216 + 0.999986 * 600.00) // 599.987143224553
      TEST_REAL_SIMILAR(data[2], -0.00428216 + 0.999986 * 699.97) // 699.955714551266
      TEST_REAL_SIMILAR(data[3], -0.00428216 + 0.999986 * 800.02) // 800.004284734697
  }

  {
      // Test with PASEF flag on, correction should only occur based on the first spectrum and thus should acheive the same results as above
      auto p = mc.getDefaults();
      p.setValue("mz_correction_function", "unweighted_regression");
      p.setValue("mz_extraction_window", 0.05);
      mc.setParameters(p);
      bool pasef = true;

      std::vector< OpenSwath::SwathMap > swath_maps_pasef;
      swath_maps_pasef.push_back(map);
      swath_maps_pasef.push_back(mapPASEF);

      mc.correctMZ(transition_group_map, targ_exp, swath_maps_pasef, pasef);
      data = swath_maps_pasef[0].sptr->getSpectrumById(0)->getMZArray()->data;
      std::cout << data[0] << ";" << data[1] << ";" <<data[2] << ";" << data[3] << std::endl;
      TEST_REAL_SIMILAR(data[0], -0.00428216 + 0.999986 * 500.02) // 500.00857204075
      TEST_REAL_SIMILAR(data[1], -0.00428216 + 0.999986 * 600.00) // 599.987143224553
      TEST_REAL_SIMILAR(data[2], -0.00428216 + 0.999986 * 699.97) // 699.955714551266
      TEST_REAL_SIMILAR(data[3], -0.00428216 + 0.999986 * 800.02) // 800.004284734697
  }

}
END_SECTION

START_SECTION( void correctIM(const std::map<String, OpenMS::MRMFeatureFinderScoring::MRMTransitionGroupType *> & transition_group_map, const std::vector< OpenSwath::SwathMap > & swath_maps, const bool pasef, TransformationDescription& im_trafo))
{

  // m/z targets for correction are : 500.00, 600.00, 700.00, 800.00, 900.00, 950.00
  // mobility targets               :  11.00,  11.00,  11.00,  11.00,  15.00,  20.00
  // "measured data" as input       :  22.00,  21.50,  20.50,  21.00,  24.00,  31.00

  MRMFeature feature;
  feature.setRT(3120);
  OpenMS::MRMFeatureFinderScoring::MRMTransitionGroupType gr1, gr2, gr3;
  gr1.addFeature(feature);
  gr2.addFeature(feature);
  gr3.addFeature(feature);
  OpenSwath::LightTargetedExperiment targ_exp = addTransitions(gr1);
  addTransitionsPep2(gr2, targ_exp);
  addTransitionsPep3(gr3, targ_exp);

  bool pasef = false;

  // Add one group to the map
  std::map<String, OpenMS::MRMFeatureFinderScoring::MRMTransitionGroupType *> transition_group_map;
  transition_group_map["group1"] = &gr1;
  transition_group_map["group2"] = &gr2;
  transition_group_map["group3"] = &gr3;

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
    p.setMZ(900.02);
    p.setIntensity(150);
    spec.push_back(p);
    p.setMZ(950.02);
    p.setIntensity(150);
    spec.push_back(p);
    spec.setRT(3121); // 3120 is the feature
    DataArrays::FloatDataArray ion_mobility;
    ion_mobility.push_back(22.0);
    ion_mobility.push_back(21.5);
    ion_mobility.push_back(20.5);
    ion_mobility.push_back(21.0);
    ion_mobility.push_back(24.0);
    ion_mobility.push_back(31.0);
    IMDataConverter::setIMUnit(ion_mobility, DriftTimeUnit::MILLISECOND);
    ion_mobility.setName("Ion Mobility");
    auto& fda = spec.getFloatDataArrays();
    fda.push_back(ion_mobility);

    spec.setFloatDataArrays(fda);
    exp->addSpectrum(spec);
  }

  // Create a mock pasef spectrum, should not be used
  boost::shared_ptr<PeakMap > exp2(new PeakMap);
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
    p.setMZ(900.02);
    p.setIntensity(150);
    spec.push_back(p);
    p.setMZ(950.02);
    p.setIntensity(150);
    spec.push_back(p);
    spec.setRT(3121); // 3120 is the feature
    DataArrays::FloatDataArray ion_mobility;
    ion_mobility.push_back(22.3);
    ion_mobility.push_back(21.5);
    ion_mobility.push_back(20.5);
    ion_mobility.push_back(21.5);
    ion_mobility.push_back(24.4);
    ion_mobility.push_back(31.0);
    IMDataConverter::setIMUnit(ion_mobility, DriftTimeUnit::MILLISECOND);
    ion_mobility.setName("Ion Mobility");
    auto& fda = spec.getFloatDataArrays();
    fda.push_back(ion_mobility);

    spec.setFloatDataArrays(fda);
    exp->addSpectrum(spec);
  }

  boost::shared_ptr<PeakMap > exp_ms1(new PeakMap);
  {
    MSSpectrum spec;
    Peak1D p;

    p.setMZ(412.02);
    p.setIntensity(150);
    spec.push_back(p);
    p.setMZ(500.02);
    p.setIntensity(150);
    spec.push_back(p);
    p.setMZ(600.01);
    p.setIntensity(150.0);
    spec.push_back(p);
    spec.setRT(3121); // 3120 is the feature
    DataArrays::FloatDataArray ion_mobility;
    ion_mobility.push_back(22.0);
    ion_mobility.push_back(24.0);
    ion_mobility.push_back(31.0);
    IMDataConverter::setIMUnit(ion_mobility, DriftTimeUnit::MILLISECOND);
    ion_mobility.setName("Ion Mobility");
    auto& fda = spec.getFloatDataArrays();
    fda.push_back(ion_mobility);
    spec.setFloatDataArrays(fda);
    exp_ms1->addSpectrum(spec);
  }
  OpenSwath::SpectrumAccessPtr sptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp);
  OpenSwath::SpectrumAccessPtr sptr2 = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp2);
  OpenSwath::SpectrumAccessPtr sptr_ms1 = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp_ms1);

  OpenSwath::SwathMap map;
  map.sptr = sptr;
  map.lower = 400;
  map.upper = 800;
  map.center = 412.5;
  map.imLower = 10;
  map.imUpper = 30;
  map.ms1 = false;

  OpenSwath::SwathMap mapPASEF;
  mapPASEF.sptr = sptr2;
  mapPASEF.lower = 400;
  mapPASEF.upper = 800;
  mapPASEF.center = 412.5;
  mapPASEF.ms1 = false;
  mapPASEF.imLower=100; // although this does not make complete sense with spectrum (since spectrum has values in 20s) this is ok because do not want to use this window
  mapPASEF.imUpper=200;

  OpenSwath::SwathMap ms1_map;
  ms1_map.sptr = sptr_ms1;
  ms1_map.ms1 = true;

  SwathMapMassCorrection mc;

  // should work with empty maps
  std::vector< OpenSwath::SwathMap > empty_swath_maps;
  TransformationDescription im_trafo;
  mc.correctIM(transition_group_map, targ_exp, empty_swath_maps, pasef, im_trafo);
  TEST_REAL_SIMILAR(im_trafo.apply(10), 10)
  TEST_REAL_SIMILAR(im_trafo.apply(100), 100)

  auto p = mc.getDefaults();
  p.setValue("mz_correction_function", "unweighted_regression");
  mc.setParameters(p);
  mc.correctIM(transition_group_map, targ_exp, empty_swath_maps, pasef, im_trafo);
  TEST_REAL_SIMILAR(im_trafo.apply(10), 10)
  TEST_REAL_SIMILAR(im_trafo.apply(100), 100)

  std::vector<double> data;
  // test MS2-based ion mobility alignment
  {
      auto p = mc.getDefaults();
      p.setValue("mz_correction_function", "none");
      p.setValue("mz_extraction_window", 1.0);
      p.setValue("im_extraction_window", 100.0);
      mc.setParameters(p);

      std::vector< OpenSwath::SwathMap > swath_maps;
      swath_maps.push_back(map);
      swath_maps.push_back(ms1_map);
      TransformationDescription trafo_result;
      mc.correctIM(transition_group_map, targ_exp, swath_maps, pasef, trafo_result);
      data = swath_maps[0].sptr->getSpectrumById(0)->getMZArray()->data;


      TEST_REAL_SIMILAR(trafo_result.apply(10), 0.889721627408994)
      TEST_REAL_SIMILAR(trafo_result.apply(20), 10.0974304068522)
      TEST_REAL_SIMILAR(trafo_result.apply(30), 19.3051391862955)
  }

  // test MS2-based ion mobility alignment without MS1 map
  {
      auto p = mc.getDefaults();
      p.setValue("mz_correction_function", "none");
      p.setValue("mz_extraction_window", 1.0);
      p.setValue("im_extraction_window", 100.0);
      mc.setParameters(p);

      std::vector< OpenSwath::SwathMap > swath_maps;
      swath_maps.push_back(map);
      // swath_maps.push_back(ms1_map);
      TransformationDescription trafo_result;
      mc.correctIM(transition_group_map, targ_exp, swath_maps, pasef, trafo_result);
      data = swath_maps[0].sptr->getSpectrumById(0)->getMZArray()->data;


      TEST_REAL_SIMILAR(trafo_result.apply(10), 0.889721627408994)
      TEST_REAL_SIMILAR(trafo_result.apply(20), 10.0974304068522)
      TEST_REAL_SIMILAR(trafo_result.apply(30), 19.3051391862955)
  }

  // test MS2-based ion mobility from a single peptide
  {
      auto p = mc.getDefaults();
      p.setValue("mz_correction_function", "none");
      p.setValue("mz_extraction_window", 1.0);
      p.setValue("im_extraction_window", 100.0);
      mc.setParameters(p);

      std::vector< OpenSwath::SwathMap > swath_maps;
      OpenSwath::SwathMap map_single = map;
      map_single.upper = 425;

      swath_maps.push_back(map_single);
      // swath_maps.push_back(ms1_map);
      TransformationDescription trafo_result;
      mc.correctIM(transition_group_map, targ_exp, swath_maps, pasef, trafo_result);
      data = swath_maps[0].sptr->getSpectrumById(0)->getMZArray()->data;

      // only got a single peptide, so regression is only intercept
      TEST_REAL_SIMILAR(trafo_result.apply(10), 11)
      TEST_REAL_SIMILAR(trafo_result.apply(20), 11)
      TEST_REAL_SIMILAR(trafo_result.apply(30), 11)
  }

  // test MS2-based ion mobility alignment for PASEF data (should exclude second spectrum)
  {
      auto p = mc.getDefaults();
      bool pasef = true;
      p.setValue("mz_correction_function", "none");
      p.setValue("mz_extraction_window", 1.0);
      p.setValue("im_extraction_window", 100.0);
      mc.setParameters(p);

      std::vector< OpenSwath::SwathMap > swath_maps;
      swath_maps.push_back(map);
      swath_maps.push_back(mapPASEF);
      swath_maps.push_back(ms1_map);
      TransformationDescription trafo_result;
      mc.correctIM(transition_group_map, targ_exp, swath_maps, pasef, trafo_result);
      data = swath_maps[0].sptr->getSpectrumById(0)->getMZArray()->data;

      TEST_REAL_SIMILAR(trafo_result.apply(10), 0.889721627408994)
      TEST_REAL_SIMILAR(trafo_result.apply(20), 10.0974304068522)
      TEST_REAL_SIMILAR(trafo_result.apply(30), 19.3051391862955)
  }

  /*
  // test MS1 map when no MS1 is present
  {
      auto p = mc.getDefaults();
      p.setValue("mz_correction_function", "none");
      p.setValue("mz_extraction_window", 1.0);
      p.setValue("im_extraction_window", 100.0);
      p.setValue("ms1_im_calibration", "true");
      mc.setParameters(p);

      map.upper = 800;
      std::vector< OpenSwath::SwathMap > swath_maps;
      swath_maps.push_back(map);
      // swath_maps.push_back(ms1_map);
      TransformationDescription trafo_result;
      TEST_EXCEPTION(OpenMS::Exception::UnableToFit, mc.correctIM(transition_group_map, targ_exp, swath_maps, pasef, trafo_result));
  }
  */

  // test MS1 map when no MS2 is present
  // This test does not work due to preconditions
  {
      auto p = mc.getDefaults();
      p.setValue("mz_correction_function", "none");
      p.setValue("mz_extraction_window", 1.0);
      p.setValue("im_extraction_window", 100.0);
      p.setValue("ms1_im_calibration", "true");
      mc.setParameters(p);

      map.upper = 800;
      std::vector< OpenSwath::SwathMap > swath_maps;
      // swath_maps.push_back(map);
      swath_maps.push_back(ms1_map);
      TransformationDescription trafo_result;
      TEST_EXCEPTION(OpenMS::Exception::UnableToFit, mc.correctIM(transition_group_map, targ_exp, swath_maps, pasef, trafo_result));
      // this could work in principle but in practice this just fails as an MS2 is expected
  }

  // test MS1 ion mobility alignment
  {
      auto p = mc.getDefaults();
      p.setValue("mz_correction_function", "none");
      p.setValue("mz_extraction_window", 1.0);
      p.setValue("im_extraction_window", 100.0);
      p.setValue("ms1_im_calibration", "true");
      mc.setParameters(p);

      map.upper = 800;
      std::vector< OpenSwath::SwathMap > swath_maps;
      swath_maps.push_back(map);
      swath_maps.push_back(ms1_map);
      TransformationDescription trafo_result;
      mc.correctIM(transition_group_map, targ_exp, swath_maps, pasef, trafo_result);

      TEST_REAL_SIMILAR(trafo_result.apply(10), 0.835820895522389)
      TEST_REAL_SIMILAR(trafo_result.apply(20), 10.089552238806)
      TEST_REAL_SIMILAR(trafo_result.apply(30), 19.3432835820896)
  }

}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
