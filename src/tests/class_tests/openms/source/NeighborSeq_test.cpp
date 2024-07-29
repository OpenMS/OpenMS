// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/ANALYSIS/ID/NeighborSeq.h>
#include <vector>


using namespace OpenMS;
using namespace std;

START_TEST(NeighborSeq, "$Id$")

// NS()=delete;


// Test section for the generateSpectrum function
// The spectra were generated via TOPPView and contained b-and y-ion
START_SECTION(MSSpectrum generateSpectrum(const String& peptide_sequence))
{
  MSSpectrum spec_1 = NeighborSeq::generateSpectrum("PEPT");
  MSSpectrum spec_2 = NeighborSeq::generateSpectrum("AR");
  MSSpectrum spec_3 = NeighborSeq::generateSpectrum("VGLPINQR");

  // Test "PEPT"
  TEST_REAL_SIMILAR(spec_1[0].getPosition()[0], 120.0655);
  TEST_REAL_SIMILAR(spec_1[1].getPosition()[0], 217.1182);
  TEST_REAL_SIMILAR(spec_1[2].getPosition()[0], 227.1026);
  TEST_REAL_SIMILAR(spec_1[3].getPosition()[0], 324.1553);
  TEST_REAL_SIMILAR(spec_1[4].getPosition()[0], 346.1608);
  

  // Test "AR"
  TEST_REAL_SIMILAR(spec_2[0].getPosition()[0], 175.1189);

  // Test "VGLPINQR"
  TEST_REAL_SIMILAR(spec_3[0].getPosition()[0], 157.0971);
  TEST_REAL_SIMILAR(spec_3[1].getPosition()[0], 175.1189);
  TEST_REAL_SIMILAR(spec_3[2].getPosition()[0], 270.1812);
  TEST_REAL_SIMILAR(spec_3[3].getPosition()[0], 303.1775);
  TEST_REAL_SIMILAR(spec_3[4].getPosition()[0], 367.2339);
  TEST_REAL_SIMILAR(spec_3[5].getPosition()[0], 417.2204);
  TEST_REAL_SIMILAR(spec_3[6].getPosition()[0], 480.3180);
  TEST_REAL_SIMILAR(spec_3[7].getPosition()[0], 530.3045);
  TEST_REAL_SIMILAR(spec_3[8].getPosition()[0], 594.3609);
  TEST_REAL_SIMILAR(spec_3[9].getPosition()[0], 627.3578);
  TEST_REAL_SIMILAR(spec_3[10].getPosition()[0], 722.4195);
  TEST_REAL_SIMILAR(spec_3[11].getPosition()[0], 740.4413);
  TEST_REAL_SIMILAR(spec_3[12].getPosition()[0], 797.4628);

}
END_SECTION

// Test section for the compareSpectra function
START_SECTION(bool compareSpectra(const MSSpectrum& spec1, const MSSpectrum& spec2, const double& min_shared_ion_fraction, const double& mz_bin_size))
{
  MSSpectrum spec1;
  spec1.push_back(Peak1D(100.0, 1.0));
  spec1.push_back(Peak1D(200.0, 1.0));
  spec1.push_back(Peak1D(300.0, 1.0));

  MSSpectrum spec2;
  spec2.push_back(Peak1D(100.05, 1.0));
  spec2.push_back(Peak1D(200.05, 1.0));
  spec2.push_back(Peak1D(300.05, 1.0));

  MSSpectrum spec3;
  spec3.push_back(Peak1D(101.00, 1.0));
  spec3.push_back(Peak1D(201.00, 1.0));
  spec3.push_back(Peak1D(301.00, 1.0));

  MSSpectrum spec4;
  spec4.push_back(Peak1D(100.05, 1.0));
  spec4.push_back(Peak1D(201.00, 1.0));
  spec4.push_back(Peak1D(300.05, 1.0));
  spec4.push_back(Peak1D(301.00, 1.0));


  double min_shared_ion_fraction = 0.5; 
  

  //bin interval is from [a,b[
  bool compare_1_2_low = NeighborSeq::compareSpectra(spec1, spec2, min_shared_ion_fraction, 1.0);
  TEST_TRUE(compare_1_2_low)
  bool compare_1_3_low = NeighborSeq::compareSpectra(spec1, spec3, min_shared_ion_fraction, 1.0);
  TEST_FALSE(compare_1_3_low)
  bool compare_1_4_low = NeighborSeq::compareSpectra(spec1, spec4, min_shared_ion_fraction, 1.0);
  TEST_TRUE(compare_1_4_low)
  bool compare_2_3_low = NeighborSeq::compareSpectra(spec2, spec3, min_shared_ion_fraction, 1.0);
  TEST_FALSE(compare_2_3_low)
  bool compare_2_4_low = NeighborSeq::compareSpectra(spec2, spec4, min_shared_ion_fraction, 1.0);
  TEST_TRUE(compare_2_4_low)
  bool compare_3_4_low = NeighborSeq::compareSpectra(spec3, spec4, min_shared_ion_fraction, 1.0);
  TEST_TRUE(compare_3_4_low)

  bool compare_1_2_high = NeighborSeq::compareSpectra(spec1, spec2, min_shared_ion_fraction, 0.05);
  TEST_FALSE(compare_1_2_high)  
  bool compare_1_3_high = NeighborSeq::compareSpectra(spec1, spec3, min_shared_ion_fraction, 0.05);
  TEST_FALSE(compare_1_3_high)
  bool compare_1_4_high = NeighborSeq::compareSpectra(spec1, spec4, min_shared_ion_fraction, 0.05);
  TEST_FALSE(compare_1_4_high)
  bool compare_2_3_high = NeighborSeq::compareSpectra(spec2, spec3, min_shared_ion_fraction, 0.05);
  TEST_FALSE(compare_2_3_high)
  bool compare_2_4_high = NeighborSeq::compareSpectra(spec2, spec4, min_shared_ion_fraction, 0.05);
  TEST_TRUE(compare_2_4_high) 
  bool compare_3_4_high = NeighborSeq::compareSpectra(spec3, spec4, min_shared_ion_fraction, 0.05);
  TEST_TRUE(compare_3_4_high)
}
END_SECTION


// Test section for the findCandidatePositions function
START_SECTION(vector<int> findCandidatePositions(const map<double, vector<int>>& mass_position_map,
                                                 const double& mono_weight,
                                                 const double& mass_tolerance))
{
  map<double, vector<int>> mass_position_map_1 = {{100.0, {1, 2}}, {101.5, {3}}, {102.0, {4, 5}}};
  double mono_weight_1 = 101.0;
  double mass_tolerance_1 = 1.0;

  vector<int> expected_1 = {1, 2, 3, 4 , 5};
  vector<int> result_1 = NeighborSeq::findCandidatePositions(mass_position_map_1, mono_weight_1, mass_tolerance_1);
  TEST_EQUAL(result_1,expected_1)

  map<double, vector<int>> mass_position_map_2 = {{100.0, {1, 2}}, {101.5, {3}}, {102.0, {4, 5}}};
  double mono_weight_2 = 100.5;
  double mass_tolerance_2 = 1.0;

  vector<int> expected_2 = {1, 2, 3};
  vector<int> result_2 = NeighborSeq::findCandidatePositions(mass_position_map_2, mono_weight_2, mass_tolerance_2);
  TEST_EQUAL(result_2,expected_2)

  map<double, vector<int>> mass_position_map_3 = {{100.0, {1, 2}}, {101.5, {3}}, {102.0, {4, 5}}};
  double mono_weight_3 = 105.0;
  double mass_tolerance_3 = 1.0;

  vector<int> expected_3 = {};
  vector<int> result_3 = NeighborSeq::findCandidatePositions(mass_position_map_3, mono_weight_3, mass_tolerance_3);
  TEST_EQUAL(result_3, expected_3)

  map<double, vector<int>> mass_position_map_4 = {{100.0, {1, 2}}, {101.5, {3}}, {102.0, {4, 5}}};
  double mono_weight_4 = 100.5;
  double mass_tolerance_4 = -1.0;

  vector<int> expected_4 = {1, 2, 3};
  vector<int> result_4 = NeighborSeq::findCandidatePositions(mass_position_map_4, mono_weight_4, mass_tolerance_4);
  TEST_EQUAL(result_4, expected_4)
}
END_SECTION

/*
START_SECTION(map<double, vector<int>> NeighborSeq::createMassPositionMap(const vector<AASequence>& candidates))
{
  vector<AASequence> candidates_1;
  AASequence seq1 = AASequence::fromString("PEPTIDE"); // 728.371 
  AASequence seq2 = AASequence::fromString("PROTEIN"); // 799.349
  AASequence seq3 = AASequence::fromString("SEQUENCE");// 837.270


  candidates_1.push_back(seq1); // 1 position
  candidates_1.push_back(seq2); // 2 position
  candidates_1.push_back(seq3); // 3 position
  map<double, vector<int>> expected_1 = {{728.371, {2}}, {799.349, {1}}, {837.270, {3}}};
  map<double, vector<int>> result_1 = NeighborSeq::createMassPositionMap(candidates_1);
 // TEST_EQUAL(result_1,expected_1)

  candidates_1.push_back(seq1); // 4 position
  map<double, vector<int>> expected_2 = {{728.371, {2}}, {799.349, {1,4}}, {837.270, {3}}};
  map<double, vector<int>> result_2 = NeighborSeq::createMassPositionMap(candidates_1);
  //TEST_EQUAL(result_2,expected_2)
  //TEST_TRY(result_2, expected_2);
}
END_SECTION
*/


// Test section for the findNeighborPeptides_ function
// The spectra were generated via TOPPView
START_SECTION(vector<int> findNeighborPeptides(const AASequence& peptides,
                                               const vector<AASequence>& neighbor_candidate,
                                               const vector<int>& candidate_position,
                                               const double& min_shared_ion_fraction,
                                               const double& mz_bin_size))
{

  AASequence seq1 = AASequence::fromString("RIPLANGR");
  AASequence seq2 = AASequence::fromString("PPPPP");
  AASequence seq3 = AASequence::fromString("SEQUENCE");
  AASequence seq4 = AASequence::fromString("VGLPINQR");
  AASequence seq5 = AASequence::fromString("A");
  AASequence seq6 = AASequence::fromString("X");
  vector<AASequence> neighbor_candidate;
  neighbor_candidate.push_back(seq1);
  neighbor_candidate.push_back(seq2);
  neighbor_candidate.push_back(seq3);
  neighbor_candidate.push_back(seq1);
  neighbor_candidate.push_back(seq4);
  neighbor_candidate.push_back(seq5);
  neighbor_candidate.push_back(seq6);
  vector<int> candidate_position = {0, 1, 2, 3, 4, 5};

  vector<int> expected = {0,3, 4};
  double min_shared_ion_fraction = 0.25;
  double mz_bin_size = 0.05;

  vector<int> result = NeighborSeq::findNeighborPeptides(seq1, neighbor_candidate, candidate_position, min_shared_ion_fraction, mz_bin_size);
  
   TEST_EQUAL(expected, result)
}
END_SECTION


// Test section for the compareSpectra function
START_SECTION(bool compareShareSpectra(const MSSpectrum& spec1, const MSSpectrum& spec2, const double& mz_bin_size))
{
  MSSpectrum spec1;
  spec1.push_back(Peak1D(100.0, 1.0));
  spec1.push_back(Peak1D(200.0, 1.0));
  spec1.push_back(Peak1D(300.0, 1.0));

  MSSpectrum spec2;
  spec2.push_back(Peak1D(100.06, 1.0));
  spec2.push_back(Peak1D(200.06, 1.0));
  spec2.push_back(Peak1D(300.06, 1.0));

  MSSpectrum spec3;
  spec3.push_back(Peak1D(101.00, 1.0));
  spec3.push_back(Peak1D(201.00, 1.0));
  spec3.push_back(Peak1D(301.00, 1.0));

  MSSpectrum spec4;
  spec4.push_back(Peak1D(100.06, 1.0));
  spec4.push_back(Peak1D(201.00, 1.0));
  spec4.push_back(Peak1D(300.06, 1.0));
  spec4.push_back(Peak1D(301.00, 1.0));


  double min_shared_ion_fraction = 0.5;


  // bin interval is from [a,b[
  bool compare_1_2_low = NeighborSeq::compareShareSpectra(spec1, spec2, 1.0);
  TEST_EQUAL(compare_1_2_low,3)
  bool compare_1_3_low = NeighborSeq::compareShareSpectra(spec1, spec3, 1.0);
  TEST_EQUAL(compare_1_3_low,0)
  bool compare_1_4_low = NeighborSeq::compareShareSpectra(spec1, spec4, 1.0);
  TEST_EQUAL(compare_1_4_low,2)
  bool compare_2_3_low = NeighborSeq::compareShareSpectra(spec2, spec3, 1.0);
  TEST_EQUAL(compare_2_3_low,0)
  bool compare_2_4_low = NeighborSeq::compareShareSpectra(spec2, spec4, 1.0);
  TEST_EQUAL(compare_2_4_low,3)
  bool compare_3_4_low = NeighborSeq::compareShareSpectra(spec3, spec4, 1.0);
  TEST_EQUAL(compare_3_4_low,2)

  bool compare_1_2_high = NeighborSeq::compareShareSpectra(spec1, spec2, 0.05);
  TEST_EQUAL(compare_1_2_high,0)
  bool compare_1_3_high = NeighborSeq::compareShareSpectra(spec1, spec3, 0.05);
  TEST_EQUAL(compare_1_3_high,0)
  bool compare_1_4_high = NeighborSeq::compareShareSpectra(spec1, spec4, 0.05);
  TEST_EQUAL(compare_1_4_high,0)
  bool compare_2_3_high = NeighborSeq::compareShareSpectra(spec2, spec3, 0.05);
  TEST_EQUAL(compare_2_3_high,0)
  bool compare_2_4_high = NeighborSeq::compareShareSpectra(spec2, spec4, 0.05);
  TEST_EQUAL(compare_2_4_high,2)
  bool compare_3_4_high = NeighborSeq::compareShareSpectra(spec3, spec4, 0.05);
  TEST_EQUAL(compare_3_4_high,2)
  
}
END_SECTION

END_TEST