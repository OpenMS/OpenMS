// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/ANALYSIS/ID/NeighborSeq.h>

using namespace OpenMS;
using namespace std;

START_TEST(NeighborSeq, "$Id$")

//NS()=delete;


START_SECTION(double calculateMass(const AASequence& peptide))
{
  double mass_1 = NeighborSeq::calculateMass("PEPTIDE");
  double mass_2 = NeighborSeq::calculateMass("PPPPPPP");
  double mass_3 = NeighborSeq::calculateMass("PEPEP");
  // Compare the calculated masses with expected values
  TEST_REAL_SIMILAR(mass_1, 799.815);
  TEST_REAL_SIMILAR(mass_2, 697.814);
  TEST_REAL_SIMILAR(mass_3, 567.579);
}
END_SECTION

// Test section for the generateSpectrum function
START_SECTION(MSSpectrum generateSpectrum(const String& peptide_sequence))
{
  MSSpectrum spec_1 = NeighborSeq::generateSpectrum("PEPT");
  MSSpectrum spec_2 = NeighborSeq::generateSpectrum("PPP");
  MSSpectrum spec_3 = NeighborSeq::generateSpectrum("AR");
  
  // Define the expected spectrum for "PEPT"
  MSSpectrum spec_1_test;
 
  Peak1D peak1;
  peak1.setMZ(120.06552075);
  peak1.setIntensity(1);
  spec_1_test.push_back(peak1);

  Peak1D peak2;
  peak2.setMZ(217.11828498);
  peak2.setIntensity(1);
  spec_1_test.push_back(peak2);

  Peak1D peak3;
  peak3.setMZ(227.10263491);
  peak3.setIntensity(1);
  spec_1_test.push_back(peak3);

  Peak1D peak4;
  peak4.setMZ(324.15539914);
  peak4.setIntensity(1);
  spec_1_test.push_back(peak4);

  Peak1D peak5;
  peak5.setMZ(346.16087920);
  peak5.setIntensity(1);
  spec_1_test.push_back(peak5);

// Define the expected spectrum for "PPP"
  MSSpectrum spec_2_test;
  Peak1D peak1;
  peak1.setMZ(116.07060575);
  peak1.setIntensity(1);
  spec_2_test.push_back(peak1);

  Peak1D peak2;
  peak2.setMZ(195.11280491);
  peak2.setIntensity(1);
  spec_2_test.push_back(peak2);

  Peak1D peak3;
  peak3.setMZ(213.12336998);
  peak3.setIntensity(1);
  spec_3_test.push_back(peak3);

// Define the expected spectrum for "AR"
  MSSpectrum spec_3_test;
  Peak1D peak1;
  peak1.setMZ(175.11895291);
  peak1.setIntensity(1);
  spec_3_test.push_back(peak1);

// Compare the generated spectra with the expected spectra
  TEST_EQUAL(spec_1, spec_1_test);
  TEST_EQUAL(spec_1, spec_2_test);
  TEST_EQUAL(spec_1, spec_3_test);
  TEST_EQUAL(spec_2, spec_1_test);
  TEST_EQUAL(spec_2, spec_2_test);
  TEST_EQUAL(spec_2, spec_3_test);
  TEST_EQUAL(spec_3, spec_1_test);
  TEST_EQUAL(spec_3, spec_2_test);
  TEST_EQUAL(spec_3, spec_3_test);
}
END_SECTION

// Test section for the compareSpectra function
START_SECTION(bool compareSpectra(const MSSpectrum& spec1, const MSSpectrum& spec2, const double& ion_tolerance, const string& resolution))
{
  MSSpectrum spec1;
  spec1.push_back(Peak1D(100.0, 1.0));
  spec1.push_back(Peak1D(200.0, 1.0));
  spec1.push_back(Peak1D(300.0, 1.0));
  spec1.push_back(Peak1D(400.0, 1.0));


  MSSpectrum spec2;
  spec2.push_back(Peak1D(100.05, 1.0));
  spec2.push_back(Peak1D(200.05, 1.0));
  spec2.push_back(Peak1D(300.05, 1.0));
  spec2.push_back(Peak1D(400.06, 1.0));


  MSSpectrum spec3;
  spec3.push_back(Peak1D(101.05, 1.0));
  spec3.push_back(Peak1D(201.05, 1.0));
  spec3.push_back(Peak1D(301.05, 1.0));
  spec3.push_back(Peak1D(301.06, 1.0));


  double ion_tolerance = 0.25;
  
  bool compare_1_2_low = NeighborSeq::compareSpectra(spec1, spec2, ion_tolerance, "low");
  bool compare_1_3_low = NeighborSeq::compareSpectra(spec1, spec3, ion_tolerance, "low");
  bool compare_2_3_low = NeighborSeq::compareSpectra(spec2, spec3, ion_tolerance, "low");

  bool compare_1_2_high = NeighborSeq::compareSpectra(spec1, spec2, ion_tolerance, "high");
  bool compare_1_3_high = NeighborSeq::compareSpectra(spec1, spec3, ion_tolerance, "high");
  bool compare_2_3_high = NeighborSeq::compareSpectra(spec2, spec3, ion_tolerance, "high");

// Verify the results of the comparisons
  TEST_EQUAL(compare_1_2_low, true)
  TEST_EQUAL(compare_1_3_low, false)
  TEST_EQUAL(compare_1_3_low, false)
  TEST_EQUAL(compare_1_2_high, true)
  TEST_EQUAL(compare_1_3_high, false)
  TEST_EQUAL(compare_2_3_high, true)
}
END_SECTION


// Test section for the findNeighborPeptides_ function
START_SECTION(vector<FASTAFile::FASTAEntry> findNeighborPeptides(const AASequence& peptides,
                             const vector<FASTAFile::FASTAEntry>& neighbor_file,
                             const double& mass_tolerance,
                             const double& ion_tolerance,
                             const String& resolution))
{
  AASequence peptide = AASequence::fromString("VGLPINQR")

  vector<FASTAFile::FASTAEntry> neighbor_file;
  neighbor_file.push_back(FASTAFile::FASTAEntry("protein1", "RIPLANGR"));
  neighbor_file.push_back(FASTAFile::FASTAEntry("protein2", "PPPPPPPP"));
  neighbor_file.push_back(FASTAFile::FASTAEntry("protein3", "VGLPINQR"));
  neighbor_file.push_back(FASTAFile::FASTAEntry("protein4", "RIPLANGR"));

  double mass_tolerance = 0.00004;
  double ion_tolerance = 0.25;
  String resolution = "low";

  vector<FASTAFile::FASTAEntry> neighbors = findNeighborPeptides(peptide, neighbor_file, mass_tolerance, ion_tolerance, resolution);
  vector<FASTAFile::FASTAEntry> neighbors_test;
  neighbor_test.push_back(FASTAFile::FASTAEntry("RIPLANGR"));
  neighbor_test.push_back(FASTAFile::FASTAEntry("VGLPINQRgit clone git@github.com:_YOURUSERNAME_/OpenMS.git"));
  
  // Verify the found neighbors match the expected neighbors
  TEST_EQUAL(neighbor, neighbor_test);

}
END_SECTION

END_TEST
