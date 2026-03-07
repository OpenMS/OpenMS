// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FEATUREFINDER/Biosaur2Algorithm.h>
///////////////////////////

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/Feature.h>

using namespace OpenMS;
using namespace std;

START_TEST(Biosaur2Algorithm, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Biosaur2Algorithm* ptr = nullptr;
Biosaur2Algorithm* nullPointer = nullptr;
START_SECTION(Biosaur2Algorithm())
{
	ptr = new Biosaur2Algorithm();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~Biosaur2Algorithm())
{
	delete ptr;
}
END_SECTION

START_SECTION(void run(FeatureMap& feature_map))
{
  // Test case 1: Isotope calibration bias check
  // Create a synthetic dataset with a known isotope pattern
  // We want to test if the calibration correctly picks the highest intensity peak
  // when multiple peaks are within the tolerance window.

  Biosaur2Algorithm algo;
  MSExperiment exp;
  
  // Create a few spectra
  for (int i = 0; i < 10; ++i)
  {
    MSSpectrum spec;
    spec.setRT(i * 1.0);
    spec.setMSLevel(1);
    
    // Monoisotopic peak at 1000.0
    Peak1D p1;
    p1.setMZ(1000.0);
    p1.setIntensity(10000.0);
    spec.push_back(p1);

    // Isotope 1 at 1001.003355 (approx +1.003355 Da)
    // We add a "noise" peak very close to it, but with lower intensity
    // to see if the algorithm picks the correct one (higher intensity)
    // or just the first one it encounters.
    
    // Correct isotope peak
    Peak1D p2;
    p2.setMZ(1001.003355);
    p2.setIntensity(6000.0); // High intensity
    spec.push_back(p2);

    // Distracting peak (noise) slightly shifted but within tolerance
    // Let's say tolerance is 20ppm. 
    // 20 ppm at 1000 m/z is 0.02 Da.
    // So we put a peak at 1001.003355 - 0.01 = 1000.993355
    Peak1D p3;
    p3.setMZ(1000.993355);
    p3.setIntensity(100.0); // Low intensity
    spec.push_back(p3);

    exp.addSpectrum(spec);
  }

  algo.setMSData(exp);
  
  Param p = algo.getParameters();
  p.setValue("itol", 20.0); // 20 ppm tolerance
  p.setValue("minmz", 900.0);
  p.setValue("maxmz", 1100.0);
  algo.setParameters(p);

  FeatureMap fmap;
  algo.run(fmap);

  // We expect one feature
  TEST_EQUAL(fmap.size(), 1);
  
  if (!fmap.empty())
  {
     const Feature& f = fmap[0];
     TEST_REAL_SIMILAR(f.getMZ(), 1000.0);
     TEST_EQUAL(f.getCharge(), 1);
     // Check if we have isotopes
     // Biosaur2 stores convex hulls differently depending on "convex_hulls" parameter.
     // Default is "bounding_box" which stores 1 hull.
     // We can check n_isotopes meta value or switch to mass_traces mode.
     
     TEST_EQUAL(f.getMetaValue("n_isotopes"), 2);
  }
}
END_SECTION

START_SECTION(void setMSData(const MSExperiment& ms_data))
{
  Biosaur2Algorithm algo;
  MSExperiment exp;
  MSSpectrum spec;
  spec.setMSLevel(1);
  exp.addSpectrum(spec);
  
  algo.setMSData(exp);
  TEST_EQUAL(algo.getMSData().size(), 1);
}
END_SECTION

START_SECTION(void setMSData(MSExperiment&& ms_data))
{
  Biosaur2Algorithm algo;
  MSExperiment exp;
  MSSpectrum spec;
  spec.setMSLevel(1);
  exp.addSpectrum(spec);
  
  algo.setMSData(std::move(exp));
  TEST_EQUAL(algo.getMSData().size(), 1);
  TEST_EQUAL(exp.size(), 0); // Should be moved
}
END_SECTION

START_SECTION(MSExperiment& getMSData())
{
  Biosaur2Algorithm algo;
  MSExperiment& exp = algo.getMSData();
  TEST_EQUAL(exp.size(), 0);
  
  MSSpectrum spec;
  spec.setMSLevel(1);
  exp.addSpectrum(spec);
  
  TEST_EQUAL(algo.getMSData().size(), 1);
}
END_SECTION

START_SECTION(const MSExperiment& getMSData() const)
{
  Biosaur2Algorithm algo;
  const Biosaur2Algorithm& const_algo = algo;
  TEST_EQUAL(const_algo.getMSData().size(), 0);
}
END_SECTION

START_SECTION(void run(FeatureMap& feature_map, std::vector<Hill>& hills, std::vector<PeptideFeature>& peptide_features))
{
  Biosaur2Algorithm algo;
  MSExperiment exp;
  
  // Create a simple hill
  for (int i = 0; i < 5; ++i)
  {
    MSSpectrum spec;
    spec.setRT(i * 1.0);
    spec.setMSLevel(1);
    Peak1D p;
    p.setMZ(500.0);
    p.setIntensity(1000.0);
    spec.push_back(p);
    exp.addSpectrum(spec);
  }
  
  algo.setMSData(exp);
  
  Param p = algo.getParameters();
  p.setValue("minmz", 400.0);
  p.setValue("maxmz", 600.0);
  p.setValue("mini", 100.0);
  algo.setParameters(p);
  
  FeatureMap fmap;
  std::vector<Biosaur2Algorithm::Hill> hills;
  std::vector<Biosaur2Algorithm::PeptideFeature> features;
  
  algo.run(fmap, hills, features);
  
  TEST_EQUAL(hills.size() >= 1, true);
  // We might not get a feature if it doesn't have isotopes or charge state, 
  // but we should at least get a hill.
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST