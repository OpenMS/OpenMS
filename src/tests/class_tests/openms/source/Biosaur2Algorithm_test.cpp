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

#include <OpenMS/IONMOBILITY/FAIMSHelper.h>
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

START_SECTION([EXTRA] run() preserves the stored MS data on the FAIMS path)
{
  // The FAIMS path splits ms_data_ by compensation voltage. It has to put the
  // spectra back afterwards, otherwise getMSData() hands the caller an emptied
  // experiment and a downstream algorithm silently works on nothing
  // (https://github.com/OpenMS/OpenMS/issues/9980).
  Biosaur2Algorithm algo;
  MSExperiment exp;

  // two CV groups, interleaved in acquisition order
  const double cvs[2] = {-45.0, -65.0};
  for (int i = 0; i < 10; ++i)
  {
    for (int g = 0; g < 2; ++g)
    {
      MSSpectrum spec;
      spec.setRT(i * 2.0 + g);
      spec.setMSLevel(1);
      spec.setDriftTime(cvs[g]);
      spec.setDriftTimeUnit(DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE);

      Peak1D p1;
      p1.setMZ(500.0 + g);
      p1.setIntensity(10000.0);
      spec.push_back(p1);

      Peak1D p2;
      p2.setMZ(500.0 + g + 1.003355);
      p2.setIntensity(5000.0);
      spec.push_back(p2);

      exp.addSpectrum(spec);
    }
  }
  // an MS2 spectrum and a chromatogram, neither of which should reappear as MS1
  MSSpectrum ms2;
  ms2.setRT(1.5);
  ms2.setMSLevel(2);
  ms2.setDriftTime(cvs[0]);
  ms2.setDriftTimeUnit(DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE);
  ms2.push_back(Peak1D(300.0, 100.0));
  exp.addSpectrum(ms2);

  MSChromatogram chrom;
  chrom.setNativeID("TIC");
  exp.addChromatogram(chrom);

  TEST_EQUAL(FAIMSHelper::getCompensationVoltages(exp).size(), 2)

  algo.setMSData(exp);

  Param p = algo.getParameters();
  p.setValue("minmz", 400.0);
  p.setValue("maxmz", 600.0);
  p.setValue("mini", 100.0);
  algo.setParameters(p);

  FeatureMap fmap;
  algo.run(fmap);

  // the 20 MS1 spectra are back (the MS2 spectrum is dropped by run(), as on the
  // non-FAIMS path), and so is the chromatogram
  const MSExperiment& back = algo.getMSData();
  TEST_EQUAL(back.size(), 20)
  TEST_EQUAL(back.getChromatograms().size(), 1)

  // both CV groups are represented, and acquisition (RT) order is restored
  // rather than left grouped by CV
  TEST_EQUAL(FAIMSHelper::getCompensationVoltages(back).size(), 2)
  bool rt_sorted = true;
  for (Size i = 1; i < back.size(); ++i)
  {
    if (back[i].getRT() < back[i - 1].getRT()) { rt_sorted = false; }
  }
  TEST_EQUAL(rt_sorted, true)
  TEST_REAL_SIMILAR(back[0].getRT(), 0.0)
  TEST_REAL_SIMILAR(back[back.size() - 1].getRT(), 19.0)
}
END_SECTION

START_SECTION([EXTRA] run() preserves the stored MS data on the non-FAIMS path)
{
  // the counterpart of the section above: both paths must leave getMSData() usable
  Biosaur2Algorithm algo;
  MSExperiment exp;
  for (int i = 0; i < 10; ++i)
  {
    MSSpectrum spec;
    spec.setRT(i * 1.0);
    spec.setMSLevel(1);
    Peak1D p1;
    p1.setMZ(500.0);
    p1.setIntensity(10000.0);
    spec.push_back(p1);
    exp.addSpectrum(spec);
  }
  TEST_EQUAL(FAIMSHelper::getCompensationVoltages(exp).empty(), true)

  algo.setMSData(exp);
  Param p = algo.getParameters();
  p.setValue("minmz", 400.0);
  p.setValue("maxmz", 600.0);
  p.setValue("mini", 100.0);
  algo.setParameters(p);

  FeatureMap fmap;
  algo.run(fmap);

  TEST_EQUAL(algo.getMSData().size(), 10)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST