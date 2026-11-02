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

#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/IONMOBILITY/FAIMSHelper.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/Feature.h>

#include <OpenMS/IONMOBILITY/IMDataArrayUtils.h>

#include <cmath>
#include <random>
#include <tuple>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace OpenMS;
using namespace std;

namespace
{
  // Synthetic MS1 run for the parallel paths: 600 scans (crosses the 256-scan block of the hill detection), enough
  // hills for several parallel chunks, isotope envelopes of charge 1-10, gaps, equal-intensity twins within htol,
  // empty scans and noise. With `im`, every peak carries an ion mobility (1/K0) in a float data array; with `faims`, the
  // scans alternate between two FAIMS CVs (two groups, processed in parallel).
  MSExperiment syntheticRun(bool im, bool faims = false)
  {
    std::mt19937 rng(42);
    auto uni = [&](double lo, double hi) { return lo + (hi - lo) * (rng() / 4294967296.0); };
    struct Env { double mz; int z; int start; int len; double height; double im; bool twin; };
    vector<Env> envs;
    for (int e = 0; e < 1200; ++e)
    {
      const int z = 1 + static_cast<int>(rng() % 10);
      envs.push_back({uni(400.0, 1400.0), z, static_cast<int>(rng() % 560), 8 + static_cast<int>(rng() % 50),
                      uni(1e4, 1e6), uni(0.7, 1.3), rng() % 10 == 0});
    }
    const double ratio[] = {1.0, 0.8, 0.5, 0.25, 0.1};
    MSExperiment exp;
    for (int s = 0; s < 600; ++s)
    {
      MSSpectrum spec;
      spec.setRT(s * 1.0);
      spec.setMSLevel(1);
      if (faims)
      {
        spec.setDriftTime(s % 2 ? -60.0 : -45.0);
        spec.setDriftTimeUnit(DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE);
      }
      vector<float> ims;
      auto add = [&](double mz, double it, double mob) { spec.push_back(Peak1D(mz, static_cast<float>(it))); ims.push_back(static_cast<float>(mob)); };
      if (s % 97 != 50)   // a few empty scans, some of them next to the block boundaries
      {
        for (const Env& e : envs)
        {
          if (s < e.start || s >= e.start + e.len || rng() % 20 == 0) continue;   // outside the elution, or a gap
          const double x = (s - e.start - 0.5 * e.len) / (0.25 * e.len);
          const double h = e.height * std::exp(-0.5 * x * x);
          for (int k = 0; k < 5; ++k)
          {
            add(e.mz + k * 1.003355 / e.z, h * ratio[k], e.im + uni(-0.005, 0.005));
            if (e.twin && k == 0) add(e.mz * (1.0 + 2e-6), h, e.im);   // same intensity, 2 ppm away
          }
        }
        for (int n = 0; n < 150; ++n) add(uni(350.0, 1500.0), uni(50.0, 5000.0), uni(0.6, 1.4));
      }
      if (im)
      {
        spec.getFloatDataArrays().resize(1);
        auto& fda = spec.getFloatDataArrays()[0];
        fda.assign(ims.begin(), ims.end());
        IMDataArrayUtils::setIMUnit(fda, DriftTimeUnit::VSSC);
        spec.setIMFormat(IMFormat::IM_PEAK);
      }
      spec.sortByPosition();
      exp.addSpectrum(spec);
    }
    return exp;
  }

  void runAt(int threads, const MSExperiment& exp, vector<Biosaur2Algorithm::Hill>& hills,
             vector<Biosaur2Algorithm::PeptideFeature>& features)
  {
#ifdef _OPENMP
    const int before = omp_get_max_threads();
    omp_set_num_threads(threads);
#endif
    Biosaur2Algorithm algo;
    algo.setMSData(exp);
    Param p = algo.getParameters();
    p.setValue("minmz", 350.0);
    p.setValue("maxmz", 1500.0);
    p.setValue("mini", 100.0);
    p.setValue("cmax", 10);
    algo.setParameters(p);
    FeatureMap fmap;
    algo.run(fmap, hills, features);
#ifdef _OPENMP
    omp_set_num_threads(before);
#else
    (void)threads;
#endif
  }

  auto hillKey(const Biosaur2Algorithm::Hill& h)
  {
    return std::tie(h.scan_indices, h.peak_indices, h.mz_values, h.intensities, h.rt_values, h.drift_times, h.ion_mobilities,
                    h.mz_weighted_mean, h.rt_start, h.rt_end, h.rt_apex, h.intensity_apex, h.intensity_sum,
                    h.drift_time_median, h.ion_mobility_median, h.length, h.hill_idx);
  }

  auto featureKey(const Biosaur2Algorithm::PeptideFeature& f)
  {
    return std::tie(f.mz, f.rt_start, f.rt_end, f.rt_apex, f.intensity_apex, f.intensity_sum, f.charge, f.n_isotopes,
                    f.n_scans, f.mass_calib, f.drift_time, f.ion_mobility, f.mono_hill_idx);
  }

  bool sameHills(const vector<Biosaur2Algorithm::Hill>& a, const vector<Biosaur2Algorithm::Hill>& b)
  {
    return std::equal(a.begin(), a.end(), b.begin(), b.end(), [](const auto& x, const auto& y) { return hillKey(x) == hillKey(y); });
  }

  bool sameFeatures(const vector<Biosaur2Algorithm::PeptideFeature>& a, const vector<Biosaur2Algorithm::PeptideFeature>& b)
  {
    auto isoKey = [](const Biosaur2Algorithm::IsotopeCandidate& i) { return std::tie(i.hill_idx, i.isotope_number, i.mass_diff_ppm, i.cos_corr); };
    return std::equal(a.begin(), a.end(), b.begin(), b.end(), [&](const auto& x, const auto& y)
    {
      return featureKey(x) == featureKey(y) &&
             std::equal(x.isotopes.begin(), x.isotopes.end(), y.isotopes.begin(), y.isotopes.end(),
                        [&](const auto& u, const auto& v) { return isoKey(u) == isoKey(v); });
    });
  }
}

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
  // The ranges have to describe the restored data: splitByFAIMSCV() ends in
  // exp.clear(true), which clears them, so they are only correct if recomputed.
  // getMinMZ() is deliberately not asserted: MSExperiment::updateRanges() folds each
  // chromatogram's m/z into the combined range, and the TIC above carries no precursor,
  // so the minimum is that chromatogram's 0 rather than the lowest peak.
  TEST_REAL_SIMILAR(back.getMinRT(), 0.0)
  TEST_REAL_SIMILAR(back.getMaxRT(), 19.0)
  TEST_REAL_SIMILAR(back.getMaxMZ(), 502.003355)
}
END_SECTION

START_SECTION([EXTRA] run() preserves chromatograms and settings on the FAIMS profile-mode path)
{
  // profile_mode centroids into a fresh experiment; replacing the whole experiment
  // would drop the chromatograms and experimental settings before the FAIMS split
  // ever saw them, so the reassembly above would have nothing to put back.
  Biosaur2Algorithm algo;
  MSExperiment exp;

  const double cvs[2] = {-45.0, -65.0};
  for (int i = 0; i < 10; ++i)
  {
    for (int g = 0; g < 2; ++g)
    {
      MSSpectrum spec;
      spec.setRT(i * 2.0 + g);
      spec.setMSLevel(1);
      spec.setType(SpectrumSettings::SpectrumType::PROFILE);
      spec.setDriftTime(cvs[g]);
      spec.setDriftTimeUnit(DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE);
      // a small profile peak shape around 500 + g
      for (int k = -3; k <= 3; ++k)
      {
        Peak1D p;
        p.setMZ(500.0 + g + k * 0.01);
        p.setIntensity(10000.0f * static_cast<float>(std::exp(-0.5 * k * k)));
        spec.push_back(p);
      }
      exp.addSpectrum(spec);
    }
  }

  MSChromatogram chrom;
  chrom.setNativeID("TIC");
  exp.addChromatogram(chrom);
  exp.getExperimentalSettings().setDateTime(DateTime::fromString("2019-09-07T09:40:04"));

  algo.setMSData(exp);

  Param p = algo.getParameters();
  p.setValue("profile", "true");
  p.setValue("minmz", 400.0);
  p.setValue("maxmz", 600.0);
  p.setValue("mini", 100.0);
  algo.setParameters(p);

  FeatureMap fmap;
  algo.run(fmap);

  const MSExperiment& back = algo.getMSData();
  TEST_EQUAL(back.size(), 20)
  TEST_EQUAL(back.getChromatograms().size(), 1)
  TEST_EQUAL(back.getExperimentalSettings().getDateTime().toString(), "2019-09-07T09:40:04")
  // centroiding must not strip the FAIMS annotation, or the split never happens
  TEST_EQUAL(FAIMSHelper::getCompensationVoltages(back).size(), 2)
}
END_SECTION

START_SECTION([EXTRA] run() refreshes ranges on the non-FAIMS profile-mode path)
{
  // The in-place path rewrites the spectra too: centroidProfileSpectra_ replaces every
  // peak, so ranges cached before run() describe data that no longer exists. This is the
  // non-FAIMS counterpart to the range assertions in the FAIMS section above.
  Biosaur2Algorithm algo;
  MSExperiment exp;

  for (int i = 0; i < 10; ++i)
  {
    MSSpectrum spec;
    spec.setRT(i * 1.0);
    spec.setMSLevel(1);
    spec.setType(SpectrumSettings::SpectrumType::PROFILE);
    for (int k = -3; k <= 3; ++k)
    {
      Peak1D p;
      p.setMZ(500.0 + k * 0.01);
      p.setIntensity(10000.0f * static_cast<float>(std::exp(-0.5 * k * k)));
      spec.push_back(p);
    }
    exp.addSpectrum(spec);
  }
  exp.updateRanges(); // ranges describing the *profile* peaks, stale once centroided

  algo.setMSData(exp);

  Param p = algo.getParameters();
  p.setValue("profile", "true");
  p.setValue("minmz", 400.0);
  p.setValue("maxmz", 600.0);
  p.setValue("mini", 100.0);
  algo.setParameters(p);

  FeatureMap fmap;
  algo.run(fmap);

  const MSExperiment& back = algo.getMSData();
  TEST_EQUAL(back.size(), 10)
  TEST_EQUAL(FAIMSHelper::getCompensationVoltages(back).empty(), true)
  TEST_REAL_SIMILAR(back.getMinRT(), 0.0)
  TEST_REAL_SIMILAR(back.getMaxRT(), 9.0)
  // Compare the cached ranges against ones recomputed from the very same spectra: this
  // fails whenever the cache is stale, whatever the values happen to be. Asserting a
  // window around the peak would not -- the pre-centroiding profile span lies inside any
  // window wide enough to hold the centroid, so a stale cache would satisfy it.
  MSExperiment recomputed = back;
  recomputed.updateRanges();
  TEST_REAL_SIMILAR(back.getMinMZ(), recomputed.getMinMZ())
  TEST_REAL_SIMILAR(back.getMaxMZ(), recomputed.getMaxMZ())
  TEST_REAL_SIMILAR(back.getMinRT(), recomputed.getMinRT())
  TEST_REAL_SIMILAR(back.getMaxRT(), recomputed.getMaxRT())
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
  // ranges describe the returned spectra here too, not just on the FAIMS path
  TEST_REAL_SIMILAR(algo.getMSData().getMinRT(), 0.0)
  TEST_REAL_SIMILAR(algo.getMSData().getMaxRT(), 9.0)
}
END_SECTION

START_SECTION([EXTRA] run() gives identical hills and features at 1 and 4 threads: plain / ion mobility / two FAIMS CVs at cmax 10)
{
  // pins the parallel hill detection, processing, splitting, isotope search and feature build to the serial result;
  // cmax 10 also covers charges above 8, which once aborted the candidate search
  for (int mode : {0, 1, 2})
  {
    const MSExperiment run = syntheticRun(mode == 1, mode == 2);
    vector<Biosaur2Algorithm::Hill> h1, h4;
    vector<Biosaur2Algorithm::PeptideFeature> f1, f4;
    runAt(1, run, h1, f1);
    runAt(4, run, h4, f4);
    TEST_EQUAL(h1.size() > 5000, true)   // several parallel chunks
    TEST_EQUAL(f1.size() > 100, true)
    bool high_charge = false;
    for (const auto& f : f1) high_charge |= f.charge >= 9;
    TEST_EQUAL(high_charge, true)
    TEST_EQUAL(sameHills(h1, h4), true)
    TEST_EQUAL(sameFeatures(f1, f4), true)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST