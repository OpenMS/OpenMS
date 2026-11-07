// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FEATUREFINDER/MultiplexDeltaMasses.h>
#include <OpenMS/FEATUREFINDER/MultiplexFilteringProfile.h>
#include <OpenMS/FEATUREFINDER/MultiplexFilteredMSExperiment.h>
#include <OpenMS/FEATUREFINDER/MultiplexClustering.h>

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/ML/CLUSTERING/GridBasedClustering.h>

#include <random>

using namespace OpenMS;

START_TEST(MultiplexFilteringProfile, "$Id$")

// read data
MSExperiment exp;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MultiplexClustering.mzML"), exp);
exp.updateRanges();

// pick data
PeakPickerHiRes picker;
Param param = picker.getParameters();
param.setValue("ms_levels", ListUtils::create<Int>("1"));
param.setValue("signal_to_noise", 0.0);
picker.setParameters(param);
std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > boundaries_exp_s;
std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > boundaries_exp_c;
MSExperiment exp_picked;
picker.pickExperiment(exp, exp_picked, boundaries_exp_s, boundaries_exp_c);

// set parameters
int charge_min = 1;
int charge_max = 4;
int isotopes_per_peptide_min = 3;
int isotopes_per_peptide_max = 6;
double intensity_cutoff = 10.0;
double rt_band = 3;
double rt_typical = 90;
double mz_tolerance = 40;
bool mz_tolerance_unit = true;    // ppm (true), Da (false)
double peptide_similarity = 0.8;
double averagine_similarity = 0.75;
double averagine_similarity_scaling = 0.75;
std::string averagine_type="peptide";

// construct list of peak patterns
MultiplexDeltaMasses shifts1;
shifts1.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(0, "no_label"));
shifts1.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(8.0443702794, "Arg8"));
MultiplexDeltaMasses shifts2;
shifts2.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(0, "no_label"));
MultiplexDeltaMasses::LabelSet label_set;
label_set.insert("Arg8");
label_set.insert("Arg8");
shifts2.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(2*8.0443702794, label_set));
std::vector<MultiplexIsotopicPeakPattern> patterns;
for (int c = charge_max; c >= charge_min; --c)
{
    MultiplexIsotopicPeakPattern pattern1(c, isotopes_per_peptide_max, shifts1, 0);
    patterns.push_back(pattern1);
    MultiplexIsotopicPeakPattern pattern2(c, isotopes_per_peptide_max, shifts2, 1);
    patterns.push_back(pattern2);
}

MultiplexFilteringProfile filtering(exp, exp_picked, boundaries_exp_s, patterns, isotopes_per_peptide_min, isotopes_per_peptide_max, intensity_cutoff, rt_band, mz_tolerance, mz_tolerance_unit, peptide_similarity, averagine_similarity, averagine_similarity_scaling, averagine_type);
std::vector<MultiplexFilteredMSExperiment> filter_results = filtering.filter();

MultiplexClustering* nullPointer = nullptr;
MultiplexClustering* ptr;

START_SECTION(MultiplexClustering(const MSExperiment& exp_profile, const MSExperiment& exp_picked, const std::vector<std::vector<PeakPickerHiRes::PeakBoundary> >& boundaries, double rt_typical))
    // constructor-only test; the expensive cluster() run is exercised (and its
    // result asserted on) in the cluster() section below, so do not run it here.
    ptr = new MultiplexClustering(exp, exp_picked, boundaries_exp_s, rt_typical);
    TEST_NOT_EQUAL(ptr, nullPointer);
    delete ptr;
END_SECTION

MultiplexClustering clustering(exp, exp_picked, boundaries_exp_s, rt_typical);

START_SECTION(cluster(const std::vector<MultiplexFilteredMSExperiment>& filter_results))
    std::vector<std::map<int,GridBasedCluster> > cluster_results = clustering.cluster(filter_results);
    TEST_EQUAL(cluster_results[0].size(), 0);
    TEST_EQUAL(cluster_results[1].size(), 0);
    TEST_EQUAL(cluster_results[2].size(), 0);
    TEST_EQUAL(cluster_results[3].size(), 0);
    TEST_EQUAL(cluster_results[4].size(), 2);
    TEST_EQUAL(cluster_results[5].size(), 0);
    TEST_EQUAL(cluster_results[6].size(), 0);
    TEST_EQUAL(cluster_results[7].size(), 0);
END_SECTION

// cluster() clusters groups of points separated by empty m/z columns of the grid on their own. The clusters must be those of
// the clustering of all points of a filter result at once (as before): same indices, points, centres and bounding boxes.
START_SECTION(([EXTRA] cluster() of groups of points equals the clustering of all points))
{
  // an experiment that spans the grid (centroided constructor, m/z tolerance in Da)
  MSExperiment exp_range;
  for (int s = 0; s <= 300; s += 50)
  {
    MSSpectrum spectrum;
    spectrum.setMSLevel(1);
    spectrum.setRT(s);
    spectrum.push_back(Peak1D(400.0, 1.0));
    spectrum.push_back(Peak1D(410.0, 1.0));
    exp_range.addSpectrum(spectrum);
  }
  exp_range.updateRanges();
  const double tolerance = 0.01; // Da, i.e. the width of the m/z columns of the grid
  const double rt_typical_range = 30.0;
  MultiplexClustering clustering_range(exp_range, tolerance, false, rt_typical_range);

  // the grid of MultiplexClustering(exp, mz_tolerance, false, rt_typical)
  std::vector<double> grid_mz, grid_rt;
  const double mz_min = exp_range.spectrumRanges().byMSLevel(1).getMinMZ() - 1e-2, mz_max = exp_range.spectrumRanges().byMSLevel(1).getMaxMZ() + 1e-2;
  const double rt_min = exp_range.spectrumRanges().byMSLevel(1).getMinRT() - 1e-2, rt_max = exp_range.spectrumRanges().byMSLevel(1).getMaxRT() + 1e-2;
  for (double mz = mz_min; mz < mz_max; mz = mz + tolerance) { grid_mz.push_back(mz); }
  grid_mz.push_back(mz_max);
  for (double rt = rt_min; rt < rt_max; rt = rt + rt_typical_range) { grid_rt.push_back(rt); }
  grid_rt.push_back(rt_max);

  // Islands of points (elution profiles of a few m/z traces each), separated by 1 to 6 empty m/z columns, and one long
  // chain of points over many columns. Points of neighbouring islands can be nearest neighbours of each other.
  std::vector<MultiplexFilteredMSExperiment> results(2);
  std::mt19937 rng(7);
  std::uniform_real_distribution<double> jitter(-0.3, 0.3);
  double mz = 401.0;
  size_t index = 0;
  for (int island = 0; island < 60; ++island)
  {
    const int traces = 1 + island % 3;
    for (int trace = 0; trace < traces; ++trace)
    {
      for (double rt = 10.0 + (island * 7) % 50; rt < 200.0 + (island * 13) % 90; rt += 1.7)
      {
        const double point_mz = mz + (trace + 0.5 + jitter(rng)) * tolerance;
        results[0].addPeak(MultiplexFilteredPeak(point_mz, rt, index, index));
        ++index;
      }
    }
    mz += (traces + 1 + island % 6) * tolerance; // traces, then 1 to 6 empty columns
  }
  for (int k = 0; k < 400; ++k)
  {
    results[1].addPeak(MultiplexFilteredPeak(405.0 + k * 0.004, 5.0 + (k % 37) * 7.5, k, k));
  }

  std::vector<std::map<int, GridBasedCluster>> clusters = clustering_range.cluster(results);
  TEST_EQUAL(clusters.size(), 2)
  for (size_t r = 0; r < results.size(); ++r)
  {
    GridBasedClustering<MultiplexClustering::MultiplexDistance> all(MultiplexClustering::MultiplexDistance(tolerance / rt_typical_range), results[r].getMZ(), results[r].getRT(), grid_mz, grid_rt);
    all.cluster();
    const std::map<int, GridBasedCluster> expected = all.getResults();
    TEST_EQUAL(expected.size() > 10, true)
    TEST_EQUAL(clusters[r].size(), expected.size())
    Size differences = 0;
    for (auto it = clusters[r].cbegin(), it_expected = expected.cbegin(); it != clusters[r].cend() && it_expected != expected.cend(); ++it, ++it_expected)
    {
      const GridBasedCluster& c = it->second;
      const GridBasedCluster& e = it_expected->second;
      differences += !(it->first == it_expected->first && c.getPoints() == e.getPoints() && c.getCentre() == e.getCentre()
                       && c.getBoundingBox() == e.getBoundingBox() && c.getPropertiesB() == e.getPropertiesB());
    }
    TEST_EQUAL(differences, 0)
  }

  // a point outside of the grid is reported as before
  std::vector<MultiplexFilteredMSExperiment> outside(1);
  outside[0].addPeak(MultiplexFilteredPeak(1000.0, 100.0, 0, 0));
  TEST_EXCEPTION(Exception::IllegalArgument, clustering_range.cluster(outside))
}
END_SECTION

START_SECTION(([EXTRA] constructors reject a non-positive rt_typical instead of looping forever))
{
  // rt_typical <= 0 used to leave the "rt = rt + rt_typical" grid-building loop stuck forever,
  // growing its vector until std::bad_alloc (OpenMS issue #10770, seen when a chromatographic
  // FWHM estimate of 0 was passed down as rt_typical). Both constructors must now reject it
  // with a clear exception instead.
  TEST_EXCEPTION(Exception::IllegalArgument, MultiplexClustering(exp, exp_picked, boundaries_exp_s, 0.0))
  TEST_EXCEPTION(Exception::IllegalArgument, MultiplexClustering(exp, exp_picked, boundaries_exp_s, -1.0))
  TEST_EXCEPTION(Exception::IllegalArgument, MultiplexClustering(exp, 40.0, true, 0.0))
  TEST_EXCEPTION(Exception::IllegalArgument, MultiplexClustering(exp, 40.0, true, -1.0))
}
END_SECTION

END_TEST
