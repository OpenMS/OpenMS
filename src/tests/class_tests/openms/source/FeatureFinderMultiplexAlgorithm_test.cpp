// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>

#include <OpenMS/FEATUREFINDER/FeatureFinderMultiplexAlgorithm.h>
#include <OpenMS/FEATUREFINDER/MultiplexDeltaMasses.h>
#include <OpenMS/FEATUREFINDER/MultiplexFilteredMSExperiment.h>
#include <OpenMS/FEATUREFINDER/MultiplexFilteredPeak.h>
#include <OpenMS/FEATUREFINDER/MultiplexIsotopicPeakPattern.h>
#include <OpenMS/FEATUREFINDER/MultiplexSatelliteCentroided.h>
#include <OpenMS/ML/CLUSTERING/GridBasedCluster.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/KERNEL/FeatureMap.h>

#include <cmath>
#include <map>
#include <vector>

using namespace OpenMS;
using namespace std;

// Gives the tests access to the feature generation of single multiplets from given satellites.
class FeatureFinderMultiplexAlgorithmTester : public FeatureFinderMultiplexAlgorithm
{
public:
  void setInput(const MSExperiment& exp, unsigned isotopes_per_peptide_max)
  {
    setParameters(getDefaults()); // e.g. 'algorithm:rt_min' = 2 s
    exp_centroid_ = exp;
    isotopes_per_peptide_max_ = isotopes_per_peptide_max;
  }
  std::vector<double> peptideIntensities(const MultiplexIsotopicPeakPattern& pattern, const std::multimap<size_t, MultiplexSatelliteCentroided>& satellites)
  {
    return determinePeptideIntensitiesCentroided_(pattern, satellites);
  }
  void generateMaps(const std::vector<MultiplexIsotopicPeakPattern>& patterns, const std::vector<MultiplexFilteredMSExperiment>& filter_results,
                    std::vector<std::map<int, GridBasedCluster>>& cluster_results)
  {
    generateMapsCentroided_(patterns, filter_results, cluster_results);
  }
};

// Four centroided MS1 spectra (RT 100, 101, 102 and 103 s) of a SILAC doublet (Arg10, 2+) with three isotopic peaks per
// peptide: peak index 0-2 are the light, 3-5 the heavy isotopes.
MSExperiment doubletSpectra()
{
  const double arg10 = 10.008268600;
  const std::vector<double> isotopes = {1.0, 0.6, 0.25};
  MSExperiment exp;
  for (int scan = 0; scan < 4; ++scan)
  {
    MSSpectrum spectrum;
    spectrum.setMSLevel(1);
    spectrum.setRT(100.0 + scan);
    for (const double mono : {500.0, 500.0 + arg10 / 2})
    {
      for (Size k = 0; k < isotopes.size(); ++k)
      {
        spectrum.push_back(Peak1D(mono + k * Constants::C13C12_MASSDIFF_U / 2, 1e5 * isotopes[k]));
      }
    }
    exp.addSpectrum(spectrum);
  }
  exp.updateRanges();
  return exp;
}

MultiplexIsotopicPeakPattern doubletPattern()
{
  std::vector<MultiplexDeltaMasses::DeltaMass> shifts = {MultiplexDeltaMasses::DeltaMass(0.0, "no_label"),
                                                         MultiplexDeltaMasses::DeltaMass(10.008268600, "Arg10")};
  return MultiplexIsotopicPeakPattern(2, 3, MultiplexDeltaMasses(shifts), 0);
}

// satellites of the given mass traces (peptide * 3 + isotope; the peak index is the same) in the given scans
void addSatellites(std::multimap<size_t, MultiplexSatelliteCentroided>& satellites, const std::vector<size_t>& traces, const std::vector<size_t>& scans)
{
  for (size_t trace : traces)
  {
    for (size_t scan : scans) { satellites.insert(std::make_pair(trace, MultiplexSatelliteCentroided(scan, trace))); }
  }
}
void addSatellites(MultiplexFilteredPeak& peak, const std::vector<size_t>& traces, const std::vector<size_t>& scans)
{
  for (size_t trace : traces)
  {
    for (size_t scan : scans) { peak.addSatellite(scan, trace, trace); }
  }
}

START_TEST(FeatureFinderMultiplexAlgorithm, "$Id$")

FeatureFinderMultiplexAlgorithm* ptr = 0;
FeatureFinderMultiplexAlgorithm* null_ptr = 0;
START_SECTION(FeatureFinderMultiplexAlgorithm())
{
  ptr = new FeatureFinderMultiplexAlgorithm();
  TEST_NOT_EQUAL(ptr, null_ptr);
}
END_SECTION

START_SECTION(~FeatureFinderMultiplexAlgorithm())
{
  delete ptr;
}
END_SECTION

START_SECTION((virtual void run()))
{
  MzMLFile mzml_file;
  MSExperiment exp;
  ConsensusMap result;
  
  mzml_file.getOptions().addMSLevel(1);
  mzml_file.load(OPENMS_GET_TEST_DATA_PATH("FeatureFinderMultiplex_1_input.mzML"), exp);
  exp.updateRanges();
  
  Param param;
  ParamXMLFile paramFile;
  paramFile.load(OPENMS_GET_TEST_DATA_PATH("FeatureFinderMultiplex_1_parameters.ini"), param);
  param = param.copy("FeatureFinderMultiplex:1:",true);
  param.remove("in");
  param.remove("out");
  param.remove("out_multiplets");
  param.remove("log");
  param.remove("debug");
  param.remove("threads");
  param.remove("no_progress");
  param.remove("force");
  param.remove("test");
  
  FeatureFinderMultiplexAlgorithm algorithm;
  algorithm.setParameters(param);
  algorithm.run(exp, true);
  result = algorithm.getConsensusMap();
  
  TEST_EQUAL(result.size(), 2);
  
  double L = result[0].getFeatures().begin()->getIntensity();
  double H = (++(result[0].getFeatures().begin()))->getIntensity();

  // Check that the HEAVY:LIGHT ratio is close to the expected 3:1 ratio
  TOLERANCE_ABSOLUTE(0.2);
  TEST_REAL_SIMILAR(H/L, 3.0);
  TOLERANCE_ABSOLUTE(1e-5);

  // every handle refers to the feature it was built from (its elution region), by unique id
  const FeatureMap& features = algorithm.getFeatureMap();
  TEST_EQUAL(features.size(), 4);
  std::map<UInt64, const Feature*> feature_of;
  for (const Feature& f : features) { feature_of[f.getUniqueId()] = &f; }
  for (const ConsensusFeature& cf : result)
  {
    for (const FeatureHandle& handle : cf.getFeatures())
    {
      TEST_EQUAL(feature_of.count(handle.getUniqueId()), 1)
      if (feature_of.count(handle.getUniqueId()) == 1)
      {
        const Feature& f = *feature_of[handle.getUniqueId()];
        TEST_REAL_SIMILAR(f.getMZ(), handle.getMZ())
        TEST_REAL_SIMILAR(f.getRT(), handle.getRT())
        TEST_EQUAL(f.getConvexHulls().empty(), false)
      }
    }
  }

  // the peaks of the detected patterns and the peaks the detection worked on are available
  TEST_EQUAL(algorithm.getBlacklist().empty(), false)
  TEST_EQUAL(algorithm.getCentroidedExperiment().empty(), false)
}
END_SECTION

// A triple SILAC run (centroided) with one peptide whose light form is absent: medium and heavy only.
// With knock-outs the two have to be found as one doublet, i.e. the doublet pattern must not
// require a peak at the absent light position.
START_SECTION((knock-out detection of a multiplet without its lightest channel))
{
  const double mz_medium = 600.0, charge = 2.0;
  const double lys4 = 4.0251069836, lys8 = 8.0141988132;
  const double mz_heavy = mz_medium + (lys8 - lys4) / charge;
  const std::vector<double> isotopes = {1.0, 0.6, 0.25, 0.08};

  MSExperiment exp;
  for (int i = 0; i < 60; ++i)
  {
    const double rt = 100.0 + i * 1.0;
    const double elution = std::exp(-0.5 * std::pow((rt - 130.0) / 5.0, 2));
    MSSpectrum spectrum;
    spectrum.setMSLevel(1);
    spectrum.setRT(rt);
    spectrum.setType(SpectrumSettings::SpectrumType::CENTROID);
    for (const double mz : {mz_medium, mz_heavy})
    {
      for (Size k = 0; k < isotopes.size(); ++k)
      {
        if (1e6 * elution * isotopes[k] > 100.0)
        {
          spectrum.push_back(Peak1D(mz + k * Constants::C13C12_MASSDIFF_U / charge, 1e6 * elution * isotopes[k]));
        }
      }
    }
    spectrum.sortByPosition();
    exp.addSpectrum(spectrum);
  }
  exp.updateRanges();

  FeatureFinderMultiplexAlgorithm algorithm;
  Param param = algorithm.getParameters();
  param.setValue("algorithm:labels", "[][Lys4,Arg6][Lys8,Arg10]");
  param.setValue("algorithm:charge", "2:2");
  param.setValue("algorithm:isotopes_per_peptide", "3:4");
  param.setValue("algorithm:knock_out", "true");
  algorithm.setParameters(param);
  MSExperiment input = exp;
  algorithm.run(input, false);

  const ConsensusMap& multiplets = algorithm.getConsensusMap();
  TEST_EQUAL(multiplets.size(), 1)
  if (multiplets.size() == 1)
  {
    TEST_EQUAL(multiplets[0].size(), 2)
    TOLERANCE_ABSOLUTE(0.01);
    TEST_REAL_SIMILAR(multiplets[0].getFeatures().begin()->getMZ(), mz_medium)
    TEST_REAL_SIMILAR(multiplets[0].getFeatures().rbegin()->getMZ(), mz_heavy)
  }
  // every configured channel has a column header with its own channel id, also the third one, which no multiplet uses (with
  // knock-outs the map index is the position in the detected pattern, so the doublet fills columns 0 and 1)
  TEST_EQUAL(multiplets.getColumnHeaders().size(), 3)
  for (const auto& [index, header] : multiplets.getColumnHeaders())
  {
    TEST_EQUAL(header.metaValueExists("channel_id"), true)
    if (header.metaValueExists("channel_id")) { TEST_EQUAL(static_cast<Size>(header.getMetaValue("channel_id")), index) }
  }
  TEST_EQUAL(multiplets.getColumnHeaders().at(2).size, 0)

  // centroided input has a blacklist, too (the resolver tells absent from not quantifiable channels with it)
  Size blacklisted = 0;
  for (const MSSpectrum& spectrum : algorithm.getBlacklist()) { blacklisted += spectrum.size(); }
  TEST_NOT_EQUAL(blacklisted, 0)
}
END_SECTION

// 'algorithm:rt_min' (2 s) applies to the RT range of each peptide's own satellites
START_SECTION(([EXTRA] peptide intensities: minimum elution time per peptide of a multiplet))
{
  FeatureFinderMultiplexAlgorithmTester algorithm;
  algorithm.setInput(doubletSpectra(), 3);
  const MultiplexIsotopicPeakPattern pattern = doubletPattern();

  // light elutes for 3 s, heavy is seen in one scan only: the heavy intensity is not reliable
  {
    std::multimap<size_t, MultiplexSatelliteCentroided> satellites;
    addSatellites(satellites, {0, 1, 2}, {0, 1, 2, 3});
    addSatellites(satellites, {3, 4, 5}, {1});
    const std::vector<double> intensities = algorithm.peptideIntensities(pattern, satellites);
    TEST_EQUAL(intensities.size(), 2)
    TEST_EQUAL(intensities[0] > 0.0, true)
    TEST_REAL_SIMILAR(intensities[1], -1.0)
  }

  // the light mono-isotopic trace is seen at 100, 101 and 103 s, the other light isotopes at 102 and 103 s: light elutes for
  // 3 s (the range of all its satellites, not only from the last mono-isotopic one on)
  {
    std::multimap<size_t, MultiplexSatelliteCentroided> satellites;
    addSatellites(satellites, {0}, {0, 1, 3});
    addSatellites(satellites, {1, 2}, {2, 3});
    addSatellites(satellites, {3, 4, 5}, {0, 1, 2, 3});
    const std::vector<double> intensities = algorithm.peptideIntensities(pattern, satellites);
    TEST_EQUAL(intensities.size(), 2)
    TEST_EQUAL(intensities[0] > 0.0, true)
    TEST_EQUAL(intensities[1] > 0.0, true)
  }
}
END_SECTION

START_SECTION(([EXTRA] multiplets from clusters: union of satellites and channel counts))
{
  FeatureFinderMultiplexAlgorithmTester algorithm;
  algorithm.setInput(doubletSpectra(), 3);
  const std::vector<MultiplexIsotopicPeakPattern> patterns = {doubletPattern()};
  const std::vector<size_t> all = {0, 1, 2, 3, 4, 5};

  MultiplexFilteredMSExperiment filtered;
  // cluster 0: two filtered peaks whose satellites overlap in scan 1 (as with neighbouring RT bands); the duplicate comes
  // first among the satellites of the second peak, which still adds scans 2 and 3, so the doublet elutes for 3 s
  MultiplexFilteredPeak first(500.0, 100.0f, 0, 0);
  addSatellites(first, all, {0, 1});
  filtered.addPeak(first);
  MultiplexFilteredPeak second(500.0, 102.0f, 0, 2);
  second.addSatellite(1, 0, 0);
  addSatellites(second, all, {2, 3});
  filtered.addPeak(second);
  // cluster 1: the heavy mono-isotopic trace is seen in a single scan, so the doublet is not reported (after its light
  // channel was built already)
  MultiplexFilteredPeak short_heavy(500.0, 100.0f, 0, 0);
  addSatellites(short_heavy, {0, 1, 2, 4, 5}, {0, 1, 2, 3});
  addSatellites(short_heavy, {3}, {1});
  filtered.addPeak(short_heavy);

  std::vector<std::map<int, GridBasedCluster>> clusters(1);
  GridBasedCluster::Rectangle box;
  box.enlarge(GridBasedCluster::Point(500.0, 100.0));
  box.enlarge(GridBasedCluster::Point(500.0, 103.0));
  clusters[0].insert(std::make_pair(0, GridBasedCluster(GridBasedCluster::Point(500.0, 101.5), box, std::vector<int>{0, 1})));
  clusters[0].insert(std::make_pair(1, GridBasedCluster(GridBasedCluster::Point(500.0, 101.5), box, std::vector<int>{2})));

  algorithm.generateMaps(patterns, {filtered}, clusters);
  const ConsensusMap& multiplets = algorithm.getConsensusMap();
  TEST_EQUAL(multiplets.size(), 1)
  if (multiplets.size() == 1)
  {
    TEST_EQUAL(multiplets[0].size(), 2)
    TEST_REAL_SIMILAR(multiplets[0].getFeatures().begin()->getMZ(), 500.0)
  }
  // only the reported doublet counts in the column headers
  TEST_EQUAL(multiplets.getColumnHeaders().at(0).size, 1)
  TEST_EQUAL(multiplets.getColumnHeaders().at(1).size, 1)
  TEST_EQUAL(algorithm.getFeatureMap().size(), 2)
}
END_SECTION

END_TEST
