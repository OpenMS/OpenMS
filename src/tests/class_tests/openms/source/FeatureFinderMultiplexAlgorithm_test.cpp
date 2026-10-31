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
#include <OpenMS/FEATUREFINDER/MultiplexResolverAlgorithm.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
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
                    std::vector<std::map<int, GridBasedCluster>>& cluster_results, std::vector<std::map<int, std::vector<UInt64>>>* multiplet_ids = nullptr)
  {
    generateMapsCentroided_(patterns, filter_results, cluster_results, multiplet_ids);
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

  // every blacklisted peak records the feature it was reserved for (profile input), a reported feature or none
  TEST_TRUE(MultiplexResolverAlgorithm::hasBlacklistOwners(algorithm.getBlacklist()))
  Size owned = 0, unknown = 0;
  for (const MSSpectrum& spectrum : algorithm.getBlacklist())
  {
    for (Size k = 0; k < spectrum.size(); ++k)
    {
      const UInt64 owner = MultiplexResolverAlgorithm::blacklistOwner(spectrum, k);
      if (owner == UniqueIdInterface::INVALID) { continue; }
      ++owned;
      unknown += (feature_of.count(owner) == 0);
    }
  }
  TEST_NOT_EQUAL(owned, 0)
  TEST_EQUAL(unknown, 0)
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

  // cluster 2: the light channel has no mono-isotopic satellites, so it gets no feature: the multiplet is reported with its
  // heavy channel only, which is not the first feature built for it by position
  MultiplexFilteredPeak no_light_mono(500.0, 100.0f, 0, 0);
  addSatellites(no_light_mono, {1, 2, 3, 4, 5}, {0, 1, 2, 3});
  filtered.addPeak(no_light_mono);
  clusters[0].insert(std::make_pair(2, GridBasedCluster(GridBasedCluster::Point(500.0, 101.5), box, std::vector<int>{3})));

  std::vector<std::map<int, std::vector<UInt64>>> multiplet_ids;
  algorithm.generateMaps(patterns, {filtered}, clusters, &multiplet_ids);
  const ConsensusMap& multiplets = algorithm.getConsensusMap();
  TEST_EQUAL(multiplets.size(), 2)
  if (multiplets.size() == 2)
  {
    TEST_EQUAL(multiplets[0].size(), 2)
    TEST_REAL_SIMILAR(multiplets[0].getFeatures().begin()->getMZ(), 500.0)
    TEST_EQUAL(multiplets[1].size(), 1)
  }
  // only the reported multiplets count in the column headers
  TEST_EQUAL(multiplets.getColumnHeaders().at(0).size, 1)
  TEST_EQUAL(multiplets.getColumnHeaders().at(1).size, 2)
  TEST_EQUAL(algorithm.getFeatureMap().size(), 3)

  // the features of every reported multiplet per (pattern, cluster, channel), by provisional id; the aborted multiplet
  // (cluster 1) is not recorded, although its light feature had an id (which the next multiplet reuses)
  TEST_EQUAL(multiplet_ids.size(), 1)
  ABORT_IF(multiplet_ids.size() != 1)
  TEST_EQUAL(multiplet_ids[0].size(), 2)
  TEST_EQUAL(multiplet_ids[0].count(1), 0)
  ABORT_IF(multiplet_ids[0].count(0) != 1 || multiplet_ids[0].count(2) != 1)
  const FeatureMap& features = algorithm.getFeatureMap();
  TEST_EQUAL(multiplet_ids[0][0].size(), 2)
  TEST_EQUAL(multiplet_ids[0][0][0], features[0].getUniqueId())
  TEST_EQUAL(multiplet_ids[0][0][1], features[1].getUniqueId())
  TEST_EQUAL(multiplet_ids[0][2].size(), 2)
  TEST_EQUAL(multiplet_ids[0][2][0], 0)
  TEST_EQUAL(multiplet_ids[0][2][1], features[2].getUniqueId())
  TEST_REAL_SIMILAR(features[2].getMZ(), 500.0 + 10.008268600 / 2)
}
END_SECTION

// Every reserved peak records the unique id of the feature it was reserved for: the feature of its channel in the reported
// multiplet, or none if the detection dropped the multiplet after the peak was reserved (here: elutes shorter than rt_min).
START_SECTION(([EXTRA] blacklist owners))
{
  const double lys8 = 8.0141988132, charge = 2.0;
  const double mz_light = 600.0, mz_heavy = mz_light + lys8 / charge, mz_short = 700.0;
  CoarseIsotopePatternGenerator generator(4);
  auto isotopes = [&](double mz)
  {
    std::vector<double> intensities;
    for (const Peak1D& peak : generator.estimateFromPeptideWeight(mz * charge)) { intensities.push_back(peak.getIntensity()); }
    return intensities;
  };
  const std::vector<double> doublet_isotopes = isotopes(mz_light), short_isotopes = isotopes(mz_short);

  MSExperiment exp;
  for (int i = 0; i < 60; ++i)
  {
    const double rt = 100.0 + i;
    const double elution = std::exp(-0.5 * std::pow((rt - 130.0) / 5.0, 2));
    MSSpectrum spectrum;
    spectrum.setMSLevel(1);
    spectrum.setRT(rt);
    spectrum.setType(SpectrumSettings::SpectrumType::CENTROID);
    for (const double mz : {mz_light, mz_heavy})
    {
      for (Size k = 0; k < doublet_isotopes.size(); ++k)
      {
        if (1e6 * elution * doublet_isotopes[k] > 2000.0) { spectrum.push_back(Peak1D(mz + k * Constants::C13C12_MASSDIFF_U / charge, 1e6 * elution * doublet_isotopes[k])); }
      }
    }
    // a peptide seen in two scans only: it passes the filters there, but its multiplet is dropped for eluting less than rt_min
    if (rt == 150.0 || rt == 151.0)
    {
      for (Size k = 0; k < short_isotopes.size(); ++k) { spectrum.push_back(Peak1D(mz_short + k * Constants::C13C12_MASSDIFF_U / charge, 1e6 * short_isotopes[k])); }
    }
    spectrum.sortByPosition();
    exp.addSpectrum(spectrum);
  }
  exp.updateRanges();

  FeatureFinderMultiplexAlgorithm algorithm;
  Param param = algorithm.getParameters();
  param.setValue("algorithm:labels", "[][Lys8]");
  param.setValue("algorithm:charge", "2:2");
  param.setValue("algorithm:isotopes_per_peptide", "3:4");
  param.setValue("algorithm:knock_out", "true");
  algorithm.setParameters(param);
  MSExperiment input = exp;
  algorithm.run(input, false);

  const ConsensusMap& multiplets = algorithm.getConsensusMap();
  TEST_EQUAL(multiplets.size(), 1)
  ABORT_IF(multiplets.size() != 1)
  TEST_EQUAL(multiplets[0].size(), 2)
  ABORT_IF(multiplets[0].size() != 2)
  const UInt64 light_id = multiplets[0].getFeatures().begin()->getUniqueId();
  const UInt64 heavy_id = multiplets[0].getFeatures().rbegin()->getUniqueId();
  TEST_NOT_EQUAL(light_id, heavy_id)

  const MSExperiment& blacklist = algorithm.getBlacklist();
  TEST_EQUAL(blacklist.size(), exp.size())
  TEST_TRUE(MultiplexResolverAlgorithm::hasBlacklistOwners(blacklist))
  Size light = 0, heavy = 0, orphans = 0, wrong = 0;
  for (const MSSpectrum& spectrum : blacklist)
  {
    for (Size k = 0; k < spectrum.size(); ++k)
    {
      const double mz = spectrum[k].getMZ();
      const UInt64 owner = MultiplexResolverAlgorithm::blacklistOwner(spectrum, k);
      if (mz > mz_short - 0.1) { orphans += (owner == UniqueIdInterface::INVALID); wrong += (owner != UniqueIdInterface::INVALID); }
      else if (mz > mz_heavy - 0.1) { heavy += (owner == heavy_id); wrong += (owner != heavy_id); }
      else { light += (owner == light_id); wrong += (owner != light_id); }
    }
  }
  TEST_NOT_EQUAL(light, 0)
  TEST_NOT_EQUAL(heavy, 0)
  TEST_NOT_EQUAL(orphans, 0)
  TEST_EQUAL(wrong, 0)
}
END_SECTION

// With an RT band, a pattern reserves the peaks of its mass traces over their RT range extended by twice the band, also peaks
// that an earlier pattern reserved. Here the light channel of a doublet outlasts the heavy one: the doublet reserves the light
// tail up to 2 * rt_band after its last light satellite. A singlet then finds the two spectra of the tail beyond that and is
// dropped (shorter than rt_min), but its extended traces reach back over the doublet's reserved tail. Those peaks stay the
// doublet's; only the peaks the singlet alone reserved belong to no reported feature.
START_SECTION(([EXTRA] blacklist owners: a dropped pattern does not take over the peaks of a reported multiplet))
{
  const double lys8 = 8.0141988132, charge = 2.0, rt_band = 5.0;
  const double mz_light = 600.0, mz_heavy = mz_light + lys8 / charge;
  const double heavy_end = 125.0, light_end = 141.0;
  CoarseIsotopePatternGenerator generator(4);
  std::vector<double> isotopes;
  for (const Peak1D& peak : generator.estimateFromPeptideWeight(mz_light * charge)) { isotopes.push_back(peak.getIntensity()); }

  MSExperiment exp;
  for (int i = 0; i < 70; ++i)
  {
    const double rt = 90.0 + i;
    const double elution = std::exp(-0.5 * std::pow((rt - 115.0) / 4.0, 2));
    MSSpectrum spectrum;
    spectrum.setMSLevel(1);
    spectrum.setRT(rt);
    spectrum.setType(SpectrumSettings::SpectrumType::CENTROID);
    for (Size k = 0; k < isotopes.size(); ++k)
    {
      const double light = rt <= light_end ? 1e6 * std::max(elution, rt > 115.0 ? 0.3 : 0.0) * isotopes[k] : 0.0;
      const double heavy = rt <= heavy_end ? 1e6 * elution * isotopes[k] : 0.0;
      if (light > 2000.0) { spectrum.push_back(Peak1D(mz_light + k * Constants::C13C12_MASSDIFF_U / charge, light)); }
      if (heavy > 2000.0) { spectrum.push_back(Peak1D(mz_heavy + k * Constants::C13C12_MASSDIFF_U / charge, heavy)); }
    }
    spectrum.sortByPosition();
    exp.addSpectrum(spectrum);
  }
  exp.updateRanges();

  FeatureFinderMultiplexAlgorithm algorithm;
  Param param = algorithm.getParameters();
  param.setValue("algorithm:labels", "[][Lys8]");
  param.setValue("algorithm:charge", "2:2");
  param.setValue("algorithm:isotopes_per_peptide", "3:4");
  param.setValue("algorithm:rt_band", rt_band);
  param.setValue("algorithm:knock_out", "true");
  algorithm.setParameters(param);
  MSExperiment input = exp;
  algorithm.run(input, false);

  // the doublet only (the singlet on the tail is dropped)
  const ConsensusMap& multiplets = algorithm.getConsensusMap();
  TEST_EQUAL(multiplets.size(), 1)
  ABORT_IF(multiplets.size() != 1)
  TEST_EQUAL(multiplets[0].size(), 2)
  ABORT_IF(multiplets[0].size() != 2)
  const UInt64 light_id = multiplets[0].getFeatures().begin()->getUniqueId();
  double last_satellite = 0.0; // of the light channel
  for (const Feature& feature : algorithm.getFeatureMap())
  {
    if (feature.getUniqueId() == light_id) { last_satellite = feature.getConvexHull().getBoundingBox().maxX(); }
  }
  TEST_EQUAL(last_satellite > heavy_end && last_satellite + 2 * rt_band < light_end, true)

  // light peaks up to 2 * rt_band after the last satellite are the doublet's, the rest (the singlet's own) are orphans
  Size doublet = 0, orphans = 0, wrong = 0;
  const MSExperiment& blacklist = algorithm.getBlacklist();
  for (const MSSpectrum& spectrum : blacklist)
  {
    for (Size k = 0; k < spectrum.size(); ++k)
    {
      if (spectrum[k].getMZ() > mz_heavy - 0.1) { continue; }
      const UInt64 owner = MultiplexResolverAlgorithm::blacklistOwner(spectrum, k);
      if (spectrum.getRT() <= last_satellite + 2 * rt_band) { doublet += (owner == light_id); wrong += (owner != light_id); }
      else { orphans += (owner == UniqueIdInterface::INVALID); wrong += (owner != UniqueIdInterface::INVALID); }
    }
  }
  TEST_NOT_EQUAL(doublet, 0)
  TEST_NOT_EQUAL(orphans, 0)
  TEST_EQUAL(wrong, 0)
}
END_SECTION

END_TEST
