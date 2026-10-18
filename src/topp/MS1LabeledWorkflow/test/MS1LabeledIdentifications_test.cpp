// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include "MS1LabeledIdentifications.h"

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FEATUREFINDER/MultiplexResolverAlgorithm.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>

#include <map>
#include <optional>
#include <string>
#include <vector>

using namespace OpenMS;

namespace
{
// SILAC duplex ([][Lys8,Arg10], the resolver's default labels): channel 0 = light, channel 1 = heavy
const double LYS8 = 8.0141988132;

PeptideIdentification makeId(const std::string& sequence, const std::string& spectrum, double rt, std::optional<Size> map_index,
                             double score = 0.01, int charge = 2)
{
  const AASequence peptide = AASequence::fromString(sequence);
  PeptideIdentification id;
  id.setIdentifier("run0");
  id.setSpectrumReference(spectrum);
  id.setRT(rt);
  id.setMZ(peptide.getMZ(charge));
  id.setScoreType("q-value");
  id.setHigherScoreBetter(false);
  PeptideHit hit;
  hit.setSequence(peptide);
  hit.setCharge(charge);
  hit.setScore(score);
  id.insertHit(hit);
  if (map_index) { id.setMetaValue("map_index", *map_index); }
  return id;
}

FeatureHandle makeHandle(Size map_index, double mz, double rt)
{
  FeatureHandle handle;
  handle.setUniqueId();
  handle.setMapIndex(map_index);
  handle.setMZ(mz);
  handle.setRT(rt);
  handle.setCharge(2);
  handle.setIntensity(100.0);
  return handle;
}

/// A multiplet at @p rt of a peptide of charge 2 whose light channel is at @p light_mz: its channels @p channels
/// (0 = light, 1 = heavy), numbered from 0 like a detected (partial) pattern
ConsensusFeature makeMultiplet(double light_mz, double rt, const std::vector<Size>& channels)
{
  ConsensusFeature multiplet;
  multiplet.setUniqueId();
  multiplet.setRT(rt);
  multiplet.setCharge(2);
  Size index = 0;
  for (const Size channel : channels) { multiplet.insert(makeHandle(index++, light_mz + channel * LYS8 / 2, rt)); }
  multiplet.setMZ(multiplet.getFeatures().begin()->getMZ());
  return multiplet;
}

/// A duplex map of @p multiplets, which link the identifications @p linked (by multiplet), with the unassigned ones @p unassigned
ConsensusMap makeMap(const std::vector<ConsensusFeature>& multiplets, const std::map<Size, PeptideIdentificationList>& linked,
                     const PeptideIdentificationList& unassigned = {})
{
  ConsensusMap map;
  map.setExperimentType("labeled_MS1");
  for (Size channel = 0; channel < 2; ++channel)
  {
    auto& header = map.getColumnHeaders()[channel];
    header.filename = "run0.mzML";
    header.label = channel == 0 ? "no_label" : "Lys8Arg10";
  }
  for (const auto& multiplet : multiplets) { map.push_back(multiplet); }
  MapPeptideIdentifications ids;
  ProteinIdentification protein;
  protein.setIdentifier("run0");
  protein.setPrimaryMSRunPath({"run0.mzML"});
  ids.proteins = {protein};
  for (const auto& [m, peptides] : linked) { ids.features[{m}] = peptides; }
  ids.unassigned = unassigned;
  IdentificationDataConverter::importMapIdentifications(map, ids);
  return map;
}

std::string join(const std::vector<IdentificationData::QueryMatches>& entries)
{
  std::string spectra;
  for (const auto& entry : entries) { spectra += (spectra.empty() ? "" : ",") + entry.query->data_id; }
  return spectra;
}

/// The spectra of the identifications that multiplet @p m links, in order
std::string linked(const ConsensusMap& map, Size m)
{
  return join(map[m].getLinkedIdentifications(map.getIdentificationData()));
}

/// The spectra of the unassigned identifications, in order
std::string unassigned(const ConsensusMap& map)
{
  return join(map.getUnassignedIdentifications());
}

/// The channel (meta value 'map_index') of the identification of spectrum @p spectrum that multiplet @p m links, -1 without
int channelOf(const ConsensusMap& map, Size m, const std::string& spectrum)
{
  for (const auto& entry : map[m].getLinkedIdentifications(map.getIdentificationData()))
  {
    if (entry.query->data_id == spectrum) { return entry.query->metaValueExists("map_index") ? int(entry.query->getMetaValue("map_index")) : -1; }
  }
  return -2;
}

const double light_mz = AASequence::fromString("PEPTIDEK").getMZ(2);
} // namespace

START_TEST(MS1LabeledIdentifications, "$Id$")

START_SECTION((static void moveIdentifications(ConsensusMap&, const std::vector<Move>&)))
{
  ConsensusMap map = makeMap({makeMultiplet(light_mz, 100.0, {0, 1}), makeMultiplet(light_mz + 1.0, 200.0, {0, 1})},
                             {{0, {makeId("PEPTIDEK", "scan=1", 100.0, 0, 0.02), makeId("PEPTIDEK", "scan=2", 100.0, 0)}},
                              {1, {makeId("PEPTIDEK", "scan=3", 200.0, 1)}}},
                             {makeId("PEPTIDEK", "scan=4", 300.0, std::nullopt)});
  const auto entries = map[0].getLinkedIdentifications(map.getIdentificationData());
  ABORT_IF(entries.size() != 2)
  const IdentificationData::QueryReference first {entries[0].run->getUuid(), entries[0].query->getId()};
  const IdentificationData::QueryReference second {entries[1].run->getUuid(), entries[1].query->getId()};

  // nothing to move, nothing changes
  ConsensusMap unchanged = map;
  MS1LabeledIdentifications::moveIdentifications(map, {});
  TEST_TRUE(map == unchanged)

  // an unknown identification, one moved twice, a multiplet that does not exist: nothing changes
  TEST_EXCEPTION(Exception::InvalidValue, MS1LabeledIdentifications::moveIdentifications(map, {{{"unknown", first.query}, 1, 1}}))
  TEST_EXCEPTION(Exception::InvalidValue, MS1LabeledIdentifications::moveIdentifications(map, {{first, 1, 1}, {first, std::nullopt, std::nullopt}}))
  TEST_EXCEPTION(Exception::InvalidValue, MS1LabeledIdentifications::moveIdentifications(map, {{second, 1, 1}, {first, 2, 0}}))
  TEST_TRUE(map == unchanged)

  // the moved identifications come last, with their matches and scores, on the multiplet with the channel given or unassigned
  // without channel
  MS1LabeledIdentifications::moveIdentifications(map, {{second, 1, 1}, {first, std::nullopt, std::nullopt}});
  TEST_EQUAL(linked(map, 0), "")
  TEST_EQUAL(linked(map, 1), "scan=3,scan=2")
  TEST_EQUAL(channelOf(map, 1, "scan=2"), 1)
  TEST_EQUAL(unassigned(map), "scan=4,scan=1")
  const auto moved = map.getUnassignedIdentifications();
  ABORT_IF(moved.size() != 2 || moved[1].matches.size() != 1)
  TEST_FALSE(moved[1].query->metaValueExists("map_index"))
  TEST_EQUAL(moved[1].matches[0]->representation, "PEPTIDEK")
  TEST_REAL_SIMILAR(*moved[1].run->getScore(moved[1].matches[0]->getId(), *moved[1].run->getPrimaryScore()), 0.02)
  TEST_EQUAL(map.getIdentificationData().getRuns().front().getNumberOfIdentifications(), 4)
}
END_SECTION

START_SECTION((static Size keepConsistentIdentifications(ConsensusMap&, const MultiplexResolverAlgorithm&)))
{
  MultiplexResolverAlgorithm resolver;
  const std::string light = "PEPTIDEK", heavy = "PEPTIDEK(Label:13C(6)15N(2))", arg10 = "PEPTIDER(Label:13C(6)15N(4))";
  // 0: a Lys8 doublet with the light peptide on the light channel and the heavy one on the heavy channel (both fit), and a
  //    better-scoring Arg10-labelled peptide on the light channel (its heavy partner would sit at the Arg10 shift)
  // 1: a Lys8 doublet that none of its identifications fits
  // 2: a Lys8 doublet with a single identification, which does not fit
  ConsensusMap map = makeMap({makeMultiplet(light_mz, 100.0, {0, 1}), makeMultiplet(light_mz + 1.0, 200.0, {0, 1}),
                              makeMultiplet(light_mz + 2.0, 300.0, {0, 1})},
                             {{0, {makeId(light, "scan=1", 100.0, 0), makeId(arg10, "scan=2", 100.0, 0, 0.001), makeId(heavy, "scan=3", 100.0, 1)}},
                              {1, {makeId(arg10, "scan=4", 200.0, 0), makeId(arg10, "scan=5", 200.0, 1)}},
                              {2, {makeId(arg10, "scan=6", 300.0, 0)}}},
                             {makeId(light, "scan=7", 400.0, std::nullopt)});
  TEST_EQUAL(MS1LabeledIdentifications::keepConsistentIdentifications(map, resolver), 1)
  TEST_EQUAL(linked(map, 0), "scan=1,scan=3")
  TEST_EQUAL(linked(map, 1), "scan=4,scan=5")
  TEST_EQUAL(linked(map, 2), "scan=6")
  // the contradicting identification stays a PSM, unassigned and without channel, after the other unassigned ones
  TEST_EQUAL(unassigned(map), "scan=7,scan=2")
  const auto moved = map.getUnassignedIdentifications();
  ABORT_IF(moved.size() != 2)
  TEST_FALSE(moved[1].query->metaValueExists("map_index"))

  // nothing left to prefer
  const ConsensusMap resolved = map;
  TEST_EQUAL(MS1LabeledIdentifications::keepConsistentIdentifications(map, resolver), 0)
  TEST_TRUE(map == resolved)
}
END_SECTION

START_SECTION((static Size keepMultiAssignedOnClosest(ConsensusMap&)))
{
  const std::string heavy = "PEPTIDEK(Label:13C(6)15N(2))";
  // the same spectrum match (scan=1) was mapped onto two doublets of PEPTIDEK, the second one 20 ppm off; distinct spectra
  // of the peptide (scan=2, scan=3) sit on one doublet each; a heavy spectrum (scan=4) whose precursor is the light m/z (an
  // isotope off) was mapped with the light channel
  PeptideIdentification precursor_off = makeId(heavy, "scan=4", 100.0, 0);
  precursor_off.setMZ(light_mz);
  const double off = light_mz * 20e-6;
  ConsensusMap map = makeMap({makeMultiplet(light_mz, 100.0, {0, 1}), makeMultiplet(light_mz + off, 104.0, {0, 1})},
                             {{0, {makeId("PEPTIDEK", "scan=1", 101.0, 1), makeId("PEPTIDEK", "scan=2", 100.0, 1), precursor_off}},
                              {1, {makeId("PEPTIDEK", "scan=1", 101.0, 1), makeId("PEPTIDEK", "scan=3", 104.0, 0)}}});
  TEST_EQUAL(MS1LabeledIdentifications::keepMultiAssignedOnClosest(map), 1)
  TEST_EQUAL(linked(map, 0), "scan=1,scan=2,scan=4")
  TEST_EQUAL(linked(map, 1), "scan=3")
  // the channel is the one closest to the identified peptidoform
  TEST_EQUAL(channelOf(map, 0, "scan=1"), 0)
  TEST_EQUAL(channelOf(map, 0, "scan=2"), 0)
  TEST_EQUAL(channelOf(map, 0, "scan=4"), 1)
  TEST_EQUAL(channelOf(map, 1, "scan=3"), 0)
  TEST_EQUAL(unassigned(map), "")
}
END_SECTION

START_SECTION((static Size attachToUndetectedChannels(ConsensusMap&, const FeatureMap&, const MultiplexResolverAlgorithm&, Size, const Tolerances&)))
{
  MultiplexResolverAlgorithm resolver;
  const double heavy_mz = light_mz + LYS8 / 2;
  const double other_mz = AASequence::fromString("SAMPLEK").getMZ(2);
  // 0: a partial multiplet of PEPTIDEK: only its heavy channel was detected, with its heavy spectrum
  // 1: a complete doublet of SAMPLEK
  ConsensusFeature partial = makeMultiplet(heavy_mz, 100.0, {0});
  ConsensusMap map = makeMap({partial, makeMultiplet(other_mz, 100.0, {0, 1})},
                             {{0, {makeId("PEPTIDEK(Label:13C(6)15N(2))", "scan=1", 100.0, 0)}}},
                             {makeId("PEPTIDEK", "scan=2", 101.0, std::nullopt),               // light spectrum of the partial multiplet
                              makeId("PEPTIDEK", "scan=3", 150.0, std::nullopt),               // elutes elsewhere
                              makeId("PEPTIDEK", "scan=4", 101.0, std::nullopt, 0.01, 3),      // another charge
                              makeId("SAMPLEK(Label:13C(6)15N(2))", "scan=5", 100.0, std::nullopt), // its channel was detected
                              makeId("PEPTIDEK", "scan=6", 102.0, std::nullopt)});             // another light spectrum
  // the channels with their elution regions (95 - 105 s)
  FeatureMap channels;
  for (const auto& multiplet : map)
  {
    for (const auto& handle : multiplet.getFeatures())
    {
      Feature channel;
      channel.setUniqueId(handle.getUniqueId());
      channel.setMZ(handle.getMZ());
      channel.setRT(handle.getRT());
      channel.setCharge(handle.getCharge());
      ConvexHull2D hull;
      hull.addPoint(DPosition<2>(95.0, handle.getMZ() - 0.001));
      hull.addPoint(DPosition<2>(105.0, handle.getMZ() + 0.001));
      channel.getConvexHulls().push_back(hull);
      channels.push_back(channel);
    }
  }
  const MS1LabeledIdentifications::Tolerances tolerances {5.0, 10.0, true};

  TEST_EQUAL(MS1LabeledIdentifications::attachToUndetectedChannels(map, channels, resolver, 2, tolerances), 2)
  // the light spectra join the partial multiplet after its own identification, with the light channel
  TEST_EQUAL(linked(map, 0), "scan=1,scan=2,scan=6")
  TEST_EQUAL(channelOf(map, 0, "scan=2"), 0)
  TEST_EQUAL(unassigned(map), "scan=3,scan=4,scan=5")
  TEST_EQUAL(linked(map, 1), "")
  // the light channel joins the multiplet once, as a feature of intensity 0 without unique id, and the channels are numbered
  // in m/z order
  ABORT_IF(map[0].size() != 2)
  const FeatureHandle& added = *map[0].getFeatures().begin();
  const FeatureHandle& detected = *map[0].getFeatures().rbegin();
  TEST_EQUAL(added.getMapIndex(), 0)
  TEST_REAL_SIMILAR(added.getMZ(), light_mz)
  TEST_REAL_SIMILAR(added.getIntensity(), 0.0)
  TEST_FALSE(added.hasValidUniqueId())
  TEST_EQUAL(detected.getMapIndex(), 1)
  TEST_EQUAL(detected.getUniqueId(), partial.getFeatures().begin()->getUniqueId())
  TEST_REAL_SIMILAR(map[0].getMZ(), light_mz)
  TEST_EQUAL(map[1].size(), 2)

  // a multiplet with as many detected channels as the labels have is complete: nothing is attached
  ConsensusMap complete = makeMap({partial}, {}, {makeId("PEPTIDEK", "scan=2", 101.0, std::nullopt)});
  TEST_EQUAL(MS1LabeledIdentifications::attachToUndetectedChannels(complete, channels, resolver, 1, tolerances), 0)
  TEST_EQUAL(unassigned(complete), "scan=2")
}
END_SECTION

END_TEST
