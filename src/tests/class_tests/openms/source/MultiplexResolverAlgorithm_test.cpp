// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/ClassTest.h>
#include <set>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FEATUREFINDER/MultiplexResolverAlgorithm.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/NativeIdentificationTest.h>
#include <OpenMS/test_config.h>
#include <algorithm>
#include <cmath>

using namespace OpenMS;

// A SILAC duplex ([][Lys8,Arg10]) test map, built the way FeatureFinderMultiplexAlgorithm + IDMapper
// (annotate_ids_with_subelements) would leave it: column 0 = light, column 1 = heavy, and the
// identification of a feature carries the map index of the handle it was mapped to.
const double LYS8 = 8.0141988132;
const double ARG10 = 10.008268600;
const double LYS4 = 4.0251069836;
const double ARG6 = 6.0201290268;
const double DIMETHYL0 = 28.031300;
const double DIMETHYL6 = 34.063117;

FeatureHandle makeHandle(Size map_index, double rt, double mz, int charge, double intensity)
{
  FeatureHandle h;
  h.setMapIndex(map_index);
  h.setRT(rt);
  h.setMZ(mz);
  h.setCharge(charge);
  h.setIntensity(intensity);
  return h;
}

PeptideIdentification makeId(const std::string& sequence, int charge, Size map_index)
{
  PeptideHit hit;
  hit.setSequence(AASequence::fromString(sequence));
  hit.setCharge(charge);
  hit.setScore(0.01);
  PeptideIdentification id;
  id.setScoreType("Posterior Error Probability");
  id.setHigherScoreBetter(false);
  id.insertHit(hit);
  id.setMetaValue("map_index", map_index);
  return id;
}

/// Give the consensus features of @p map the identifications @p ids (index of the feature, identification) as its
/// native identification data
void identify(ConsensusMap& map, const std::vector<std::pair<Size, PeptideIdentification>>& ids)
{
  MapPeptideIdentifications mapped;
  for (const auto& [index, id] : ids)
  {
    mapped.features[{index}].push_back(id);
  }
  Internal::ClassTest::addSearchRun(mapped);
  Internal::ClassTest::setIdentifications(map, mapped);
}

ConsensusMap makeDuplexMap()
{
  ConsensusMap map;
  map.setExperimentType("labeled_MS1");
  map.getColumnHeaders()[0].filename = "run.mzML";
  map.getColumnHeaders()[0].label = "no_label";
  map.getColumnHeaders()[0].setMetaValue("channel_id", 0);
  map.getColumnHeaders()[1].filename = "run.mzML";
  map.getColumnHeaders()[1].label = "Lys8Arg10";
  map.getColumnHeaders()[1].setMetaValue("channel_id", 1);
  MapPeptideIdentifications ids;

  // A: complete doublet (Lys8 shift), heavy peptide identified
  {
    ConsensusFeature cf;
    cf.setRT(100.0);
    cf.setMZ(500.0);
    cf.setCharge(2);
    cf.setIntensity(300.0);
    cf.insert(makeHandle(0, 100.0, 500.0, 2, 100.0));
    cf.insert(makeHandle(1, 100.0, 500.0 + LYS8 / 2, 2, 200.0));
    ids.features[{map.size()}].push_back(makeId("PEPTIDEK(Label:13C(6)15N(2))", 2, 1));
    map.push_back(cf);
  }
  // B: only the heavy partner was detected; nothing blacklisted around the light position
  {
    ConsensusFeature cf;
    cf.setRT(200.0);
    cf.setMZ(600.0);
    cf.setCharge(2);
    cf.setIntensity(50.0);
    cf.insert(makeHandle(0, 200.0, 600.0, 2, 50.0));
    ids.features[{map.size()}].push_back(makeId("AAAAK(Label:13C(6)15N(2))", 2, 0));
    map.push_back(cf);
  }
  // C: doublet with an Arg10 shift, but the sequence carries a Lys8 label -> conflict
  {
    ConsensusFeature cf;
    cf.setRT(300.0);
    cf.setMZ(700.0);
    cf.setCharge(2);
    cf.setIntensity(300.0);
    cf.insert(makeHandle(0, 300.0, 700.0, 2, 100.0));
    cf.insert(makeHandle(1, 300.0, 700.0 + ARG10 / 2, 2, 200.0));
    ids.features[{map.size()}].push_back(makeId("AAAAK(Label:13C(6)15N(2))", 2, 1));
    map.push_back(cf);
  }
  // D: no identification -> conflict output, unchanged
  {
    ConsensusFeature cf;
    cf.setRT(400.0);
    cf.setMZ(800.0);
    cf.setCharge(2);
    cf.setIntensity(300.0);
    cf.insert(makeHandle(0, 400.0, 800.0, 2, 100.0));
    cf.insert(makeHandle(1, 400.0, 800.0 + LYS8 / 2, 2, 200.0));
    map.push_back(cf);
  }
  // E: only the heavy partner was detected, and the light position is blacklisted
  {
    ConsensusFeature cf;
    cf.setRT(500.0);
    cf.setMZ(900.0);
    cf.setCharge(2);
    cf.setIntensity(50.0);
    cf.insert(makeHandle(0, 500.0, 900.0, 2, 50.0));
    ids.features[{map.size()}].push_back(makeId("AAAAR(Label:13C(6)15N(4))", 2, 0));
    map.push_back(cf);
  }
  map.applyMemberFunction(&UniqueIdInterface::setUniqueId);
  Internal::ClassTest::addSearchRun(ids);
  Internal::ClassTest::setIdentifications(map, ids);
  return map;
}

START_TEST(MultiplexResolverAlgorithm, "$Id$")

/////////////////////////////////////////////////////////////

MultiplexResolverAlgorithm* ptr = nullptr;
MultiplexResolverAlgorithm* null_ptr = nullptr;

START_SECTION(MultiplexResolverAlgorithm())
  ptr = new MultiplexResolverAlgorithm();
  TEST_NOT_EQUAL(ptr, null_ptr)
  TEST_EQUAL(ptr->getParameters().getValue("algorithm:labels").toString(), "[][Lys8,Arg10]")
  TEST_EQUAL((int)ptr->getParameters().getValue("algorithm:max_nr_labelled_aas"), 0)
  TEST_REAL_SIMILAR((double)ptr->getParameters().getValue("labels:Lys8"), LYS8)
  delete ptr;
END_SECTION

START_SECTION((void resolve(const ConsensusMap& map_in, ConsensusMap& map_out, ConsensusMap& map_conflicts, const MSExperiment& blacklist = MSExperiment(), bool keep_conflicting_identifications = false) const))
{
  MultiplexResolverAlgorithm resolver;

  // blacklist: one MS1 spectrum at the RT of feature E with a peak at its light position
  MSExperiment blacklist;
  MSSpectrum spec;
  spec.setRT(500.0);
  spec.setMSLevel(1);
  spec.push_back(Peak1D(900.0 - ARG10 / 2, 1000.0));
  blacklist.addSpectrum(spec);
  blacklist.updateRanges();

  ConsensusMap in = makeDuplexMap();
  ConsensusMap out, conflicts;
  resolver.resolve(in, out, conflicts, blacklist);

  // meta data is inherited by both outputs
  TEST_EQUAL(out.getExperimentType(), "labeled_MS1")
  TEST_EQUAL(out.getColumnHeaders().size(), 2)
  TEST_EQUAL(conflicts.getColumnHeaders().size(), 2)
  TEST_EQUAL(out.getColumnHeaders()[0].label, "no_label")
  TEST_EQUAL(out.getColumnHeaders()[1].label, "Lys8Arg10")

  TEST_EQUAL(out.size(), 3)
  TEST_EQUAL(conflicts.size(), 2)
  TEST_EQUAL(out.getColumnHeaders()[0].size, 3)
  TEST_EQUAL(out.getColumnHeaders()[1].size, 3)

  // A: unchanged
  TEST_EQUAL(out[0].size(), 2)
  TEST_REAL_SIMILAR(out[0].getMZ(), 500.0)
  TEST_EQUAL(Internal::ClassTest::peptidesOf(out, {0}).size(), 1)

  // B: completed with a zero-intensity light dummy; the detected feature moved to column 1
  {
    const ConsensusFeature& cf = out[1];
    TEST_EQUAL(cf.size(), 2)
    TEST_REAL_SIMILAR(cf.getMZ(), 600.0 - LYS8 / 2)
    TEST_EQUAL(cf.getCharge(), 2)
    auto it = cf.getFeatures().begin();
    TEST_EQUAL(it->getMapIndex(), 0)
    TEST_REAL_SIMILAR(it->getMZ(), 600.0 - LYS8 / 2)
    TEST_REAL_SIMILAR(it->getIntensity(), 0.0)
    TEST_EQUAL(it->getCharge(), 2)
    ++it;
    TEST_EQUAL(it->getMapIndex(), 1)
    TEST_REAL_SIMILAR(it->getMZ(), 600.0)
    TEST_REAL_SIMILAR(it->getIntensity(), 50.0)
    // the new map index of the identified feature is recorded on the hit
    TEST_EQUAL((int)Internal::ClassTest::peptidesOf(out, {1})[0].getHits()[0].getMetaValue("map_index"), 1)
  }

  // E: completed, but the light dummy is not quantifiable because its region was blacklisted
  {
    const ConsensusFeature& cf = out[2];
    TEST_EQUAL(cf.size(), 2)
    auto it = cf.getFeatures().begin();
    TEST_EQUAL(it->getMapIndex(), 0)
    TEST_REAL_SIMILAR(it->getMZ(), 900.0 - ARG10 / 2)
    TEST_TRUE(std::isnan(it->getIntensity()))
    ++it;
    TEST_EQUAL(it->getMapIndex(), 1)
    TEST_REAL_SIMILAR(it->getIntensity(), 50.0)
  }

  // C and D went to the conflicts, unchanged
  TEST_REAL_SIMILAR(conflicts[0].getMZ(), 700.0)
  TEST_EQUAL(Internal::ClassTest::peptidesOf(conflicts, {0}).size(), 1)
  TEST_REAL_SIMILAR(conflicts[1].getMZ(), 800.0)
  TEST_EQUAL(Internal::ClassTest::peptidesOf(conflicts, {1}).size(), 0)

  // without a blacklist, E's dummy is reported as absent (0) instead of not quantifiable (NaN)
  ConsensusMap out2, conflicts2;
  resolver.resolve(in, out2, conflicts2);
  TEST_EQUAL(out2.size(), 3)
  TEST_REAL_SIMILAR(out2[2].getFeatures().begin()->getIntensity(), 0.0)

  // with one allowed missed cleavage, a doubly labelled peptide (Lys8 + Arg10) resolves as well
  {
    Param p = resolver.getParameters();
    p.setValue("algorithm:max_nr_labelled_aas", 1);
    resolver.setParameters(p);

    ConsensusMap in2;
    in2.setExperimentType("labeled_MS1");
    in2.getColumnHeaders() = in.getColumnHeaders();
    ConsensusFeature cf;
    cf.setRT(100.0);
    cf.setMZ(500.0);
    cf.setCharge(2);
    cf.insert(makeHandle(0, 100.0, 500.0, 2, 100.0));
    cf.insert(makeHandle(1, 100.0, 500.0 + (LYS8 + ARG10) / 2, 2, 200.0));
    MapPeptideIdentifications ids2;
    ids2.features[{0}].push_back(makeId("PEPTIDEK(Label:13C(6)15N(2))AAR(Label:13C(6)15N(4))", 2, 1));
    in2.push_back(cf);
    Internal::ClassTest::addSearchRun(ids2);
    Internal::ClassTest::setIdentifications(in2, ids2);

    ConsensusMap out3, conflicts3;
    resolver.resolve(in2, out3, conflicts3);
    TEST_EQUAL(out3.size(), 1)
    TEST_EQUAL(conflicts3.size(), 0)

    // ... whereas it is a conflict when only one labelled amino acid per peptide is allowed
    p.setValue("algorithm:max_nr_labelled_aas", 0);
    resolver.setParameters(p);
    resolver.resolve(in2, out3, conflicts3);
    TEST_EQUAL(out3.size(), 0)
    TEST_EQUAL(conflicts3.size(), 1)
  }

  // a light peptide with K and R: its heavy channel carries Lys8 and Arg10
  {
    MultiplexResolverAlgorithm resolver_kr;
    Param p = resolver_kr.getParameters();
    p.setValue("algorithm:max_nr_labelled_aas", 1);
    resolver_kr.setParameters(p);

    ConsensusMap in_kr;
    in_kr.setExperimentType("labeled_MS1");
    in_kr.getColumnHeaders() = in.getColumnHeaders();
    // a light singlet is completed at the combined shift, not at the Lys8 shift alone
    ConsensusFeature singlet;
    singlet.setRT(100.0);
    singlet.setMZ(500.0);
    singlet.setCharge(2);
    singlet.insert(makeHandle(0, 100.0, 500.0, 2, 100.0));
    in_kr.push_back(singlet);
    // a doublet at the Lys8 shift cannot be this peptide
    ConsensusFeature lys8_doublet = singlet;
    lys8_doublet.setRT(200.0);
    lys8_doublet.insert(makeHandle(1, 200.0, 500.0 + LYS8 / 2, 2, 100.0));
    in_kr.push_back(lys8_doublet);
    identify(in_kr, {{0, makeId("PEPTIDEKAAR", 2, 0)}, {1, makeId("PEPTIDEKAAR", 2, 0)}});

    ConsensusMap out_kr, conflicts_kr;
    resolver_kr.resolve(in_kr, out_kr, conflicts_kr);
    TEST_EQUAL(out_kr.size(), 1)
    TEST_EQUAL(conflicts_kr.size(), 1)
    ABORT_IF(out_kr.size() != 1)
    TEST_EQUAL(out_kr[0].size(), 2)
    auto it = out_kr[0].getFeatures().begin();
    TEST_EQUAL(it->getMapIndex(), 0)
    TEST_REAL_SIMILAR(it->getIntensity(), 100.0)
    ++it;
    TEST_EQUAL(it->getMapIndex(), 1)
    TEST_REAL_SIMILAR(it->getMZ(), 500.0 + (LYS8 + ARG10) / 2)
    TEST_REAL_SIMILAR(it->getIntensity(), 0.0)
    TEST_REAL_SIMILAR(conflicts_kr[0].getRT(), 200.0)
  }

  // an identification without the 'map_index' annotation cannot be resolved
  {
    ConsensusMap in3 = makeDuplexMap();
    auto ids3 = Internal::ClassTest::identificationsOf(in3);
    ids3.features[{0}][0].removeMetaValue("map_index");
    Internal::ClassTest::setIdentifications(in3, ids3);
    ConsensusMap out4, conflicts4;
    TEST_EXCEPTION(Exception::MissingInformation, resolver.resolve(in3, out4, conflicts4))
  }
}
END_SECTION

START_SECTION(([EXTRA] resolve() on identification data))
{
  // the outputs get the identification data of the input
  ConsensusMap in = makeDuplexMap();
  auto ids = Internal::ClassTest::identificationsOf(in);
  ids.unassigned.push_back(makeId("AAAAK", 2, 0));
  Internal::ClassTest::addSearchRun(ids);
  Internal::ClassTest::setIdentifications(in, ids);
  MultiplexResolverAlgorithm resolver;
  ConsensusMap out, conflicts;
  resolver.resolve(in, out, conflicts);
  ABORT_IF(out.size() != 3 || conflicts.size() != 2)

  // each output keeps the identifications of its consensus features and the unassigned one of the input
  TEST_EQUAL(out.getIdentificationData().getRuns()[0].getNumberOfIdentifications(), 4)
  TEST_EQUAL(conflicts.getIdentificationData().getRuns()[0].getNumberOfIdentifications(), 2)
  TEST_EQUAL(out.getUnassignedIdentifications().size(), 1)
  TEST_EQUAL(conflicts.getUnassignedIdentifications().size(), 1)

  // the new map index of the identified feature of a completed multiplet is recorded on its match
  const auto completed = out[1].getLinkedIdentifications(out.getIdentificationData());
  ABORT_IF(completed.size() != 1 || completed[0].matches.empty())
  TEST_EQUAL((int)completed[0].matches[0]->getMetaValue("map_index"), 1)
  const auto unchanged = out[0].getLinkedIdentifications(out.getIdentificationData());
  ABORT_IF(unchanged.size() != 1 || unchanged[0].matches.empty())
  TEST_EQUAL(unchanged[0].matches[0]->metaValueExists("map_index"), false)
  // the input is not changed
  TEST_EQUAL(in[1].getLinkedIdentifications(in.getIdentificationData())[0].matches[0]->metaValueExists("map_index"), false)

  // with keep_conflicting_identifications, the resolved output keeps those of the conflicting multiplets as unassigned ones
  ConsensusMap kept_out, kept_conflicts;
  resolver.resolve(in, kept_out, kept_conflicts, MSExperiment(), true);
  TEST_EQUAL(kept_out.size(), 3)
  TEST_EQUAL(kept_out.getIdentificationData().getRuns()[0].getNumberOfIdentifications(), 5)
  TEST_EQUAL(kept_out.getUnassignedIdentifications().size(), 2)
  TEST_EQUAL(kept_conflicts.getIdentificationData().getRuns()[0].getNumberOfIdentifications(), 2)
  // the identifications of the conflicting multiplets are those that the conflict output links
  std::set<IdentificationData::QueryReference> conflicting;
  for (const auto& cf : kept_conflicts)
  {
    for (const auto& entry : cf.getLinkedIdentifications(kept_conflicts.getIdentificationData()))
    {
      conflicting.insert({entry.run->getUuid(), entry.query->getId()});
    }
  }
  TEST_EQUAL(conflicting.size(), 1)
  Size found = 0;
  for (const auto& entry : kept_out.getUnassignedIdentifications())
  {
    found += conflicting.count({entry.run->getUuid(), entry.query->getId()});
  }
  TEST_EQUAL(found, 1)
}
END_SECTION

START_SECTION(([EXTRA] resolve() keeps the identifications of conflicting multiplets on request (as peptide identifications)))
{
  using namespace Internal::ClassTest;
  ConsensusMap in = makeDuplexMap();
  MultiplexResolverAlgorithm resolver;
  ConsensusMap out, conflicts;
  resolver.resolve(in, out, conflicts);
  const Size unassigned = unassignedOf(out).size();
  Size conflicting = 0;
  for (Size i = 0; i < conflicts.size(); ++i) { conflicting += peptidesOf(conflicts, {i}).size(); }
  TEST_EQUAL(conflicting > 0, true)
  resolver.resolve(in, out, conflicts, MSExperiment(), true);
  TEST_EQUAL(unassignedOf(out).size(), unassigned + conflicting)
  TEST_EQUAL(out.size(), 3)
}
END_SECTION

START_SECTION((std::vector<bool> checkIdentifications(const ConsensusFeature& consensus, const IdentificationData& data) const))
{
  MultiplexResolverAlgorithm resolver;
  // a Lys8 doublet carrying: the light peptide on the light feature (consistent), a better-scoring
  // Arg10-labelled peptide on the light feature (its heavy partner would sit at the Arg10 shift),
  // the heavy peptide on the heavy feature (consistent), an identification naming no handle,
  // one without map index and one without matches
  ConsensusFeature cf;
  cf.setRT(100.0);
  cf.setMZ(500.0);
  cf.setCharge(2);
  cf.insert(makeHandle(0, 100.0, 500.0, 2, 100.0));
  cf.insert(makeHandle(1, 100.0, 500.0 + LYS8 / 2, 2, 200.0));
  // a duplex map with the multiplet, which links @p peptides
  const auto map_with = [&cf](const PeptideIdentificationList& peptides)
  {
    ConsensusMap map;
    map.setExperimentType("labeled_MS1");
    map.getColumnHeaders() = makeDuplexMap().getColumnHeaders();
    map.push_back(cf);
    map.applyMemberFunction(&UniqueIdInterface::setUniqueId);
    MapPeptideIdentifications ids;
    ids.features[{0}] = peptides;
    Internal::ClassTest::addSearchRun(ids);
    Internal::ClassTest::setIdentifications(map, ids);
    return map;
  };
  PeptideIdentificationList peptides;
  peptides.push_back(makeId("PEPTIDEK", 2, 0));
  PeptideIdentification wrong = makeId("PEPTIDER(Label:13C(6)15N(4))", 2, 0);
  wrong.getHits()[0].setScore(0.001);
  peptides.push_back(wrong);
  peptides.push_back(makeId("PEPTIDEK(Label:13C(6)15N(2))", 2, 1));
  peptides.push_back(makeId("PEPTIDEK", 2, 5));
  PeptideIdentification no_index = makeId("PEPTIDEK", 2, 0);
  no_index.removeMetaValue("map_index");
  peptides.push_back(no_index);
  PeptideIdentification no_hits = makeId("PEPTIDEK", 2, 0);
  no_hits.setHits({});
  peptides.push_back(no_hits);
  const ConsensusMap map = map_with(peptides);
  // in the order of the peptide identifications
  ABORT_IF(Internal::ClassTest::peptidesOf(map, {0}).size() != 6)
  TEST_EQUAL(Internal::ClassTest::peptidesOf(map, {0})[1].getHits()[0].getSequence().toString(), "PEPTIDER(Label:13C(6)15N(4))")

  const std::vector<bool> consistent = resolver.checkIdentifications(map[0], map.getIdentificationData());
  TEST_EQUAL(consistent.size(), 6)
  TEST_TRUE(consistent[0])
  TEST_FALSE(consistent[1])
  TEST_TRUE(consistent[2])
  TEST_FALSE(consistent[3])
  TEST_FALSE(consistent[4])
  TEST_FALSE(consistent[5])

  // the same verdict as resolve() reaches for the first identification
  ConsensusMap out, conflicts;
  resolver.resolve(map_with({wrong}), out, conflicts);
  TEST_EQUAL(out.size(), 0)
  TEST_EQUAL(conflicts.size(), 1)

  // a multiplet without features has nothing to agree with
  TEST_EQUAL(resolver.checkIdentifications(ConsensusFeature(), map.getIdentificationData()).size(), 0)

  // a light peptide with K and R fits the doublet at the Lys8 + Arg10 shift, not the one at the Lys8 shift
  {
    Param p = resolver.getParameters();
    p.setValue("algorithm:max_nr_labelled_aas", 1);
    resolver.setParameters(p);
    ConsensusMap doublets;
    ConsensusFeature lys8_doublet;
    lys8_doublet.setRT(100.0);
    lys8_doublet.setMZ(500.0);
    lys8_doublet.setCharge(2);
    lys8_doublet.insert(makeHandle(0, 100.0, 500.0, 2, 100.0));
    lys8_doublet.insert(makeHandle(1, 100.0, 500.0 + LYS8 / 2, 2, 200.0));
    doublets.push_back(lys8_doublet);

    ConsensusFeature kr_doublet;
    kr_doublet.setRT(100.0);
    kr_doublet.setMZ(500.0);
    kr_doublet.setCharge(2);
    kr_doublet.insert(makeHandle(0, 100.0, 500.0, 2, 100.0));
    kr_doublet.insert(makeHandle(1, 100.0, 500.0 + (LYS8 + ARG10) / 2, 2, 200.0));
    doublets.push_back(kr_doublet);
    identify(doublets, {{0, makeId("PEPTIDEKAAR", 2, 0)}, {1, makeId("PEPTIDEKAAR", 2, 0)}});
    TEST_FALSE(resolver.checkIdentifications(doublets[0], doublets.getIdentificationData())[0])
    TEST_TRUE(resolver.checkIdentifications(doublets[1], doublets.getIdentificationData())[0])
  }
}
END_SECTION

START_SECTION((std::vector<double> expectedChannelMZs(const AASequence& sequence, int charge, Size& position) const))
{
  MultiplexResolverAlgorithm resolver;
  Size position = 99;
  // the light peptide of a K-terminal peptide: light, then Lys8
  const AASequence light = AASequence::fromString("PEPTIDEK");
  std::vector<double> mzs = resolver.expectedChannelMZs(light, 2, position);
  TEST_EQUAL(mzs.size(), 2)
  TEST_EQUAL(position, 0)
  TEST_REAL_SIMILAR(mzs[0], light.getMZ(2))
  TEST_REAL_SIMILAR(mzs[1], light.getMZ(2) + LYS8 / 2)

  // its heavy form names the heavy channel; an R-terminal light peptide gets the Arg10 shift
  mzs = resolver.expectedChannelMZs(AASequence::fromString("PEPTIDEK(Label:13C(6)15N(2))"), 2, position);
  TEST_EQUAL(mzs.size(), 2)
  TEST_EQUAL(position, 1)
  TEST_REAL_SIMILAR(mzs[0], light.getMZ(2))
  mzs = resolver.expectedChannelMZs(AASequence::fromString("PEPTIDER"), 3, position);
  TEST_EQUAL(mzs.size(), 2)
  TEST_REAL_SIMILAR(mzs[1] - mzs[0], ARG10 / 3)

  // no pattern: two labelled residues with max_nr_labelled_aas = 0 (two K, or K and R), or no charge
  TEST_EQUAL(resolver.expectedChannelMZs(AASequence::fromString("PEKPTIDEK"), 2, position).size(), 0)
  const AASequence kr = AASequence::fromString("PEPTIDEKAAR");
  TEST_EQUAL(resolver.expectedChannelMZs(kr, 2, position).size(), 0)
  TEST_EQUAL(resolver.expectedChannelMZs(light, 0, position).size(), 0)

  // a peptide with K and R is shifted by both labels in the heavy channel
  Param p = resolver.getParameters();
  p.setValue("algorithm:max_nr_labelled_aas", 1);
  resolver.setParameters(p);
  mzs = resolver.expectedChannelMZs(kr, 2, position);
  TEST_EQUAL(mzs.size(), 2)
  TEST_EQUAL(position, 0)
  TEST_REAL_SIMILAR(mzs[0], kr.getMZ(2))
  TEST_REAL_SIMILAR(mzs[1] - mzs[0], (LYS8 + ARG10) / 2)
  mzs = resolver.expectedChannelMZs(AASequence::fromString("PEPTIDEK(Label:13C(6)15N(2))AAR(Label:13C(6)15N(4))"), 2, position);
  TEST_EQUAL(mzs.size(), 2)
  TEST_EQUAL(position, 1)
  TEST_REAL_SIMILAR(mzs[0], kr.getMZ(2))
  // K-only and R-only peptides keep their single shift
  mzs = resolver.expectedChannelMZs(light, 2, position);
  TEST_EQUAL(mzs.size(), 2)
  TEST_REAL_SIMILAR(mzs[1] - mzs[0], LYS8 / 2)
  mzs = resolver.expectedChannelMZs(AASequence::fromString("PEPTIDER"), 2, position);
  TEST_EQUAL(mzs.size(), 2)
  TEST_REAL_SIMILAR(mzs[1] - mzs[0], ARG10 / 2)
  // no pattern: two K and one R exceed two labelled residues (no Arg10-only shift), no K or R, or a partially
  // labelled peptide (heavy K, light R) that no sample produces
  TEST_EQUAL(resolver.expectedChannelMZs(AASequence::fromString("PEKPTIDEKR"), 2, position).size(), 0)
  TEST_EQUAL(resolver.expectedChannelMZs(AASequence::fromString("PEPTIDE"), 2, position).size(), 0)
  TEST_EQUAL(resolver.expectedChannelMZs(AASequence::fromString("PEPTIDEK(Label:13C(6)15N(2))AAR"), 2, position).size(), 0)

  // triple SILAC: both labelled channels carry both labels
  p.setValue("algorithm:labels", "[][Lys4,Arg6][Lys8,Arg10]");
  resolver.setParameters(p);
  mzs = resolver.expectedChannelMZs(kr, 2, position);
  TEST_EQUAL(mzs.size(), 3)
  TEST_EQUAL(position, 0)
  TEST_REAL_SIMILAR(mzs[1] - mzs[0], (LYS4 + ARG6) / 2)
  TEST_REAL_SIMILAR(mzs[2] - mzs[0], (LYS8 + ARG10) / 2)

  // the overload also gives the labels of every channel, one entry per labelled residue
  std::vector<MultiplexDeltaMasses::LabelSet> label_sets;
  mzs = resolver.expectedChannelMZs(AASequence::fromString("PEK(Label:2H(4))PTIDEK(Label:2H(4))"), 2, position, label_sets);
  TEST_EQUAL(mzs.size(), 3)
  TEST_EQUAL(position, 1)
  TEST_EQUAL(label_sets.size(), 3)
  TEST_EQUAL(MultiplexDeltaMasses::labelSetToString(label_sets[0]), "no_label")
  TEST_EQUAL(label_sets[1].count("Lys4"), 2)
  TEST_EQUAL(label_sets[1].size(), 2)
  TEST_EQUAL(label_sets[2].count("Lys8"), 2)
  TEST_REAL_SIMILAR(mzs[1] - mzs[0], 2 * LYS4 / 2)
  // no pattern: no label sets either
  TEST_EQUAL(resolver.expectedChannelMZs(AASequence::fromString("PEPTIDE"), 2, position, label_sets).size(), 0)
  TEST_EQUAL(label_sets.size(), 0)

  // the sample of a channel is found by its labels, not by its position (implicit light sample)
  p.setValue("algorithm:labels", "[Lys8,Arg10]");
  resolver.setParameters(p);
  mzs = resolver.expectedChannelMZs(kr, 2, position);
  TEST_EQUAL(mzs.size(), 2)
  TEST_REAL_SIMILAR(mzs[1] - mzs[0], (LYS8 + ARG10) / 2)

  // a residue kind the sample does not label imposes nothing: with Lys8 labelling only, R is not shifted
  p.setValue("algorithm:labels", "[][Lys8]");
  resolver.setParameters(p);
  mzs = resolver.expectedChannelMZs(kr, 2, position);
  TEST_EQUAL(mzs.size(), 2)
  TEST_REAL_SIMILAR(mzs[1] - mzs[0], LYS8 / 2)

  // labels that are not counted per residue (Dimethyl: N-terminus and K) follow the labels of the sequence
  p.setValue("algorithm:labels", "[Dimethyl0][Dimethyl6]");
  resolver.setParameters(p);
  mzs = resolver.expectedChannelMZs(AASequence::fromString(".(Dimethyl)PEPTIDEK(Dimethyl)"), 2, position);
  TEST_EQUAL(mzs.size(), 2)
  TEST_EQUAL(position, 0)
  TEST_REAL_SIMILAR(mzs[1] - mzs[0], 2 * (DIMETHYL6 - DIMETHYL0) / 2)
}
END_SECTION

START_SECTION((FAIMS metadata and CV - specific blacklisting))
{
  auto in = makeDuplexMap();
  in[4].setMetaValue(Constants::UserParam::FAIMS_CV, -65.0);
  MSExperiment blacklist;
  MSSpectrum spec;
  spec.setMSLevel(1);
  spec.setRT(500.0);
  spec.setMetaValue(Constants::UserParam::FAIMS_CV, -45.0);
  spec.push_back(Peak1D(900.0 - ARG10 / 2, 1000.0));
  blacklist.addSpectrum(spec);
  blacklist.updateRanges();
  MultiplexResolverAlgorithm resolver;
  ConsensusMap out, conflicts;
  resolver.resolve(in, out, conflicts, blacklist);
  TEST_REAL_SIMILAR(out[2].getFeatures().begin()->getIntensity(), 0.0)
  TEST_REAL_SIMILAR((double)out[2].getMetaValue(Constants::UserParam::FAIMS_CV), -65.0)

  // A peak at the same CV must still mask the missing partner.
  blacklist[0].setMetaValue(Constants::UserParam::FAIMS_CV, -65.0);
  resolver.resolve(in, out, conflicts, blacklist);
  TEST_TRUE(std::isnan(out[2].getFeatures().begin()->getIntensity()))
  // Also accept CV metadata in its native mzML representation.
  blacklist[0].removeMetaValue(Constants::UserParam::FAIMS_CV);
  blacklist[0].setDriftTimeUnit(DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE);
  blacklist[0].setDriftTime(-65.0);
  resolver.resolve(in, out, conflicts, blacklist);
  TEST_TRUE(std::isnan(out[2].getFeatures().begin()->getIntensity()))
}
END_SECTION

START_SECTION((static UInt64 blacklistOwner(const MSSpectrum& spectrum, Size peak_index)))
{
  // 64-bit unique ids, also with the upper bit of either half set, are split into two 32-bit integer arrays and back
  const std::vector<UInt64> owners = {0xF123456789ABCDEFULL, UniqueIdInterface::INVALID, 0x00000001FFFFFFFFULL, 42};
  MSSpectrum spectrum;
  spectrum.setRT(10.0);
  for (Size i = 0; i < owners.size(); ++i) { spectrum.push_back(Peak1D(500.0 + i, 1000.0)); }
  TEST_EQUAL(MultiplexResolverAlgorithm::blacklistOwner(spectrum, 0), UniqueIdInterface::INVALID) // nothing recorded
  MultiplexResolverAlgorithm::setBlacklistOwners(spectrum, owners);
  TEST_EQUAL(spectrum.getIntegerDataArrays().size(), 2)
  for (Size i = 0; i < owners.size(); ++i) { TEST_EQUAL(MultiplexResolverAlgorithm::blacklistOwner(spectrum, i), owners[i]) }
  TEST_EQUAL(MultiplexResolverAlgorithm::blacklistOwner(spectrum, owners.size()), UniqueIdInterface::INVALID)
  // setting them again replaces them
  MultiplexResolverAlgorithm::setBlacklistOwners(spectrum, std::vector<UInt64>(owners.size(), 7));
  TEST_EQUAL(spectrum.getIntegerDataArrays().size(), 2)
  TEST_EQUAL(MultiplexResolverAlgorithm::blacklistOwner(spectrum, 0), 7)
  MultiplexResolverAlgorithm::setBlacklistOwners(spectrum, owners);

  // the owners stay aligned with the peaks when these are sorted, and survive an mzML file
  MSExperiment blacklist;
  MSSpectrum empty;
  empty.setRT(5.0);
  blacklist.addSpectrum(empty); // no peaks, no owners
  MSSpectrum reversed = spectrum;
  std::reverse(reversed.begin(), reversed.end());
  for (auto& array : reversed.getIntegerDataArrays()) { std::reverse(array.begin(), array.end()); }
  reversed.sortByPosition();
  blacklist.addSpectrum(reversed);
  std::string file;
  NEW_TMP_FILE(file)
  MzMLFile().store(file, blacklist);
  MSExperiment loaded;
  MzMLFile().load(file, loaded);
  TEST_EQUAL(loaded.size(), 2)
  ABORT_IF(loaded.size() != 2)
  TEST_EQUAL(loaded[1].size(), owners.size())
  for (Size i = 0; i < owners.size() && i < loaded[1].size(); ++i)
  {
    TEST_REAL_SIMILAR(loaded[1][i].getMZ(), 500.0 + i)
    TEST_EQUAL(MultiplexResolverAlgorithm::blacklistOwner(loaded[1], i), owners[i])
  }
  TEST_TRUE(MultiplexResolverAlgorithm::hasBlacklistOwners(loaded))
}
END_SECTION

START_SECTION((static bool hasBlacklistOwners(const MSExperiment& blacklist)))
{
  MSExperiment blacklist;
  TEST_FALSE(MultiplexResolverAlgorithm::hasBlacklistOwners(blacklist)) // no peaks
  MSSpectrum spectrum;
  spectrum.setRT(10.0);
  spectrum.push_back(Peak1D(500.0, 1000.0));
  spectrum.push_back(Peak1D(501.0, 1000.0));
  blacklist.addSpectrum(spectrum);
  TEST_FALSE(MultiplexResolverAlgorithm::hasBlacklistOwners(blacklist)) // anonymous
  MultiplexResolverAlgorithm::setBlacklistOwners(blacklist[0], {1, 2});
  TEST_TRUE(MultiplexResolverAlgorithm::hasBlacklistOwners(blacklist))
  MSSpectrum empty;
  empty.setRT(20.0);
  blacklist.addSpectrum(empty);
  TEST_TRUE(MultiplexResolverAlgorithm::hasBlacklistOwners(blacklist)) // spectra without peaks need no owners
  blacklist.addSpectrum(spectrum);
  TEST_FALSE(MultiplexResolverAlgorithm::hasBlacklistOwners(blacklist)) // a spectrum with peaks but without owners
  MultiplexResolverAlgorithm::setBlacklistOwners(blacklist[2], {1});
  TEST_FALSE(MultiplexResolverAlgorithm::hasBlacklistOwners(blacklist)) // not aligned with the peaks
}
END_SECTION

START_SECTION((static OwnerRelation classifyBlacklistOwner(UInt64 owner, const ConsensusFeature& consensus, const OwnerPeptides& peptides)))
{
  using Relation = MultiplexResolverAlgorithm::OwnerRelation;
  ConsensusMap map;
  std::vector<std::pair<Size, PeptideIdentification>> ids;
  // a multiplet, to be added to the map next
  auto multiplet = [&](UInt64 id, int charge, const std::string& sequence)
  {
    ConsensusFeature cf;
    cf.setCharge(charge);
    FeatureHandle h = makeHandle(0, 100.0, 500.0, charge, 100.0);
    h.setUniqueId(id);
    cf.insert(h);
    if (!sequence.empty()) { ids.emplace_back(map.size(), makeId(sequence, charge, 0)); }
    return cf;
  };
  ConsensusFeature own = multiplet(11, 2, "PEPTIDEK");
  own.insert(makeHandle(1, 100.0, 500.0 + LYS8 / 2, 2, 0.0)); // a dummy (no unique id) changes nothing
  map.push_back(own);
  map.push_back(multiplet(12, 2, "PEPTIDEK(Label:13C(6)15N(2))")); // same peptide, other label state
  map.push_back(multiplet(13, 3, "PEPTIDEK"));                     // other charge
  map.push_back(multiplet(14, 2, "EPPTIDEK"));                     // same composition, other sequence
  map.push_back(multiplet(15, 2, "PEPTM(Oxidation)IDEK"));         // other modification state
  map.push_back(multiplet(16, 2, ""));                             // unidentified
  identify(map, ids);

  MultiplexResolverAlgorithm resolver;
  const MultiplexResolverAlgorithm::OwnerPeptides peptides = resolver.ownerPeptides(map);
  TEST_EQUAL(peptides.size(), 5)
  TEST_EQUAL(peptides.at(11), peptides.at(12))

  auto classify = [&](UInt64 owner) { return MultiplexResolverAlgorithm::classifyBlacklistOwner(owner, own, peptides); };
  TEST_TRUE(classify(UniqueIdInterface::INVALID) == Relation::ORPHAN)
  TEST_TRUE(classify(11) == Relation::OWN)
  TEST_TRUE(classify(12) == Relation::SAME_PEPTIDE)
  TEST_TRUE(classify(13) == Relation::OTHER_PEPTIDE)
  TEST_TRUE(classify(14) == Relation::OTHER_PEPTIDE)
  TEST_TRUE(classify(15) == Relation::OTHER_PEPTIDE)
  TEST_TRUE(classify(16) == Relation::UNIDENTIFIED)
  TEST_TRUE(classify(99) == Relation::UNIDENTIFIED) // not in the map at all

  TEST_FALSE(MultiplexResolverAlgorithm::isForeign(Relation::ORPHAN))
  TEST_FALSE(MultiplexResolverAlgorithm::isForeign(Relation::OWN))
  TEST_TRUE(MultiplexResolverAlgorithm::isForeign(Relation::SAME_PEPTIDE)) // quantified by its own multiplet
  TEST_TRUE(MultiplexResolverAlgorithm::isForeign(Relation::UNIDENTIFIED))
  TEST_TRUE(MultiplexResolverAlgorithm::isForeign(Relation::OTHER_PEPTIDE))
}
END_SECTION

START_SECTION((std::string peptideKey(const std::string& sequence, int charge) const))
{
  MultiplexResolverAlgorithm resolver; // [][Lys8,Arg10]
  auto key = [&](const std::string& sequence, int charge) { return resolver.peptideKey(AASequence::fromString(sequence).toString(), charge); };
  TEST_EQUAL(key("PEPTIDEK", 2), "PEPTIDEK/2")
  TEST_EQUAL(key("PEPTIDEK(Label:13C(6)15N(2))", 2), "PEPTIDEK/2")              // the label states of a peptide share it
  TEST_EQUAL(key("PEPTIDEKR(Label:13C(6)15N(4))", 2), key("PEPTIDEKR", 2))
  TEST_EQUAL(key("PEPTIDEK(Label:13C(6)15N(2))R(Label:13C(6)15N(4))", 2), key("PEPTIDEKR", 2))
  TEST_NOT_EQUAL(key("PEPTIDEK", 3), key("PEPTIDEK", 2))                      // charge
  TEST_NOT_EQUAL(key("EPPTIDEK", 2), key("PEPTIDEK", 2))                      // same composition, other sequence
  TEST_EQUAL(key("PEPTM(Oxidation)IDEK(Label:13C(6)15N(2))", 2), "PEPTM(Oxidation)IDEK/2") // other modifications are kept

  // the key of ownerPeptides()
  ConsensusMap map;
  ConsensusFeature cf;
  cf.setCharge(2);
  FeatureHandle h = makeHandle(0, 100.0, 500.0, 2, 100.0);
  h.setUniqueId(21);
  cf.insert(h);
  map.push_back(cf);
  identify(map, {{0, makeId("PEPTIDEK(Label:13C(6)15N(2))", 2, 0)}});
  TEST_EQUAL(resolver.ownerPeptides(map).at(21), key("PEPTIDEK", 2))

  // labels of other configurations: Dimethyl sits on the N-terminus as well
  MultiplexResolverAlgorithm dimethyl;
  Param p = dimethyl.getParameters();
  p.setValue("algorithm:labels", "[Dimethyl0][Dimethyl6]");
  dimethyl.setParameters(p);
  const std::string light = AASequence::fromString(".(Dimethyl)PEPTIDEK(Dimethyl)").toString();
  const std::string heavy = AASequence::fromString(".(Dimethyl:2H(4)13C(2))PEPTIDEK(Dimethyl:2H(4)13C(2))").toString();
  const std::string other = AASequence::fromString(".(Dimethyl)EPPTIDEK(Dimethyl)").toString();
  TEST_EQUAL(dimethyl.peptideKey(light, 2), dimethyl.peptideKey(heavy, 2))
  TEST_NOT_EQUAL(dimethyl.peptideKey(light, 2), dimethyl.peptideKey(other, 2))
}
END_SECTION

START_SECTION((static OwnerRelation classifyBlacklistOwner(UInt64 owner, const std::string& peptide, const OwnerPeptides& peptides)))
{
  using Relation = MultiplexResolverAlgorithm::OwnerRelation;
  MultiplexResolverAlgorithm resolver;
  ConsensusMap map;
  std::vector<std::pair<Size, PeptideIdentification>> ids;
  auto add = [&](UInt64 id, const std::string& sequence)
  {
    ConsensusFeature cf;
    cf.setCharge(2);
    FeatureHandle h = makeHandle(0, 100.0, 500.0, 2, 100.0);
    h.setUniqueId(id);
    cf.insert(h);
    if (!sequence.empty()) { ids.emplace_back(map.size(), makeId(sequence, 2, 0)); }
    map.push_back(cf);
  };
  add(31, "PEPTIDEK(Label:13C(6)15N(2))"); // the peptide, in another label state
  add(32, "EPPTIDEK");                     // another peptide
  add(33, "");                             // unidentified
  identify(map, ids);
  const MultiplexResolverAlgorithm::OwnerPeptides peptides = resolver.ownerPeptides(map);
  // a peptide without a multiplet of its own, e.g. an identification that mapped onto none
  const std::string peptide = resolver.peptideKey("PEPTIDEK", 2);
  auto classify = [&](UInt64 owner) { return MultiplexResolverAlgorithm::classifyBlacklistOwner(owner, peptide, peptides); };
  TEST_TRUE(classify(UniqueIdInterface::INVALID) == Relation::ORPHAN)
  TEST_TRUE(classify(31) == Relation::SAME_PEPTIDE)
  TEST_TRUE(classify(32) == Relation::OTHER_PEPTIDE)
  TEST_TRUE(classify(33) == Relation::UNIDENTIFIED)
  TEST_TRUE(classify(99) == Relation::UNIDENTIFIED)
}
END_SECTION

START_SECTION(([EXTRA] resolve with a blacklist that records owners))
{
  // E: only the heavy channel of AAAAR (Arg10) detected at RT 500; its light channel is missing at 900 - ARG10/2.
  // S: the light channel of the same peptide, detected on its own; F: another peptide; U: an unidentified multiplet.
  ConsensusMap in = makeDuplexMap();
  in.clear(false);
  std::vector<std::pair<Size, PeptideIdentification>> ids;
  // a multiplet, to be added to the map next
  auto multiplet = [&](UInt64 id, double rt, double mz, const std::string& sequence)
  {
    ConsensusFeature cf;
    cf.setRT(rt);
    cf.setMZ(mz);
    cf.setCharge(2);
    cf.setIntensity(50.0);
    FeatureHandle h = makeHandle(0, rt, mz, 2, 50.0);
    h.setUniqueId(id);
    cf.insert(h);
    if (!sequence.empty()) { ids.emplace_back(in.size(), makeId(sequence, 2, 0)); }
    return cf;
  };
  const UInt64 E = 101, S = 201, F = 301, U = 401;
  const double dummy_mz = 900.0 - ARG10 / 2;
  in.push_back(multiplet(E, 500.0, 900.0, "AAAAR(Label:13C(6)15N(4))"));
  in.push_back(multiplet(S, 520.0, dummy_mz, "AAAAR"));
  in.push_back(multiplet(F, 300.0, 700.0, "CCCCK"));
  in.push_back(multiplet(U, 400.0, 800.0, ""));
  identify(in, ids);

  auto blacklistAt = [](double rt, double mz, UInt64 owner)
  {
    MSExperiment blacklist;
    MSSpectrum spec;
    spec.setRT(rt);
    spec.setMSLevel(1);
    spec.push_back(Peak1D(mz, 1000.0));
    MultiplexResolverAlgorithm::setBlacklistOwners(spec, {owner});
    blacklist.addSpectrum(spec);
    blacklist.updateRanges();
    return blacklist;
  };
  // the owners are used only on request (by a caller that measures the channels left at 0, like MS1LabeledWorkflow)
  MultiplexResolverAlgorithm ignoring_owners;
  MultiplexResolverAlgorithm resolver;
  {
    Param p = resolver.getParameters();
    TEST_EQUAL(p.getValue("algorithm:use_blacklist_owners").toString(), "false")
    p.setValue("algorithm:use_blacklist_owners", "true");
    resolver.setParameters(p);
  }
  // the intensity of E's light dummy, given the blacklist
  auto dummyOf = [&](const MultiplexResolverAlgorithm& r, const MSExperiment& blacklist) -> double
  {
    ConsensusMap out, conflicts;
    r.resolve(in, out, conflicts, blacklist);
    for (const ConsensusFeature& cf : out)
    {
      if (cf.getRT() == 500.0) { return cf.getFeatures().begin()->getIntensity(); }
    }
    return -1.0;
  };
  auto dummyOfE = [&](const MSExperiment& blacklist) { return dummyOf(resolver, blacklist); };

  // a peak at the dummy position reserved by no reported feature or by E itself: 0 (not occupied), to be re-quantified
  TEST_REAL_SIMILAR(dummyOfE(blacklistAt(500.0, dummy_mz, UniqueIdInterface::INVALID)), 0.0)
  TEST_REAL_SIMILAR(dummyOfE(blacklistAt(500.0, dummy_mz, E)), 0.0)
  // ... reserved by another multiplet of the same peptide (S quantifies it), another peptide or an unidentified multiplet:
  // not quantifiable
  TEST_TRUE(std::isnan(dummyOfE(blacklistAt(500.0, dummy_mz, S))))
  TEST_TRUE(std::isnan(dummyOfE(blacklistAt(500.0, dummy_mz, F))))
  TEST_TRUE(std::isnan(dummyOfE(blacklistAt(500.0, dummy_mz, U))))
  // ... the same peak without owners, or with owners that are not used: not quantifiable (as before), whoever reserved it
  {
    MSExperiment anonymous = blacklistAt(500.0, dummy_mz, S);
    anonymous[0].getIntegerDataArrays().clear();
    TEST_TRUE(std::isnan(dummyOfE(anonymous)))
    TEST_TRUE(std::isnan(dummyOf(ignoring_owners, blacklistAt(500.0, dummy_mz, UniqueIdInterface::INVALID))))
    TEST_TRUE(std::isnan(dummyOf(ignoring_owners, blacklistAt(500.0, dummy_mz, E))))
    TEST_TRUE(std::isnan(dummyOf(ignoring_owners, blacklistAt(500.0, dummy_mz, S))))
    TEST_REAL_SIMILAR(dummyOf(ignoring_owners, blacklistAt(500.0, dummy_mz - 0.002, F)), 0.0) // one-sided, as without owners
  }

  // with owners, peaks below the expected position count as well (the anonymous check looks at the first peak above it only)
  {
    MSExperiment below = blacklistAt(500.0, dummy_mz - 0.002, F); // 2.2 ppm
    TEST_TRUE(std::isnan(dummyOfE(below)))
    below[0].getIntegerDataArrays().clear();
    TEST_REAL_SIMILAR(dummyOfE(below), 0.0)
    // outside the tolerance (10 ppm) or the RT window (5 s): nothing
    TEST_REAL_SIMILAR(dummyOfE(blacklistAt(500.0, dummy_mz - 0.01, F)), 0.0)
    TEST_REAL_SIMILAR(dummyOfE(blacklistAt(506.0, dummy_mz, F)), 0.0)
    // the M+1 position
    TEST_TRUE(std::isnan(dummyOfE(blacklistAt(500.0, dummy_mz + Constants::C13C12_MASSDIFF_U / 2, F))))
  }

  // a channel added for an identification (intensity 0, no unique id) is checked with owners only, also in a complete multiplet
  {
    // E with its heavy channel detected and the light one added (the identification stays linked)
    ConsensusMap in2 = in;
    ConsensusFeature& added = in2[0];
    added.setMZ(dummy_mz);
    added.setFeatures(ConsensusFeature::HandleSetType());
    FeatureHandle heavy = makeHandle(1, 500.0, 900.0, 2, 50.0);
    heavy.setUniqueId(E);
    added.insert(heavy);
    added.insert(makeHandle(0, 500.0, dummy_mz, 2, 0.0));
    auto ids2 = Internal::ClassTest::identificationsOf(in2);
    ids2.features[{0}][0].setMetaValue("map_index", 1);
    Internal::ClassTest::setIdentifications(in2, ids2);
    auto addedChannel = [&](const MSExperiment& blacklist) -> double
    {
      ConsensusMap out, conflicts;
      resolver.resolve(in2, out, conflicts, blacklist);
      TEST_EQUAL(out[0].size(), 2)
      return out[0].getFeatures().begin()->getIntensity();
    };
    TEST_TRUE(std::isnan(addedChannel(blacklistAt(500.0, dummy_mz, F))))
    TEST_TRUE(std::isnan(addedChannel(blacklistAt(500.0, dummy_mz, S))))
    TEST_REAL_SIMILAR(addedChannel(blacklistAt(500.0, dummy_mz, UniqueIdInterface::INVALID)), 0.0)
    MSExperiment anonymous = blacklistAt(500.0, dummy_mz, F);
    anonymous[0].getIntegerDataArrays().clear();
    TEST_REAL_SIMILAR(addedChannel(anonymous), 0.0)
    {
      ConsensusMap out, conflicts;
      ignoring_owners.resolve(in2, out, conflicts, blacklistAt(500.0, dummy_mz, F));
      TEST_REAL_SIMILAR(out[0].getFeatures().begin()->getIntensity(), 0.0) // not checked
    }

    // triple SILAC: an added light channel next to the detected heavy one, the medium one completed with a dummy
    MultiplexResolverAlgorithm triple;
    Param p = triple.getParameters();
    p.setValue("algorithm:labels", "[][Lys4,Arg6][Lys8,Arg10]");
    p.setValue("algorithm:use_blacklist_owners", "true");
    triple.setParameters(p);
    ConsensusMap in3 = in2;
    in3.getColumnHeaders()[2] = in3.getColumnHeaders()[1];
    ConsensusMap out3, conflicts3;
    triple.resolve(in3, out3, conflicts3, blacklistAt(500.0, dummy_mz, F));
    TEST_EQUAL(out3[0].size(), 3)
    ABORT_IF(out3[0].size() != 3)
    auto it = out3[0].getFeatures().begin();
    TEST_TRUE(std::isnan(it->getIntensity())) // the added light channel
    ++it;
    TEST_REAL_SIMILAR(it->getMZ(), dummy_mz + ARG6 / 2)
    TEST_REAL_SIMILAR(it->getIntensity(), 0.0) // the medium dummy
    ++it;
    TEST_REAL_SIMILAR(it->getIntensity(), 50.0)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
END_TEST
