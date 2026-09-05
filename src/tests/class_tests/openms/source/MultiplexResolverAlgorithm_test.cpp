// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/FEATUREFINDER/MultiplexResolverAlgorithm.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <cmath>

using namespace OpenMS;

// A SILAC duplex ([][Lys8,Arg10]) test map, built the way FeatureFinderMultiplexAlgorithm + IDMapper
// (annotate_ids_with_subelements) would leave it: column 0 = light, column 1 = heavy, and the
// identification of a feature carries the map index of the handle it was mapped to.
const double LYS8 = 8.0141988132;
const double ARG10 = 10.008268600;

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

  // A: complete doublet (Lys8 shift), heavy peptide identified
  {
    ConsensusFeature cf;
    cf.setRT(100.0);
    cf.setMZ(500.0);
    cf.setCharge(2);
    cf.setIntensity(300.0);
    cf.insert(makeHandle(0, 100.0, 500.0, 2, 100.0));
    cf.insert(makeHandle(1, 100.0, 500.0 + LYS8 / 2, 2, 200.0));
    cf.getPeptideIdentifications().push_back(makeId("PEPTIDEK(Label:13C(6)15N(2))", 2, 1));
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
    cf.getPeptideIdentifications().push_back(makeId("AAAAK(Label:13C(6)15N(2))", 2, 0));
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
    cf.getPeptideIdentifications().push_back(makeId("AAAAK(Label:13C(6)15N(2))", 2, 1));
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
    cf.getPeptideIdentifications().push_back(makeId("AAAAR(Label:13C(6)15N(4))", 2, 0));
    map.push_back(cf);
  }
  map.applyMemberFunction(&UniqueIdInterface::setUniqueId);
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

START_SECTION((void resolve(const ConsensusMap& map_in, ConsensusMap& map_out, ConsensusMap& map_conflicts, const MSExperiment& blacklist = MSExperiment()) const))
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
  TEST_EQUAL(out[0].getPeptideIdentifications().size(), 1)

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
    TEST_EQUAL((int)cf.getPeptideIdentifications()[0].getHits()[0].getMetaValue("map_index"), 1)
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
  TEST_EQUAL(conflicts[0].getPeptideIdentifications().size(), 1)
  TEST_REAL_SIMILAR(conflicts[1].getMZ(), 800.0)
  TEST_EQUAL(conflicts[1].getPeptideIdentifications().size(), 0)

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
    cf.getPeptideIdentifications().push_back(makeId("PEPTIDEK(Label:13C(6)15N(2))AAR(Label:13C(6)15N(4))", 2, 1));
    in2.push_back(cf);

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

  // an identification without the 'map_index' annotation cannot be resolved
  {
    ConsensusMap in3 = makeDuplexMap();
    in3[0].getPeptideIdentifications()[0].removeMetaValue("map_index");
    ConsensusMap out4, conflicts4;
    TEST_EXCEPTION(Exception::MissingInformation, resolver.resolve(in3, out4, conflicts4))
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
END_TEST
