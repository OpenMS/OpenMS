// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: David Voigt $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/METADATA/AnnotatedMSRun.h>
#include <unordered_set>
#include <unordered_map>

using namespace OpenMS;

namespace
{
  using ID = IdentificationData;

  /// Two MS2 spectra (RT 10 and 20), each linking an identification of a search with one match (PEPTIDER, PEPTIDAR)
  AnnotatedMSRun makeRun()
  {
    MSExperiment exp;
    MSSpectrum spec1, spec2;
    spec1.setRT(10.0);
    spec1.setMSLevel(2);
    spec2.setRT(20.0);
    spec2.setMSLevel(2);
    exp.addSpectrum(spec1);
    exp.addSpectrum(spec2);
    AnnotatedMSRun annotated_data(std::move(exp));

    auto& run = annotated_data.getIdentificationData().addRun("search");
    ID::ScoreDefinition score;
    score.name = "score";
    run.setPrimaryScore(run.addScore(score));
    ID::SourceFile file;
    file.path = "test.mzML";
    const auto source = run.addSource(file);
    for (Size index = 0; index < 2; ++index)
    {
      ID::Observation observation;
      observation.rt = annotated_data.getMSExperiment()[index].getRT();
      const auto query = run.addIdentification(source, observation);
      ID::MatchData match;
      match.representation = index == 0 ? "PEPTIDER" : "PEPTIDAR";
      match.charge = 2;
      run.addMatch(query, match, {10.0 * (index + 1)});
      annotated_data.addSpectrumIDQuery(index, {run.getUuid(), query});
    }
    return annotated_data;
  }
}

START_TEST(AnnotatedMSRun, "$Id$")

// Default constructor
AnnotatedMSRun* ptr = nullptr;
AnnotatedMSRun* nullPointer = nullptr;

START_SECTION((AnnotatedMSRun()))
  ptr = new AnnotatedMSRun();
  TEST_NOT_EQUAL(ptr, nullPointer)
  TEST_EQUAL(ptr->getIdentificationData().empty(), true)
  TEST_EQUAL(ptr->hasSpectrumIdentifications(), false)
END_SECTION

START_SECTION((~AnnotatedMSRun()))
  delete ptr;
END_SECTION

START_SECTION((explicit AnnotatedMSRun(MSExperiment&& experiment)))
  MSExperiment exp;
  MSSpectrum spec;
  spec.setRT(42.0);
  spec.setMSLevel(2);
  exp.addSpectrum(spec);

  AnnotatedMSRun annotated_data(std::move(exp));
  TEST_EQUAL(annotated_data.getMSExperiment().size(), 1)
  TEST_REAL_SIMILAR(annotated_data.getMSExperiment()[0].getRT(), 42.0)
END_SECTION

START_SECTION((IdentificationData& getIdentificationData()))
  AnnotatedMSRun annotated_data;
  annotated_data.getIdentificationData().addRun("Test");
  TEST_EQUAL(annotated_data.getIdentificationData().getRuns().size(), 1)
  TEST_EQUAL(annotated_data.getIdentificationData().getRuns()[0].getIdentifier(), "Test")
END_SECTION

START_SECTION((const IdentificationData& getIdentificationData() const))
  AnnotatedMSRun annotated_data;
  annotated_data.getIdentificationData().addRun("Test");
  const AnnotatedMSRun& const_data = annotated_data;
  TEST_EQUAL(const_data.getIdentificationData().getRuns()[0].getIdentifier(), "Test")
END_SECTION

START_SECTION((const std::set<IdentificationData::QueryReference>& getSpectrumIDQueries(Size index) const))
  AnnotatedMSRun annotated_data = makeRun();
  const auto& run = annotated_data.getIdentificationData().getRuns()[0];
  TEST_EQUAL(annotated_data.getSpectrumIDQueries(0).size(), 1)
  TEST_EQUAL(annotated_data.getSpectrumIDQueries(1).size(), 1)
  TEST_EQUAL(annotated_data.getSpectrumIDQueries(0).begin()->run_uuid, run.getUuid())
  TEST_EQUAL(annotated_data.getSpectrumIDQueries(0) == annotated_data.getSpectrumIDQueries(1), false)
  // a spectrum without links (also beyond the spectra) has none
  TEST_EQUAL(annotated_data.getSpectrumIDQueries(5).empty(), true)
END_SECTION

START_SECTION((void addSpectrumIDQuery(Size index, const IdentificationData::QueryReference& query)))
  AnnotatedMSRun annotated_data = makeRun();
  // an identification can be linked by several spectra
  const auto query = *annotated_data.getSpectrumIDQueries(0).begin();
  annotated_data.addSpectrumIDQuery(1, query);
  TEST_EQUAL(annotated_data.getSpectrumIDQueries(1).size(), 2)
  annotated_data.addSpectrumIDQuery(1, query); // links are a set
  TEST_EQUAL(annotated_data.getSpectrumIDQueries(1).size(), 2)
  TEST_EXCEPTION(Exception::IndexOverflow, annotated_data.addSpectrumIDQuery(2, query))
END_SECTION

START_SECTION((std::vector<IdentificationData::QueryMatches> getSpectrumIdentifications(Size index) const))
  AnnotatedMSRun annotated_data = makeRun();
  auto identifications = annotated_data.getSpectrumIdentifications(0);
  TEST_EQUAL(identifications.size(), 1)
  TEST_EQUAL(identifications[0].matches.size(), 1)
  TEST_EQUAL(identifications[0].matches[0]->representation, "PEPTIDER")
  TEST_REAL_SIMILAR(*identifications[0].query->rt, 10.0)

  // several identifications, in the order of export, each with all its matches
  annotated_data.addSpectrumIDQuery(1, *annotated_data.getSpectrumIDQueries(0).begin());
  identifications = annotated_data.getSpectrumIdentifications(1);
  TEST_EQUAL(identifications.size(), 2)
  TEST_EQUAL(identifications[0].matches[0]->representation, "PEPTIDER")
  TEST_EQUAL(identifications[1].matches[0]->representation, "PEPTIDAR")

  TEST_EQUAL(annotated_data.getSpectrumIdentifications(5).empty(), true)

  // links to identifications that are not there are rejected
  AnnotatedMSRun dangling = annotated_data;
  dangling.getIdentificationData().clear();
  TEST_EXCEPTION(Exception::MissingInformation, dangling.getSpectrumIdentifications(0))
END_SECTION

START_SECTION((bool hasSpectrumIdentifications() const))
  AnnotatedMSRun annotated_data = makeRun();
  TEST_EQUAL(annotated_data.hasSpectrumIdentifications(), true)
  annotated_data.clearSpectrumIDQueries();
  TEST_EQUAL(annotated_data.hasSpectrumIdentifications(), false)
END_SECTION

START_SECTION((void clearSpectrumIDQueries()))
  AnnotatedMSRun annotated_data = makeRun();
  annotated_data.clearSpectrumIDQueries();
  TEST_EQUAL(annotated_data.getSpectrumIDQueries(0).empty(), true)
  TEST_EQUAL(annotated_data.getSpectrumIDQueries(1).empty(), true)
  // the identifications stay
  TEST_EQUAL(annotated_data.getIdentificationData().getRuns()[0].getNumberOfMatches(), 2)
END_SECTION

START_SECTION((void clearIdentifications()))
  AnnotatedMSRun annotated_data = makeRun();
  annotated_data.clearIdentifications();
  TEST_EQUAL(annotated_data.getIdentificationData().empty(), true)
  TEST_EQUAL(annotated_data.hasSpectrumIdentifications(), false)
  TEST_EQUAL(annotated_data.getMSExperiment().size(), 2)
END_SECTION

START_SECTION((Size eraseMatches(const std::function<bool(const IdentificationData::Run&, const IdentificationData::Identification&, const IdentificationData::Match&)>& remove)))
  AnnotatedMSRun annotated_data = makeRun();
  TEST_EQUAL(annotated_data.eraseMatches([](const ID::Run&, const ID::Identification&, const ID::Match& match) {
    return match.representation == "PEPTIDER";
  }), 1)
  // the identification stays, without matches, and so does the link
  const auto identifications = annotated_data.getSpectrumIdentifications(0);
  TEST_EQUAL(identifications.size(), 1)
  TEST_EQUAL(identifications[0].matches.empty(), true)
  TEST_EQUAL(annotated_data.getSpectrumIdentifications(1)[0].matches.size(), 1)
END_SECTION

START_SECTION((Size eraseIdentifications(const std::function<bool(const IdentificationData::Run&, const IdentificationData::Identification&)>& remove)))
  AnnotatedMSRun annotated_data = makeRun();
  annotated_data.addSpectrumIDQuery(1, *annotated_data.getSpectrumIDQueries(0).begin());
  TEST_EQUAL(annotated_data.eraseIdentifications([](const ID::Run&, const ID::Identification& query) {
    return *query.rt < 15.0;
  }), 1)
  // the links of all spectra to it are gone
  TEST_EQUAL(annotated_data.getSpectrumIDQueries(0).empty(), true)
  TEST_EQUAL(annotated_data.getSpectrumIDQueries(1).size(), 1)
  TEST_EQUAL(annotated_data.getSpectrumIdentifications(1)[0].matches[0]->representation, "PEPTIDAR")
END_SECTION

START_SECTION((MSExperiment& getMSExperiment()))
  AnnotatedMSRun annotated_data;
  MSExperiment exp;
  MSSpectrum spec;
  spec.setRT(42.0);
  spec.setMSLevel(2);
  exp.addSpectrum(spec);

  annotated_data.getMSExperiment() = std::move(exp);
  TEST_EQUAL(annotated_data.getMSExperiment().size(), 1)
  TEST_REAL_SIMILAR(annotated_data.getMSExperiment()[0].getRT(), 42.0)
END_SECTION

START_SECTION((const MSExperiment& getMSExperiment() const))
  AnnotatedMSRun annotated_data;
  MSExperiment exp;
  MSSpectrum spec;
  spec.setRT(42.0);
  spec.setMSLevel(2);
  exp.addSpectrum(spec);

  annotated_data.getMSExperiment() = std::move(exp);

  const AnnotatedMSRun& const_data = annotated_data;
  TEST_EQUAL(const_data.getMSExperiment().size(), 1)
  TEST_REAL_SIMILAR(const_data.getMSExperiment()[0].getRT(), 42.0)
END_SECTION

START_SECTION((bool operator==(const AnnotatedMSRun& rhs) const))
  // a copy is equal (it has the same identification runs and links)
  AnnotatedMSRun run1 = makeRun();
  AnnotatedMSRun run2 = run1;

  // Test equality
  TEST_EQUAL(run1 == run2, true)
  TEST_EQUAL(run1 != run2, false)

  // Test inequality - modify run2
  run2.getMSExperiment()[0].setRT(100.0);
  TEST_EQUAL(run1 == run2, false)
  TEST_EQUAL(run1 != run2, true)

  // other links
  run2 = run1;
  run2.addSpectrumIDQuery(1, *run1.getSpectrumIDQueries(0).begin());
  TEST_EQUAL(run1 == run2, false)

  // other identifications
  run2 = run1;
  run2.getIdentificationData().addRun("other");
  TEST_EQUAL(run1 == run2, false)

  // no links equals links that were removed
  const auto first = *run1.getSpectrumIDQueries(0).begin();
  const auto is_first = [&](const ID::Run& run, const ID::Identification& query) { return ID::QueryReference {run.getUuid(), query.getId()} == first; };
  run1.clearSpectrumIDQueries();
  run2 = run1;
  run2.addSpectrumIDQuery(1, first);
  run1.eraseIdentifications(is_first);
  run2.eraseIdentifications(is_first); // leaves spectra with empty links
  TEST_EQUAL(run1 == run2, true)
END_SECTION

START_SECTION(([EXTRA] std::hash<AnnotatedMSRun>))
  // Test 1: Equal objects have equal hashes
  AnnotatedMSRun run1 = makeRun();
  AnnotatedMSRun run2 = run1;

  TEST_EQUAL(run1 == run2, true)
  TEST_EQUAL(std::hash<AnnotatedMSRun>{}(run1), std::hash<AnnotatedMSRun>{}(run2))

  // Test 2: Different objects have different hashes (likely, not guaranteed)
  AnnotatedMSRun run3;
  MSExperiment exp3;
  MSSpectrum spec3;
  spec3.setRT(100.0);  // Different RT
  spec3.setMSLevel(1); // Different MS level
  exp3.addSpectrum(spec3);
  run3.setMSExperiment(exp3);

  TEST_EQUAL(run1 == run3, false)
  TEST_NOT_EQUAL(std::hash<AnnotatedMSRun>{}(run1), std::hash<AnnotatedMSRun>{}(run3))

  // Test 3: Test use in unordered_set
  std::unordered_set<AnnotatedMSRun> run_set;
  run_set.insert(run1);
  run_set.insert(run2);  // Should not increase size (equal to run1)
  run_set.insert(run3);

  TEST_EQUAL(run_set.size(), 2)
  TEST_EQUAL(run_set.count(run1), 1)
  TEST_EQUAL(run_set.count(run3), 1)

  // Test 4: Test use in unordered_map
  std::unordered_map<AnnotatedMSRun, std::string> run_map;
  run_map[run1] = "first";
  run_map[run3] = "third";

  TEST_EQUAL(run_map.size(), 2)
  TEST_EQUAL(run_map[run1], "first")
  TEST_EQUAL(run_map[run2], "first")  // run2 == run1, so should get same value
  TEST_EQUAL(run_map[run3], "third")

  // Test 5: Empty AnnotatedMSRun has consistent hash
  AnnotatedMSRun empty1, empty2;
  TEST_EQUAL(empty1 == empty2, true)
  TEST_EQUAL(std::hash<AnnotatedMSRun>{}(empty1), std::hash<AnnotatedMSRun>{}(empty2))
END_SECTION

END_TEST
