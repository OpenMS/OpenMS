// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/SYSTEM/SysInfo.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>

#include <type_traits> // to check if movable

///////////////////////////

START_TEST(IdentificationData, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

using ID = IdentificationData;

IdentificationData* ptr = nullptr;
IdentificationData* null = nullptr;
START_SECTION((IdentificationData()))
  ptr = new IdentificationData();
  TEST_NOT_EQUAL(ptr, null);
END_SECTION

START_SECTION((movable))
  TEST_TRUE(std::is_nothrow_move_constructible_v<IdentificationData>);
  TEST_TRUE(std::is_nothrow_move_assignable_v<IdentificationData>);  
END_SECTION


START_SECTION((~IdentificationData()))
  delete ptr;
END_SECTION

IdentificationData data;
ID::InputFileRef file_ref;
ID::ProcessingSoftwareRef sw_ref;
ID::SearchParamRef param_ref;
ID::ProcessingStepRef step_ref;
ID::ScoreTypeRef score_ref;
ID::ObservationRef obs_ref;
ID::ParentSequenceRef protein_ref, rna_ref;
ID::IdentifiedPeptideRef peptide_ref;
ID::IdentifiedOligoRef oligo_ref;
ID::IdentifiedCompoundRef compound_ref;
ID::AdductRef adduct_ref;
ID::ObservationMatchRef match_ref1, match_ref2, match_ref3;

START_SECTION((const InputFiles& getInputFiles() const))
{
  TEST_EQUAL(data.getInputFiles().empty(), true);
  // tested further below
}
END_SECTION

START_SECTION((InputFileRef registerInputFile(const InputFile& file)))
{
  ID::InputFile file("test.mzML");
  file_ref = data.registerInputFile(file);
  TEST_EQUAL(data.getInputFiles().size(), 1);
  TEST_STRING_EQUAL(file_ref->name, file.name);
  // re-registering doesn't lead to redundant entries:
  data.registerInputFile(file);
  TEST_EQUAL(data.getInputFiles().size(), 1);
}
END_SECTION

START_SECTION((const ProcessingSoftwares& getProcessingSoftwares() const))
{
  TEST_EQUAL(data.getProcessingSoftwares().empty(), true);
  // tested further below
}
END_SECTION

START_SECTION((ProcessingSoftwareRef registerProcessingSoftware(const Software& software)))
{
  ID::ProcessingSoftware sw("Tool", "1.0");
  sw_ref = data.registerProcessingSoftware(sw);
  TEST_EQUAL(data.getProcessingSoftwares().size(), 1);
  TEST_EQUAL(*sw_ref == sw, true); // "TEST_EQUAL(*sw_ref, sw)" doesn't compile - same below
  // re-registering doesn't lead to redundant entries:
  data.registerProcessingSoftware(sw);
  TEST_EQUAL(data.getProcessingSoftwares().size(), 1);
}
END_SECTION

START_SECTION((const DBSearchParams& getDBSearchParams() const))
{
  TEST_EQUAL(data.getDBSearchParams().empty(), true);
  // tested further below
}
END_SECTION

START_SECTION((SearchParamRef registerDBSearchParam(const DBSearchParam& param)))
{
  ID::DBSearchParam param;
  param.database = "test-db.fasta";
  param.precursor_mass_tolerance = 1;
  param.fragment_mass_tolerance = 2;
  param_ref = data.registerDBSearchParam(param);
  TEST_EQUAL(data.getDBSearchParams().size(), 1);
  TEST_EQUAL(*param_ref == param, true);
  // re-registering doesn't lead to redundant entries:
  data.registerDBSearchParam(param);
  TEST_EQUAL(data.getDBSearchParams().size(), 1);
}
END_SECTION

START_SECTION((const ProcessingSteps& getProcessingSteps() const))
{
  TEST_EQUAL(data.getProcessingSteps().empty(), true);
  // tested further below
}
END_SECTION

START_SECTION((ProcessingStepRef registerProcessingStep(const ProcessingStep& step)))
{
  vector<ID::InputFileRef> file_refs(1, file_ref);
  ID::ProcessingStep step(sw_ref, file_refs);
  step_ref = data.registerProcessingStep(step);
  TEST_EQUAL(data.getProcessingSteps().size(), 1);
  TEST_EQUAL(*step_ref == step, true);
  // re-registering doesn't lead to redundant entries:
  data.registerProcessingStep(step);
  TEST_EQUAL(data.getProcessingSteps().size(), 1);
}
END_SECTION

START_SECTION((const ProcessingSteps& getDBSearchSteps() const))
{
  TEST_EQUAL(data.getDBSearchSteps().empty(), true);
  // tested further below
}
END_SECTION

START_SECTION((ProcessingStepRef registerProcessingStep(const ProcessingStep& step, SearchParamRef search_ref)))
{
  ID::ProcessingStep step(sw_ref);
  step_ref = data.registerProcessingStep(step, param_ref);
  TEST_EQUAL(data.getProcessingSteps().size(), 2);
  TEST_EQUAL(*step_ref == step, true);
  TEST_EQUAL(data.getDBSearchSteps().size(), 1);
  TEST_EQUAL(data.getDBSearchSteps().at(step_ref), param_ref);
  // re-registering doesn't lead to redundant entries:
  data.registerProcessingStep(step, param_ref);
  TEST_EQUAL(data.getProcessingSteps().size(), 2);
  TEST_EQUAL(data.getDBSearchSteps().size(), 1);
}
END_SECTION

START_SECTION((const ScoreTypes& getScoreTypes() const))
{
  TEST_EQUAL(data.getScoreTypes().empty(), true);
  // tested further below
}
END_SECTION

START_SECTION((ScoreTypeRef registerScoreType(const ScoreType& score)))
{
  ID::ScoreType score("test_score", true);
  score_ref = data.registerScoreType(score);
  TEST_EQUAL(data.getScoreTypes().size(), 1);
  TEST_EQUAL(*score_ref == score, true);
  // re-registering doesn't lead to redundant entries:
  data.registerScoreType(score);
  TEST_EQUAL(data.getScoreTypes().size(), 1);
}
END_SECTION

START_SECTION((const Observations& getObservations() const))
{
  TEST_EQUAL(data.getObservations().empty(), true);
  // tested further below
}
END_SECTION

START_SECTION((ObservationRef registerObservation(const Observation& obs)))
{
  ID::Observation obs("spectrum_1", file_ref, 100.0, 1000.0);
  obs_ref = data.registerObservation(obs);
  TEST_EQUAL(data.getObservations().size(), 1);
  TEST_EQUAL(*obs_ref == obs, true);
  // re-registering doesn't lead to redundant entries:
  data.registerObservation(obs);
  TEST_EQUAL(data.getObservations().size(), 1);
}
END_SECTION

START_SECTION((const ParentSequences& getParentSequences() const))
{
  TEST_EQUAL(data.getParentSequences().empty(), true);
  // tested further below
}
END_SECTION

START_SECTION((ParentSequenceRef registerParentSequence(const ParentSequence& parent)))
{
  ID::ParentSequence protein("");
  // can't register a parent sequence without accession:
  TEST_EXCEPTION(Exception::IllegalArgument,
                 data.registerParentSequence(protein));
  TEST_EQUAL(data.getParentSequences().empty(), true);

  protein.accession = "protein_1";
  protein.sequence = "TESTPEPTIDEAAA";
  protein_ref = data.registerParentSequence(protein);
  TEST_EQUAL(data.getParentSequences().size(), 1);
  TEST_EQUAL(*protein_ref == protein, true);

  ID::ParentSequence rna("rna_1", ID::MoleculeType::RNA);
  rna_ref = data.registerParentSequence(rna);
  TEST_EQUAL(data.getParentSequences().size(), 2);
  TEST_EQUAL(*rna_ref == rna, true);
  // re-registering doesn't lead to redundant entries:
  data.registerParentSequence(rna);
  TEST_EQUAL(data.getParentSequences().size(), 2);
}
END_SECTION

START_SECTION((const ParentGroupSets& getParentGroupSets() const))
{
  TEST_EQUAL(data.getParentGroupSets().empty(), true);
  // tested further below
}
END_SECTION

START_SECTION((void registerParentGroupSet(const ParentGroupSet& groups)))
{
  ID::ParentGroup group;
  group.parent_refs.insert(protein_ref);
  group.parent_refs.insert(rna_ref);
  ID::ParentGroupSet groups;
  groups.label = "test_grouping";
  groups.groups.insert(group);
  data.registerParentGroupSet(groups);
  TEST_EQUAL(data.getParentGroupSets().size(), 1);
  TEST_EQUAL(data.getParentGroupSets()[0].groups.size(), 1);
  TEST_EQUAL(data.getParentGroupSets()[0].groups.begin()->parent_refs.size(), 2);
}
END_SECTION

START_SECTION((const IdentifiedPeptides& getIdentifiedPeptides() const))
{
  TEST_EQUAL(data.getIdentifiedPeptides().empty(), true);
  // tested further below
}
END_SECTION

START_SECTION((IdentifiedPeptideRef registerIdentifiedPeptide(const IdentifiedPeptide& peptide)))
{
  ID::IdentifiedPeptide peptide(AASequence::fromString(""));
  // can't register a peptide without a sequence:
  TEST_EXCEPTION(Exception::IllegalArgument,
                 data.registerIdentifiedPeptide(peptide));
  TEST_EQUAL(data.getIdentifiedPeptides().empty(), true);

  // peptide without protein reference:
  peptide.sequence = AASequence::fromString("TEST");
  peptide_ref = data.registerIdentifiedPeptide(peptide);
  TEST_EQUAL(data.getIdentifiedPeptides().size(), 1);
  TEST_EQUAL(*peptide_ref == peptide, true);

  // peptide with protein reference:
  peptide.sequence = AASequence::fromString("PEPTIDE");
  peptide.parent_matches[protein_ref].insert(ID::
                                             ParentMatch(4, 10));
  peptide_ref = data.registerIdentifiedPeptide(peptide);
  TEST_EQUAL(data.getIdentifiedPeptides().size(), 2);
  TEST_EQUAL(*peptide_ref == peptide, true);

  // re-registering doesn't lead to redundant entries:
  data.registerIdentifiedPeptide(peptide);
  TEST_EQUAL(data.getIdentifiedPeptides().size(), 2);

  // registering a peptide with RNA reference doesn't work:
  peptide.parent_matches[rna_ref];
  TEST_EXCEPTION(Exception::IllegalArgument,
                 data.registerIdentifiedPeptide(peptide));
}
END_SECTION

START_SECTION((const IdentifiedOligos& getIdentifiedOligos() const))
{
  TEST_EQUAL(data.getIdentifiedOligos().empty(), true);
  // tested further below
}
END_SECTION

START_SECTION((IdentifiedOligoRef registerIdentifiedOligo(const IdentifiedOligo& oligo)))
{
  ID::IdentifiedOligo oligo(NASequence::fromString(""));
  // can't register an oligo without a sequence:
  TEST_EXCEPTION(Exception::IllegalArgument,
                 data.registerIdentifiedOligo(oligo));
  TEST_EQUAL(data.getIdentifiedOligos().empty(), true);

  // oligo without RNA reference:
  oligo.sequence = NASequence::fromString("ACGU");
  oligo_ref = data.registerIdentifiedOligo(oligo);
  TEST_EQUAL(data.getIdentifiedOligos().size(), 1);
  TEST_EQUAL(*oligo_ref == oligo, true);

  // oligo with RNA reference:
  oligo.sequence = NASequence::fromString("UGCA");
  oligo.parent_matches[rna_ref];
  oligo_ref = data.registerIdentifiedOligo(oligo);
  TEST_EQUAL(data.getIdentifiedOligos().size(), 2);
  TEST_EQUAL(*oligo_ref == oligo, true);

  // re-registering doesn't lead to redundant entries:
  data.registerIdentifiedOligo(oligo);
  TEST_EQUAL(data.getIdentifiedOligos().size(), 2);

  // registering an oligo with protein reference doesn't work:
  oligo.parent_matches[protein_ref];
  TEST_EXCEPTION(Exception::IllegalArgument,
                 data.registerIdentifiedOligo(oligo));
}
END_SECTION

START_SECTION((const IdentifiedCompounds& getIdentifiedCompounds() const))
{
  TEST_EQUAL(data.getIdentifiedCompounds().empty(), true);
  // tested further below
}
END_SECTION

START_SECTION((IdentifiedCompoundRef registerIdentifiedCompound(const IdentifiedCompound& compound)))
{
  ID::IdentifiedCompound compound("");
  // can't register a compound without identifier:
  TEST_EXCEPTION(Exception::IllegalArgument,
                 data.registerIdentifiedCompound(compound));
  TEST_EQUAL(data.getIdentifiedCompounds().empty(), true);

  compound = ID::IdentifiedCompound("compound_1", EmpiricalFormula("C2H5OH"),
                                    "ethanol");
  compound_ref = data.registerIdentifiedCompound(compound);
  TEST_EQUAL(data.getIdentifiedCompounds().size(), 1);
  TEST_EQUAL(*compound_ref == compound, true);

  // re-registering doesn't lead to redundant entries:
  data.registerIdentifiedCompound(compound);
  TEST_EQUAL(data.getIdentifiedCompounds().size(), 1);
}
END_SECTION

START_SECTION((const Adducts& getAdducts() const))
{
  TEST_EQUAL(data.getAdducts().empty(), true);
  // tested further below
}
END_SECTION

START_SECTION((AdductRef registerAdduct(const AdductInfo& adduct)))
{
  AdductInfo adduct("Na+", EmpiricalFormula("Na"), 1);
  adduct_ref = data.registerAdduct(adduct);
  TEST_EQUAL(data.getAdducts().size(), 1);
  TEST_EQUAL(*adduct_ref == adduct, true);
}
END_SECTION

START_SECTION((const ObservationMatches& getObservationMatches() const))
{
  TEST_EQUAL(data.getObservationMatches().empty(), true);
  // tested further below
}
END_SECTION

START_SECTION((ObservationMatchRef registerObservationMatch(const ObservationMatch& match)))
{
  // match with a peptide:
  ID::ObservationMatch match(peptide_ref, obs_ref, 3);
  match_ref1 = data.registerObservationMatch(match);
  TEST_EQUAL(data.getObservationMatches().size(), 1);
  TEST_EQUAL(*match_ref1 == match, true);

  // match with an oligo (+ adduct):
  match = ID::ObservationMatch(oligo_ref, obs_ref, 2, adduct_ref);
  match_ref2 = data.registerObservationMatch(match);
  TEST_EQUAL(data.getObservationMatches().size(), 2);
  TEST_EQUAL(*match_ref2 == match, true);
  TEST_EQUAL((*match_ref2->adduct_opt)->getName(), "Na+");

  // match with a compound:
  match = ID::ObservationMatch(compound_ref, obs_ref, 1);
  match_ref3 = data.registerObservationMatch(match);
  TEST_EQUAL(data.getObservationMatches().size(), 3);
  TEST_EQUAL(*match_ref3 == match, true);

  // re-registering doesn't lead to redundant entries:
  data.registerObservationMatch(match);
  TEST_EQUAL(data.getObservationMatches().size(), 3);
}
END_SECTION

START_SECTION((const ObservationMatchGroups& getObservationMatchGroups() const))
{
  TEST_EQUAL(data.getObservationMatchGroups().empty(), true);
  // tested further below
}
END_SECTION

START_SECTION((MatchGroupRef registerObservationMatchGroup(const ObservationMatchGroup& group)))
{
  ID::ObservationMatchGroup group;
  group.observation_match_refs.insert(match_ref1);
  group.observation_match_refs.insert(match_ref2);
  group.observation_match_refs.insert(match_ref3);

  data.registerObservationMatchGroup(group);
  TEST_EQUAL(data.getObservationMatchGroups().size(), 1);
  TEST_EQUAL(*data.getObservationMatchGroups().begin() == group, true);
}
END_SECTION

START_SECTION((void addScore(ObservationMatchRef match_ref, ScoreTypeRef score_ref, double value)))
{
  TEST_EQUAL(match_ref1->steps_and_scores.empty(), true);
  data.addScore(match_ref1, score_ref, 100.0);
  TEST_EQUAL(match_ref1->steps_and_scores.size(), 1);
  TEST_EQUAL(match_ref1->steps_and_scores.back().scores.begin()->first,
             score_ref);
  TEST_EQUAL(match_ref1->steps_and_scores.back().scores.begin()->second, 100.0);
  TEST_EQUAL(match_ref2->steps_and_scores.empty(), true);
  data.addScore(match_ref2, score_ref, 200.0);
  TEST_EQUAL(match_ref2->steps_and_scores.size(), 1);
  TEST_EQUAL(match_ref2->steps_and_scores.back().scores.begin()->first,
             score_ref);
  TEST_EQUAL(match_ref2->steps_and_scores.back().scores.begin()->second, 200.0);
}
END_SECTION

START_SECTION((ProcessingStepRef getCurrentProcessingStep()))
{
  TEST_EQUAL(data.getCurrentProcessingStep() == data.getProcessingSteps().end(), true);
  // tested further below
}
END_SECTION

START_SECTION((void setCurrentProcessingStep(ProcessingStepRef step_ref)))
{
  data.setCurrentProcessingStep(step_ref);
  TEST_EQUAL(data.getCurrentProcessingStep() == step_ref, true);
  // registering new data automatically adds the processing step:
  ID::IdentifiedPeptide peptide(AASequence::fromString("EDIT"));
  peptide.parent_matches[protein_ref];
  peptide_ref = data.registerIdentifiedPeptide(peptide);
  TEST_EQUAL(peptide_ref->steps_and_scores.size(), 1);
  TEST_EQUAL(peptide_ref->steps_and_scores.front().processing_step_opt ==
             step_ref, true);
}
END_SECTION

START_SECTION((void clearCurrentProcessingStep()))
{
  data.clearCurrentProcessingStep();
  TEST_EQUAL(data.getCurrentProcessingStep() == data.getProcessingSteps().end(), true);
}
END_SECTION

START_SECTION((pair<ID::ObservationMatchRef, ID::ObservationMatchRef> getMatchesForObservation(ObservationRef obs_ref) const))
{
  pair<ID::ObservationMatchRef, ID::ObservationMatchRef> result =
    data.getMatchesForObservation(obs_ref);
  TEST_EQUAL(distance(result.first, result.second), 3);
  for (; result.first != result.second; ++result.first)
  {
    TEST_EQUAL((result.first == match_ref1) || (result.first == match_ref2) ||
               (result.first == match_ref3), true);
  }
}
END_SECTION

START_SECTION((ScoreTypeRef findScoreType(const String& score_name) const))
{
  // non-existent score:
  TEST_EQUAL(data.findScoreType("fake_score") == data.getScoreTypes().end(), true);
  // registered score:
  TEST_EQUAL(data.findScoreType("test_score") == score_ref, true);
}
END_SECTION

START_SECTION((void calculateCoverages(bool check_molecule_length = false)))
{
  TEST_EQUAL(protein_ref->coverage, 0.0);
  data.calculateCoverages();
  TEST_REAL_SIMILAR(protein_ref->coverage, 0.5);
  // partially overlapping peptide:
  ID::IdentifiedPeptide peptide(AASequence::fromString("TESTPEP"));
  peptide.parent_matches[protein_ref].insert(ID::ParentMatch(0, 6));
  data.registerIdentifiedPeptide(peptide);
  data.calculateCoverages();
  TEST_REAL_SIMILAR(protein_ref->coverage, 11.0/14.0);
}
END_SECTION

START_SECTION((void cleanup(bool require_observation_match = true, bool require_identified_sequence = true, bool require_parent_match = true, bool require_parent_group = false, bool require_match_group = false)))
{
  TEST_EQUAL(data.getIdentifiedPeptides().size(), 4);
  TEST_EQUAL(data.getIdentifiedOligos().size(), 2);
  data.cleanup(false);
  // identified peptide/oligo without parent match is removed:
  TEST_EQUAL(data.getIdentifiedPeptides().size(), 3);
  TEST_EQUAL(data.getIdentifiedOligos().size(), 1);
  data.cleanup();
  // identified peptides without matches are removed:
  TEST_EQUAL(data.getIdentifiedPeptides().size(), 1);
  TEST_EQUAL(data.getIdentifiedOligos().size(), 1);
}
END_SECTION

START_SECTION((ProcessingStepRef merge(const IdentificationData& other)))
{
  TEST_EQUAL(data.getIdentifiedPeptides().size(), 1);
  TEST_EQUAL(data.getIdentifiedOligos().size(), 1);
  TEST_EQUAL(data.getParentSequences().size(), 2);
  data.merge(data); // self-merge shouldn't change anything
  TEST_EQUAL(data.getIdentifiedPeptides().size(), 1);
  TEST_EQUAL(data.getIdentifiedOligos().size(), 1);
  TEST_EQUAL(data.getParentSequences().size(), 2);
  IdentificationData other;
  ID::IdentifiedPeptide peptide(AASequence::fromString("MASSSPEC"));
  other.registerIdentifiedPeptide(peptide);
  data.merge(other);
  TEST_EQUAL(data.getIdentifiedPeptides().size(), 2);
  TEST_EQUAL(data.getIdentifiedOligos().size(), 1);
  TEST_EQUAL(data.getParentSequences().size(), 2);
}
END_SECTION

START_SECTION((IdentificationData(const IdentificationData& other)))
{
  IdentificationData copy(data);
  TEST_EQUAL(copy.getIdentifiedPeptides().size(), 2);
  TEST_EQUAL(copy.getIdentifiedOligos().size(), 1);
  TEST_EQUAL(copy.getParentSequences().size(), 2);
  TEST_EQUAL(copy.getObservationMatches().size(), 3);
  // focus on processing steps and scores for observation matches:
  IdentificationData data2;
  ID::InputFile file("test.mzML");
  auto file_ref = data2.registerInputFile(file);
  ID::ProcessingSoftware sw("Tool", "1.0");
  auto sw_ref = data2.registerProcessingSoftware(sw);
  ID::ProcessingStep step(sw_ref, {file_ref});
  auto step_ref = data2.registerProcessingStep(step);
  data2.setCurrentProcessingStep(step_ref);
  ID::Observation obs("spectrum_1", file_ref, 100.0, 1000.0);
  auto obs_ref = data2.registerObservation(obs);
  ID::IdentifiedPeptide peptide(AASequence::fromString("PEPTIDE"));
  auto pep_ref = data2.registerIdentifiedPeptide(peptide);
  ID::ObservationMatch match(pep_ref, obs_ref, 2);
  ID::ScoreType score("score1", true);
  auto score_ref1 = data2.registerScoreType(score);
  score = ID::ScoreType("score2", false);
  auto score_ref2 = data2.registerScoreType(score);
  // add first score, not connected to a processing step:
  match.addScore(score_ref1, 1.0);
  auto match_ref = data2.registerObservationMatch(match);
  // add second score, automatically connected to last processing step:
  data2.addScore(match_ref, score_ref2, 2.0);
  TEST_EQUAL(data2.getObservationMatches().begin()->steps_and_scores.size(), 2);
  TEST_EQUAL(data2.getObservationMatches().begin()->getNumberOfScores(), 2);
  // look up scores by score type:
  TEST_EQUAL(data2.getObservationMatches().begin()->getScore(score_ref1).first, 1.0);
  TEST_EQUAL(data2.getObservationMatches().begin()->getScore(score_ref2).first, 2.0);
  // look up score by score type and (wrong) processing step -> fails:
  TEST_EQUAL(data2.getObservationMatches().begin()->getScore(score_ref1, step_ref).second, false);
  // look up score by score type and (correct) processing step -> succeeds:
  TEST_EQUAL(data2.getObservationMatches().begin()->getScore(score_ref2, step_ref).first, 2.0);
  auto triple = data2.getObservationMatches().begin()->getMostRecentScore();
  TEST_EQUAL(std::get<0>(triple), 2.0);
  TEST_EQUAL(std::get<1>(triple) == score_ref2, true);
  // after copying:
  IdentificationData copy2(data2);
  TEST_EQUAL(copy2.getObservationMatches().begin()->steps_and_scores.size(), 2);
  TEST_EQUAL(copy2.getObservationMatches().begin()->getNumberOfScores(), 2);
  score_ref1 = copy2.findScoreType("score1");
  TEST_EQUAL(copy2.getObservationMatches().begin()->getScore(score_ref1).first, 1.0);
  score_ref2 = copy2.findScoreType("score2");
  TEST_EQUAL(copy2.getObservationMatches().begin()->getScore(score_ref2).first, 2.0);
  step_ref = copy2.getCurrentProcessingStep();
  TEST_EQUAL(copy2.getObservationMatches().begin()->getScore(score_ref1, step_ref).second, false);
  TEST_EQUAL(copy2.getObservationMatches().begin()->getScore(score_ref2, step_ref).first, 2.0);
  triple = copy2.getObservationMatches().begin()->getMostRecentScore();
  TEST_EQUAL(std::get<0>(triple), 2.0);
  TEST_EQUAL(std::get<1>(triple) == score_ref2, true);
}
END_SECTION

START_SECTION((vector<ObservationMatchRef> getBestMatchPerObservation(ScoreTypeRef score_ref) const))
{
  // add a second observation and match (without score):
  ID::Observation obs("spectrum_2", file_ref, 200.0, 2000.0);
  ID::ObservationRef obs_ref2 = data.registerObservation(obs);
  ID::ObservationMatch match(oligo_ref, obs_ref2, 2);
  ID::ObservationMatchRef match_ref4 = data.registerObservationMatch(match);
  TEST_EQUAL(data.getObservationMatches().size(), 4);
  // best matches, requiring score:
  vector<ID::ObservationMatchRef> results = data.getBestMatchPerObservation(score_ref, true);
  TEST_EQUAL(results.size(), 1);
  TEST_EQUAL(results[0] == match_ref2, true);
  // best matches, no score required:
  results = data.getBestMatchPerObservation(score_ref, false);
  TEST_EQUAL(results.size(), 2);
  ABORT_IF(results.size() != 2);
  if (results[0] == match_ref2) // can't be sure about the order
  {
    TEST_EQUAL(results[1] == match_ref4, true);
  }
  else
  {
    TEST_EQUAL(results[0] == match_ref4, true);
    TEST_EQUAL(results[1] == match_ref2, true);
  }
}
END_SECTION

START_SECTION(([EXTRA] UseCaseBuildBottomUpProteomicsID()))
{
  IdentificationData id;

  ID::InputFile file("file://ROOT/FOLDER/SPECTRA.mzML");
  auto file_ref = id.registerInputFile(file);

  // register a score type
  ID::ScoreType score("MySearchEngineScore", true);
  auto score_ref = id.registerScoreType(score);

  // register software (connected to score)
  ID::ProcessingSoftware sw("MySearchEngineTool", "1.0");
  sw.assigned_scores.push_back(score_ref);
  auto sw_ref = id.registerProcessingSoftware(sw);

  // all supported search settings
  ID::DBSearchParam search_param;
  search_param.database = "file://ROOT/FOLDER/DATABASE.fasta";
  search_param.database_version = "nextprot1234";
  search_param.taxonomy = "Homo Sapiens";
  search_param.charges = {2,3,4,5};
  search_param.precursor_mass_tolerance = 8.0;
  search_param.precursor_tolerance_ppm = true;
  search_param.fixed_mods = {"Carbamidomethyl (C)"};
  search_param.variable_mods = {"Oxidation (M)"};
  search_param.digestion_enzyme = ProteaseDB::getInstance()->getEnzyme("Trypsin");
  search_param.enzyme_term_specificity = EnzymaticDigestion::SPEC_SEMI;
  search_param.missed_cleavages = 2;
  search_param.min_length = 6;
  search_param.max_length = 40;
  search_param.fragment_mass_tolerance = 0.3;
  search_param.fragment_tolerance_ppm = true;
  auto search_param_ref = id.registerDBSearchParam(search_param);

  // file has been processed by software
  ID::ProcessingStep step(sw_ref);
  step.input_file_refs.push_back(file_ref);
  auto step_ref = id.registerProcessingStep(step, search_param_ref);
  // all further data comes from this processing step
  id.setCurrentProcessingStep(step_ref);

  // register spectrum
  ID::Observation obs("spectrum_1", file_ref, 100.0, 1000.0);
  auto obs_ref = id.registerObservation(obs);

  // peptide without protein reference (yet)
  ID::IdentifiedPeptide peptide(AASequence::fromString("TESTPEPTIDR")); // seq. is required
  auto peptide_ref = id.registerIdentifiedPeptide(peptide);
  TEST_EQUAL(peptide_ref->parent_matches.size(), 0);

  // peptide-spectrum match
  ID::ObservationMatch match(peptide_ref, obs_ref); // both refs. are required
  match.addScore(score_ref, 123, step_ref);
  id.registerObservationMatch(match);

  // some calculations, inference etc. could take place ...
  ID::ParentSequence protein("protein_1"); // accession is required
  protein.sequence = "PRTTESTPEPTIDRPRT";
  protein.description = "Human Random Protein 1";
  auto protein_ref = id.registerParentSequence(protein);

  // add reference to parent (protein) and update peptide
  ID::IdentifiedPeptide augmented_pep = *peptide_ref;
  // @TODO: wrap this in a convenience function (like "match.addScore" above)
  augmented_pep.parent_matches[protein_ref].insert(ID::ParentMatch(3, 13));
  id.registerIdentifiedPeptide(augmented_pep); // protein reference will be added
  // peptide_ref should still be valid and now contain link to protein
  TEST_EQUAL(peptide_ref->sequence, augmented_pep.sequence);
  TEST_EQUAL(peptide_ref->parent_matches.size(), 1);

  // and now update protein coverage of all proteins
  id.calculateCoverages();
  TEST_NOT_EQUAL(protein_ref->coverage, 0.0);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
