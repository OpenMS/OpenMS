// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

///////////////////////////

START_TEST(IdentificationData, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

IdentificationData* ptr = 0;
IdentificationData* null = 0;
START_SECTION((IdentificationData()))
  ptr = new IdentificationData();
  TEST_NOT_EQUAL(ptr, null);
END_SECTION

IdentificationData data;
IdentificationData::InputFileRef file_ref;
IdentificationData::ProcessingSoftwareRef sw_ref;
IdentificationData::SearchParamRef param_ref;
IdentificationData::ProcessingStepRef step_ref;
IdentificationData::ScoreTypeRef score_ref;
IdentificationData::DataQueryRef query_ref;
IdentificationData::ParentMoleculeRef protein_ref, rna_ref;
IdentificationData::IdentifiedPeptideRef peptide_ref;
IdentificationData::IdentifiedOligoRef oligo_ref;
IdentificationData::IdentifiedCompoundRef compound_ref;
IdentificationData::QueryMatchRef match_ref1, match_ref2, match_ref3;

START_SECTION((const InputFiles& getInputFiles() const))
{
  TEST_EQUAL(data.getInputFiles().empty(), true);
  // tested further below
}
END_SECTION

START_SECTION((InputFileRef registerInputFile(const InputFile& file)))
{
  IdentificationData::InputFile file("test.mzML");
  file_ref = data.registerInputFile(file);
  TEST_EQUAL(data.getInputFiles().size(), 1);
  TEST_STRING_EQUAL(file_ref->name, file.name);
  // re-registering doesn't lead to redundant entries:
  data.registerInputFile(file);
  TEST_EQUAL(data.getInputFiles().size(), 1);
}
END_SECTION

START_SECTION((const DataProcessingSoftwares& getDataProcessingSoftwares() const))
{
  TEST_EQUAL(data.getDataProcessingSoftwares().empty(), true);
  // tested further below
}
END_SECTION

START_SECTION((ProcessingSoftwareRef registerDataProcessingSoftware(const Software& software)))
{
  IdentificationData::DataProcessingSoftware sw("Tool", "1.0");
  sw_ref = data.registerDataProcessingSoftware(sw);
  TEST_EQUAL(data.getDataProcessingSoftwares().size(), 1);
  TEST_EQUAL(*sw_ref == sw, true); // "TEST_EQUAL(*sw_ref, sw)" doesn't compile - same below
  // re-registering doesn't lead to redundant entries:
  data.registerDataProcessingSoftware(sw);
  TEST_EQUAL(data.getDataProcessingSoftwares().size(), 1);
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
  IdentificationData::DBSearchParam param;
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

START_SECTION((const DataProcessingSteps& getDataProcessingSteps() const))
{
  TEST_EQUAL(data.getDataProcessingSteps().empty(), true);
  // tested further below
}
END_SECTION

START_SECTION((ProcessingStepRef registerDataProcessingStep(const DataProcessingStep& step)))
{
  vector<IdentificationData::InputFileRef> file_refs(1, file_ref);
  IdentificationData::DataProcessingStep step(sw_ref, file_refs);
  step_ref = data.registerDataProcessingStep(step);
  TEST_EQUAL(data.getDataProcessingSteps().size(), 1);
  TEST_EQUAL(*step_ref == step, true);
  // re-registering doesn't lead to redundant entries:
  data.registerDataProcessingStep(step);
  TEST_EQUAL(data.getDataProcessingSteps().size(), 1);
}
END_SECTION

START_SECTION((const DataProcessingSteps& getDBSearchSteps() const))
{
  TEST_EQUAL(data.getDBSearchSteps().empty(), true);
  // tested further below
}
END_SECTION

START_SECTION((ProcessingStepRef registerDataProcessingStep(const DataProcessingStep& step, SearchParamRef search_ref)))
{
  IdentificationData::DataProcessingStep step(sw_ref);
  step_ref = data.registerDataProcessingStep(step, param_ref);
  TEST_EQUAL(data.getDataProcessingSteps().size(), 2);
  TEST_EQUAL(*step_ref == step, true);
  TEST_EQUAL(data.getDBSearchSteps().size(), 1);
  TEST_EQUAL(data.getDBSearchSteps().at(step_ref), param_ref);
  // re-registering doesn't lead to redundant entries:
  data.registerDataProcessingStep(step, param_ref);
  TEST_EQUAL(data.getDataProcessingSteps().size(), 2);
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
  IdentificationData::ScoreType score("test_score", true);
  score_ref = data.registerScoreType(score);
  TEST_EQUAL(data.getScoreTypes().size(), 1);
  TEST_EQUAL(*score_ref == score, true);
  // re-registering doesn't lead to redundant entries:
  data.registerScoreType(score);
  TEST_EQUAL(data.getScoreTypes().size(), 1);
}
END_SECTION

START_SECTION((const DataQueries& getDataQueries() const))
{
  TEST_EQUAL(data.getDataQueries().empty(), true);
  // tested further below
}
END_SECTION

START_SECTION((DataQueryRef registerDataQuery(const DataQuery& query)))
{
  IdentificationData::DataQuery query("spectrum_1", file_ref, 100.0, 1000.0);
  query_ref = data.registerDataQuery(query);
  TEST_EQUAL(data.getDataQueries().size(), 1);
  TEST_EQUAL(*query_ref == query, true);
  // re-registering doesn't lead to redundant entries:
  data.registerDataQuery(query);
  TEST_EQUAL(data.getDataQueries().size(), 1);
}
END_SECTION

START_SECTION((const ParentMolecules& getParentMolecules() const))
{
  TEST_EQUAL(data.getParentMolecules().empty(), true);
  // tested further below
}
END_SECTION

START_SECTION((ParentMoleculeRef registerParentMolecule(const ParentMolecule& parent)))
{
  IdentificationData::ParentMolecule protein("");
  // can't register a parent molecule without accession:
  TEST_EXCEPTION(Exception::IllegalArgument,
                 data.registerParentMolecule(protein));
  TEST_EQUAL(data.getParentMolecules().empty(), true);

  protein.accession = "protein_1";
  protein.sequence = "TESTPEPTIDEAAA";
  protein_ref = data.registerParentMolecule(protein);
  TEST_EQUAL(data.getParentMolecules().size(), 1);
  TEST_EQUAL(*protein_ref == protein, true);

  IdentificationData::ParentMolecule rna("rna_1",
                                         IdentificationData::MoleculeType::RNA);
  rna_ref = data.registerParentMolecule(rna);
  TEST_EQUAL(data.getParentMolecules().size(), 2);
  TEST_EQUAL(*rna_ref == rna, true);
  // re-registering doesn't lead to redundant entries:
  data.registerParentMolecule(rna);
  TEST_EQUAL(data.getParentMolecules().size(), 2);
}
END_SECTION

START_SECTION((const ParentMoleculeGroupings& getParentMoleculeGroupings() const))
{
  TEST_EQUAL(data.getParentMoleculeGroupings().empty(), true);
  // tested further below
}
END_SECTION

START_SECTION((void registerParentMoleculeGrouping(const ParentMoleculeGrouping& grouping)))
{
  IdentificationData::ParentMoleculeGroup group;
  group.parent_molecule_refs.insert(protein_ref);
  group.parent_molecule_refs.insert(rna_ref);
  IdentificationData::ParentMoleculeGrouping grouping;
  grouping.label = "test_grouping";
  grouping.groups.insert(group);
  data.registerParentMoleculeGrouping(grouping);
  TEST_EQUAL(data.getParentMoleculeGroupings().size(), 1);
  TEST_EQUAL(data.getParentMoleculeGroupings()[0].groups.size(), 1);
  TEST_EQUAL(data.getParentMoleculeGroupings()[0].groups.begin()->parent_molecule_refs.size(), 2);
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
  IdentificationData::IdentifiedPeptide peptide(AASequence::fromString(""));
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
  peptide.parent_matches[protein_ref].insert(IdentificationData::
                                             MoleculeParentMatch(4, 10));
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
  IdentificationData::IdentifiedOligo oligo(NASequence::fromString(""));
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
  IdentificationData::IdentifiedCompound compound("");
  // can't register a compound without identifier:
  TEST_EXCEPTION(Exception::IllegalArgument,
                 data.registerIdentifiedCompound(compound));
  TEST_EQUAL(data.getIdentifiedCompounds().empty(), true);

  compound = IdentificationData::IdentifiedCompound("compound_1",
                                                    EmpiricalFormula("C2H5OH"),
                                                    "ethanol");
  compound_ref = data.registerIdentifiedCompound(compound);
  TEST_EQUAL(data.getIdentifiedCompounds().size(), 1);
  TEST_EQUAL(*compound_ref == compound, true);

  // re-registering doesn't lead to redundant entries:
  data.registerIdentifiedCompound(compound);
  TEST_EQUAL(data.getIdentifiedCompounds().size(), 1);
}
END_SECTION

START_SECTION((const MoleculeQueryMatches& getMoleculeQueryMatches() const))
{
  TEST_EQUAL(data.getMoleculeQueryMatches().empty(), true);
  // tested further below
}
END_SECTION

START_SECTION((QueryMatchRef registerMoleculeQueryMatch(const MoleculeQueryMatch& match)))
{
  // match with a peptide:
  IdentificationData::MoleculeQueryMatch match(peptide_ref, query_ref, 3);
  match_ref1 = data.registerMoleculeQueryMatch(match);
  TEST_EQUAL(data.getMoleculeQueryMatches().size(), 1);
  TEST_EQUAL(*match_ref1 == match, true);

  // match with an oligo:
  match = IdentificationData::MoleculeQueryMatch(oligo_ref, query_ref, 2);
  match_ref2 = data.registerMoleculeQueryMatch(match);
  TEST_EQUAL(data.getMoleculeQueryMatches().size(), 2);
  TEST_EQUAL(*match_ref2 == match, true);

  // match with a compound:
  match = IdentificationData::MoleculeQueryMatch(compound_ref, query_ref, 1);
  match_ref3 = data.registerMoleculeQueryMatch(match);
  TEST_EQUAL(data.getMoleculeQueryMatches().size(), 3);
  TEST_EQUAL(*match_ref3 == match, true);

  // re-registering doesn't lead to redundant entries:
  data.registerMoleculeQueryMatch(match);
  TEST_EQUAL(data.getMoleculeQueryMatches().size(), 3);
}
END_SECTION

START_SECTION((const QueryMatchGroups& getQueryMatchGroups() const))
{
  TEST_EQUAL(data.getQueryMatchGroups().empty(), true);
  // tested further below
}
END_SECTION

START_SECTION((MatchGroupRef registerQueryMatchGroup(const QueryMatchGroup& group)))
{
  IdentificationData::QueryMatchGroup group;
  group.query_match_refs.insert(match_ref1);
  group.query_match_refs.insert(match_ref2);
  group.query_match_refs.insert(match_ref3);

  data.registerQueryMatchGroup(group);
  TEST_EQUAL(data.getQueryMatchGroups().size(), 1);
  TEST_EQUAL(*data.getQueryMatchGroups().begin() == group, true);
}
END_SECTION

START_SECTION((void addScore(QueryMatchRef match_ref, ScoreTypeRef score_ref, double value)))
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
  TEST_EQUAL(data.getCurrentProcessingStep() == data.getDataProcessingSteps().end(), true);
  // tested further below
}
END_SECTION

START_SECTION((void setCurrentProcessingStep(ProcessingStepRef step_ref)))
{
  data.setCurrentProcessingStep(step_ref);
  TEST_EQUAL(data.getCurrentProcessingStep() == step_ref, true);
  // registering new data automatically adds the processing step:
  IdentificationData::IdentifiedPeptide peptide(AASequence::fromString("EDIT"));
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
  TEST_EQUAL(data.getCurrentProcessingStep() == data.getDataProcessingSteps().end(), true);
}
END_SECTION

START_SECTION((vector<QueryMatchRef> getBestMatchPerQuery(ScoreTypeRef score_ref) const))
{
  vector<IdentificationData::QueryMatchRef> result = data.getBestMatchPerQuery(score_ref);
  TEST_EQUAL(result.size(), 1);
  TEST_EQUAL(result[0] == match_ref2, true);
}
END_SECTION

START_SECTION((pair<ScoreTypeRef, bool> findScoreType(const String& score_name) const))
{
  // non-existent score:
  TEST_EQUAL(data.findScoreType("fake_score").second, false);
  // registered score:
  auto result = data.findScoreType("test_score");
  TEST_EQUAL(result.first == score_ref, true);
  TEST_EQUAL(result.second, true);
}
END_SECTION

START_SECTION((void calculateCoverages(bool check_molecule_length = false)))
{
  TEST_EQUAL(protein_ref->coverage, 0.0);
  data.calculateCoverages();
  TEST_REAL_SIMILAR(protein_ref->coverage, 0.5);
  // partially overlapping peptide:
  IdentificationData::IdentifiedPeptide peptide(AASequence::
                                                fromString("TESTPEP"));
  peptide.parent_matches[protein_ref].insert(IdentificationData::
                                             MoleculeParentMatch(0, 6));
  data.registerIdentifiedPeptide(peptide);
  data.calculateCoverages();
  TEST_REAL_SIMILAR(protein_ref->coverage, 11.0/14.0);
}
END_SECTION

START_SECTION((void cleanup(bool require_query_match = true, bool require_identified_sequence = true, bool require_parent_match = true, bool require_parent_group = false, bool require_match_group = false)))
{
  TEST_EQUAL(data.getIdentifiedPeptides().size(), 4);
  TEST_EQUAL(data.getIdentifiedOligos().size(), 2);
  data.cleanup(false);
  // identified peptide/oligo without parent match is removed:
  TEST_EQUAL(data.getIdentifiedPeptides().size(), 3);
  TEST_EQUAL(data.getIdentifiedOligos().size(), 1);
  data.cleanup();
  // identified peptides without query matches are removed:
  TEST_EQUAL(data.getIdentifiedPeptides().size(), 1);
  TEST_EQUAL(data.getIdentifiedOligos().size(), 1);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
