// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

#include <OpenMS/METADATA/IdentificationData.h>

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

START_SECTION((const InputFiles& getInputFiles() const))
{
  TEST_EQUAL(data.getInputFiles().empty(), true);
  // tested further below
}
END_SECTION

START_SECTION((InputFileRef registerInputFile(const String& file)))
{
  String file = "test.mzML";
  file_ref = data.registerInputFile(file);
  TEST_EQUAL(data.getInputFiles().size(), 1);
  TEST_STRING_EQUAL(*file_ref, file);
}
END_SECTION

START_SECTION((const DataProcessingSoftware& getDataProcessingSoftware() const))
{
  TEST_EQUAL(data.getDataProcessingSoftware().empty(), true);
  // tested further below
}
END_SECTION

START_SECTION((ProcessingSoftwareRef registerDataProcessingSoftware(const Software& software)))
{
  Software sw("Tool", "1.0");
  sw_ref = data.registerDataProcessingSoftware(sw);
  TEST_EQUAL(data.getDataProcessingSoftware().size(), 1);
  TEST_EQUAL(*sw_ref == sw, true); // "TEST_EQUAL(*sw_ref, sw)" doesn't compile
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
  TEST_EQUAL(*param_ref == param, true); // "TEST_EQUAL(*param_ref, param)" doesn't compile
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
  TEST_EQUAL(*step_ref == step, true); // "TEST_EQUAL(*step_ref, step)" doesn't compile
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
  TEST_EQUAL(*step_ref == step, true); // "TEST_EQUAL(*step_ref, step)" doesn't compile
  TEST_EQUAL(data.getDBSearchSteps().size(), 1);
  TEST_EQUAL(data.getDBSearchSteps().at(step_ref), param_ref);
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
  IdentificationData::ScoreType score("test_score", true, sw_ref);
  score_ref = data.registerScoreType(score);
  TEST_EQUAL(data.getScoreTypes().size(), 1);
  TEST_EQUAL(*score_ref == score, true); // "TEST_EQUAL(*score_ref, score)" doesn't compile
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
  TEST_EQUAL(*query_ref == query, true); // "TEST_EQUAL(*query_ref, query)" doesn't compile
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
  IdentificationData::ParentMolecule protein("protein_1");
  protein_ref = data.registerParentMolecule(protein);
  TEST_EQUAL(data.getParentMolecules().size(), 1);
  TEST_EQUAL(*protein_ref == protein, true); // "TEST_EQUAL(*parent_ref, parent)" doesn't compile
  IdentificationData::ParentMolecule rna("rna_1",
                                         IdentificationData::MoleculeType::RNA);
  rna_ref = data.registerParentMolecule(rna);
  TEST_EQUAL(data.getParentMolecules().size(), 2);
  TEST_EQUAL(*rna_ref == rna, true); // "TEST_EQUAL(*parent_ref, parent)" doesn't compile
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

/*
START_SECTION((IdentifiedPeptideRef registerIdentifiedPeptide(const IdentifiedPeptide& peptide)))
{
  IdentificationData::IdentifiedPeptide peptide;
  peptide_ref = data.registerIdentifiedPeptide(peptide);
  TEST_EQUAL(data.getIdentifiedPeptides().size(), 1);
  TEST_EQUAL(*peptide_ref == peptide, true); // "TEST_EQUAL(*peptide_ref, peptide)" doesn't compile
}
END_SECTION
*/

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
