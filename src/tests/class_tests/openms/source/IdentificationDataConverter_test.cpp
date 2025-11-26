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
#include <functional>

///////////////////////////

#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/PepXMLFile.h>

///////////////////////////

using namespace OpenMS;
using namespace std;
using namespace std::placeholders;

struct ComparePIdSize
{
      bool operator()(const ProteinIdentification& lhs, const ProteinIdentification& rhs) const
      {
        return lhs.getHits().size() < rhs.getHits().size();
      }
};

START_TEST(IdentificationDataConverter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((void importIDs(IdentificationData&, const vector<ProteinIdentification>&, const PeptideIdentificationList&)))
{
  vector<ProteinIdentification> proteins_in;
  PeptideIdentificationList peptides_in;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IdXMLFile_whole.idXML"), proteins_in, peptides_in);
  // IdentificationData doesn't allow score types with the same name, but different orientations:
  peptides_in[0].setHigherScoreBetter(true);

  IdentificationData ids;
  IdentificationDataConverter::importIDs(ids, proteins_in, peptides_in);

  vector<ProteinIdentification> proteins_out;
  PeptideIdentificationList peptides_out;
  IdentificationDataConverter::exportIDs(ids, proteins_out, peptides_out);

  TEST_EQUAL(peptides_in.size(), peptides_out.size());
  vector<PeptideHit> hits_in, hits_out;
  for (const auto& pep : peptides_in)
  {
    hits_in.insert(hits_in.end(), pep.getHits().begin(), pep.getHits().end());
  }
  for (const auto& pep : peptides_out)
  {
    hits_out.insert(hits_out.end(), pep.getHits().begin(), pep.getHits().end());
  }
  TEST_EQUAL(hits_in.size(), hits_out.size());
  // order of hits is different, check that every output one is in the input:
  for (const auto& hit : hits_out)
  {
    TEST_EQUAL(find(hits_in.begin(), hits_in.end(), hit) != hits_in.end(),
               true);
  }

  std::sort(proteins_in.begin(), proteins_in.end(), ComparePIdSize());
  std::sort(proteins_out.begin(), proteins_out.end(), ComparePIdSize());
  TEST_EQUAL(proteins_in.size(), proteins_out.size());
  TEST_EQUAL(proteins_in[0].getHits().size(), 1) // is sorted
  TEST_EQUAL(proteins_in[1].getHits().size(), 2) // is sorted

  // the exporter adds target/decoy information (default: target):
  for (auto& hit : proteins_in[0].getHits()) hit.setMetaValue("target_decoy", "target");
  for (auto& hit : proteins_in[1].getHits()) hit.setMetaValue("target_decoy", "target");

  // TEST_EQUAL(proteins_in[0].getIdentifier(), proteins_out[0].getIdentifier() ) // identifiers are not equal
  // TEST_EQUAL(proteins_in[1].getIdentifier(), proteins_out[1].getIdentifier() ) // identifiers are not equal

  TEST_EQUAL(proteins_in[0].getHits().size(), proteins_out[0].getHits().size());
  TEST_EQUAL(proteins_in[1].getHits().size(), proteins_out[1].getHits().size());
  TEST_EQUAL(proteins_in[0].getHits() == proteins_out[0].getHits(), true);
  TEST_EQUAL(proteins_in[1].getHits() == proteins_out[1].getHits(), true);

  TEST_EQUAL(proteins_in[0].getDateTime().get(),
             proteins_out[0].getDateTime().get());
  TEST_EQUAL(proteins_in[1].getDateTime().get(),
             proteins_out[1].getDateTime().get());

  TEST_EQUAL(proteins_in[0].getSearchParameters() == proteins_out[0].getSearchParameters(), true);
  TEST_EQUAL(proteins_in[1].getSearchParameters() == proteins_out[1].getSearchParameters(), true);
  // if something breaks and the search parameters don't match, find where the difference is:
  /*
  for (Size i = 0; i <= 1; ++i)
  {
    TEST_EQUAL(static_cast<MetaInfoInterface>(proteins_in[i].getSearchParameters()) ==
               static_cast<MetaInfoInterface>(proteins_out[i].getSearchParameters()), true);
    TEST_EQUAL(proteins_in[i].getSearchParameters().db, proteins_out[i].getSearchParameters().db);
    TEST_EQUAL(proteins_in[i].getSearchParameters().db_version, proteins_out[i].getSearchParameters().db_version);
    TEST_EQUAL(proteins_in[i].getSearchParameters().taxonomy, proteins_out[i].getSearchParameters().taxonomy);
    TEST_EQUAL(proteins_in[i].getSearchParameters().charges, proteins_out[i].getSearchParameters().charges);
    TEST_EQUAL(proteins_in[i].getSearchParameters().mass_type, proteins_out[i].getSearchParameters().mass_type);
    TEST_EQUAL(proteins_in[i].getSearchParameters().fixed_modifications ==
               proteins_out[i].getSearchParameters().fixed_modifications, true);
    TEST_EQUAL(proteins_in[i].getSearchParameters().variable_modifications ==
               proteins_out[i].getSearchParameters().variable_modifications, true);
    TEST_EQUAL(proteins_in[i].getSearchParameters().missed_cleavages, proteins_out[i].getSearchParameters().missed_cleavages);
    TEST_EQUAL(proteins_in[i].getSearchParameters().fragment_mass_tolerance, proteins_out[i].getSearchParameters().fragment_mass_tolerance);
    TEST_EQUAL(proteins_in[i].getSearchParameters().fragment_mass_tolerance_ppm, proteins_out[i].getSearchParameters().fragment_mass_tolerance_ppm);
    TEST_EQUAL(proteins_in[i].getSearchParameters().precursor_mass_tolerance, proteins_out[i].getSearchParameters().precursor_mass_tolerance);
    TEST_EQUAL(proteins_in[i].getSearchParameters().precursor_mass_tolerance_ppm, proteins_out[i].getSearchParameters().precursor_mass_tolerance_ppm);
    TEST_EQUAL(proteins_in[i].getSearchParameters().digestion_enzyme ==  proteins_out[i].getSearchParameters().digestion_enzyme, true);
  }
  */
  // String filename = OPENMS_GET_TEST_DATA_PATH("IdentificationDataConverter_out.idXML");
  // IdXMLFile().store(filename, proteins_out, peptides_out);
}
END_SECTION

START_SECTION((void importSequences(IdentificationData&, const vector<FASTAFile::FASTAEntry>&, IdentificationData::MoleculeType, const String&)))
{
  vector<FASTAFile::FASTAEntry> fasta;
  FASTAFile().load(OPENMS_GET_TEST_DATA_PATH("FASTAFile_test.fasta"), fasta);
  IdentificationData ids;
  IdentificationDataConverter::importSequences(ids, fasta);
  TEST_EQUAL(ids.getParentSequences().size(), 5);
}
END_SECTION

START_SECTION((void exportIDs(const IdentificationData&, vector<ProteinIdentification>&, PeptideIdentificationList&)))
{
  vector<ProteinIdentification> proteins_in;
  PeptideIdentificationList peptides_in;

  String filename = OPENMS_GET_TEST_DATA_PATH("../../../topp/THIRDPARTY/FidoAdapter_4_output.idXML");
  //String filename = OPENMS_GET_TEST_DATA_PATH("debug_fraction_1_IDs_after_transfer.idXML");
  IdXMLFile().load(filename, proteins_in, peptides_in);

  IdentificationData ids;
  IdentificationDataConverter::importIDs(ids, proteins_in, peptides_in);

  vector<ProteinIdentification> proteins_out;
  PeptideIdentificationList peptides_out;
  IdentificationDataConverter::exportIDs(ids, proteins_out, peptides_out);

  TEST_EQUAL(proteins_in.size(), proteins_out.size());
  TEST_EQUAL(proteins_in[0].getHits().size(),
             proteins_out[0].getHits().size());
  TEST_EQUAL(proteins_in[0].getHits() == proteins_out[0].getHits(), true);

  TEST_EQUAL(proteins_in[0].getIndistinguishableProteins() ==
             proteins_out[0].getIndistinguishableProteins(), true);
  TEST_EQUAL(proteins_in[0].getProteinGroups() ==
             proteins_out[0].getProteinGroups(), true);

  TEST_EQUAL(peptides_in.size(), peptides_out.size());
  // no "operator<" for PeptideHit, otherwise we could use a set:
  vector<PeptideHit> hits_in, hits_out;
  for (const auto& pep : peptides_in)
  {
    hits_in.insert(hits_in.end(), pep.getHits().begin(), pep.getHits().end());
  }
  for (const auto& pep : peptides_out)
  {
    hits_out.insert(hits_out.end(), pep.getHits().begin(), pep.getHits().end());
  }
  for (auto& hit : hits_in)
  {
    // "target+decoy" is counted as "target" in IdentificationData:
    if (hit.getMetaValue("target_decoy") == "target+decoy")
    {
      hit.setMetaValue("target_decoy", "target");
    }
  }
  TEST_EQUAL(hits_in.size(), hits_out.size());
  // order of hits is different, check that every output one is in the input:
  TEST_EQUAL(all_of(hits_out.begin(), hits_out.end(), [&hits_in](const PeptideHit& hit)
  {
    return find(hits_in.begin(), hits_in.end(), hit) != hits_in.end();
  }), true);

  // and the other way round!
  TEST_EQUAL(all_of(hits_in.begin(), hits_in.end(), [&hits_out](const PeptideHit& hit)
  {
    return find(hits_out.begin(), hits_out.end(), hit) != hits_out.end();
  }), true);

  auto mzrtcomp = [](const PeptideIdentification& p1, const PeptideIdentification& p2)
      {return p1.getMZ() == p2.getMZ() && p1.getRT() == p2.getRT();};

  TEST_EQUAL(peptides_in.size(), peptides_out.size());
  // order of ids is different, check that every output one is in the input:
  TEST_EQUAL(all_of(peptides_out.begin(), peptides_out.end(), [&peptides_in, &mzrtcomp](const PeptideIdentification& hit) -> bool
  {
    return std::find_if(peptides_in.begin(), peptides_in.end(), std::bind(mzrtcomp, hit, std::placeholders::_1)) != peptides_in.end();
  }), true);

  // and the other way round!
  TEST_EQUAL(all_of(peptides_in.begin(), peptides_in.end(), [&peptides_out, &mzrtcomp](const PeptideIdentification& hit) -> bool
  {
    return std::find_if(peptides_out.begin(), peptides_out.end(), std::bind(mzrtcomp, hit, std::placeholders::_1)) != peptides_out.end();
  }), true);

  // filename = OPENMS_GET_TEST_DATA_PATH("IdentificationDataConverter_out2.idXML");
  // IdXMLFile().store(filename, proteins_out, peptides_out);
}
END_SECTION

START_SECTION((MzTab exportMzTab(const IdentificationData& id_data)))
{
  vector<ProteinIdentification> proteins_in;
  PeptideIdentificationList peptides_in;
  String filename = OPENMS_GET_TEST_DATA_PATH("../../../topp/THIRDPARTY/FidoAdapter_4_output.idXML");
  IdXMLFile().load(filename, proteins_in, peptides_in);

  IdentificationData ids;
  IdentificationDataConverter::importIDs(ids, proteins_in, peptides_in);

  MzTab mztab = IdentificationDataConverter::exportMzTab(ids);
  NEW_TMP_FILE(filename);
  MzTabFile().store(filename, mztab);

  TEST_FILE_SIMILAR(filename, OPENMS_GET_TEST_DATA_PATH("IdentificationDataConverter_out1.mzTab"));

  // RNA data, oligonucleotide that matches several times in the same RNA:
  IdentificationData rna_ids;
  IdentificationData::ParentSequence rna("test", IdentificationData::MoleculeType::RNA, "AUCGAUCG");
  IdentificationData::ParentSequenceRef ref = rna_ids.registerParentSequence(rna);
  IdentificationData::IdentifiedOligo oli(NASequence::fromString("AUCG"));
  IdentificationData::ParentMatch match1(0, 3), match2(4, 7);
  oli.parent_matches[ref].insert(match1);
  oli.parent_matches[ref].insert(match2);
  rna_ids.registerIdentifiedOligo(oli);

  mztab = IdentificationDataConverter::exportMzTab(rna_ids);
  NEW_TMP_FILE(filename);
  MzTabFile().store(filename, mztab);

  TEST_FILE_SIMILAR(filename, OPENMS_GET_TEST_DATA_PATH("IdentificationDataConverter_out2.mzTab"));

  // Test average mass export for RNA OSMs:
  // Create RNA identification data with observation match and average mass search param
  IdentificationData rna_avg_ids;

  // Register input file
  IdentificationData::InputFile input("test_input.mzML");
  IdentificationData::InputFileRef input_ref = rna_avg_ids.registerInputFile(input);

  // Register processing software
  IdentificationData::ProcessingSoftware sw("TestEngine", "1.0");
  IdentificationData::ProcessingSoftwareRef sw_ref = rna_avg_ids.registerProcessingSoftware(sw);

  // Register search parameters with AVERAGE mass type
  IdentificationData::DBSearchParam search_param;
  search_param.molecule_type = IdentificationData::MoleculeType::RNA;
  search_param.mass_type = IdentificationData::MassType::AVERAGE;
  IdentificationData::SearchParamRef search_ref = rna_avg_ids.registerDBSearchParam(search_param);

  // Register processing step linked to search parameters
  IdentificationData::ProcessingStep step(sw_ref);
  step.input_file_refs.push_back(input_ref);
  IdentificationData::ProcessingStepRef step_ref = rna_avg_ids.registerProcessingStep(step, search_ref);

  // Register parent sequence (RNA)
  IdentificationData::ParentSequence rna_parent("test_rna", IdentificationData::MoleculeType::RNA, "AUCGAUCG");
  IdentificationData::ParentSequenceRef rna_ref = rna_avg_ids.registerParentSequence(rna_parent);

  // Register identified oligonucleotide
  NASequence oligo_seq = NASequence::fromString("AUCG");
  IdentificationData::IdentifiedOligo oligo_id(oligo_seq);
  IdentificationData::ParentMatch pmatch(0, 3);
  oligo_id.parent_matches[rna_ref].insert(pmatch);
  IdentificationData::IdentifiedOligoRef oligo_ref = rna_avg_ids.registerIdentifiedOligo(oligo_id);

  // Register observation
  IdentificationData::Observation obs("scan_1", input_ref, 100.0, 500.0);
  IdentificationData::ObservationRef obs_ref = rna_avg_ids.registerObservation(obs);

  // Register observation match with the processing step
  IdentificationData::ObservationMatch obs_match(oligo_ref, obs_ref, 2);
  IdentificationData::AppliedProcessingStep applied_step(step_ref);
  obs_match.addProcessingStep(applied_step);
  rna_avg_ids.registerObservationMatch(obs_match);

  // Export to MzTab and check OSM section
  mztab = IdentificationDataConverter::exportMzTab(rna_avg_ids);
  const MzTabOSMSectionRows& osm_rows = mztab.getOSMSectionRows();

  TEST_EQUAL(osm_rows.size(), 1);
  if (!osm_rows.empty())
  {
    // The calculated mass should be the average weight, not monoisotopic
    double expected_avg_mass = oligo_seq.getAverageWeight(NASequence::Full, 2) / 2; // charge 2
    double expected_mono_mass = oligo_seq.getMonoWeight(NASequence::Full, 2) / 2; // charge 2
    double actual_mass = osm_rows[0].calc_mass_to_charge.get();
    // Check that actual mass is closer to average mass than monoisotopic mass
    TEST_REAL_SIMILAR(actual_mass, expected_avg_mass);
    // Verify it's different from monoisotopic (sanity check)
    TEST_NOT_EQUAL(fabs(expected_avg_mass - expected_mono_mass) < 0.001, true);
  }

  // Also test with MONOISOTOPIC mass type for comparison
  IdentificationData rna_mono_ids;
  input_ref = rna_mono_ids.registerInputFile(input);
  sw_ref = rna_mono_ids.registerProcessingSoftware(sw);

  // Register search parameters with MONOISOTOPIC mass type
  IdentificationData::DBSearchParam mono_search_param;
  mono_search_param.molecule_type = IdentificationData::MoleculeType::RNA;
  mono_search_param.mass_type = IdentificationData::MassType::MONOISOTOPIC;
  search_ref = rna_mono_ids.registerDBSearchParam(mono_search_param);

  step = IdentificationData::ProcessingStep(sw_ref);
  step.input_file_refs.push_back(input_ref);
  step_ref = rna_mono_ids.registerProcessingStep(step, search_ref);

  rna_ref = rna_mono_ids.registerParentSequence(rna_parent);
  oligo_id = IdentificationData::IdentifiedOligo(oligo_seq);
  oligo_id.parent_matches[rna_ref].insert(pmatch);
  oligo_ref = rna_mono_ids.registerIdentifiedOligo(oligo_id);
  obs_ref = rna_mono_ids.registerObservation(obs);
  obs_match = IdentificationData::ObservationMatch(oligo_ref, obs_ref, 2);
  applied_step = IdentificationData::AppliedProcessingStep(step_ref);
  obs_match.addProcessingStep(applied_step);
  rna_mono_ids.registerObservationMatch(obs_match);

  mztab = IdentificationDataConverter::exportMzTab(rna_mono_ids);
  const MzTabOSMSectionRows& mono_osm_rows = mztab.getOSMSectionRows();

  TEST_EQUAL(mono_osm_rows.size(), 1);
  if (!mono_osm_rows.empty())
  {
    // The calculated mass should be the monoisotopic weight
    double expected_mono_mass = oligo_seq.getMonoWeight(NASequence::Full, 2) / 2; // charge 2
    double actual_mass = mono_osm_rows[0].calc_mass_to_charge.get();
    TEST_REAL_SIMILAR(actual_mass, expected_mono_mass);
  }
}
END_SECTION

/*
// performance test on a large file:
START_SECTION(([[EXTRA]] void importIDs(IdentificationData&, const vector<ProteinIdentification>&, const PeptideIdentificationList&)))
{
  SysInfo::MemUsage mem_usage;
  vector<ProteinIdentification> proteins_in;
  PeptideIdentificationList peptides_in;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("large_test.idXML"), proteins_in, peptides_in);
  STATUS(mem_usage.delta("PeptideIdentification/ProteinIdentification"));

  TEST_EQUAL(proteins_in.size(), 1);
  TEST_EQUAL(proteins_in[0].getHits().size(), 11098);
  TEST_EQUAL(peptides_in.size(), 328591);
  TEST_EQUAL(proteins_in[0].getIndistinguishableProteins().size(), 10853);
  TEST_EQUAL(proteins_in[0].getProteinGroups().size(), 9092);

  mem_usage.reset();
  mem_usage.before();
  IdentificationData ids;
  IdentificationDataConverter::importIDs(ids, proteins_in, peptides_in);
  STATUS(mem_usage.delta("IdentificationData"));

  TEST_EQUAL(ids.getParentSequences().size(), 11098);
  // problem: input data comes from multiple files, spectra with matching names
  // in different files get merged together -> lower number of input items:
  TEST_EQUAL(ids.getObservations().size(), 55522);
  TEST_EQUAL(ids.getIdentifiedPeptides().size(), 73950);
  // according to "grep" on the input file, there should be 335250 peptide hits
  // in total - maybe some duplicates?:
  TEST_EQUAL(ids.getObservationMatches().size(), 332778);

  TEST_EQUAL(ids.getParentGroupSets().size(), 2);
  TEST_EQUAL(ids.getParentGroupSets()[0].groups.size(), 10853);
  TEST_EQUAL(ids.getParentGroupSets()[1].groups.size(), 9092);
}
END_SECTION
*/

FeatureMap features; // persist through sections

START_SECTION((void importFeatureIDs(FeatureMap& features, bool clear_original)))
{
  FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_1.featureXML"), features);
  // protein and peptide IDs use same score type (name) with different orientations;
  // IdentificationData doesn't allow this, so change it here:
  for (auto& run : features.getProteinIdentifications())
  {
    run.setScoreType(run.getScoreType() + "_protein");
  }
  IdentificationDataConverter::importFeatureIDs(features);
  TEST_EQUAL(features.getIdentificationData().getObservations().size(), 5);
  TEST_EQUAL(features.getIdentificationData().getObservationMatches().size(), 7);
  TEST_EQUAL(features.getIdentificationData().getIdentifiedPeptides().size(), 7);
  TEST_EQUAL(features.getIdentificationData().getParentSequences().size(), 3);
  TEST_EQUAL(features[0].getIDMatches().size(), 3);
  TEST_EQUAL(features[1].getIDMatches().size(), 1);
  TEST_EQUAL(features.getUnassignedIDMatches().size(), 3);
  // check that original IDs were cleared:
  TEST_EQUAL(features.getProteinIdentifications().size(), 0);
  TEST_EQUAL(features.getUnassignedPeptideIdentifications().size(), 0);
  TEST_EQUAL(features[0].getPeptideIdentifications().size(), 0);
  TEST_EQUAL(features[1].getPeptideIdentifications().size(), 0);
}
END_SECTION

START_SECTION((void exportFeatureIDs(FeatureMap& features, bool clear_original)))
{
  // convert IDs from previous test back:
  IdentificationDataConverter::exportFeatureIDs(features);
  TEST_EQUAL(features.getProteinIdentifications().size(), 2);
  TEST_EQUAL(features.getUnassignedPeptideIdentifications().size(), 2);
  TEST_EQUAL(features[0].getPeptideIdentifications().size(), 2);
  TEST_EQUAL(features[1].getPeptideIdentifications().size(), 1);
  // check that "original" IDs were cleared:
  TEST_EQUAL(features.getIdentificationData().empty(), true);
}
END_SECTION

ConsensusMap consensus; // persist through sections

START_SECTION((void importConsensusIDs(ConsensusMap& consensus, bool clear_original)))
{
  ConsensusXMLFile().load(OPENMS_GET_TEST_DATA_PATH("ConsensusXMLFile_1.consensusXML"), consensus);
  // protein and peptide IDs use same score type (name) with different orientations;
  // IdentificationData doesn't allow this, so change it here:
  for (auto& run : consensus.getProteinIdentifications())
  {
    run.setScoreType(run.getScoreType() + "_protein");
  }
  IdentificationDataConverter::importConsensusIDs(consensus);
  TEST_EQUAL(consensus.getIdentificationData().getObservations().size(), 5);
  TEST_EQUAL(consensus.getIdentificationData().getObservationMatches().size(), 7);
  TEST_EQUAL(consensus.getIdentificationData().getIdentifiedPeptides().size(), 7);
  TEST_EQUAL(consensus.getIdentificationData().getParentSequences().size(), 3);
  TEST_EQUAL(consensus[0].getIDMatches().size(), 3);
  TEST_EQUAL(consensus[1].getIDMatches().size(), 1);
  TEST_EQUAL(consensus.getUnassignedIDMatches().size(), 3);
  // check that original IDs were cleared:
  TEST_EQUAL(consensus.getProteinIdentifications().size(), 0);
  TEST_EQUAL(consensus.getUnassignedPeptideIdentifications().size(), 0);
  TEST_EQUAL(consensus[0].getPeptideIdentifications().size(), 0);
  TEST_EQUAL(consensus[1].getPeptideIdentifications().size(), 0);
}
END_SECTION

START_SECTION((void exportConsensusIDs(ConsensusMap& consensus, bool clear_original)))
{
  // convert IDs from previous test back:
  IdentificationDataConverter::exportConsensusIDs(consensus);
  TEST_EQUAL(consensus.getProteinIdentifications().size(), 2);
  TEST_EQUAL(consensus.getUnassignedPeptideIdentifications().size(), 2);
  TEST_EQUAL(consensus[0].getPeptideIdentifications().size(), 2);
  TEST_EQUAL(consensus[1].getPeptideIdentifications().size(), 1);
  // check that "original" IDs were cleared:
  TEST_EQUAL(consensus.getIdentificationData().empty(), true);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
