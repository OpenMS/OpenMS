// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/CONCEPT/FuzzyStringComparator.h>
#include <OpenMS/CHEMISTRY/CrossLinksDB.h>
#include <OpenMS/CONCEPT/Constants.h>

#include <OpenMS/FORMAT/HANDLERS/MzIdentMLDOMHandler.h>

#include <set>
#include <regex>
#include <fstream>
#include <sstream>


using namespace OpenMS;
using namespace std;

// helpers for the protein group tests (issue #1105)

// accessions of a protein group, for comparisons
static std::vector<std::string> accessionsOf(const ProteinIdentification::ProteinGroup& group)
{
  return group.accessions;
}

// A protein inference run as written by Epifany: indistinguishable proteins A and B, protein C, decoy protein
// DECOY_D, and a general protein group (cluster) of A, B and C. PSMs: PEPTIDEK (A, B, C), ELVISK (A, B), the
// decoy GLYCINEK (DECOY_D), and a second-ranked hit SAMPLER (C) that is not used as evidence.
static void makeInferenceRun(std::vector<ProteinIdentification>& proteins, PeptideIdentificationList& peptides)
{
  ProteinIdentification run;
  run.setIdentifier("inference_run");
  run.setSearchEngine("MS-GF+");
  run.setSearchEngineVersion("2024.03.26");
  ProteinIdentification::SearchParameters search_parameters;
  search_parameters.db = "db.fasta";
  run.setSearchParameters(search_parameters);
  run.setInferenceEngine("Epifany"); // stored in the search parameters
  run.setInferenceEngineVersion("3.6.0");
  run.setMetaValue("spectra_data", StringList{"run.mzML"});
  run.setScoreType("Posterior Probability");
  run.setHigherScoreBetter(true);
  auto add_hit = [&run](const std::string& acc, double score, const std::string& seq, bool decoy)
  {
    ProteinHit hit(score, 1, acc, seq);
    hit.setTargetDecoyType(decoy ? ProteinHit::TargetDecoyType::DECOY : ProteinHit::TargetDecoyType::TARGET);
    run.insertHit(hit);
  };
  add_hit("A", 0.98, "MPEPTIDEKELVISK", false);
  add_hit("B", 0.98, "MPEPTIDEKELVISKR", false);
  add_hit("C", 0.5, "MPEPTIDEKSAMPLER", false);
  add_hit("DECOY_D", 0.3, "MGLYCINEK", true);
  run.getHits()[0].setCoverage(93.3);
  ProteinIdentification::ProteinGroup ab, c, d, abc;
  ab.probability = 0.98;
  ab.accessions = {"A", "B"};
  c.probability = 0.5;
  c.accessions = {"C"};
  d.probability = 0.3;
  d.accessions = {"DECOY_D"};
  abc.probability = 0.99;
  abc.accessions = {"A", "B", "C"};
  run.insertIndistinguishableProteins(ab);
  run.insertIndistinguishableProteins(c);
  run.insertIndistinguishableProteins(d);
  run.insertProteinGroup(abc);
  proteins.push_back(run);

  auto add_psm = [&peptides](const std::string& spectrum, const std::vector<std::pair<std::string, std::vector<std::string>>>& hits, bool decoy)
  {
    PeptideIdentification id;
    id.setIdentifier("inference_run");
    id.setScoreType("q-value");
    id.setHigherScoreBetter(false);
    id.setSpectrumReference(spectrum);
    id.setRT(100.0);
    id.setMZ(500.0);
    double score = 0.001;
    UInt rank = 0;
    for (const auto& [sequence, accessions] : hits)
    {
      PeptideHit hit(score, rank++, 2, AASequence::fromString(sequence));
      score *= 10;
      for (const std::string& acc : accessions)
      {
        PeptideEvidence evidence;
        evidence.setProteinAccession(acc);
        hit.addPeptideEvidence(evidence);
      }
      hit.setMetaValue(Constants::UserParam::TARGET_DECOY, decoy ? "decoy" : "target");
      id.insertHit(hit);
    }
    peptides.push_back(id);
  };
  add_psm("index=1", {{"PEPTIDEK", {"A", "B", "C"}}, {"SAMPLER", {"C"}}}, false);
  add_psm("index=2", {{"ELVISK", {"A", "B"}}}, false);
  add_psm("index=3", {{"GLYCINEK", {"DECOY_D"}}}, true);
}

static Size countOf(const std::string& text, const std::string& what)
{
  Size n = 0;
  for (std::string::size_type pos = text.find(what); pos != std::string::npos; pos = text.find(what, pos + what.size()))
  {
    ++n;
  }
  return n;
}

static std::string readFile(const std::string& filename)
{
  std::ifstream in(filename);
  std::stringstream buffer;
  buffer << in.rdbuf();
  return buffer.str();
}

START_TEST(MzIdentMLFile, "$Id")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


MzIdentMLFile* ptr = nullptr;
MzIdentMLFile* nullPointer = nullptr;

START_SECTION((MzIdentMLFile()))
{
  ptr = new MzIdentMLFile;
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((~MzIdentMLFile()))
{
  delete ptr;
}
END_SECTION

START_SECTION(void load(const std::string& filename, std::vector<ProteinIdentification>& protein_ids, PeptideIdentificationList& peptide_ids) )
{
  std::vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  std::vector<std::string> fm{"Carbamidomethyl (C)", "Xlink:DTSSP[88] (Protein N-term)"};
  MzIdentMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_msgf_mini.mzid"), protein_ids, peptide_ids);

  TEST_EQUAL(protein_ids.size(),2)
  TEST_EQUAL(protein_ids[0].getHits().size(),2)
  TEST_EQUAL(protein_ids[1].getHits().size(),1)
  TEST_EQUAL(peptide_ids.size(),5)

  // loading again into the same containers must replace, not accumulate
  // (regression: the DOM handler only appends; load() now clears first, matching IdXMLFile::load)
  MzIdentMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_msgf_mini.mzid"), protein_ids, peptide_ids);
  TEST_EQUAL(protein_ids.size(),2)
  TEST_EQUAL(peptide_ids.size(),5)
  TEST_EQUAL(peptide_ids[0].getHits().size(),1)
  TEST_EQUAL(peptide_ids[1].getHits().size(),1)
  TEST_EQUAL(peptide_ids[2].getHits().size(),1)
  TEST_EQUAL(peptide_ids[3].getHits().size(),1)
  TEST_EQUAL(peptide_ids[4].getHits().size(),1)

  /////////////// protein id 1 //////////////////
  TEST_EQUAL(protein_ids[0].getSearchEngine(),"MS-GF+")
  TEST_EQUAL(protein_ids[0].getSearchEngineVersion(),"Beta (v9979)")
  TEST_NOT_EQUAL(protein_ids[0].getDateTime().getDate(),"0000-00-00")
  TEST_NOT_EQUAL(protein_ids[0].getDateTime().getTime(),"00:00:00")
  TEST_EQUAL(protein_ids[0].getSearchParameters().db,"database.fasta")
  TEST_EQUAL(protein_ids[0].getSearchParameters().missed_cleavages, 1000)
  ABORT_IF(protein_ids[0].getSearchParameters().fixed_modifications.size() != fm.size())
  for (size_t i = 0; i < fm.size(); ++i)
  {
    TEST_EQUAL(protein_ids[0].getSearchParameters().fixed_modifications[i], fm[i]);
  }
  TEST_REAL_SIMILAR(protein_ids[0].getSearchParameters().fragment_mass_tolerance,0)
  TEST_REAL_SIMILAR(protein_ids[0].getSearchParameters().precursor_mass_tolerance,20)

  //ProteinGroups not nupported yet, also no ProteinDetection, too few input here
  //  TEST_EQUAL(protein_ids[0].getProteinGroups().size(), 0);
  //  TEST_EQUAL(protein_ids[0].getIndistinguishableProteins().size(), 0);

  //protein hit 1
  TEST_EQUAL(protein_ids[0].getHits()[0].getAccession(),"sp|P0A9K9|SLYD_ECOLI")
  TEST_EQUAL(protein_ids[0].getHits()[0].getSequence(),"")
  //protein hit 2
  TEST_EQUAL(protein_ids[0].getHits()[1].getAccession(),"sp|P0A786|PYRB_ECOLI")
  TEST_EQUAL(protein_ids[0].getHits()[1].getSequence(),"")

  //peptide id s
  TEST_EQUAL(peptide_ids[0].getScoreType(),"MS-GF:RawScore")
  TEST_REAL_SIMILAR(peptide_ids[0].getHits()[0].getScore(),195)
  TEST_EQUAL(peptide_ids[0].getHits()[0].getSequence().toString(),"LATEFSGNVPVLNAGDGSNQHPTQTLLDLFTIQETQGR")
  TEST_EQUAL(peptide_ids[0].getSpectrumReference(),"controllerType=0 controllerNumber=1 scan=32805")
  TEST_EQUAL(peptide_ids[1].getScoreType(),"MS-GF:RawScore")
  TEST_REAL_SIMILAR(peptide_ids[1].getHits()[0].getScore(),182)
  TEST_EQUAL(peptide_ids[1].getHits()[0].getSequence().toString(),"FLAETDQGPVPVEITAVEDDHVVVDGNHMLAGQNLK")
  TEST_EQUAL(peptide_ids[1].getSpectrumReference(),"controllerType=0 controllerNumber=1 scan=26090")
  TEST_EQUAL(peptide_ids[2].getScoreType(),"MS-GF:RawScore")
  TEST_REAL_SIMILAR(peptide_ids[2].getHits()[0].getScore(),191)
  TEST_EQUAL(peptide_ids[2].getHits()[0].getSequence().toString(),"FLAETDQGPVPVEITAVEDDHVVVDGNHMLAGQNLK")
  TEST_EQUAL(peptide_ids[2].getSpectrumReference(),"controllerType=0 controllerNumber=1 scan=26157")
  TEST_EQUAL(peptide_ids[3].getScoreType(),"MS-GF:RawScore")
  TEST_REAL_SIMILAR(peptide_ids[3].getHits()[0].getScore(),211)
  TEST_EQUAL(peptide_ids[3].getHits()[0].getSequence().toString(),"VGAGPFPTELFDETGEFLC(Carbamidomethyl)K")
  TEST_EQUAL(peptide_ids[3].getSpectrumReference(),"controllerType=0 controllerNumber=1 scan=15094")
}
END_SECTION

START_SECTION(([EXTRA] read mzIdentML Modification without the optional location attribute (issue #5443)))
{
  // The 'location' attribute of <Modification> is optional in mzIdentML; some search engines omit it for
  // terminal modifications. Reading such a file used to crash with an IndexOverflow exception. The reader
  // must now infer a terminal position (or skip the modification) instead of applying a negative index.
  std::vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  MzIdentMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_missing_mod_location.mzid"), protein_ids, peptide_ids);

  ABORT_IF(peptide_ids.size() != 5)
  for (Size i = 0; i < peptide_ids.size(); ++i)
  {
    ABORT_IF(peptide_ids[i].getHits().empty())
  }

  // 1) N-terminal Acetyl without 'location' -> inferred as N-terminal
  const AASequence& acetyl_seq = peptide_ids[0].getHits()[0].getSequence();
  TEST_EQUAL(acetyl_seq.toUnmodifiedString(), "PEPTIDEK")
  TEST_TRUE(acetyl_seq.hasNTerminalModification())
  TEST_FALSE(acetyl_seq.hasCTerminalModification())
  TEST_EQUAL(acetyl_seq.getNTerminalModificationName(), "Acetyl")

  // 2) Oxidation without 'location' -> position not uniquely terminal -> skipped (but no crash and no misplacement)
  const AASequence& oxidation_seq = peptide_ids[1].getHits()[0].getSequence();
  TEST_EQUAL(oxidation_seq.toUnmodifiedString(), "PEPTIDEM")
  TEST_FALSE(oxidation_seq.isModified())
  TEST_EQUAL(oxidation_seq.toString(), "PEPTIDEM")

  // 3) Amidated without 'location' -> exclusively C-terminal -> inferred as C-terminal
  const AASequence& amidated_seq = peptide_ids[2].getHits()[0].getSequence();
  TEST_EQUAL(amidated_seq.toUnmodifiedString(), "PEPTIDER")
  TEST_TRUE(amidated_seq.hasCTerminalModification())
  TEST_FALSE(amidated_seq.hasNTerminalModification())
  TEST_EQUAL(amidated_seq.getCTerminalModificationName(), "Amidated")

  // 4) Control: a regular internal modification WITH a valid 'location' is still parsed correctly
  const AASequence& carbamidomethyl_seq = peptide_ids[3].getHits()[0].getSequence();
  TEST_EQUAL(carbamidomethyl_seq.toString(), "PEPTIDEC(Carbamidomethyl)K")

  // 5) Acetyl without 'location' but with residues="K" -> residue-specific (internal), position unknown,
  //    so it is skipped rather than forced onto the N-terminus even though Acetyl supports the N-terminus
  const AASequence& acetyl_on_k_seq = peptide_ids[4].getHits()[0].getSequence();
  TEST_EQUAL(acetyl_on_k_seq.toUnmodifiedString(), "PEPTIKDE")
  TEST_FALSE(acetyl_on_k_seq.isModified())
}
END_SECTION

START_SECTION(void store(std::string filename, const std::vector<ProteinIdentification>& protein_ids, const PeptideIdentificationList& peptide_ids) )
{
  //store and load data from various sources, starting with idxml, contents already checked above, so checking integrity of the data over repeated r/w
  std::vector<ProteinIdentification> protein_ids, protein_ids2;
  PeptideIdentificationList peptide_ids, peptide_ids2;
  std::string input_path = OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_whole.mzid");
  MzIdentMLFile().load(input_path, protein_ids2, peptide_ids2);
  std::string filename;
  NEW_TMP_FILE(filename)
  MzIdentMLFile().store(filename, protein_ids2, peptide_ids2);

  MzIdentMLFile().load(filename, protein_ids, peptide_ids);
  TEST_EQUAL(protein_ids.size(),protein_ids2.size())
  TEST_EQUAL(protein_ids[0].getHits().size(),protein_ids2[0].getHits().size())
  TEST_EQUAL(peptide_ids.size(),peptide_ids2.size())
  TEST_EQUAL(peptide_ids[0].getHits().size(),peptide_ids2[0].getHits().size())
  TEST_EQUAL(peptide_ids[1].getHits().size(),peptide_ids2[1].getHits().size())
  TEST_EQUAL(peptide_ids[2].getHits().size(),peptide_ids2[2].getHits().size())

  /////////////// protein id 1 //////////////////
  TEST_EQUAL(protein_ids[0].getSearchEngine(),protein_ids2[0].getSearchEngine())
  TEST_EQUAL(protein_ids[0].getSearchEngineVersion(),protein_ids2[0].getSearchEngineVersion())
  TEST_EQUAL(protein_ids[0].getDateTime().getDate(),protein_ids2[0].getDateTime().getDate())
  TEST_EQUAL(protein_ids[0].getDateTime().getTime(),protein_ids2[0].getDateTime().getTime())
  TEST_EQUAL(protein_ids[0].getSearchParameters().db,protein_ids2[0].getSearchParameters().db)
  TEST_EQUAL(protein_ids[0].getSearchParameters().db_version,protein_ids2[0].getSearchParameters().db_version)
  TEST_EQUAL(protein_ids[0].getSearchParameters().digestion_enzyme.getName(),protein_ids2[0].getSearchParameters().digestion_enzyme.getName())
  TEST_EQUAL(protein_ids[0].getSearchParameters().charges,protein_ids2[0].getSearchParameters().charges)
  TEST_EQUAL(protein_ids[0].getSearchParameters().mass_type,protein_ids2[0].getSearchParameters().mass_type)
  TEST_REAL_SIMILAR(protein_ids[0].getSearchParameters().fragment_mass_tolerance,protein_ids2[0].getSearchParameters().fragment_mass_tolerance)
  TEST_REAL_SIMILAR(protein_ids[0].getSearchParameters().precursor_mass_tolerance,protein_ids2[0].getSearchParameters().precursor_mass_tolerance)

  TEST_EQUAL(protein_ids[0].getSearchParameters().variable_modifications.size(),protein_ids2[0].getSearchParameters().variable_modifications.size())
  for (size_t i = 0; i < protein_ids[0].getSearchParameters().variable_modifications.size(); ++i)
  {
    TEST_STRING_EQUAL(protein_ids[0].getSearchParameters().variable_modifications[i],protein_ids2[0].getSearchParameters().variable_modifications[i])
  }
  TEST_STRING_EQUAL(protein_ids[0].getSearchParameters().variable_modifications.back(),"Acetyl (N-term)")
  TEST_EQUAL(protein_ids[0].getSearchParameters().fixed_modifications.size(),protein_ids2[0].getSearchParameters().fixed_modifications.size())
  for (size_t i = 0; i < protein_ids[0].getSearchParameters().fixed_modifications.size(); ++i)
  {
    TEST_STRING_EQUAL(protein_ids[0].getSearchParameters().fixed_modifications[i],protein_ids2[0].getSearchParameters().fixed_modifications[i])
  }

  //ProteinGroups not nupported yet, also no ProteinDetection, too few input here
//  TEST_EQUAL(protein_ids[0].getProteinGroups().size(), 0);
//  TEST_EQUAL(protein_ids[0].getIndistinguishableProteins().size(), 0);

  //protein hit 1
  TEST_EQUAL(protein_ids[0].getHits()[0].getAccession(),protein_ids2[0].getHits()[0].getAccession())
  TEST_EQUAL(protein_ids[0].getHits()[0].getSequence(),protein_ids2[0].getHits()[0].getSequence())
  //protein hit 2
  TEST_EQUAL(protein_ids[0].getHits()[1].getAccession(),protein_ids2[0].getHits()[1].getAccession())
  TEST_EQUAL(protein_ids[0].getHits()[1].getSequence(),protein_ids2[0].getHits()[1].getSequence())

  //peptide id 1
  TEST_EQUAL(peptide_ids[0].getScoreType(),peptide_ids2[0].getScoreType())
  TEST_EQUAL(peptide_ids[0].isHigherScoreBetter(),peptide_ids2[0].isHigherScoreBetter())
  TEST_REAL_SIMILAR(peptide_ids[0].getMZ(),peptide_ids2[0].getMZ())
  TEST_REAL_SIMILAR(peptide_ids[0].getRT(),peptide_ids2[0].getRT())
  TEST_EQUAL(peptide_ids[0].getSpectrumReference(),peptide_ids2[0].getSpectrumReference())
  //peptide hit 1
  TEST_REAL_SIMILAR(peptide_ids[0].getHits()[0].getScore(),peptide_ids2[0].getHits()[0].getScore())
  TEST_EQUAL(peptide_ids[0].getHits()[0].getSequence(),peptide_ids2[0].getHits()[0].getSequence())
  TEST_EQUAL(peptide_ids[0].getHits()[0].getCharge(),peptide_ids2[0].getHits()[0].getCharge())
  for (size_t i = 0; i < peptide_ids[0].getHits()[0].getPeptideEvidences().size(); ++i){
    //AA before/after tested by peptide evidences vector equality check - not working if the order of proteins is pertubated
    //TEST_EQUAL(peptide_ids[0].getHits()[0].getPeptideEvidences()[i]==peptide_ids2[0].getHits()[0].getPeptideEvidences()[i],true)
    TEST_EQUAL(peptide_ids[0].getHits()[0].getPeptideEvidences()[i].getStart(),peptide_ids2[0].getHits()[0].getPeptideEvidences()[i].getStart())
    TEST_EQUAL(peptide_ids[0].getHits()[0].getPeptideEvidences()[i].getEnd(),peptide_ids2[0].getHits()[0].getPeptideEvidences()[i].getEnd())
    TEST_EQUAL(peptide_ids[0].getHits()[0].getPeptideEvidences()[i].getAABefore(),peptide_ids2[0].getHits()[0].getPeptideEvidences()[i].getAABefore())
    TEST_EQUAL(peptide_ids[0].getHits()[0].getPeptideEvidences()[i].getAAAfter(),peptide_ids2[0].getHits()[0].getPeptideEvidences()[i].getAAAfter())
//    TEST_EQUAL(peptide_ids[0].getHits()[0].getPeptideEvidences()[i].getProteinAccession(),peptide_ids2[0].getHits()[0].getPeptideEvidences()[i].getProteinAccession())
  }
  //peptide hit 2
  TEST_REAL_SIMILAR(peptide_ids[0].getHits()[1].getScore(),peptide_ids2[0].getHits()[1].getScore())
  TEST_EQUAL(peptide_ids[0].getHits()[1].getSequence(),peptide_ids2[0].getHits()[1].getSequence())
  TEST_EQUAL(peptide_ids[0].getHits()[1].getCharge(),peptide_ids2[0].getHits()[1].getCharge())
  for (size_t i = 0; i < peptide_ids[0].getHits()[1].getPeptideEvidences().size(); ++i){
//    TEST_EQUAL(peptide_ids[0].getHits()[1].getPeptideEvidences()[i]==peptide_ids2[0].getHits()[1].getPeptideEvidences()[i],true)
    TEST_EQUAL(peptide_ids[0].getHits()[1].getPeptideEvidences()[i].getStart(),peptide_ids2[0].getHits()[1].getPeptideEvidences()[i].getStart())
    TEST_EQUAL(peptide_ids[0].getHits()[1].getPeptideEvidences()[i].getEnd(),peptide_ids2[0].getHits()[1].getPeptideEvidences()[i].getEnd())
    TEST_EQUAL(peptide_ids[0].getHits()[1].getPeptideEvidences()[i].getAABefore(),peptide_ids2[0].getHits()[1].getPeptideEvidences()[i].getAABefore())
    TEST_EQUAL(peptide_ids[0].getHits()[1].getPeptideEvidences()[i].getAAAfter(),peptide_ids2[0].getHits()[1].getPeptideEvidences()[i].getAAAfter())
  }

  //peptide id 2
  TEST_EQUAL(peptide_ids[1].getScoreType(),peptide_ids2[1].getScoreType())
  TEST_EQUAL(peptide_ids[1].isHigherScoreBetter(),peptide_ids2[1].isHigherScoreBetter())
  TEST_REAL_SIMILAR(peptide_ids[1].getMZ(),peptide_ids2[1].getMZ())
  TEST_REAL_SIMILAR(peptide_ids[1].getRT(),peptide_ids2[1].getRT())
  //peptide hit 1
  TEST_REAL_SIMILAR(peptide_ids[1].getHits()[0].getScore(),peptide_ids2[1].getHits()[0].getScore())
  TEST_EQUAL(peptide_ids[1].getHits()[0].getSequence(),peptide_ids2[1].getHits()[0].getSequence())
  TEST_EQUAL(peptide_ids[1].getHits()[0].getCharge(),peptide_ids2[1].getHits()[0].getCharge())
  for (size_t i = 0; i < peptide_ids[1].getHits()[0].getPeptideEvidences().size(); ++i)
    TEST_EQUAL(peptide_ids[1].getHits()[0].getPeptideEvidences()[i]==peptide_ids2[1].getHits()[0].getPeptideEvidences()[i],true)
  //peptide hit 2
  TEST_REAL_SIMILAR(peptide_ids[1].getHits()[1].getScore(),peptide_ids2[1].getHits()[1].getScore())
  TEST_EQUAL(peptide_ids[1].getHits()[1].getSequence(),peptide_ids2[1].getHits()[1].getSequence())
  TEST_EQUAL(peptide_ids[1].getHits()[1].getCharge(),peptide_ids2[1].getHits()[1].getCharge())
  for (size_t i = 0; i < peptide_ids[1].getHits()[1].getPeptideEvidences().size(); ++i)
    TEST_EQUAL(peptide_ids[1].getHits()[1].getPeptideEvidences()[i]==peptide_ids2[1].getHits()[1].getPeptideEvidences()[i],true)
  //peptide id 3
  TEST_EQUAL(peptide_ids[2].getScoreType(),peptide_ids2[2].getScoreType())
  TEST_EQUAL(peptide_ids[2].isHigherScoreBetter(),peptide_ids2[2].isHigherScoreBetter())
  TEST_REAL_SIMILAR(peptide_ids[2].getMZ(),peptide_ids2[2].getMZ())
  TEST_REAL_SIMILAR(peptide_ids[2].getRT(),peptide_ids2[2].getRT())
  //peptide hit 1
  TEST_REAL_SIMILAR(peptide_ids[2].getHits()[0].getScore(),peptide_ids2[2].getHits()[0].getScore())
  TEST_EQUAL(peptide_ids[2].getHits()[0].getSequence(),peptide_ids2[2].getHits()[0].getSequence())
  TEST_EQUAL(peptide_ids[2].getHits()[0].getCharge(),peptide_ids2[2].getHits()[0].getCharge())
  for (size_t i = 0; i < peptide_ids[2].getHits()[0].getPeptideEvidences().size(); ++i)
    TEST_EQUAL(peptide_ids[2].getHits()[0].getPeptideEvidences()[i]==peptide_ids2[2].getHits()[0].getPeptideEvidences()[i],true)
  //peptide hit 2
  TEST_REAL_SIMILAR(peptide_ids[1].getHits()[1].getScore(),peptide_ids2[1].getHits()[1].getScore())
  TEST_EQUAL(peptide_ids[2].getHits()[1].getSequence(),peptide_ids2[2].getHits()[1].getSequence())
  TEST_EQUAL(peptide_ids[2].getHits()[1].getCharge(),peptide_ids2[2].getHits()[1].getCharge())
  for (size_t i = 0; i < peptide_ids[2].getHits()[1].getPeptideEvidences().size(); ++i)
    TEST_EQUAL(peptide_ids[2].getHits()[1].getPeptideEvidences()[i]==peptide_ids2[2].getHits()[1].getPeptideEvidences()[i],true)
}
END_SECTION

START_SECTION(([EXTRA] multiple runs))
{
  std::vector<ProteinIdentification> protein_ids, protein_ids2;
  PeptideIdentificationList peptide_ids, peptide_ids2;
  std::string input_path = OPENMS_GET_TEST_DATA_PATH("MzIdentML_3runs.mzid");
  MzIdentMLFile().load(input_path, protein_ids2, peptide_ids2);
  std::string filename;
  NEW_TMP_FILE(filename)
  MzIdentMLFile().store(filename, protein_ids2, peptide_ids2);

  MzIdentMLFile().load(filename, protein_ids, peptide_ids);

  TEST_EQUAL(protein_ids.size(),protein_ids2.size())

  TEST_EQUAL(protein_ids[0].getHits().size(),protein_ids2[0].getHits().size())
  TEST_EQUAL(protein_ids[1].getHits().size(),protein_ids2[1].getHits().size())
  TEST_EQUAL(protein_ids[2].getHits().size(),protein_ids2[2].getHits().size())

  TEST_EQUAL(protein_ids[0].getSearchParameters().precursor_mass_tolerance_ppm, true)
}
END_SECTION

START_SECTION(([EXTRA] thresholds))
{
  std::vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  std::string input_path = OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_whole.mzid");
  MzIdentMLFile().load(input_path, protein_ids, peptide_ids);

  TEST_EQUAL(protein_ids.size(),1)
  TEST_EQUAL(protein_ids[0].getSignificanceThreshold(),0.5)

  TEST_EQUAL(peptide_ids.size(),3)
  for (size_t i = 0; i < peptide_ids.size(); ++i)
  {
    if (peptide_ids[i].getSpectrumReference() == "17")
    {
      TEST_EQUAL(peptide_ids[i].getHits().size(),2)
      for (size_t j = 0; j < peptide_ids[i].getHits().size(); ++j)
      {
        TEST_EQUAL(peptide_ids[i].getHits()[j].getMetaValue("pass_threshold"),false)
      }
      PeptideHit x = peptide_ids[i].getHits().back();
      x.removeMetaValue("pass_threshold");
      x.setSequence(AASequence::fromString("TESTER"));
      x.setScore(0.4);
      peptide_ids[i].insertHit(x);
    }
  }

  std::string filename;
  NEW_TMP_FILE(filename)
  MzIdentMLFile().store(filename, protein_ids, peptide_ids);
  protein_ids.clear();
  peptide_ids.clear();
  MzIdentMLFile().load(filename, protein_ids, peptide_ids);

  TEST_EQUAL(peptide_ids.size(),3)
  for (size_t i = 0; i < peptide_ids.size(); ++i)
  {
    if (peptide_ids[i].getSpectrumReference() == "17")
    {
      TEST_EQUAL(peptide_ids[i].getHits().size(),3)
      for (size_t j = 0; j < peptide_ids[i].getHits().size(); ++j)
      {
        if (peptide_ids[i].getHits()[j].getScore() > protein_ids[0].getSignificanceThreshold())
        {
          TEST_EQUAL(peptide_ids[i].getHits()[j].getMetaValue("pass_threshold"),false)
        }
        else
          TEST_EQUAL(peptide_ids[i].getHits()[j].getMetaValue("pass_threshold"),true)
      }
    }
  }

}
END_SECTION

START_SECTION(([EXTRA] regression test for file loading on example files))
{
  std::vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  std::string input_path = OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_whole.mzid");
  MzIdentMLFile().load(input_path, protein_ids, peptide_ids);
//  input_path = OPENMS_GET_TEST_DATA_PATH("Mascot_MSMS_example.mzid");
//  MzIdentMLFile().load(input_path, protein_ids, peptide_ids);
  input_path = OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_msgf_mini.mzid");
  MzIdentMLFile().load(input_path, protein_ids, peptide_ids);
  input_path = OPENMS_GET_TEST_DATA_PATH("MzIdentML_3runs.mzid");
  MzIdentMLFile().load(input_path, protein_ids, peptide_ids);
}
END_SECTION


START_SECTION(([EXTRA] compability issues))
{
//  MzIdentMLFile mzidfile;
//  vector<ProteinIdentification> protein_ids;
//  PeptideIdentificationList peptide_ids;
//  mzidfile.load(OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_no_proteinhits.mzid"), protein_ids, peptide_ids);

//  TEST_EQUAL(protein_ids.size(), 1)
//  TEST_EQUAL(protein_ids[0].getHits().size(), 0)
//  TEST_EQUAL(peptide_ids.size(), 10)
//  TEST_EQUAL(peptide_ids[0].getHits().size(), 1)

//  std::string filename;
//  NEW_TMP_FILE(filename)
//  mzidfile.store(filename , protein_ids, peptide_ids);

//  vector<ProteinIdentification> protein_ids2;
//  PeptideIdentificationList peptide_ids2;
//  mzidfile.load(filename, protein_ids2, peptide_ids2);

//  TEST_TRUE(protein_ids == protein_ids2)
//  TEST_TRUE(peptide_ids == peptide_ids2)

//  Misplaced Elements ignored in ParamGroup

//  Converting unknown score type to search engine specific score CV. #should not occur  scoretype is whatever

//  PSM without peptide evidences registered in the given search database found. This will cause an invalid MzIdentML file (which OpenMS still can consume). #might occur when reading idxml. no protein reference accession

//  No RT #might occurr when reading idxml. no rt to peptidehit
//  No MZ #might occurr when reading idxml. no mz to peptidehit

//  PeptideEvidence without reference to the positional in originating sequence found. #will always occur when reading idxml  no start end positional arguments
}
END_SECTION

START_SECTION(([EXTRA] XLMS data labeled cross-linker))
{
  vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  vector<ProteinIdentification> protein_ids2;
  PeptideIdentificationList peptide_ids2;

  std::string input_file= OPENMS_GET_TEST_DATA_PATH("MzIdentML_XLMS_labelled.mzid");
  MzIdentMLFile().load(input_file, protein_ids, peptide_ids);
  //
  TEST_EQUAL(peptide_ids[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1), 3)
  TEST_EQUAL(peptide_ids[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2), 4)
  TEST_EQUAL(peptide_ids[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_ALPHA), "ANYWHERE")
  TEST_EQUAL(peptide_ids[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_BETA_SEQUENCE), "SAVIKTSTR")
  TEST_EQUAL(peptide_ids[1].getHits()[0].getSequence().toString(), "FIVKASSGPR")

  // Reading and writing
  std::string filename;
  NEW_TMP_FILE(filename)
  MzIdentMLFile().store(filename, protein_ids, peptide_ids);
  MzIdentMLFile().load(filename, protein_ids2, peptide_ids2);

  // parameters from written and reloaded file
  // ProteinIdentification
  TEST_EQUAL(protein_ids2[0].getSearchParameters().fragment_mass_tolerance_ppm, false)
  TEST_EQUAL(protein_ids2[0].getSearchParameters().precursor_mass_tolerance_ppm, true)
  TEST_EQUAL(protein_ids2[0].getSearchParameters().getMetaValue("cross_link:residue1"), "[K, N-term]")
  TEST_EQUAL(protein_ids2[0].getSearchParameters().getMetaValue("cross_link:residue2"), "[K, N-term]")
  TEST_REAL_SIMILAR(StringUtils::toDouble(StringUtils::toStr(protein_ids2[0].getSearchParameters().getMetaValue("cross_link:mass"))), 138.0680796)
  TEST_REAL_SIMILAR(StringUtils::toDouble(StringUtils::toStr(protein_ids2[0].getSearchParameters().getMetaValue("cross_link:mass_isoshift"))), 12.075321)
  TEST_EQUAL(protein_ids2[0].getSearchParameters().getMetaValue("extra_features"), "precursor_mz_error_ppm,\
OpenPepXL:score,isotope_error,OpenPepXL:xquest_score,OpenPepXL:xcorr xlink,\
OpenPepXL:xcorr common,OpenPepXL:match-odds,OpenPepXL:intsum,OpenPepXL:wTIC,OpenPepXL:TIC,OpenPepXL:prescore,OpenPepXL:log_occupancy,\
OpenPepXL:log_occupancy_alpha,OpenPepXL:log_occupancy_beta,matched_xlink_alpha,matched_xlink_beta,matched_linear_alpha,\
matched_linear_beta,ppm_error_abs_sum_linear_alpha,ppm_error_abs_sum_linear_beta,ppm_error_abs_sum_xlinks_alpha,\
ppm_error_abs_sum_xlinks_beta,ppm_error_abs_sum_linear,ppm_error_abs_sum_xlinks,ppm_error_abs_sum_alpha,ppm_error_abs_sum_beta,\
ppm_error_abs_sum,precursor_total_intensity,precursor_target_intensity,precursor_signal_proportion,precursor_target_peak_count,\
precursor_residual_peak_count")
  TEST_EQUAL(protein_ids[0].getMetaValue("SpectrumIdentificationProtocol"), "MS:1002494") // crosslinking search


  // PeptideIdentification (Indices may change, without making the reading/writing invalid, if e.g. more is added to the test file)
  TEST_EQUAL(peptide_ids2.size(), 10)
  TEST_EQUAL(peptide_ids2[1].getRT(), peptide_ids2[2].getRT())
  TEST_REAL_SIMILAR(peptide_ids2[1].getRT(), 2132.4757)
  TEST_REAL_SIMILAR(peptide_ids2[1].getMZ(), 721.0845)
  TEST_EQUAL(peptide_ids2[1].getMetaValue(Constants::UserParam::SPECTRUM_REFERENCE), "spectrum=131,spectrum=113")

  // PeptideHit
  TEST_EQUAL(peptide_ids2[0].getHits().size(), 1)
  TEST_EQUAL(peptide_ids2[3].getHits().size(), 1)
  TEST_EQUAL(peptide_ids2[1].getHits().size(), 1)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE), "cross-link")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_RT), 2125.5966796875)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_MZ), 725.109252929687841)
  TEST_REAL_SIMILAR(peptide_ids2[1].getHits()[0].getScore(), -0.190406834856118)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getSequence().toString(), "FIVKASSGPR")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_BETA_SEQUENCE), "SAVIKTSTR")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1), 3)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2), 4)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_ALPHA), "ANYWHERE")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_BETA), "ANYWHERE")
  TEST_REAL_SIMILAR(StringUtils::toDouble(StringUtils::toStr(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_MASS))), 138.0680796)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_MOD), "DSS")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[0].annotation, "[alpha|ci$b2]")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[0].charge, 1)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[1].annotation, "[beta|ci$y2]")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[8].annotation, "[alpha|xi$b4]")
  TEST_EQUAL(peptide_ids2[0].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE), "mono-link")
  TEST_EQUAL(peptide_ids2[0].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1), 5)
  TEST_EQUAL(peptide_ids2[0].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2), "-")
}
END_SECTION

START_SECTION(([EXTRA] XLMS data unlabeled cross-linker))
{
  vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  vector<ProteinIdentification> protein_ids2;
  PeptideIdentificationList peptide_ids2;

  std::string input_file = OPENMS_GET_TEST_DATA_PATH("MzIdentML_XLMS_unlabelled.mzid");
  MzIdentMLFile().load(input_file, protein_ids, peptide_ids);

  // Reading and writing
  std::string filename;
  NEW_TMP_FILE(filename)
  MzIdentMLFile().store(filename, protein_ids, peptide_ids);
  MzIdentMLFile().load(filename, protein_ids2, peptide_ids2);

  // ProteinIdentification
  TEST_EQUAL(protein_ids2[0].getSearchParameters().fragment_mass_tolerance_ppm, true)
  TEST_EQUAL(protein_ids2[0].getSearchParameters().precursor_mass_tolerance_ppm, true)
  TEST_EQUAL(protein_ids2[0].getSearchParameters().getMetaValue("cross_link:residue1"), "[K, N-term]")
  TEST_EQUAL(protein_ids2[0].getSearchParameters().getMetaValue("cross_link:residue2"), "[K, N-term]")
  TEST_EQUAL(StringUtils::toDouble(StringUtils::toStr(protein_ids2[0].getSearchParameters().getMetaValue("cross_link:mass"))), 138.0680796)
  TEST_EQUAL(protein_ids[0].getMetaValue("SpectrumIdentificationProtocol"), "MS:1002494") // crosslinking search

  // PeptideIdentification (Indices may change, without making the reading/writing invalid, if e.g. more is added to the test file)
  TEST_EQUAL(peptide_ids2.size(), 3)
  TEST_REAL_SIMILAR(peptide_ids2[0].getRT(), 2175.3003)
  TEST_REAL_SIMILAR(peptide_ids2[0].getMZ(), 787.740356445313)
  TEST_EQUAL(peptide_ids2[0].getMetaValue(Constants::UserParam::SPECTRUM_REFERENCE), "controllerType=0 controllerNumber=1 scan=2395")

  // PeptideHit
  TEST_EQUAL(peptide_ids2[0].getHits().size(), 1)
  TEST_EQUAL(peptide_ids2[1].getHits().size(), 1)
  TEST_EQUAL(peptide_ids2[2].getHits().size(), 1)
  TEST_EQUAL(peptide_ids2[0].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE), "mono-link")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE), "cross-link")
  TEST_EQUAL(peptide_ids2[2].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE), "mono-link")

  TEST_EQUAL(peptide_ids2[0].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1), 5)
  TEST_EQUAL(peptide_ids2[0].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2), "-")
  TEST_EQUAL(peptide_ids2[0].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_ALPHA), "ANYWHERE")
  TEST_EQUAL(peptide_ids2[0].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_BETA), "ANYWHERE")

  TEST_EQUAL(peptide_ids2[1].getHits()[0].getSequence().toString(), "KNVPIEFPVIDR")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_BETA_SEQUENCE), "LGCKALHVLFER")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1), 0)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2), 3)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_MASS), 138.0680796)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_MOD), "DSS")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_ALPHA), "ANYWHERE")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_BETA), "ANYWHERE")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations().size(), 5)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[0].annotation, "[alpha|ci$y5]")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[0].charge, 1)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[1].annotation, "[alpha|ci$y7]")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[2].annotation, "[beta|ci$y7]")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[3].annotation, "[alpha|ci$y8]")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[2].charge, 1)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[4].charge, 2)

  TEST_EQUAL(peptide_ids2[2].getHits()[0].getSequence().toString(), "VEPSWLGPLFPDK(Xlink:DSS[156])TSNLR")
  TEST_EQUAL(peptide_ids2[2].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_BETA_SEQUENCE), "-")
  TEST_EQUAL(peptide_ids2[2].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1), 12)
  TEST_EQUAL(peptide_ids2[2].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2), "-")
}
END_SECTION

START_SECTION(([EXTRA] mzIdentML 1.3 crosslinking scores_and_thresholds))
{
  vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  std::string input_file = OPENMS_GET_TEST_DATA_PATH("MzIdentML_v1.3_crosslinking.mzid");
  MzIdentMLFile().load(input_file, protein_ids, peptide_ids);
  TEST_EQUAL(!protein_ids.empty(), true)
  TEST_EQUAL(!peptide_ids.empty(), true)

  // every parsed identification carries a spectrum reference; collect hits across all ids
  // (a SpectrumIdentificationResult may legitimately yield an id without hits, so guard the hit access)
  Size total_hits = 0;
  for (const PeptideIdentification& pid : peptide_ids)
  {
    TEST_EQUAL(pid.getSpectrumReference().empty(), false)
    total_hits += pid.getHits().size();
    if (!pid.getHits().empty())
    {
      TEST_EQUAL(pid.getHits()[0].getSequence().empty(), false)
    }
  }
  TEST_EQUAL(total_hits > 0, true)

  // store/load roundtrip: the writer must default to mzIdentML 1.3.0 and the reader must read it back without loss
  std::string filename;
  NEW_TMP_FILE(filename)
  MzIdentMLFile().store(filename, protein_ids, peptide_ids);

  // the default writer path emits version 1.3.0
  TextFile written(filename);
  bool has_version_130 = false;
  for (TextFile::ConstIterator it = written.begin(); it != written.end(); ++it)
  {
    if (StringUtils::hasSubstring(*it, "version=\"1.3.0\"")) { has_version_130 = true; break; }
  }
  TEST_EQUAL(has_version_130, true)

  vector<ProteinIdentification> protein_ids2;
  PeptideIdentificationList peptide_ids2;
  MzIdentMLFile().load(filename, protein_ids2, peptide_ids2);
  TEST_EQUAL(!peptide_ids2.empty(), true)
  TEST_EQUAL(!protein_ids2.empty(), true)
}
END_SECTION

START_SECTION(([EXTRA] mzIdentML 1.3 noncovalent association))
{
  vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  std::string input_file = OPENMS_GET_TEST_DATA_PATH("MzIdentML_v1.3_noncov_assoc.mzid");
  MzIdentMLFile().load(input_file, protein_ids, peptide_ids);
  TEST_EQUAL(!protein_ids.empty(), true)
  TEST_EQUAL(!peptide_ids.empty(), true)

  // parsed hits must have valid sequences (modification/cvParam fallbacks resolved correctly)
  Size total_hits = 0;
  for (const PeptideIdentification& pid : peptide_ids)
  {
    total_hits += pid.getHits().size();
    if (!pid.getHits().empty())
    {
      TEST_EQUAL(pid.getHits()[0].getSequence().empty(), false)
    }
  }
  TEST_EQUAL(total_hits > 0, true)

  // the ProteinDetectionList: each protein is listed once (it used to be added again for each
  // ProteinDetectionHypothesis) and each ProteinAmbiguityGroup becomes a protein group
  ABORT_IF(protein_ids.size() != 1)
  std::set<std::string> accessions;
  for (const ProteinHit& hit : protein_ids[0].getHits())
  {
    accessions.insert(hit.getAccession());
  }
  TEST_EQUAL(accessions.size(), protein_ids[0].getHits().size())
  TEST_EQUAL(protein_ids[0].getIndistinguishableProteins().size(), 2)
}
END_SECTION

START_SECTION(([EXTRA] mzIdentML 1.3 EDC crosslinking))
{
  vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  std::string input_file = OPENMS_GET_TEST_DATA_PATH("MzIdentML_v1.3_xlink_edc.mzid");
  // EDC files mix crosslinked and standalone (non-XL) peptides; both must parse without throwing
  MzIdentMLFile().load(input_file, protein_ids, peptide_ids);
  TEST_EQUAL(!protein_ids.empty(), true)
  TEST_EQUAL(!peptide_ids.empty(), true)

  Size total_hits = 0;
  for (const PeptideIdentification& pid : peptide_ids)
  {
    total_hits += pid.getHits().size();
  }
  TEST_EQUAL(total_hits > 0, true)
}
END_SECTION

START_SECTION(([EXTRA] mzIdentML 1.3 multiple spectra per identification))
{
  vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  std::string input_file = OPENMS_GET_TEST_DATA_PATH("MzIdentML_v1.3_multi_spectra.mzid");
  MzIdentMLFile().load(input_file, protein_ids, peptide_ids);
  TEST_EQUAL(!protein_ids.empty(), true)
  TEST_EQUAL(!peptide_ids.empty(), true)

  // multiple SpectrumIdentificationLists yield several distinct spectrum references
  set<std::string> refs;
  for (const PeptideIdentification& pid : peptide_ids)
  {
    refs.insert(pid.getSpectrumReference());
  }
  TEST_EQUAL(refs.size() > 1, true)
}
END_SECTION

START_SECTION(([EXTRA] PeptideSequence held in a CDATA section, behind a comment or padded with whitespace))
{
  // reference: the same identifications with the plain spelling '<PeptideSequence>PEPTIDER</PeptideSequence>'
  std::vector<ProteinIdentification> ref_protein_ids;
  PeptideIdentificationList ref_peptide_ids;
  MzIdentMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_whole.mzid"), ref_protein_ids, ref_peptide_ids);

  // same file, but one sequence sits in a CDATA section, one is preceded by an XML comment inside the element
  // and one is surrounded by whitespace; the reader used to discard the first two ("Non Text Node")
  std::vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  MzIdentMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_cdata_comment_peptide_sequence.mzid"), protein_ids, peptide_ids);

  TEST_EQUAL(protein_ids.size(), ref_protein_ids.size())
  ABORT_IF(peptide_ids.size() != ref_peptide_ids.size())
  std::multiset<std::string> all_sequences;
  for (Size i = 0; i < peptide_ids.size(); ++i)
  {
    std::multiset<std::string> sequences, ref_sequences;
    for (const PeptideHit& hit : peptide_ids[i].getHits())
    {
      sequences.insert(hit.getSequence().toString());
      all_sequences.insert(hit.getSequence().toString());
    }
    for (const PeptideHit& hit : ref_peptide_ids[i].getHits())
    {
      ref_sequences.insert(hit.getSequence().toString());
    }
    TEST_EQUAL(sequences == ref_sequences, true)
  }
  // the three rewritten peptides are referenced by PSMs and came back as plain residues
  TEST_EQUAL(all_sequences.count("PEPTIDER"), 1)    // CDATA section
  TEST_EQUAL(all_sequences.count("PEPTIDERR"), 1)   // padded with whitespace
  TEST_EQUAL(all_sequences.count("PEPTIDERRRR"), 2) // comment before the residues
  TEST_EQUAL(all_sequences.count(""), 0)
}
END_SECTION

START_SECTION(([EXTRA] Peptide with an empty PeptideSequence element))
{
  // '<PeptideSequence></PeptideSequence>' and '<PeptideSequence/>' are schema-valid and are what the writer
  // emits for a hit with an empty sequence; the reader used to dereference the missing first child node of
  // such an element and crash
  std::vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  MzIdentMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_empty_peptide_sequence.mzid"), protein_ids, peptide_ids);

  // the proteins reached through the PeptideEvidences are still there
  ABORT_IF(protein_ids.size() != 1)
  TEST_EQUAL(protein_ids[0].getHits().size(), 2)

  // every PSM is kept: the ones referencing an empty peptide carry an empty sequence, the others keep theirs
  ABORT_IF(peptide_ids.size() != 3)
  const std::vector<std::multiset<std::string>> expected = {
    {"", "PEPTIDERR"},  // spectrumID 17: PEPTIDER (now '<PeptideSequence></PeptideSequence>') and PEPTIDERR
    {"", "PEPTIDERRR"}, // spectrumID 45: PEPTIDERRRR (now '<PeptideSequence/>') and PEPTIDERRR
    {"", "PEPTIDERRR"}  // spectrumID 99: PEPTIDERRR and PEPTIDERRRR (now '<PeptideSequence/>')
  };
  for (Size i = 0; i < peptide_ids.size(); ++i)
  {
    std::multiset<std::string> sequences;
    Size n_empty = 0;
    for (const PeptideHit& hit : peptide_ids[i].getHits())
    {
      sequences.insert(hit.getSequence().toString());
      if (hit.getSequence().empty()) ++n_empty;
      TEST_EQUAL(hit.getPeptideEvidences().empty(), false) // the PSM was built completely
    }
    TEST_EQUAL(sequences == expected[i], true)
    TEST_EQUAL(n_empty, 1)
  }
}
END_SECTION

START_SECTION(([EXTRA] Peptide with an out-of-range SubstitutionModification location (CPP-169)))
{
  // A SubstitutionModification's "location" is a 1-based residue position;
  // location="0" (or one past the sequence length) used to write one byte
  // before (or past) the sequence's character buffer via
  // as[toInt32(location) - 1] with no bounds check. Reuse the empty-sequence
  // fixture and inject one into the "PEPTIDERR" peptide.
  std::ifstream orig_in(OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_empty_peptide_sequence.mzid"));
  std::stringstream orig_buf;
  orig_buf << orig_in.rdbuf();
  std::string content = orig_buf.str();

  std::string from = "<PeptideSequence>PEPTIDERR</PeptideSequence>";
  std::string to = "<PeptideSequence>PEPTIDERR</PeptideSequence>"
                    // A huge location (rather than "0") makes the pre-fix write land far
                    // outside any mapped memory, so this reliably crashes without the fix
                    // instead of landing one byte before the string's own small-object
                    // buffer, which in practice can go unnoticed in a non-sanitized build.
                    "<SubstitutionModification location=\"100000000\" originalResidue=\"P\" replacementResidue=\"A\"/>";
  Size pos = content.find(from);
  TEST_NOT_EQUAL(pos, std::string::npos)
  content.replace(pos, from.size(), to);

  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  std::ofstream of(tmp_filename.c_str());
  of << content;
  of.close();

  std::vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  MzIdentMLFile().load(tmp_filename, protein_ids, peptide_ids);

  // Must not crash, and must come back with the same identification/hit
  // structure as the unmodified fixture (see the empty-PeptideSequence test
  // above: 3 identifications, each with one empty-sequence PSM and one real
  // one -- the identification carrying "PEPTIDERR" is spectrumID 17, with
  // hits {"", "PEPTIDERR"}). The out-of-range SubstitutionModification is
  // rejected (caught by the same fallback as an unreadable Peptide), so that
  // identification's "PEPTIDERR" hit now also comes back empty -- {"", ""}
  // -- instead of either applying a one-byte-out-of-bounds write or silently
  // keeping "PEPTIDERR" unmodified.
  ABORT_IF(peptide_ids.size() != 3)
  Size n_empty = 0, n_hits = 0;
  for (const auto& pep_id : peptide_ids)
  {
    for (const PeptideHit& hit : pep_id.getHits())
    {
      ++n_hits;
      if (hit.getSequence().empty()) ++n_empty;
    }
  }
  TEST_EQUAL(n_hits, 6)
  TEST_EQUAL(n_empty, 4)
}
END_SECTION

START_SECTION(([EXTRA] store and load a PeptideHit with an empty sequence))
{
  // the writer emits '<PeptideSequence></PeptideSequence>' for such a hit; loading the file used to crash
  std::vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  MzIdentMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_whole.mzid"), protein_ids, peptide_ids);
  ABORT_IF(peptide_ids.empty() || peptide_ids[0].getHits().empty())

  PeptideHit empty_hit(0.5, 3, 1, AASequence());
  empty_hit.setPeptideEvidences(peptide_ids[0].getHits()[0].getPeptideEvidences()); // a protein reference is needed to be written
  peptide_ids[0].insertHit(empty_hit);
  Size n_hits = 0;
  for (const PeptideIdentification& pid : peptide_ids)
  {
    n_hits += pid.getHits().size();
  }

  std::string filename;
  NEW_TMP_FILE(filename)
  MzIdentMLFile().store(filename, protein_ids, peptide_ids);

  std::vector<ProteinIdentification> protein_ids2;
  PeptideIdentificationList peptide_ids2;
  MzIdentMLFile().load(filename, protein_ids2, peptide_ids2);

  TEST_EQUAL(peptide_ids2.size(), peptide_ids.size())
  Size n_hits2 = 0, n_empty2 = 0;
  for (const PeptideIdentification& pid : peptide_ids2)
  {
    for (const PeptideHit& hit : pid.getHits())
    {
      ++n_hits2;
      if (hit.getSequence().empty()) ++n_empty2;
    }
  }
  TEST_EQUAL(n_hits2, n_hits)
  TEST_EQUAL(n_empty2, 1)
}
END_SECTION

START_SECTION(([EXTRA] store and load a PeptideHit with a C-terminal modification))
{
  // mzIdentML puts a C-terminal modification at location length + 1; length is the last residue
  std::vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  MzIdentMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_whole.mzid"), protein_ids, peptide_ids);
  ABORT_IF(peptide_ids.empty() || peptide_ids[0].getHits().empty())

  PeptideHit cterm_hit(0.5, 3, 2, AASequence::fromString("PEPTIDEK.(Methyl)"));
  cterm_hit.setPeptideEvidences(peptide_ids[0].getHits()[0].getPeptideEvidences());
  peptide_ids[0].insertHit(cterm_hit);

  std::string filename;
  NEW_TMP_FILE(filename)
  MzIdentMLFile().store(filename, protein_ids, peptide_ids);

  TextFile stored;
  stored.load(filename, true);
  bool c_term_location = false, last_residue_location = false;
  for (const std::string& line : stored)
  {
    if (line == "<Modification location=\"9\">") c_term_location = true;
    if (line == "<Modification location=\"8\">") last_residue_location = true;
  }
  TEST_EQUAL(c_term_location, true)
  TEST_EQUAL(last_residue_location, false)

  std::vector<ProteinIdentification> protein_ids2;
  PeptideIdentificationList peptide_ids2;
  MzIdentMLFile().load(filename, protein_ids2, peptide_ids2);
  Size n_found = 0;
  for (const PeptideIdentification& pid : peptide_ids2)
  {
    for (const PeptideHit& hit : pid.getHits())
    {
      if (hit.getSequence().toUnmodifiedString() == "PEPTIDEK")
      {
        ++n_found;
        TEST_STRING_EQUAL(hit.getSequence().toString(), "PEPTIDEK.(Methyl)")
      }
    }
  }
  TEST_EQUAL(n_found, 1)
}
END_SECTION

START_SECTION(([EXTRA] score selection with Comet:expectation value and PSM-level q-value))
{
  // Comet:expectation value (MS:1002257) is lower-is-better and must not displace a PSM-level q-value.
  // Every result holds a better (PEPTIDER) and a worse (PEPTIDERR) hit.
  std::vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  MzIdentMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_comet_scores.mzid"), protein_ids, peptide_ids);
  TEST_EQUAL(peptide_ids.size(), 5)
  ABORT_IF(peptide_ids.size() != 5)

  struct Expected
  {
    std::string score_type;
    bool higher_better;
    double best_score;
  };
  const std::vector<Expected> expected{
    {"q-value", false, 0.0001},                // expectation value and q-value
    {"q-value", false, 0.0001},                // the same, q-value listed first
    {"Comet:expectation value", false, 1e-05}, // expectation value only
    {"q-value", false, 0.0001},                // q-value only
    {"Comet:xcorr", true, 3.5}                 // xcorr, expectation value and q-value: the raw score is used
  };
  for (Size i = 0; i < expected.size(); ++i)
  {
    PeptideIdentification& pid = peptide_ids[i];
    pid.sort();
    TEST_STRING_EQUAL(pid.getScoreType(), expected[i].score_type)
    TEST_EQUAL(pid.isHigherScoreBetter(), expected[i].higher_better)
    TEST_EQUAL(pid.getHits().size(), 2)
    if (pid.getHits().empty()) continue;
    TEST_STRING_EQUAL(pid.getHits()[0].getSequence().toString(), "PEPTIDER")
    TEST_REAL_SIMILAR(pid.getHits()[0].getScore(), expected[i].best_score)
  }

  // the expectation value stays available when the q-value is used
  TEST_REAL_SIMILAR(double(peptide_ids[0].getHits()[0].getMetaValue("MS:1002257")), 1e-05)
}
END_SECTION

START_SECTION(([EXTRA] score selection with X!Tandem and OMSSA scores))
{
  // X!Tandem:expect (MS:1001330), OMSSA:evalue (MS:1001328) and OMSSA:pvalue (MS:1001329) have no score order in
  // PSI-MS; they are lower-is-better. Every result holds a better (PEPTIDER) and a worse (PEPTIDERR) hit.
  std::vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  MzIdentMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_xtandem_omssa_scores.mzid"), protein_ids, peptide_ids);
  TEST_EQUAL(peptide_ids.size(), 6)
  ABORT_IF(peptide_ids.size() != 6)

  struct Expected
  {
    std::string score_type;
    bool higher_better;
    double best_score;
  };
  const std::vector<Expected> expected{
    {"X\\!Tandem:expect", false, 0.001},     // expect and hyperscore
    {"q-value", false, 0.005},               // expect, hyperscore and q-value: the expect yields to the q-value
    {"X\\!Tandem:hyperscore", true, 45.5},   // hyperscore only
    {"OMSSA:evalue", false, 0.002},          // E-value and p-value
    {"q-value", false, 0.004},               // E-value, p-value and q-value
    {"OMSSA:pvalue", false, 2e-06}           // p-value only
  };
  for (Size i = 0; i < expected.size(); ++i)
  {
    PeptideIdentification& pid = peptide_ids[i];
    pid.sort();
    TEST_STRING_EQUAL(pid.getScoreType(), expected[i].score_type)
    TEST_EQUAL(pid.isHigherScoreBetter(), expected[i].higher_better)
    TEST_EQUAL(pid.getHits().size(), 2)
    if (pid.getHits().empty()) continue;
    TEST_STRING_EQUAL(pid.getHits()[0].getSequence().toString(), "PEPTIDER")
    TEST_REAL_SIMILAR(pid.getHits()[0].getScore(), expected[i].best_score)
  }

  // the search engine scores stay available when the q-value is used
  TEST_REAL_SIMILAR(double(peptide_ids[1].getHits()[0].getMetaValue("MS:1001330")), 0.001)
  TEST_REAL_SIMILAR(double(peptide_ids[1].getHits()[0].getMetaValue("MS:1001331")), 45.5)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
START_SECTION(([EXTRA] native DOM namespace Unicode entity and malformed-input regressions))
{
  std::ifstream source(OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_whole.mzid"));
  const std::string original((std::istreambuf_iterator<char>(source)), std::istreambuf_iterator<char>());
  ABORT_IF(original.empty())
  auto save_xml = [](const std::string& path, const std::string& contents)
  {
    std::ofstream output(path, std::ios::binary);
    output << contents;
  };
  std::string filename;
  NEW_TMP_FILE(filename)
  std::vector<ProteinIdentification> proteins;
  PeptideIdentificationList peptides;

  // All element names may carry a namespace prefix; attributes stay unqualified.
  std::string prefixed = std::regex_replace(original, std::regex("<(/?)([A-Za-z][A-Za-z0-9]*)"), "<$1mzid:$2");
  StringUtils::substitute(prefixed, "xmlns=", "xmlns:mzid=");
  StringUtils::substitute(prefixed, "<mzid:DBSequence accession=", "<mzid:DBSequence xmlns:foreign=\"urn:foreign\" foreign:accession=\"ignored\" accession=");
  StringUtils::substitute(prefixed, "accession=\"PROT1\"", "accession=\"Pröt&#x1F9EC;&amp;&quot;\"");
  save_xml(filename, prefixed);
  MzIdentMLFile().load(filename, proteins, peptides);
  ABORT_IF(proteins.size() != 1)
  TEST_EQUAL(peptides.size(), 3)
  std::set<std::string> accessions;
  for (const auto& hit : proteins[0].getHits()) { accessions.insert(hit.getAccession()); }
  TEST_TRUE(accessions.contains("Pröt🧬&\""))
  TEST_TRUE(accessions.contains("PROT2"))

  // Foreign element names do not enter the mzIdentML identifier/score maps.
  const double reference_score = peptides[0].getHits()[0].getScore();
  std::string foreign = prefixed;
  foreign.insert(foreign.find("<mzid:DBSequence "),
    "<foreign:DBSequence xmlns:foreign='urn:unrelated' accession='POISON' id='5172702226708088127'>"
    "<foreign:Seq>AAAAAAAA</foreign:Seq></foreign:DBSequence>");
  const auto item_begin = foreign.find("<mzid:SpectrumIdentificationItem ");
  ABORT_IF(item_begin == std::string::npos)
  foreign.insert(foreign.find('>', item_begin) + 1,
    "<foreign:cvParam xmlns:foreign='urn:unrelated' accession='MS:1001171' cvRef='PSI-MS' name='Mascot:score' value='999999'/>");
  save_xml(filename, foreign);
  MzIdentMLFile().load(filename, proteins, peptides);
  ABORT_IF(proteins.size() != 1 || peptides.size() != 3 || peptides[0].getHits().empty())
  TEST_REAL_SIMILAR(peptides[0].getHits()[0].getScore(), reference_score)
  for (const auto& hit : proteins[0].getHits())
  {
    TEST_NOT_EQUAL(hit.getAccession(), "POISON")
    TEST_NOT_EQUAL(hit.getSequence(), "AAAAAAAA")
  }

  // An entity's declaration tree is not a second copy of its expanded content.
  std::string entity_wrapped = original;
  const auto sil_begin = entity_wrapped.find("<SpectrumIdentificationList ");
  const auto sil_end_tag = entity_wrapped.find("</SpectrumIdentificationList>", sil_begin);
  ABORT_IF(sil_begin == std::string::npos || sil_end_tag == std::string::npos)
  const auto sil_length = sil_end_tag + std::string("</SpectrumIdentificationList>").size() - sil_begin;
  std::string entity_value = entity_wrapped.substr(sil_begin, sil_length);
  StringUtils::substitute(entity_value, "'", "&#39;");
  StringUtils::substitute(entity_value, "%", "&#37;");
  entity_wrapped.replace(sil_begin, sil_length, "&identification_list;");
  entity_wrapped.insert(entity_wrapped.find("?>") + 2,
    "\n<!DOCTYPE MzIdentML [<!ENTITY identification_list '" + entity_value + "'>]>");
  save_xml(filename, entity_wrapped);
  MzIdentMLFile().load(filename, proteins, peptides);
  TEST_EQUAL(proteins.size(), 1)
  TEST_EQUAL(peptides.size(), 3)

  // A real existing local file is never loaded through an external entity.
  std::string secret_file;
  NEW_TMP_FILE(secret_file)
  save_xml(secret_file, "EXTERNAL_ENTITY_MUST_NOT_BE_READ");
  std::string external = original;
  const auto declaration = external.find("?>");
  ABORT_IF(declaration == std::string::npos)
  external.insert(declaration + 2, "\n<!DOCTYPE MzIdentML [<!ENTITY secret SYSTEM 'file://" + secret_file + "'>]>");
  StringUtils::substitute(external, "<Seq>PRGTEIN</Seq>", "<Seq>&secret;</Seq>");
  save_xml(filename, external);
  MzIdentMLFile().load(filename, proteins, peptides);
  for (const auto& protein : proteins)
  {
    for (const auto& hit : protein.getHits())
    {
      TEST_NOT_EQUAL(hit.getSequence(), "EXTERNAL_ENTITY_MUST_NOT_BE_READ")
    }
  }

  // Large ordinary text and attributes remain readable. HUGE
  // relaxes the 10 MB node limit but must retain the explicit entity budget.
  {
    const std::string large_sequence(10000001, 'A');
    std::string large = original;
    StringUtils::substitute(large, "<Seq>PRGTEIN</Seq>", "<Seq>" + large_sequence + "</Seq>");
    const auto root = large.find("<MzIdentML ");
    ABORT_IF(root == std::string::npos)
    large.insert(root + std::string("<MzIdentML").size(), " ignored=\"" + large_sequence + "\"");
    save_xml(filename, large);
    MzIdentMLFile().load(filename, proteins, peptides);
    ABORT_IF(proteins.size() != 1)
    bool found = false;
    for (const auto& hit : proteins[0].getHits())
    {
      if (hit.getAccession() == "PROT1")
      {
        found = true;
        TEST_EQUAL(hit.getSequence().size(), large_sequence.size())
        TEST_TRUE(hit.getSequence() == large_sequence)
      }
    }
    TEST_TRUE(found)
    TEST_EQUAL(peptides.size(), 3)
  }
  std::string amplified = original;
  std::string declarations = "\n<!DOCTYPE MzIdentML [<!ENTITY a '" + std::string(1000, 'A') + "'>";
  for (char name = 'b'; name <= 'e'; ++name)
  {
    declarations += "<!ENTITY " + std::string(1, name) + " '";
    for (int repeat = 0; repeat < 10; ++repeat) { declarations += "&" + std::string(1, name - 1) + ";"; }
    declarations += "'>";
  }
  declarations += "]>";
  const auto xml_declaration_end = amplified.find("?>");
  ABORT_IF(xml_declaration_end == std::string::npos)
  amplified.insert(xml_declaration_end + 2, declarations);
  // The same declaration and entity path must also work below the budget.
  std::string below_limit = amplified;
  StringUtils::substitute(below_limit, "<Seq>PRGTEIN</Seq>", "<Seq>&b;</Seq>");
  save_xml(filename, below_limit);
  MzIdentMLFile().load(filename, proteins, peptides);
  ABORT_IF(proteins.size() != 1)
  bool expanded_below_limit = false;
  for (const auto& hit : proteins[0].getHits())
  {
    if (hit.getAccession() == "PROT1")
    {
      expanded_below_limit = hit.getSequence() == std::string(10000, 'A');
    }
  }
  TEST_TRUE(expanded_below_limit)
  StringUtils::substitute(amplified, "<Seq>PRGTEIN</Seq>", "<Seq>&e;</Seq>");
  save_xml(filename, amplified);
  bool rejected_amplification = false;
  try
  {
    MzIdentMLFile().load(filename, proteins, peptides);
  }
  catch (const Exception::ParseError& error)
  {
    rejected_amplification = std::string(error.getMessage()).find("amplification") != std::string::npos;
  }
  TEST_TRUE(rejected_amplification)

  // Diagnostics keep the actual input location in the exception, not stderr.
  save_xml(filename, "<MzIdentML>\n<broken></MzIdentML>");
  bool reported_location = false;
  try
  {
    MzIdentMLFile().load(filename, proteins, peptides);
  }
  catch (const Exception::ParseError& error)
  {
    const std::string message = error.getMessage();
    reported_location = message.find(filename) != std::string::npos &&
      message.find("(line 2, column ") != std::string::npos;
  }
  TEST_TRUE(reported_location)

  // Do not traverse a partial DOM after a malformed tail, invalid reference or UTF8.
  save_xml(filename, original.substr(0, original.rfind("</MzIdentML>")));
  TEST_EXCEPTION(Exception::ParseError, MzIdentMLFile().load(filename, proteins, peptides))
  std::string malformed = original;
  StringUtils::substitute(malformed, "accession=\"PROT1\"", "accession=\"bad&-\"");
  save_xml(filename, malformed);
  TEST_EXCEPTION(Exception::ParseError, MzIdentMLFile().load(filename, proteins, peptides))
  malformed = original;
  StringUtils::substitute(malformed, "accession=\"PROT1\"", std::string("accession=\"bad") + char(0xff) + "\"");
  save_xml(filename, malformed);
  TEST_EXCEPTION(Exception::ParseError, MzIdentMLFile().load(filename, proteins, peptides))
}
END_SECTION

START_SECTION(([EXTRA] legacy DOM writer emits UTF8 namespaces and propagates output failures))
{
  ProteinIdentification protein;
  ProteinHit hit;
  hit.setAccession("Pröt🧬&\"");
  hit.setSequence("PEPTIDER");
  protein.insertHit(hit);
  protein.setSearchEngine("Mascot");
  protein.setSearchEngineVersion("2.1.0");
  const std::vector<ProteinIdentification> proteins{protein};
  const PeptideIdentificationList peptides;
  MzIdentMLFile logger;
  Internal::MzIdentMLDOMHandler writer(proteins, peptides, "1.1.0", logger);
  std::string filename;
  NEW_TMP_FILE(filename)
  writer.writeMzIdentMLFile(filename);
  std::ifstream stream(filename, std::ios::binary);
  const std::string xml((std::istreambuf_iterator<char>(stream)), std::istreambuf_iterator<char>());
  TEST_TRUE(xml.find("xmlns=\"http://psidev.info/psi/pi/mzIdentML/1.1\"") != std::string::npos)
  TEST_TRUE(xml.find("xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"") != std::string::npos)
  TEST_TRUE(xml.find("xsi:schemaLocation=") != std::string::npos)
  TEST_TRUE(xml.find("Pröt🧬&amp;&quot;") != std::string::npos)
  struct SoftwareAttributes : Internal::XMLHandler
  {
    SoftwareAttributes() : XMLHandler("memory", "1.1.0") {}
    Size count = 0;
    std::u16string engine, version;
    void onStartElement(const char16_t* name, const Internal::XMLAttributes& attributes) override
    {
      if (std::u16string_view(name) == u"AnalysisSoftware")
      {
        ++count;
        if (const auto* value = attributes.value("name")) { engine = value; }
        if (const auto* value = attributes.value("version")) { version = value; }
      }
    }
  };
  struct Reader : Internal::XMLFile
  {
    using XMLFile::parse_;
    using XMLFile::parseBuffer_;
  } reader;
  SoftwareAttributes explicit_version;
  reader.parseBuffer_(xml, &explicit_version);
  TEST_EQUAL(explicit_version.count, 1)
  TEST_TRUE(explicit_version.engine == u"Mascot")
  TEST_TRUE(explicit_version.version == u"2.1.0")

  // The reader ignores AnalysisSoftware entries with an empty version. Preserve
  // the real engine name and emit a nonempty fallback when no version is known.
  protein.setSearchEngineVersion("");
  const std::vector<ProteinIdentification> missing_version_proteins{protein};
  Internal::MzIdentMLDOMHandler missing_version_writer(missing_version_proteins, peptides, "1.1.0", logger);
  std::string missing_version_file;
  NEW_TMP_FILE(missing_version_file)
  missing_version_writer.writeMzIdentMLFile(missing_version_file);
  SoftwareAttributes missing_version;
  reader.parse_(missing_version_file, &missing_version);
  TEST_EQUAL(missing_version.count, 1)
  TEST_TRUE(missing_version.engine == u"Mascot")
  TEST_TRUE(missing_version.version == u"unknown")

  // An existing regular file cannot be used as the parent directory of an output.
  TEST_EXCEPTION(Exception::UnableToCreateFile, writer.writeMzIdentMLFile(filename + "/unwritable.mzid"))
}
END_SECTION

// ---------------------------------------------------------------------------------------------------------------------
// protein groups (ProteinDetectionList), issue #1105
// ---------------------------------------------------------------------------------------------------------------------

START_SECTION(([EXTRA] read protein groups from a ProteinDetectionList (issue #1105)))
{
  std::vector<ProteinIdentification> proteins;
  PeptideIdentificationList peptides;
  MzIdentMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzIdentML_protein_groups.mzid"), proteins, peptides);
  ABORT_IF(proteins.size() != 1)
  const ProteinIdentification& run = proteins[0];
  TEST_EQUAL(peptides.size(), 8)

  // each protein once, with the score of its ProteinDetectionHypothesis
  TEST_EQUAL(run.getHits().size(), 9)
  std::map<std::string, const ProteinHit*> hits;
  for (const ProteinHit& hit : run.getHits())
  {
    hits[hit.getAccession()] = &hit;
  }
  TEST_EQUAL(hits.size(), 9)
  ABORT_IF(!hits.contains("A") || !hits.contains("C") || !hits.contains("E") || !hits.contains("G"))
  TEST_EQUAL(run.getScoreType(), "Posterior Probability")
  TEST_EQUAL(run.isHigherScoreBetter(), true)
  TEST_REAL_SIMILAR(hits["A"]->getScore(), 0.98)
  TEST_REAL_SIMILAR(hits["C"]->getScore(), 0.2)
  TEST_REAL_SIMILAR(hits["G"]->getScore(), 0.1)
  TEST_REAL_SIMILAR(hits["A"]->getCoverage(), 93.3)
  TEST_EQUAL(run.getInferenceEngine(), "Epifany")
  TEST_EQUAL(run.getInferenceEngineVersion(), "3.6.0")

  // indistinguishable groups: the leading proteins (representative first) with the same evidence, or - without
  // annotation (PAG_3) - the proteins with the most peptides
  const auto& indist = run.getIndistinguishableProteins();
  TEST_EQUAL(indist.size(), 5)
  ABORT_IF(indist.size() != 5)
  TEST_EQUAL(ListUtils::concatenate(accessionsOf(indist[0]), ","), "B,A")
  TEST_REAL_SIMILAR(indist[0].probability, 0.97) // the group score of the PAG
  TEST_EQUAL(ListUtils::concatenate(accessionsOf(indist[1]), ","), "D")
  TEST_REAL_SIMILAR(indist[1].probability, 0.6)
  TEST_EQUAL(ListUtils::concatenate(accessionsOf(indist[2]), ","), "E,F")
  TEST_REAL_SIMILAR(indist[2].probability, 0.4) // no group score: the best protein score
  // PAG_4: leading proteins with different peptides are not indistinguishable
  TEST_EQUAL(ListUtils::concatenate(accessionsOf(indist[3]), ","), "H")
  TEST_REAL_SIMILAR(indist[3].probability, 0.9)
  TEST_EQUAL(ListUtils::concatenate(accessionsOf(indist[4]), ","), "I")
  TEST_REAL_SIMILAR(indist[4].probability, 0.88)

  // general groups: the unannotated PAG_3 (with its sub-set protein G), PAG_4 and cluster 1 (PAG_1 and PAG_2)
  const auto& general = run.getProteinGroups();
  TEST_EQUAL(general.size(), 3)
  ABORT_IF(general.size() != 3)
  TEST_EQUAL(ListUtils::concatenate(accessionsOf(general[0]), ","), "E,F,G")
  TEST_EQUAL(ListUtils::concatenate(accessionsOf(general[1]), ","), "H,I")
  TEST_REAL_SIMILAR(general[1].probability, 0.95) // the group score of the PAG
  TEST_EQUAL(ListUtils::concatenate(accessionsOf(general[2]), ","), "B,A,C,D")
  TEST_REAL_SIMILAR(general[2].probability, 0.97)

  // PAG_3 does not pass the threshold
  TEST_EQUAL(hits["E"]->getMetaValue("pass_threshold", "").toString(), "false")
  TEST_EQUAL(hits["A"]->metaValueExists("pass_threshold"), false)
  // the protein group-level FDR threshold
  TEST_REAL_SIMILAR(double(run.getMetaValue(Constants::UserParam::PROTEIN_GROUP_SCORE_THRESHOLD)), 0.05)
  TEST_EQUAL(run.getMetaValue(Constants::UserParam::PROTEIN_GROUP_SCORE_THRESHOLD_TYPE).toString(), "q-value")
}
END_SECTION

START_SECTION(([EXTRA] store and load protein groups (issue #1105)))
{
  std::vector<ProteinIdentification> proteins;
  PeptideIdentificationList peptides;
  makeInferenceRun(proteins, peptides);

  std::string filename;
  NEW_TMP_FILE_EXT(filename, ".mzid")
  MzIdentMLFile file;
  file.store(filename, proteins, peptides);

  // the file
  const std::string text = readFile(filename);
  TEST_EQUAL(countOf(text, "<ProteinDetectionList "), 1)
  TEST_EQUAL(countOf(text, "<ProteinDetection "), 1)
  TEST_EQUAL(countOf(text, "<ProteinDetectionProtocol "), 1)
  TEST_EQUAL(countOf(text, "<ProteinAmbiguityGroup "), 3)
  TEST_EQUAL(countOf(text, "<ProteinDetectionHypothesis "), 4)
  TEST_EQUAL(countOf(text, "name=\"leading protein\""), 4)
  TEST_EQUAL(countOf(text, "name=\"group representative\""), 3)
  TEST_EQUAL(countOf(text, "name=\"sequence same-set protein\""), 2)
  TEST_EQUAL(countOf(text, "name=\"EPIFANY:Protein posterior probability\""), 4)
  TEST_EQUAL(countOf(text, "name=\"protein group-level probability\""), 3)
  TEST_EQUAL(countOf(text, "name=\"EPIFANY\""), 1)
  // the decoy group is not reported as identified
  TEST_EQUAL(countOf(text, "name=\"protein group passes threshold\" value=\"true\""), 2)
  TEST_EQUAL(countOf(text, "name=\"protein group passes threshold\" value=\"false\""), 1)
  TEST_EQUAL(countOf(text, "name=\"count of identified proteins\" value=\"2\""), 1)
  TEST_EQUAL(countOf(text, "name=\"no threshold\""), 2) // PSMs and proteins
  // two clusters: A, B and C (with the score of the general group) and the decoy group
  TEST_EQUAL(countOf(text, "name=\"cluster identifier\" value=\"1\""), 2)
  TEST_EQUAL(countOf(text, "name=\"cluster identifier\" value=\"2\""), 1)
  TEST_EQUAL(countOf(text, "<userParam name=\"cluster score\""), 2)
  TEST_EQUAL(countOf(text, "name=\"count of identified clusters\" value=\"1\""), 1)
  // evidence: the top hits only (6 PeptideHypotheses: PEPTIDEK for A, B, C; ELVISK for A, B; GLYCINEK for DECOY_D),
  // not the second-ranked SAMPLER for C
  TEST_EQUAL(countOf(text, "<PeptideHypothesis "), 6)

  std::string version;
  std::stringstream schema_errors;
  TEST_EQUAL(file.isValid(filename, schema_errors, version), true)
  TEST_EQUAL(schema_errors.str(), "")
  // semantic validation (the mapping file has rules for the protein detection elements of mzIdentML 1.1+)
  StringList errors, warnings;
  TEST_EQUAL(file.isSemanticallyValid(filename, errors, warnings), true)
  TEST_EQUAL(ListUtils::concatenate(errors, "\n"), "")
  // ... which fails without the terms that mzIdentML 1.2 requires
  {
    std::string invalid;
    std::istringstream lines(text);
    for (std::string line; std::getline(lines, line);)
    {
      if (!StringUtils::hasSubstring(line, "name=\"protein group passes threshold\"") &&
          !StringUtils::hasSubstring(line, "name=\"count of identified proteins\""))
      {
        invalid += line + "\n";
      }
    }
    std::string invalid_file;
    NEW_TMP_FILE_EXT(invalid_file, ".mzid")
    std::ofstream(invalid_file) << invalid;
    StringList errors2, warnings2;
    TEST_EQUAL(file.isSemanticallyValid(invalid_file, errors2, warnings2), false)
    const std::string all_errors = ListUtils::concatenate(errors2, "\n");
    TEST_EQUAL(countOf(all_errors, "MzIdentML_ProteinAmbiguityGroupThreshold_rule"), 3) // once per group
    TEST_EQUAL(countOf(all_errors, "MzIdentML_ProteinDetectionListCount_rule"), 1)

    // proteins without a relationship to their group (leading protein etc.)
    std::string no_relationship;
    std::istringstream lines2(text);
    for (std::string line; std::getline(lines2, line);)
    {
      if (!StringUtils::hasSubstring(line, "name=\"leading protein\"") && !StringUtils::hasSubstring(line, "name=\"group representative\"") &&
          !StringUtils::hasSubstring(line, "name=\"sequence same-set protein\""))
      {
        no_relationship += line + "\n";
      }
    }
    std::ofstream(invalid_file) << no_relationship;
    StringList errors3, warnings3;
    TEST_EQUAL(file.isSemanticallyValid(invalid_file, errors3, warnings3), false)
    TEST_EQUAL(countOf(ListUtils::concatenate(errors3, "\n"), "MzIdentML_ProteinDetectionHypothesisRelationship_rule"), 4) // once per protein

    // mzIdentML 1.1 did not require these terms
    std::string version_1_1 = invalid;
    StringUtils::substitute(version_1_1, "mzIdentML/1.3", "mzIdentML/1.1");
    StringUtils::substitute(version_1_1, "version=\"1.3.0\"", "version=\"1.1.0\"");
    // the version of the search engine (in the header, after the root element) is not the mzIdentML version
    TEST_EQUAL(countOf(version_1_1, "version=\"2024.03.26\""), 1)
    StringUtils::substitute(version_1_1, "version=\"2024.03.26\"", "version=\"1.3.0\"");
    std::ofstream(invalid_file) << version_1_1;
    TEST_EQUAL(file.detectVersion(invalid_file), "1.1.0")
    StringList errors4, warnings4;
    TEST_EQUAL(file.isSemanticallyValid(invalid_file, errors4, warnings4), true)
    TEST_EQUAL(ListUtils::concatenate(errors4, "\n"), "")
  }

  // round trip
  std::vector<ProteinIdentification> proteins2;
  PeptideIdentificationList peptides2;
  file.load(filename, proteins2, peptides2);
  ABORT_IF(proteins2.size() != 1)
  const ProteinIdentification& run = proteins[0];
  const ProteinIdentification& run2 = proteins2[0];
  TEST_EQUAL(run2.getHits().size(), 4) // no duplicates
  TEST_EQUAL(run2.getScoreType(), "Posterior Probability")
  TEST_EQUAL(run2.isHigherScoreBetter(), true)
  TEST_EQUAL(run2.getInferenceEngine(), "Epifany")
  TEST_EQUAL(run2.getInferenceEngineVersion(), "3.6.0")
  std::map<std::string, const ProteinHit*> hits2;
  for (const ProteinHit& hit : run2.getHits())
  {
    hits2[hit.getAccession()] = &hit;
  }
  ABORT_IF(!hits2.contains("A") || !hits2.contains("B") || !hits2.contains("C") || !hits2.contains("DECOY_D"))
  for (const ProteinHit& hit : run.getHits())
  {
    TEST_REAL_SIMILAR(hits2[hit.getAccession()]->getScore(), hit.getScore())
  }
  TEST_REAL_SIMILAR(hits2["A"]->getCoverage(), 93.3)
  TEST_EQUAL(hits2["DECOY_D"]->getMetaValue("target_decoy", "").toString(), "decoy")
  TEST_EQUAL(hits2["DECOY_D"]->getMetaValue("pass_threshold", "").toString(), "false")
  TEST_EQUAL(run2.getIndistinguishableProteins().size(), 3)
  ABORT_IF(run2.getIndistinguishableProteins().size() != 3)
  for (Size i = 0; i < 3; ++i)
  {
    TEST_EQUAL(ListUtils::concatenate(accessionsOf(run2.getIndistinguishableProteins()[i]), ","),
               ListUtils::concatenate(accessionsOf(run.getIndistinguishableProteins()[i]), ","))
    TEST_REAL_SIMILAR(run2.getIndistinguishableProteins()[i].probability, run.getIndistinguishableProteins()[i].probability)
  }
  TEST_EQUAL(run2.getProteinGroups().size(), 1)
  ABORT_IF(run2.getProteinGroups().size() != 1)
  TEST_EQUAL(ListUtils::concatenate(accessionsOf(run2.getProteinGroups()[0]), ","), "A,B,C")
  TEST_REAL_SIMILAR(run2.getProteinGroups()[0].probability, 0.99)

  // storing again gives the same protein groups
  std::string filename2;
  NEW_TMP_FILE_EXT(filename2, ".mzid")
  file.store(filename2, proteins2, peptides2);
  const std::string text2 = readFile(filename2);
  TEST_EQUAL(countOf(text2, "<ProteinAmbiguityGroup "), 3)
  TEST_EQUAL(countOf(text2, "name=\"protein group passes threshold\" value=\"true\""), 2)
  TEST_EQUAL(countOf(text2, "<userParam name=\"cluster score\""), 2)
}
END_SECTION

START_SECTION(([EXTRA] protein groups: score type aliases and group scores only (issue #1105)))
{
  // a score type alias ("qvalue" for "q-value") is written with the q-value terms and read back
  std::vector<ProteinIdentification> proteins;
  PeptideIdentificationList peptides;
  makeInferenceRun(proteins, peptides);
  ProteinIdentification& run = proteins[0];
  run.setScoreType("qvalue");
  run.setHigherScoreBetter(false);
  for (Size i = 0; i < run.getHits().size(); ++i)
  {
    run.getHits()[i].setScore(0.001 * (i + 1));
  }
  for (Size i = 0; i < run.getIndistinguishableProteins().size(); ++i)
  {
    run.getIndistinguishableProteins()[i].probability = 0.01 * (i + 1);
  }
  run.setMetaValue(Constants::UserParam::PROTEIN_GROUP_SCORE_THRESHOLD, 0.025);
  run.setMetaValue(Constants::UserParam::PROTEIN_GROUP_SCORE_THRESHOLD_TYPE, "q-value"); // another alias
  std::string filename;
  NEW_TMP_FILE_EXT(filename, ".mzid")
  MzIdentMLFile().store(filename, proteins, peptides);
  const std::string text = readFile(filename);
  TEST_EQUAL(countOf(text, "name=\"protein-level q-value\""), 4)
  TEST_EQUAL(countOf(text, "name=\"protein group-level q-value\""), 3)
  TEST_EQUAL(countOf(text, "name=\"protein group-level global FDR\" value=\"0.025\""), 1)
  TEST_EQUAL(countOf(text, "name=\"protein group passes threshold\" value=\"true\""), 2) // the decoy group fails

  std::vector<ProteinIdentification> proteins2;
  PeptideIdentificationList peptides2;
  MzIdentMLFile().load(filename, proteins2, peptides2);
  ABORT_IF(proteins2.size() != 1 || proteins2[0].getIndistinguishableProteins().size() != 3)
  TEST_EQUAL(proteins2[0].getScoreType(), "qvalue")
  TEST_EQUAL(proteins2[0].isHigherScoreBetter(), false)
  for (const ProteinHit& hit : proteins2[0].getHits())
  {
    TEST_EQUAL(hit.getScore() > 0.0, true)
  }
  TEST_REAL_SIMILAR(proteins2[0].getIndistinguishableProteins()[1].probability, 0.02)
  TEST_EQUAL(proteins2[0].getMetaValue(Constants::UserParam::PROTEIN_GROUP_SCORE_THRESHOLD_TYPE).toString(), "qvalue")

  // a file that reports only group scores (no protein scores, no OpenMS AnalysisParams)
  std::string group_scores_only;
  std::istringstream lines(text);
  for (std::string line; std::getline(lines, line);)
  {
    if (!StringUtils::hasSubstring(line, "name=\"protein-level q-value\"") && !StringUtils::hasSubstring(line, "userParam name=\"protein score type\"") &&
        !StringUtils::hasSubstring(line, "userParam name=\"higher protein score better\""))
    {
      group_scores_only += line + "\n";
    }
  }
  std::string filename2;
  NEW_TMP_FILE_EXT(filename2, ".mzid")
  std::ofstream(filename2) << group_scores_only;
  std::vector<ProteinIdentification> proteins3;
  PeptideIdentificationList peptides3;
  MzIdentMLFile().load(filename2, proteins3, peptides3);
  ABORT_IF(proteins3.size() != 1 || proteins3[0].getIndistinguishableProteins().size() != 3)
  TEST_EQUAL(proteins3[0].getScoreType(), "q-value")
  TEST_EQUAL(proteins3[0].isHigherScoreBetter(), false)
  {
    // group statistics without an order in the CV (here: e-values) are error measures: lower is better
    std::string evalues_only = group_scores_only;
    StringUtils::substitute(evalues_only, "accession=\"MS:1002373\" cvRef=\"PSI-MS\" name=\"protein group-level q-value\"",
                            "accession=\"MS:1002372\" cvRef=\"PSI-MS\" name=\"protein group-level e-value\"");
    std::string filename3;
    NEW_TMP_FILE_EXT(filename3, ".mzid")
    std::ofstream(filename3) << evalues_only;
    std::vector<ProteinIdentification> proteins4;
    PeptideIdentificationList peptides4;
    MzIdentMLFile().load(filename3, proteins4, peptides4);
    ABORT_IF(proteins4.size() != 1 || proteins4[0].getIndistinguishableProteins().size() != 3)
    TEST_EQUAL(proteins4[0].getScoreType(), "protein group-level e-value")
    TEST_EQUAL(proteins4[0].isHigherScoreBetter(), false)
    TEST_REAL_SIMILAR(proteins4[0].getIndistinguishableProteins()[1].probability, 0.02)
  }
  TEST_EQUAL(proteins3[0].isHigherScoreBetter(), false)
  for (Size i = 0; i < 3; ++i)
  {
    TEST_REAL_SIMILAR(proteins3[0].getIndistinguishableProteins()[i].probability, 0.01 * (i + 1))
  }
}
END_SECTION

START_SECTION(([EXTRA] read protein groups: clusters and leading proteins (issue #1105)))
{
  const std::string fixture = readFile(OPENMS_GET_TEST_DATA_PATH("MzIdentML_protein_groups.mzid"));
  auto load = [](const std::string& text, std::vector<ProteinIdentification>& proteins)
  {
    std::string filename;
    NEW_TMP_FILE_EXT(filename, ".mzid")
    std::ofstream(filename) << text;
    PeptideIdentificationList peptides;
    MzIdentMLFile().load(filename, proteins, peptides);
  };

  // PAG_1 alone in its cluster: one indistinguishable group, but the cluster also holds the sub-set protein C
  {
    std::string text = fixture;
    const std::string pag2_cluster = "<cvParam cvRef=\"PSI-MS\" accession=\"MS:1002407\" name=\"cluster identifier\" value=\"1\"/>\n        </ProteinAmbiguityGroup>\n        <ProteinAmbiguityGroup id=\"PAG_3\">";
    ABORT_IF(text.find(pag2_cluster) == std::string::npos)
    std::string pag2_cluster2 = pag2_cluster;
    StringUtils::substitute(pag2_cluster2, "value=\"1\"", "value=\"2\"");
    text.replace(text.find(pag2_cluster), pag2_cluster.size(), pag2_cluster2);
    std::vector<ProteinIdentification> proteins;
    load(text, proteins);
    ABORT_IF(proteins.size() != 1)
    std::vector<std::string> general;
    for (const auto& group : proteins[0].getProteinGroups())
    {
      general.push_back(ListUtils::concatenate(accessionsOf(group), ","));
    }
    TEST_EQUAL(ListUtils::concatenate(general, ";"), "E,F,G;H,I;B,A,C") // cluster 2 (D alone) is no general group
  }

  // leading proteins annotated as same-set proteins form one indistinguishable group, even with different peptides
  {
    std::string text = fixture;
    const std::string leading_i = "<cvParam cvRef=\"PSI-MS\" accession=\"MS:1002401\" name=\"leading protein\"/>\n            <cvParam cvRef=\"PSI-MS\" accession=\"MS:1003119\" name=\"EPIFANY:Protein posterior probability\" value=\"0.88\"/>";
    ABORT_IF(text.find(leading_i) == std::string::npos)
    text.replace(text.find(leading_i), leading_i.size(),
      "<cvParam cvRef=\"PSI-MS\" accession=\"MS:1002401\" name=\"leading protein\"/>\n            <cvParam cvRef=\"PSI-MS\" accession=\"MS:1001594\" name=\"sequence same-set protein\" value=\"PDH_H\"/>\n            <cvParam cvRef=\"PSI-MS\" accession=\"MS:1003119\" name=\"EPIFANY:Protein posterior probability\" value=\"0.88\"/>");
    std::vector<ProteinIdentification> proteins;
    load(text, proteins);
    ABORT_IF(proteins.size() != 1 || proteins[0].getIndistinguishableProteins().size() != 4)
    TEST_EQUAL(ListUtils::concatenate(accessionsOf(proteins[0].getIndistinguishableProteins()[3]), ","), "H,I")
    TEST_REAL_SIMILAR(proteins[0].getIndistinguishableProteins()[3].probability, 0.95) // the group score of the PAG
  }

  // the same-set term is optional: leading proteins without it are compared by their peptides, even if another
  // leading protein of the group is annotated
  {
    std::string text = fixture;
    for (const std::string& same_set : {std::string("\n            <cvParam cvRef=\"PSI-MS\" accession=\"MS:1001594\" name=\"sequence same-set protein\" value=\"PDH_B\"/>"),
                                        std::string("\n            <cvParam cvRef=\"PSI-MS\" accession=\"MS:1001594\" name=\"sequence same-set protein\" value=\"PDH_A\"/>")})
    {
      ABORT_IF(text.find(same_set) == std::string::npos)
      text.erase(text.find(same_set), same_set.size());
    }
    // C becomes a leading protein, annotated as same-set with a protein that is not in the group
    const std::string c_roles = "<cvParam cvRef=\"PSI-MS\" accession=\"MS:1002402\" name=\"non-leading protein\"/>\n            <cvParam cvRef=\"PSI-MS\" accession=\"MS:1001596\" name=\"sequence sub-set protein\" value=\"PDH_A PDH_B\"/>";
    ABORT_IF(text.find(c_roles) == std::string::npos)
    text.replace(text.find(c_roles), c_roles.size(),
      "<cvParam cvRef=\"PSI-MS\" accession=\"MS:1002401\" name=\"leading protein\"/>\n            <cvParam cvRef=\"PSI-MS\" accession=\"MS:1001594\" name=\"sequence same-set protein\" value=\"PDH_X\"/>");
    std::vector<ProteinIdentification> proteins;
    load(text, proteins);
    ABORT_IF(proteins.size() != 1 || proteins[0].getIndistinguishableProteins().size() != 6)
    TEST_EQUAL(ListUtils::concatenate(accessionsOf(proteins[0].getIndistinguishableProteins()[0]), ","), "B,A") // same peptides
    TEST_EQUAL(ListUtils::concatenate(accessionsOf(proteins[0].getIndistinguishableProteins()[1]), ","), "C")
  }

  // an optional userParam whose value does not match its type does not stop the file from loading
  {
    std::string text = fixture;
    const std::string pag4_end = "name=\"protein group-level probability\" value=\"0.95\"/>";
    ABORT_IF(text.find(pag4_end) == std::string::npos)
    text.insert(text.find(pag4_end) + pag4_end.size(), "\n          <userParam name=\"cluster score\" type=\"xsd:double\" value=\"unavailable\"/>");
    std::vector<ProteinIdentification> proteins;
    load(text, proteins);
    ABORT_IF(proteins.size() != 1)
    TEST_EQUAL(proteins[0].getIndistinguishableProteins().size(), 5)
  }
}
END_SECTION

START_SECTION(([EXTRA] store protein groups with a recorded group score threshold (issue #1105)))
{
  std::vector<ProteinIdentification> proteins;
  PeptideIdentificationList peptides;
  makeInferenceRun(proteins, peptides);
  // as IDFilter -score:proteingroup 0.6: only the group of A and B has a better score
  proteins[0].setMetaValue(Constants::UserParam::PROTEIN_GROUP_SCORE_THRESHOLD, 0.6);
  proteins[0].setMetaValue(Constants::UserParam::PROTEIN_GROUP_SCORE_THRESHOLD_TYPE, "Posterior Probability");

  std::string filename;
  NEW_TMP_FILE_EXT(filename, ".mzid")
  MzIdentMLFile().store(filename, proteins, peptides);
  const std::string text = readFile(filename);
  TEST_EQUAL(countOf(text, "name=\"protein group-level statistical threshold\" value=\"0.6\""), 1)
  TEST_EQUAL(countOf(text, "name=\"no threshold\""), 1) // only for the PSMs
  TEST_EQUAL(countOf(text, "name=\"protein group passes threshold\" value=\"true\""), 1)
  TEST_EQUAL(countOf(text, "name=\"count of identified proteins\" value=\"1\""), 1)

  std::vector<ProteinIdentification> proteins2;
  PeptideIdentificationList peptides2;
  MzIdentMLFile().load(filename, proteins2, peptides2);
  ABORT_IF(proteins2.size() != 1)
  TEST_REAL_SIMILAR(double(proteins2[0].getMetaValue(Constants::UserParam::PROTEIN_GROUP_SCORE_THRESHOLD)), 0.6)
  TEST_EQUAL(proteins2[0].getMetaValue(Constants::UserParam::PROTEIN_GROUP_SCORE_THRESHOLD_TYPE).toString(), "Posterior Probability")

  // a threshold on another score type is not used
  proteins[0].setMetaValue(Constants::UserParam::PROTEIN_GROUP_SCORE_THRESHOLD_TYPE, "q-value");
  MzIdentMLFile().store(filename, proteins, peptides);
  const std::string text2 = readFile(filename);
  TEST_EQUAL(countOf(text2, "name=\"no threshold\""), 2)
  TEST_EQUAL(countOf(text2, "name=\"protein group passes threshold\" value=\"true\""), 2)

  // nor is one for a run without indistinguishable groups (one group per protein): no group was filtered with it,
  // and the protein scores may have been recalculated since
  proteins[0].setMetaValue(Constants::UserParam::PROTEIN_GROUP_SCORE_THRESHOLD_TYPE, "Posterior Probability");
  proteins[0].getIndistinguishableProteins().clear();
  proteins[0].getProteinGroups().clear();
  MzIdentMLFile().store(filename, proteins, peptides);
  const std::string text3 = readFile(filename);
  TEST_EQUAL(countOf(text3, "name=\"no threshold\""), 2)
  TEST_EQUAL(countOf(text3, "<ProteinAmbiguityGroup "), 4) // A, B, C, DECOY_D
  TEST_EQUAL(countOf(text3, "name=\"protein group passes threshold\" value=\"true\""), 3) // all but the decoy
}
END_SECTION

START_SECTION(([EXTRA] store protein groups: only one run may have protein inference results (issue #1105)))
{
  std::vector<ProteinIdentification> proteins;
  PeptideIdentificationList peptides;
  makeInferenceRun(proteins, peptides);
  makeInferenceRun(proteins, peptides);
  proteins[1].setIdentifier("inference_run_2");
  for (Size i = peptides.size() / 2; i < peptides.size(); ++i)
  {
    peptides[i].setIdentifier("inference_run_2");
  }
  std::string filename;
  NEW_TMP_FILE_EXT(filename, ".mzid")
  MzIdentMLFile().store(filename, proteins, peptides);
  const std::string text = readFile(filename);
  TEST_EQUAL(countOf(text, "<ProteinDetectionList "), 0) // two lists are not allowed: no protein groups are written
  TEST_EQUAL(countOf(text, "<SpectrumIdentificationList "), 2) // but the PSMs are

  // without protein inference, there are no protein results to write
  std::vector<ProteinIdentification> search_only(1, proteins[0]);
  search_only[0].getIndistinguishableProteins().clear();
  search_only[0].getProteinGroups().clear();
  search_only[0].getSearchParameters().removeMetaValue("InferenceEngine");
  PeptideIdentificationList search_peptides(std::vector<PeptideIdentification>(peptides.begin(), peptides.begin() + peptides.size() / 2));
  MzIdentMLFile().store(filename, search_only, search_peptides);
  TEST_EQUAL(countOf(readFile(filename), "<ProteinDetectionList "), 0)

  // protein inference without groups: every protein is a group of its own
  std::vector<ProteinIdentification> ungrouped(1, proteins[0]);
  ungrouped[0].getIndistinguishableProteins().clear();
  ungrouped[0].getProteinGroups().clear();
  MzIdentMLFile().store(filename, ungrouped, search_peptides);
  const std::string text2 = readFile(filename);
  TEST_EQUAL(countOf(text2, "<ProteinAmbiguityGroup "), 4)
  TEST_EQUAL(countOf(text2, "name=\"cluster identifier\""), 0)
  std::vector<ProteinIdentification> proteins2;
  PeptideIdentificationList peptides2;
  MzIdentMLFile().load(filename, proteins2, peptides2);
  ABORT_IF(proteins2.size() != 1)
  TEST_EQUAL(proteins2[0].getIndistinguishableProteins().size(), 4)
  TEST_EQUAL(proteins2[0].getProteinGroups().size(), 0)
}
END_SECTION

START_SECTION(std::string detectVersion(const std::string& filename) const)
{
  // only the root element declares the version (from its version attribute, or else from its namespace)
  const std::vector<std::pair<std::string, std::string>> cases = {
    // a comment with a root-like tag before the root
    {"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<!-- e.g. <MzIdentML version=\"1.3.0\"> -->\n"
     "<MzIdentML xmlns=\"http://psidev.info/psi/pi/mzIdentML/1.1\" version=\"1.1.0\" id=\"x\"></MzIdentML>\n", "1.1.0"},
    // '>' in an attribute value before the version
    {"<MzIdentML name=\"A > B\" xmlns=\"http://psidev.info/psi/pi/mzIdentML/1.1\" version=\"1.1.0\" id=\"x\"></MzIdentML>\n", "1.1.0"},
    // another element's version attribute
    {"<MzIdentML xmlns=\"http://psidev.info/psi/pi/mzIdentML/1.1\" version=\"1.1.0\" id=\"x\">\n<AnalysisSoftwareList>"
     "<AnalysisSoftware id=\"s\" version=\"1.3.0\"/></AnalysisSoftwareList></MzIdentML>\n", "1.1.0"},
    // namespace of a prefixed root without version attribute
    {"<mzid:MzIdentML xmlns:mzid=\"http://psidev.info/psi/pi/mzIdentML/1.2\" id=\"x\"></mzid:MzIdentML>\n", "1.2.0"},
    // mzIdentML 1.0
    {"<mzIdentML xmlns=\"http://psidev.info/psi/pi/mzIdentML/1.0\" version=\"1.0.0\" id=\"x\"></mzIdentML>\n", "1.0.0"},
    // no version: the default
    {"<MzIdentML id=\"x\"></MzIdentML>\n", "1.3.0"},
    // not XML: the default
    {"no XML\n", "1.3.0"}};
  std::string filename;
  NEW_TMP_FILE_EXT(filename, ".mzid")
  for (const auto& [content, version] : cases)
  {
    std::ofstream(filename) << content;
    TEST_EQUAL(MzIdentMLFile().detectVersion(filename), version)
  }
  TEST_EXCEPTION(Exception::FileNotFound, MzIdentMLFile().detectVersion("this_file_does_not_exist.mzid"))
}
END_SECTION

START_SECTION(([EXTRA] store PeptideEvidences per evidence, also with several runs (issue #10983)))
{
  // two runs searched against different databases. Both identify PEPTIDEK: from P1 in run 1 and from P2 in run 2. Both
  // runs also have the protein SHARED (with a different peptide), and run 1 identifies PEPTIDEK twice.
  std::vector<ProteinIdentification> proteins;
  PeptideIdentificationList peptides;
  auto add_run = [&](const std::string& identifier, const std::string& db, const std::string& protein)
  {
    ProteinIdentification run;
    run.setIdentifier(identifier);
    run.setSearchEngine("MS-GF+");
    ProteinIdentification::SearchParameters search_parameters;
    search_parameters.db = db;
    run.setSearchParameters(search_parameters);
    run.setMetaValue("spectra_data", StringList{identifier + ".mzML"});
    run.insertHit(ProteinHit(0.0, 1, protein, ""));
    run.insertHit(ProteinHit(0.0, 1, "SHARED", ""));
    proteins.push_back(run);
  };
  auto add_psm = [&](const std::string& identifier, const std::string& spectrum, const std::string& sequence, const PeptideEvidence& evidence)
  {
    PeptideIdentification id;
    id.setIdentifier(identifier);
    id.setScoreType("q-value");
    id.setHigherScoreBetter(false);
    id.setSpectrumReference(spectrum);
    id.setRT(100.0);
    id.setMZ(500.0);
    PeptideHit hit(0.01, 0, 2, AASequence::fromString(sequence));
    hit.addPeptideEvidence(evidence);
    id.insertHit(hit);
    peptides.push_back(id);
  };
  add_run("run_one", "db_one.fasta", "P1");
  add_run("run_two", "db_two.fasta", "P2");
  add_psm("run_one", "scan=1", "PEPTIDEK", PeptideEvidence("P1", 10, 17, 'K', 'A'));
  add_psm("run_one", "scan=2", "PEPTIDEK", PeptideEvidence("P1", 10, 17, 'K', 'A'));
  add_psm("run_one", "scan=3", "ELVISK", PeptideEvidence("SHARED", 0, 5, '[', 'G'));
  add_psm("run_two", "scan=7", "PEPTIDEK", PeptideEvidence("P2", 40, 47, 'R', 'G'));
  add_psm("run_two", "scan=8", "ELVISK", PeptideEvidence("SHARED", 0, 5, '[', 'G'));
  // same peptide, protein and run, but one other field of the evidence each: position, flanking residue, decoy flag
  add_psm("run_one", "scan=4", "PEPTIDEK", PeptideEvidence("P1", 20, 27, 'K', 'A'));
  add_psm("run_one", "scan=5", "PEPTIDEK", PeptideEvidence("P1", 10, 17, 'R', 'A'));
  add_psm("run_one", "scan=6", "PEPTIDEK", PeptideEvidence("P1", 10, 17, 'K', 'A'));
  peptides.back().getHits()[0].setMetaValue(Constants::UserParam::TARGET_DECOY, "decoy");

  std::string filename;
  NEW_TMP_FILE_EXT(filename, ".mzid")
  MzIdentMLFile().store(filename, proteins, peptides);
  const std::string text = readFile(filename);
  // one PeptideEvidence per evidence: P1 (shared by the first two PSMs of run 1, and three variants of it), P2, and SHARED
  // in each database
  TEST_EQUAL(countOf(text, "<PeptideEvidence "), 7)
  // the PeptideEvidence that each spectrum references, as written
  std::map<std::string, std::string> pev_attributes; // PeptideEvidence id -> its attributes
  const std::regex pev_element(R"re(<PeptideEvidence id="([^"]*)"([^>]*)/>)re");
  for (std::sregex_iterator m(text.begin(), text.end(), pev_element), end; m != end; ++m)
  {
    pev_attributes[(*m)[1]] = (*m)[2];
  }
  std::map<std::string, std::string> pev_of; // spectrum -> id of the PeptideEvidence it references
  const std::regex pev_reference(R"re(spectrumID="([^"]*)"[^>]*>\s*<SpectrumIdentificationItem[^>]*>\s*<PeptideEvidenceRef peptideEvidence_ref="([^"]*)")re");
  for (std::sregex_iterator m(text.begin(), text.end(), pev_reference), end; m != end; ++m)
  {
    pev_of[(*m)[1]] = (*m)[2];
  }
  TEST_EQUAL(pev_of.size(), 8)
  TEST_EQUAL(pev_of["scan=1"], pev_of["scan=2"])
  TEST_EQUAL(StringUtils::hasSubstring(pev_attributes[pev_of["scan=1"]], "post=\"A\" pre=\"K\" start=\"11\" end=\"18\""), true)
  TEST_EQUAL(StringUtils::hasSubstring(pev_attributes[pev_of["scan=1"]], "isDecoy"), false)
  TEST_EQUAL(StringUtils::hasSubstring(pev_attributes[pev_of["scan=4"]], "post=\"A\" pre=\"K\" start=\"21\" end=\"28\""), true)
  TEST_EQUAL(StringUtils::hasSubstring(pev_attributes[pev_of["scan=5"]], "post=\"A\" pre=\"R\" start=\"11\" end=\"18\""), true)
  TEST_EQUAL(StringUtils::hasSubstring(pev_attributes[pev_of["scan=6"]], "post=\"A\" pre=\"K\" start=\"11\" end=\"18\" isDecoy=\"1\""), true)
  TEST_EQUAL(StringUtils::hasSubstring(pev_attributes[pev_of["scan=7"]], "post=\"G\" pre=\"R\" start=\"41\" end=\"48\""), true)
  TEST_EQUAL(std::set<std::string>({pev_of["scan=1"], pev_of["scan=4"], pev_of["scan=5"], pev_of["scan=6"], pev_of["scan=7"]}).size(), 5)
  TEST_EQUAL(countOf(text, "<Peptide "), 2)
  // the protein of both runs has a DBSequence in each database
  TEST_EQUAL(countOf(text, "<DBSequence accession=\"SHARED\""), 2)
  TEST_EQUAL(countOf(text, "<DBSequence "), 4)

  std::vector<ProteinIdentification> proteins2;
  PeptideIdentificationList peptides2;
  MzIdentMLFile().load(filename, proteins2, peptides2);
  TEST_EQUAL(peptides2.size(), 8)
  // each PSM has its own evidence (before the fix: the PSM of run 2 had the evidence of P1 from run 1)
  std::map<std::string, std::string> evidence_of; // spectrum -> accession and flanking residues of its evidence
  for (const PeptideIdentification& id : peptides2)
  {
    ABORT_IF(id.getHits().size() != 1 || id.getHits()[0].getPeptideEvidences().size() != 1)
    const PeptideEvidence& evidence = id.getHits()[0].getPeptideEvidences()[0];
    evidence_of[id.getSpectrumReference()] = evidence.getProteinAccession() + " " + evidence.getAABefore() + evidence.getAAAfter();
  }
  TEST_EQUAL(evidence_of["scan=1"], "P1 KA")
  TEST_EQUAL(evidence_of["scan=2"], "P1 KA")
  TEST_EQUAL(evidence_of["scan=7"], "P2 RG")
  TEST_EQUAL(StringUtils::hasPrefix(evidence_of["scan=3"], "SHARED "), true)
  TEST_EQUAL(StringUtils::hasPrefix(evidence_of["scan=8"], "SHARED "), true)
  // ... and each run's proteins are in the database of the run
  std::map<std::string, std::string> db_of;
  for (const ProteinIdentification& run : proteins2)
  {
    for (const ProteinHit& hit : run.getHits())
    {
      db_of[run.getSearchParameters().db + " " + hit.getAccession()] = run.getIdentifier();
    }
  }
  TEST_EQUAL(db_of.count("db_one.fasta SHARED"), 1)
  TEST_EQUAL(db_of.count("db_two.fasta SHARED"), 1)
  TEST_EQUAL(db_of.count("db_two.fasta P1"), 0)
}
END_SECTION

START_SECTION(([EXTRA] store and load the units of CV-term meta values (issue #10981)))
{
  std::vector<ProteinIdentification> proteins;
  PeptideIdentificationList peptides;
  MzIdentMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_msgf_mini.mzid"), proteins, peptides);
  ABORT_IF(peptides.empty() || peptides[0].getHits().empty())
  DataValue scan_start_time(2014.221);
  scan_start_time.setUnit(10); // UO:0000010 second
  scan_start_time.setUnitType(DataValue::UnitType::UNIT_ONTOLOGY);
  peptides[0].getHits()[0].setMetaValue("MS:1000016", scan_start_time);
  // a unit of the PSI-MS CV: its CV reference is the id of PSI-MS in the cvList of the document ("PSI-MS", not "MS")
  DataValue selected_mz(500.25);
  selected_mz.setUnit(1000040); // MS:1000040 m/z
  selected_mz.setUnitType(DataValue::UnitType::MS_ONTOLOGY);
  peptides[0].getHits()[0].setMetaValue("MS:1000744", selected_mz);

  std::string filename;
  NEW_TMP_FILE_EXT(filename, ".mzid")
  MzIdentMLFile().store(filename, proteins, peptides);
  const std::string text = readFile(filename);
  TEST_EQUAL(countOf(text, "name=\"selected ion m/z\" value=\"500.25\" unitAccession=\"MS:1000040\" unitCvRef=\"PSI-MS\" unitName=\"m/z\"/>"), 1)
  TEST_EQUAL(countOf(text, "unitCvRef=\"MS\""), 0)
  std::string version;
  std::stringstream schema_errors;
  TEST_EQUAL(MzIdentMLFile().isValid(filename, schema_errors, version), true)
  TEST_EQUAL(schema_errors.str(), "")
  const std::regex with_unit(R"(name="scan start time" value="2014\.221[0-9]*" unitAccession="UO:0000010" unitCvRef="UO" unitName="second"/>)");
  TEST_EQUAL(std::distance(std::sregex_iterator(text.begin(), text.end(), with_unit), std::sregex_iterator()), 1)

  std::vector<ProteinIdentification> proteins2;
  PeptideIdentificationList peptides2;
  MzIdentMLFile().load(filename, proteins2, peptides2);
  // the PSM with the values (the order of the runs in the file follows their generated ids)
  const PeptideHit* hit = nullptr;
  for (const PeptideIdentification& id : peptides2)
  {
    if (!id.getHits().empty() && id.getHits()[0].metaValueExists("MS:1000016")) hit = &id.getHits()[0];
  }
  ABORT_IF(hit == nullptr)
  const DataValue& value = hit->getMetaValue("MS:1000016");
  TEST_REAL_SIMILAR(double(value), 2014.221)
  TEST_EQUAL(value.hasUnit(), true)
  TEST_EQUAL(value.getUnit(), 10)
  TEST_EQUAL(value.getUnitType() == DataValue::UnitType::UNIT_ONTOLOGY, true)
  const DataValue& mz = hit->getMetaValue("MS:1000744");
  TEST_REAL_SIMILAR(double(mz), 500.25)
  TEST_EQUAL(mz.getUnit(), 1000040)
  TEST_EQUAL(mz.getUnitType() == DataValue::UnitType::MS_ONTOLOGY, true)
}
END_SECTION

START_SECTION(([EXTRA] read PeptideEvidenceRefs: listed twice, or to the evidence of another peptide (issue #10983)))
{
  // a file with PEPTIDEK (scan=1), ELVISK (scan=2), a modified PEPTIDEK (scan=3) and SAMPLER (scan=4), each from its own
  // evidence; then the PeptideEvidenceRefs of the PSMs are changed as in faulty files
  std::vector<ProteinIdentification> proteins(1);
  proteins[0].setIdentifier("run");
  proteins[0].setSearchEngine("MS-GF+");
  ProteinIdentification::SearchParameters search_parameters;
  search_parameters.db = "db.fasta";
  proteins[0].setSearchParameters(search_parameters);
  proteins[0].setMetaValue("spectra_data", StringList{"run.mzML"});
  proteins[0].insertHit(ProteinHit(0.0, 1, "P1", ""));
  proteins[0].insertHit(ProteinHit(0.0, 1, "P2", ""));
  PeptideIdentificationList peptides;
  auto add_psm = [&](const std::string& spectrum, const std::string& sequence, const PeptideEvidence& evidence)
  {
    PeptideIdentification id;
    id.setIdentifier("run");
    id.setScoreType("q-value");
    id.setHigherScoreBetter(false);
    id.setSpectrumReference(spectrum);
    id.setRT(100.0);
    id.setMZ(500.0);
    PeptideHit hit(0.01, 0, 2, AASequence::fromString(sequence));
    hit.addPeptideEvidence(evidence);
    id.insertHit(hit);
    peptides.push_back(id);
  };
  add_psm("scan=1", "PEPTIDEK", PeptideEvidence("P1", 10, 17, 'K', 'A'));
  add_psm("scan=2", "ELVISK", PeptideEvidence("P2", 20, 25, 'K', 'G'));
  add_psm("scan=3", "PEPTIDEK(Label:13C(6)15N(2))", PeptideEvidence("P2", 30, 37, 'R', 'L'));
  add_psm("scan=4", "SAMPLER", PeptideEvidence("P2", 50, 56, 'K', 'A'));

  std::string filename;
  NEW_TMP_FILE_EXT(filename, ".mzid")
  MzIdentMLFile().store(filename, proteins, peptides);
  std::string text = readFile(filename);
  TEST_EQUAL(countOf(text, "<Peptide "), 4)
  std::map<std::string, std::string> pev_of; // spectrum -> id of the PeptideEvidence it references
  const std::regex pev_reference(R"re(spectrumID="([^"]*)"[^>]*>\s*<SpectrumIdentificationItem[^>]*>\s*<PeptideEvidenceRef peptideEvidence_ref="([^"]*)")re");
  for (std::sregex_iterator m(text.begin(), text.end(), pev_reference), end; m != end; ++m)
  {
    pev_of[(*m)[1]] = (*m)[2];
  }
  ABORT_IF(pev_of.size() != 4)
  auto set_references = [&](const std::string& spectrum, const std::vector<std::string>& pev_ids)
  {
    const Size start = text.find("<PeptideEvidenceRef ", text.find("spectrumID=\"" + spectrum + "\""));
    const Size end = text.find("/>", start) + 2;
    std::string references;
    for (const std::string& pev_id : pev_ids) references += "<PeptideEvidenceRef peptideEvidence_ref=\"" + pev_id + "\"/>";
    text.replace(start, end - start, references);
  };
  // its own evidence twice, the evidence of ELVISK, and a PeptideEvidence that is not in the file
  set_references("scan=1", {pev_of["scan=1"], pev_of["scan=1"], pev_of["scan=2"], "PEV_missing"});
  // only the evidence of PEPTIDEK
  set_references("scan=2", {pev_of["scan=1"]});
  // the evidence of another Peptide element with the same sequence
  set_references("scan=3", {pev_of["scan=1"]});
  // only a PeptideEvidence that is not in the file
  set_references("scan=4", {"PEV_missing"});
  NEW_TMP_FILE_EXT(filename, ".mzid")
  std::ofstream(filename) << text;

  std::vector<ProteinIdentification> proteins2;
  PeptideIdentificationList peptides2;
  MzIdentMLFile().load(filename, proteins2, peptides2);
  std::map<std::string, std::string> evidences_of; // spectrum -> accessions and flanking residues of its evidences
  for (const PeptideIdentification& id : peptides2)
  {
    ABORT_IF(id.getHits().size() != 1)
    for (const PeptideEvidence& evidence : id.getHits()[0].getPeptideEvidences())
    {
      evidences_of[id.getSpectrumReference()] += evidence.getProteinAccession() + " " + evidence.getAABefore() + evidence.getAAAfter() + ";";
    }
  }
  TEST_EQUAL(evidences_of["scan=1"], "P1 KA;") // once, without the protein of ELVISK and without an empty evidence
  TEST_EQUAL(evidences_of["scan=2"], "P2 KG;") // no usable reference: the evidence of its peptide
  TEST_EQUAL(evidences_of["scan=3"], "P1 KA;")
  TEST_EQUAL(evidences_of["scan=4"], "P2 KA;")
}
END_SECTION

END_TEST
