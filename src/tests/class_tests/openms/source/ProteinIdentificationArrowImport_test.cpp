// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/ProteinIdentificationArrowImport.h>
///////////////////////////

#include <OpenMS/config.h>

#ifdef WITH_PARQUET

#include <OpenMS/FORMAT/ProteinIdentificationArrowExport.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>

using namespace OpenMS;
using namespace std;

START_TEST(ProteinIdentificationArrowImport, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(importSearchParamsFromArrow - round trip)
{
  // Create test data
  vector<ProteinIdentification> orig_ids;

  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_search_1");
  prot_id.setSearchEngine("Comet");
  prot_id.setSearchEngineVersion("2023.01");
  prot_id.setScoreType("expect");
  prot_id.setHigherScoreBetter(false);
  prot_id.setSignificanceThreshold(0.01);
  prot_id.setPrimaryMSRunPath(StringList{"sample_1.mzML", "sample_2.mzML"});

  DateTime dt;
  dt.set("2024-06-15 14:30:00");
  prot_id.setDateTime(dt);

  ProteinIdentification::SearchParameters sp;
  sp.db = "uniprot_human";
  sp.db_version = "2024_01";
  sp.taxonomy = "Homo sapiens";
  sp.charges = "2:4";
  sp.mass_type = ProteinIdentification::PeakMassType::MONOISOTOPIC;
  sp.precursor_mass_tolerance = 10.0;
  sp.precursor_mass_tolerance_ppm = true;
  sp.fragment_mass_tolerance = 0.02;
  sp.fragment_mass_tolerance_ppm = false;
  sp.digestion_enzyme = *ProteaseDB::getInstance()->getEnzyme("Trypsin");
  sp.enzyme_term_specificity = EnzymaticDigestion::SPEC_FULL;
  sp.missed_cleavages = 2;
  sp.fixed_modifications = {"Carbamidomethyl (C)"};
  sp.variable_modifications = {"Oxidation (M)", "Acetyl (Protein N-term)"};
  sp.setMetaValue("sp_custom_param", "my_value");
  sp.setMetaValue("sp_custom_int", 123);
  prot_id.setSearchParameters(sp);

  // Set inference engine
  prot_id.setInferenceEngine("Percolator");
  prot_id.setInferenceEngineVersion("3.06");

  // Set a ProteinIdentification-level metavalue (should NOT end up on SearchParameters)
  prot_id.setMetaValue("prot_id_custom", "prot_level_value");

  orig_ids.push_back(prot_id);

  // Export to Arrow
  auto table = ProteinIdentificationArrowExport::exportSearchParamsToArrow(orig_ids);
  TEST_NOT_EQUAL(table, nullptr)

  // Import back
  vector<ProteinIdentification> imported_ids;
  TEST_EQUAL(ProteinIdentificationArrowImport::importSearchParamsFromArrow(table, imported_ids), true)
  TEST_EQUAL(imported_ids.size(), 1)

  const auto& imp = imported_ids[0];
  TEST_EQUAL(imp.getIdentifier(), "run_search_1")
  TEST_EQUAL(imp.getSearchEngine(), "Comet")
  TEST_EQUAL(imp.getSearchEngineVersion(), "2023.01")
  TEST_EQUAL(imp.getScoreType(), "expect")
  TEST_EQUAL(imp.isHigherScoreBetter(), false)
  TEST_REAL_SIMILAR(imp.getSignificanceThreshold(), 0.01)

  // Date
  TEST_EQUAL(imp.getDateTime().toString("yyyy-MM-ddThh:mm:ss"), "2024-06-15T14:30:00")

  // Primary MS run paths
  StringList ms_runs;
  imp.getPrimaryMSRunPath(ms_runs);
  TEST_EQUAL(ms_runs.size(), 2)
  TEST_EQUAL(ms_runs[0], "sample_1.mzML")
  TEST_EQUAL(ms_runs[1], "sample_2.mzML")

  // SearchParameters
  const auto& imp_sp = imp.getSearchParameters();
  TEST_EQUAL(imp_sp.db, "uniprot_human")
  TEST_EQUAL(imp_sp.db_version, "2024_01")
  TEST_EQUAL(imp_sp.taxonomy, "Homo sapiens")
  TEST_EQUAL(imp_sp.charges, "2:4")
  TEST_EQUAL(imp_sp.mass_type == ProteinIdentification::PeakMassType::MONOISOTOPIC, true)
  TEST_REAL_SIMILAR(imp_sp.precursor_mass_tolerance, 10.0)
  TEST_EQUAL(imp_sp.precursor_mass_tolerance_ppm, true)
  TEST_REAL_SIMILAR(imp_sp.fragment_mass_tolerance, 0.02)
  TEST_EQUAL(imp_sp.fragment_mass_tolerance_ppm, false)
  TEST_EQUAL(imp_sp.digestion_enzyme.getName(), "Trypsin")
  TEST_EQUAL(imp_sp.enzyme_term_specificity == EnzymaticDigestion::SPEC_FULL, true)
  TEST_EQUAL(imp_sp.missed_cleavages, 2)
  TEST_EQUAL(imp_sp.fixed_modifications.size(), 1)
  TEST_EQUAL(imp_sp.fixed_modifications[0], "Carbamidomethyl (C)")
  TEST_EQUAL(imp_sp.variable_modifications.size(), 2)
  TEST_EQUAL(imp_sp.variable_modifications[0], "Oxidation (M)")
  TEST_EQUAL(imp_sp.variable_modifications[1], "Acetyl (Protein N-term)")

  // Inference engine
  TEST_EQUAL(imp.getInferenceEngine(), "Percolator")
  TEST_EQUAL(imp.getInferenceEngineVersion(), "3.06")

  // SearchParameters metavalues should be on SearchParameters, not on ProteinIdentification
  TEST_EQUAL(String(imp_sp.getMetaValue("sp_custom_param")), "my_value")
  TEST_EQUAL(int(imp_sp.getMetaValue("sp_custom_int")), 123)

  // ProteinIdentification metavalue should be on ProteinIdentification
  TEST_EQUAL(String(imp.getMetaValue("prot_id_custom")), "prot_level_value")

  // Cross-check: sp_custom_param should NOT be on the ProteinIdentification object
  TEST_EQUAL(imp.metaValueExists("sp_custom_param"), false)
}
END_SECTION

START_SECTION(importProteinsFromArrow - round trip)
{
  // Create test data
  vector<ProteinIdentification> orig_ids;

  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_1");
  prot_id.setSearchEngine("Mascot");
  prot_id.setSearchEngineVersion("2.7");
  prot_id.setScoreType("Mascot");
  prot_id.setHigherScoreBetter(true);
  prot_id.setSignificanceThreshold(0.05);
  prot_id.setPrimaryMSRunPath(StringList{"sample_1.mzML"});

  DateTime dt;
  dt.set("2024-01-15 10:30:00");
  prot_id.setDateTime(dt);

  // Hit 1: target protein with coverage and description
  ProteinHit hit1;
  hit1.setAccession("sp|P12345|PROT1_HUMAN");
  hit1.setScore(150.5);
  hit1.setRank(1);
  hit1.setCoverage(45.2);
  hit1.setSequence("MKWVTFISLLLLFSSAYS");
  hit1.setDescription("Serum albumin");
  hit1.setMetaValue("target_decoy", "target");
  hit1.setMetaValue("custom_score", 0.99);

  // Hit 2: decoy protein with unknown coverage
  ProteinHit hit2;
  hit2.setAccession("DECOY_sp|Q67890|PROT2_HUMAN");
  hit2.setScore(25.3);
  hit2.setRank(2);
  hit2.setMetaValue("target_decoy", "decoy");

  prot_id.insertHit(hit1);
  prot_id.insertHit(hit2);
  orig_ids.push_back(prot_id);

  // Export to Arrow
  auto table = ProteinIdentificationArrowExport::exportProteinsToArrow(orig_ids);
  TEST_NOT_EQUAL(table, nullptr)

  // Import back (starting from scratch - no pre-existing shells)
  vector<ProteinIdentification> imported_ids;
  TEST_EQUAL(ProteinIdentificationArrowImport::importProteinsFromArrow(table, imported_ids), true)
  TEST_EQUAL(imported_ids.size(), 1)

  const auto& imp = imported_ids[0];
  TEST_EQUAL(imp.getIdentifier(), "run_1")
  TEST_EQUAL(imp.getSearchEngine(), "Mascot")
  TEST_EQUAL(imp.getHits().size(), 2)

  // Hit 1
  const auto& imp_hit1 = imp.getHits()[0];
  TEST_EQUAL(imp_hit1.getAccession(), "sp|P12345|PROT1_HUMAN")
  TEST_REAL_SIMILAR(imp_hit1.getScore(), 150.5)
  TEST_EQUAL(imp_hit1.getRank(), 1)
  TEST_REAL_SIMILAR(imp_hit1.getCoverage(), 45.2)
  TEST_EQUAL(imp_hit1.getSequence(), "MKWVTFISLLLLFSSAYS")
  TEST_EQUAL(imp_hit1.getDescription(), "Serum albumin")
  TEST_EQUAL(String(imp_hit1.getMetaValue("target_decoy")), "target")
  TEST_REAL_SIMILAR(double(imp_hit1.getMetaValue("custom_score")), 0.99)

  // Hit 2
  const auto& imp_hit2 = imp.getHits()[1];
  TEST_EQUAL(imp_hit2.getAccession(), "DECOY_sp|Q67890|PROT2_HUMAN")
  TEST_REAL_SIMILAR(imp_hit2.getScore(), 25.3)
  TEST_EQUAL(imp_hit2.getRank(), 2)
  TEST_REAL_SIMILAR(imp_hit2.getCoverage(), ProteinHit::COVERAGE_UNKNOWN) // was null -> default
  TEST_EQUAL(String(imp_hit2.getMetaValue("target_decoy")), "decoy")
}
END_SECTION

START_SECTION(importProteinGroupsFromArrow - round trip)
{
  // Create test data
  vector<ProteinIdentification> orig_ids;

  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_groups_1");

  // Protein group
  ProteinIdentification::ProteinGroup pg;
  pg.probability = 0.95;
  pg.accessions = {"sp|P12345|PROT1_HUMAN", "sp|P12346|PROT2_HUMAN"};
  prot_id.insertProteinGroup(pg);

  // Indistinguishable group
  ProteinIdentification::ProteinGroup ig;
  ig.probability = 0.0;
  ig.accessions = {"sp|Q11111|PROTA_HUMAN", "sp|Q22222|PROTB_HUMAN"};
  prot_id.insertIndistinguishableProteins(ig);

  orig_ids.push_back(prot_id);

  // Export to Arrow
  auto table = ProteinIdentificationArrowExport::exportProteinGroupsToArrow(orig_ids);
  TEST_NOT_EQUAL(table, nullptr)

  // Import back (create shell first to test run_id matching)
  vector<ProteinIdentification> imported_ids;
  ProteinIdentification shell;
  shell.setIdentifier("run_groups_1");
  imported_ids.push_back(shell);

  TEST_EQUAL(ProteinIdentificationArrowImport::importProteinGroupsFromArrow(table, imported_ids), true)
  TEST_EQUAL(imported_ids.size(), 1)

  // Verify protein groups
  const auto& imp = imported_ids[0];
  TEST_EQUAL(imp.getProteinGroups().size(), 1)
  TEST_REAL_SIMILAR(imp.getProteinGroups()[0].probability, 0.95)
  TEST_EQUAL(imp.getProteinGroups()[0].accessions.size(), 2)
  TEST_EQUAL(imp.getProteinGroups()[0].accessions[0], "sp|P12345|PROT1_HUMAN")
  TEST_EQUAL(imp.getProteinGroups()[0].accessions[1], "sp|P12346|PROT2_HUMAN")

  // Verify indistinguishable proteins
  TEST_EQUAL(imp.getIndistinguishableProteins().size(), 1)
  TEST_REAL_SIMILAR(imp.getIndistinguishableProteins()[0].probability, 0.0)
  TEST_EQUAL(imp.getIndistinguishableProteins()[0].accessions.size(), 2)
  TEST_EQUAL(imp.getIndistinguishableProteins()[0].accessions[0], "sp|Q11111|PROTA_HUMAN")
  TEST_EQUAL(imp.getIndistinguishableProteins()[0].accessions[1], "sp|Q22222|PROTB_HUMAN")
}
END_SECTION

START_SECTION(importSearchParamsFromParquet - file round trip)
{
  // Create test data
  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_parquet_sp");
  prot_id.setSearchEngine("Comet");
  prot_id.setScoreType("expect");
  prot_id.setHigherScoreBetter(false);

  ProteinIdentification::SearchParameters sp;
  sp.db = "uniprot_human";
  sp.precursor_mass_tolerance = 10.0;
  sp.precursor_mass_tolerance_ppm = true;
  sp.fragment_mass_tolerance = 0.02;
  sp.fragment_mass_tolerance_ppm = false;
  sp.digestion_enzyme = *ProteaseDB::getInstance()->getEnzyme("Trypsin");
  sp.missed_cleavages = 2;
  sp.fixed_modifications = {"Carbamidomethyl (C)"};
  sp.variable_modifications = {"Oxidation (M)"};
  prot_id.setSearchParameters(sp);

  vector<ProteinIdentification> orig_ids = {prot_id};

  String tmp_file;
  NEW_TMP_FILE(tmp_file)
  tmp_file += ".parquet";

  TEST_EQUAL(ProteinIdentificationArrowExport::exportSearchParamsToParquet(orig_ids, tmp_file), true)

  vector<ProteinIdentification> imported_ids;
  TEST_EQUAL(ProteinIdentificationArrowImport::importSearchParamsFromParquet(tmp_file, imported_ids), true)
  TEST_EQUAL(imported_ids.size(), 1)
  TEST_EQUAL(imported_ids[0].getIdentifier(), "run_parquet_sp")
  TEST_EQUAL(imported_ids[0].getSearchParameters().db, "uniprot_human")
  TEST_EQUAL(imported_ids[0].getSearchParameters().digestion_enzyme.getName(), "Trypsin")
}
END_SECTION

START_SECTION(importProteinsFromParquet - file round trip)
{
  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_parquet_prot");
  prot_id.setSearchEngine("MS-GF+");
  prot_id.setScoreType("q-value");
  prot_id.setHigherScoreBetter(false);

  ProteinHit hit1;
  hit1.setAccession("P12345");
  hit1.setScore(0.001);
  prot_id.setHits({hit1});

  vector<ProteinIdentification> orig_ids = {prot_id};

  String tmp_file;
  NEW_TMP_FILE(tmp_file)
  tmp_file += ".parquet";

  TEST_EQUAL(ProteinIdentificationArrowExport::exportProteinsToParquet(orig_ids, tmp_file), true)

  vector<ProteinIdentification> imported_ids;
  TEST_EQUAL(ProteinIdentificationArrowImport::importProteinsFromParquet(tmp_file, imported_ids), true)
  TEST_EQUAL(imported_ids.size(), 1)
  TEST_EQUAL(imported_ids[0].getHits().size(), 1)
  TEST_EQUAL(imported_ids[0].getHits()[0].getAccession(), "P12345")
  TEST_REAL_SIMILAR(imported_ids[0].getHits()[0].getScore(), 0.001)
}
END_SECTION

START_SECTION(importProteinGroupsFromParquet - file round trip)
{
  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_parquet_pg");

  ProteinIdentification::ProteinGroup pg;
  pg.probability = 0.95;
  pg.accessions = {"P12345", "P12346"};
  prot_id.insertProteinGroup(pg);

  vector<ProteinIdentification> orig_ids = {prot_id};

  String tmp_file;
  NEW_TMP_FILE(tmp_file)
  tmp_file += ".parquet";

  TEST_EQUAL(ProteinIdentificationArrowExport::exportProteinGroupsToParquet(orig_ids, tmp_file), true)

  vector<ProteinIdentification> imported_ids;
  ProteinIdentification shell;
  shell.setIdentifier("run_parquet_pg");
  imported_ids.push_back(shell);

  TEST_EQUAL(ProteinIdentificationArrowImport::importProteinGroupsFromParquet(tmp_file, imported_ids), true)
  TEST_EQUAL(imported_ids.size(), 1)
  TEST_EQUAL(imported_ids[0].getProteinGroups().size(), 1)
  TEST_REAL_SIMILAR(imported_ids[0].getProteinGroups()[0].probability, 0.95)
}
END_SECTION

START_SECTION(importFromParquet - full combined round trip)
{
  // Create rich test data with 2 runs
  vector<ProteinIdentification> orig_ids;

  // Run 1
  {
    ProteinIdentification prot_id;
    prot_id.setIdentifier("run_combined_1");
    prot_id.setSearchEngine("Comet");
    prot_id.setSearchEngineVersion("2023.01");
    prot_id.setScoreType("expect");
    prot_id.setHigherScoreBetter(false);
    prot_id.setSignificanceThreshold(0.01);
    prot_id.setPrimaryMSRunPath(StringList{"sample_1.mzML"});

    DateTime dt;
    dt.set("2024-06-15 14:30:00");
    prot_id.setDateTime(dt);

    ProteinIdentification::SearchParameters sp;
    sp.db = "uniprot_human";
    sp.db_version = "2024_01";
    sp.mass_type = ProteinIdentification::PeakMassType::MONOISOTOPIC;
    sp.precursor_mass_tolerance = 10.0;
    sp.precursor_mass_tolerance_ppm = true;
    sp.fragment_mass_tolerance = 0.02;
    sp.fragment_mass_tolerance_ppm = false;
    sp.digestion_enzyme = *ProteaseDB::getInstance()->getEnzyme("Trypsin");
    sp.enzyme_term_specificity = EnzymaticDigestion::SPEC_FULL;
    sp.missed_cleavages = 2;
    sp.fixed_modifications = {"Carbamidomethyl (C)"};
    sp.variable_modifications = {"Oxidation (M)"};
    prot_id.setSearchParameters(sp);

    ProteinHit hit1;
    hit1.setAccession("sp|P12345|PROT1_HUMAN");
    hit1.setScore(0.001);
    hit1.setRank(1);
    hit1.setCoverage(45.2);
    hit1.setDescription("Serum albumin");
    hit1.setMetaValue("target_decoy", "target");
    prot_id.insertHit(hit1);

    ProteinHit hit2;
    hit2.setAccession("sp|Q67890|PROT2_HUMAN");
    hit2.setScore(0.005);
    hit2.setRank(2);
    hit2.setMetaValue("target_decoy", "target");
    prot_id.insertHit(hit2);

    ProteinIdentification::ProteinGroup pg;
    pg.probability = 0.99;
    pg.accessions = {"sp|P12345|PROT1_HUMAN", "sp|Q67890|PROT2_HUMAN"};
    prot_id.insertProteinGroup(pg);

    ProteinIdentification::ProteinGroup ig;
    ig.probability = 0.0;
    ig.accessions = {"sp|P12345|PROT1_HUMAN"};
    prot_id.insertIndistinguishableProteins(ig);

    orig_ids.push_back(prot_id);
  }

  // Run 2
  {
    ProteinIdentification prot_id;
    prot_id.setIdentifier("run_combined_2");
    prot_id.setSearchEngine("MS-GF+");
    prot_id.setScoreType("SpecEValue");
    prot_id.setHigherScoreBetter(false);

    ProteinIdentification::SearchParameters sp;
    sp.db = "uniprot_mouse";
    sp.mass_type = ProteinIdentification::PeakMassType::MONOISOTOPIC;
    sp.precursor_mass_tolerance = 20.0;
    sp.precursor_mass_tolerance_ppm = true;
    sp.fragment_mass_tolerance = 0.5;
    sp.fragment_mass_tolerance_ppm = false;
    sp.digestion_enzyme = *ProteaseDB::getInstance()->getEnzyme("Trypsin/P");
    sp.enzyme_term_specificity = EnzymaticDigestion::SPEC_SEMI;
    sp.missed_cleavages = 1;
    prot_id.setSearchParameters(sp);

    ProteinHit hit;
    hit.setAccession("MOUSE_PROT1");
    hit.setScore(1e-10);
    prot_id.insertHit(hit);

    orig_ids.push_back(prot_id);
  }

  // Export all 3 tables
  String proteins_file, groups_file, params_file;
  NEW_TMP_FILE(proteins_file)
  proteins_file += ".parquet";
  NEW_TMP_FILE(groups_file)
  groups_file += ".parquet";
  NEW_TMP_FILE(params_file)
  params_file += ".parquet";

  TEST_EQUAL(ProteinIdentificationArrowExport::exportProteinsToParquet(orig_ids, proteins_file), true)
  TEST_EQUAL(ProteinIdentificationArrowExport::exportProteinGroupsToParquet(orig_ids, groups_file), true)
  TEST_EQUAL(ProteinIdentificationArrowExport::exportSearchParamsToParquet(orig_ids, params_file), true)

  // Import combined
  vector<ProteinIdentification> imported_ids;
  TEST_EQUAL(ProteinIdentificationArrowImport::importFromParquet(
    proteins_file, groups_file, params_file, imported_ids), true)

  // Should have 2 runs
  TEST_EQUAL(imported_ids.size(), 2)

  // Find runs by identifier (order may differ)
  const ProteinIdentification* imp_run1 = nullptr;
  const ProteinIdentification* imp_run2 = nullptr;
  for (const auto& pid : imported_ids)
  {
    if (pid.getIdentifier() == "run_combined_1") imp_run1 = &pid;
    else if (pid.getIdentifier() == "run_combined_2") imp_run2 = &pid;
  }
  TEST_NOT_EQUAL(imp_run1, nullptr)
  TEST_NOT_EQUAL(imp_run2, nullptr)

  // Verify run 1
  TEST_EQUAL(imp_run1->getSearchEngine(), "Comet")
  TEST_EQUAL(imp_run1->getSearchEngineVersion(), "2023.01")
  TEST_EQUAL(imp_run1->getScoreType(), "expect")
  TEST_EQUAL(imp_run1->isHigherScoreBetter(), false)
  TEST_EQUAL(imp_run1->getHits().size(), 2)
  TEST_EQUAL(imp_run1->getProteinGroups().size(), 1)
  TEST_EQUAL(imp_run1->getIndistinguishableProteins().size(), 1)

  // Verify run 1 search params
  const auto& sp1 = imp_run1->getSearchParameters();
  TEST_EQUAL(sp1.db, "uniprot_human")
  TEST_REAL_SIMILAR(sp1.precursor_mass_tolerance, 10.0)
  TEST_EQUAL(sp1.digestion_enzyme.getName(), "Trypsin")

  // Verify run 1 hits
  bool found_p12345 = false;
  for (const auto& hit : imp_run1->getHits())
  {
    if (hit.getAccession() == "sp|P12345|PROT1_HUMAN")
    {
      found_p12345 = true;
      TEST_REAL_SIMILAR(hit.getScore(), 0.001)
      TEST_REAL_SIMILAR(hit.getCoverage(), 45.2)
    }
  }
  TEST_EQUAL(found_p12345, true)

  // Verify run 1 protein groups
  TEST_REAL_SIMILAR(imp_run1->getProteinGroups()[0].probability, 0.99)
  TEST_EQUAL(imp_run1->getProteinGroups()[0].accessions.size(), 2)

  // Verify run 2
  TEST_EQUAL(imp_run2->getSearchEngine(), "MS-GF+")
  TEST_EQUAL(imp_run2->getHits().size(), 1)
  TEST_EQUAL(imp_run2->getHits()[0].getAccession(), "MOUSE_PROT1")

  const auto& sp2 = imp_run2->getSearchParameters();
  TEST_EQUAL(sp2.db, "uniprot_mouse")
  TEST_REAL_SIMILAR(sp2.precursor_mass_tolerance, 20.0)
  TEST_EQUAL(sp2.digestion_enzyme.getName(), "Trypsin/P")
}
END_SECTION

START_SECTION(empty input round trip)
{
  // Test with empty input
  vector<ProteinIdentification> empty_ids;

  auto sp_table = ProteinIdentificationArrowExport::exportSearchParamsToArrow(empty_ids);
  auto prot_table = ProteinIdentificationArrowExport::exportProteinsToArrow(empty_ids);
  auto pg_table = ProteinIdentificationArrowExport::exportProteinGroupsToArrow(empty_ids);

  vector<ProteinIdentification> imported_ids;
  TEST_EQUAL(ProteinIdentificationArrowImport::importSearchParamsFromArrow(sp_table, imported_ids), true)
  TEST_EQUAL(imported_ids.size(), 0)

  TEST_EQUAL(ProteinIdentificationArrowImport::importProteinsFromArrow(prot_table, imported_ids), true)
  TEST_EQUAL(imported_ids.size(), 0)

  TEST_EQUAL(ProteinIdentificationArrowImport::importProteinGroupsFromArrow(pg_table, imported_ids), true)
  TEST_EQUAL(imported_ids.size(), 0)
}
END_SECTION

START_SECTION(nullable fields round trip)
{
  // Test that nullable fields survive round-trip correctly
  vector<ProteinIdentification> orig_ids;

  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_nullable");
  prot_id.setSearchEngine("TestEngine");
  prot_id.setScoreType("TestScore");
  prot_id.setHigherScoreBetter(true);

  // Hit with minimal data (many nullable fields will be null)
  ProteinHit hit;
  hit.setAccession("MINIMAL_PROT");
  hit.setScore(1.0);
  // No coverage set (COVERAGE_UNKNOWN), no sequence, no description, no target_decoy
  prot_id.insertHit(hit);
  orig_ids.push_back(prot_id);

  auto table = ProteinIdentificationArrowExport::exportProteinsToArrow(orig_ids);
  TEST_NOT_EQUAL(table, nullptr)

  vector<ProteinIdentification> imported_ids;
  TEST_EQUAL(ProteinIdentificationArrowImport::importProteinsFromArrow(table, imported_ids), true)
  TEST_EQUAL(imported_ids.size(), 1)

  const auto& imp_hit = imported_ids[0].getHits()[0];
  TEST_EQUAL(imp_hit.getAccession(), "MINIMAL_PROT")
  TEST_REAL_SIMILAR(imp_hit.getScore(), 1.0)
  // Coverage should be default (COVERAGE_UNKNOWN since null was read)
  TEST_REAL_SIMILAR(imp_hit.getCoverage(), ProteinHit::COVERAGE_UNKNOWN)
  // Sequence should be empty
  TEST_EQUAL(imp_hit.getSequence(), "")
}
END_SECTION

START_SECTION(metavalue type preservation round trip)
{
  // Test that int, float, and string metavalues preserve their types
  vector<ProteinIdentification> orig_ids;

  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_mv_types");
  prot_id.setSearchEngine("TestEngine");
  prot_id.setScoreType("TestScore");
  prot_id.setHigherScoreBetter(true);

  ProteinHit hit;
  hit.setAccession("PROT_MV");
  hit.setScore(1.0);
  hit.setMetaValue("my_int", 42);
  hit.setMetaValue("my_float", 3.14);
  hit.setMetaValue("my_string", "hello_world");
  prot_id.insertHit(hit);
  orig_ids.push_back(prot_id);

  auto table = ProteinIdentificationArrowExport::exportProteinsToArrow(orig_ids);
  TEST_NOT_EQUAL(table, nullptr)

  vector<ProteinIdentification> imported_ids;
  TEST_EQUAL(ProteinIdentificationArrowImport::importProteinsFromArrow(table, imported_ids), true)
  TEST_EQUAL(imported_ids.size(), 1)

  const auto& imp_hit = imported_ids[0].getHits()[0];
  TEST_EQUAL(int(imp_hit.getMetaValue("my_int")), 42)
  TEST_REAL_SIMILAR(double(imp_hit.getMetaValue("my_float")), 3.14)
  TEST_EQUAL(String(imp_hit.getMetaValue("my_string")), "hello_world")
}
END_SECTION

START_SECTION(protein groups with data arrays round trip)
{
  // Test that float_data, string_data, and integer_data survive round-trip
  vector<ProteinIdentification> orig_ids;

  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_data_arrays");

  ProteinIdentification::ProteinGroup pg;
  pg.probability = 0.85;
  pg.accessions = {"PROT_A", "PROT_B"};

  // Add float data array
  ProteinIdentification::ProteinGroup::FloatDataArray fda;
  fda.setName("q_values");
  fda.push_back(0.01f);
  fda.push_back(0.02f);
  pg.getFloatDataArrays().push_back(fda);

  // Add string data array
  ProteinIdentification::ProteinGroup::StringDataArray sda;
  sda.setName("annotations");
  sda.push_back("group_leader");
  sda.push_back("member");
  pg.getStringDataArrays().push_back(sda);

  // Add integer data array
  ProteinIdentification::ProteinGroup::IntegerDataArray ida;
  ida.setName("counts");
  ida.push_back(10);
  ida.push_back(20);
  pg.getIntegerDataArrays().push_back(ida);

  prot_id.insertProteinGroup(pg);
  orig_ids.push_back(prot_id);

  auto table = ProteinIdentificationArrowExport::exportProteinGroupsToArrow(orig_ids);
  TEST_NOT_EQUAL(table, nullptr)

  vector<ProteinIdentification> imported_ids;
  ProteinIdentification shell;
  shell.setIdentifier("run_data_arrays");
  imported_ids.push_back(shell);

  TEST_EQUAL(ProteinIdentificationArrowImport::importProteinGroupsFromArrow(table, imported_ids), true)
  TEST_EQUAL(imported_ids[0].getProteinGroups().size(), 1)

  const auto& imp_pg = imported_ids[0].getProteinGroups()[0];
  TEST_REAL_SIMILAR(imp_pg.probability, 0.85)

  // Float data arrays
  TEST_EQUAL(imp_pg.getFloatDataArrays().size(), 1)
  TEST_EQUAL(imp_pg.getFloatDataArrays()[0].getName(), "q_values")
  TEST_EQUAL(imp_pg.getFloatDataArrays()[0].size(), 2)
  TEST_REAL_SIMILAR(imp_pg.getFloatDataArrays()[0][0], 0.01)
  TEST_REAL_SIMILAR(imp_pg.getFloatDataArrays()[0][1], 0.02)

  // String data arrays
  TEST_EQUAL(imp_pg.getStringDataArrays().size(), 1)
  TEST_EQUAL(imp_pg.getStringDataArrays()[0].getName(), "annotations")
  TEST_EQUAL(imp_pg.getStringDataArrays()[0].size(), 2)
  TEST_EQUAL(imp_pg.getStringDataArrays()[0][0], "group_leader")
  TEST_EQUAL(imp_pg.getStringDataArrays()[0][1], "member")

  // Integer data arrays
  TEST_EQUAL(imp_pg.getIntegerDataArrays().size(), 1)
  TEST_EQUAL(imp_pg.getIntegerDataArrays()[0].getName(), "counts")
  TEST_EQUAL(imp_pg.getIntegerDataArrays()[0].size(), 2)
  TEST_EQUAL(imp_pg.getIntegerDataArrays()[0][0], 10)
  TEST_EQUAL(imp_pg.getIntegerDataArrays()[0][1], 20)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST

#else // WITH_PARQUET

START_TEST(ProteinIdentificationArrowImport, "$Id$")
END_TEST

#endif // WITH_PARQUET
