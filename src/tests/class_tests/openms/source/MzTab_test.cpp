// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Timo Sachsenberg$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/MzTab.h>
///////////////////////////

#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/KERNEL/FeatureHandle.h>

START_TEST(MzTab, "$Id$")

using namespace OpenMS;
using namespace std;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MzTab* ptr = nullptr;
MzTab* null_ptr = nullptr;
START_SECTION(MzTab())
{
	ptr = new MzTab();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~MzTab())
{
	delete ptr;
}
END_SECTION

START_SECTION(std::vector<std::string> getPSMOptionalColumnNames() const)
{
  MzTab mztab;
  MzTabPSMSectionRow row;
  MzTabPSMSectionRows rows;
  MzTabString s;
  MzTabOptionalColumnEntry e;

  // row 1 //////////////////////
  row.sequence.fromCellString("NDYKAPPQPAPGK");
  row.PSM_ID.fromCellString("38");
  row.accession.fromCellString("IPI:B1");
  row.unique.fromCellString("1");
  row.database.fromCellString("null");
  row.database_version.fromCellString("null");
  row.search_engine.fromCellString("[, , Percolator, ]");
  row.search_engine_score[0].fromCellString("51.9678841193106");

  e.first = "Percolator_score";
  s.fromCellString("0.359083");
  e.second = s;
  row.opt_.push_back(e);

  e.first = "Percolator_qvalue";
  s.fromCellString("0.00649874");
  e.second = s;
  row.opt_.push_back(e);

  e.first = "Percolator_PEP";
  s.fromCellString("0.0420992");
  e.second = s;
  row.opt_.push_back(e);

  e.first = "search_engine_sequence";
  s.fromCellString("NDYKAPPQPAPGK");
  e.second = s;
  row.opt_.push_back(e);

  rows.push_back(row);

  // row 2 //////////////////////
  row.sequence.fromCellString("IRRS(Phospho)SFSSK");
  row.PSM_ID.fromCellString("39");
  row.accession.fromCellString("IPI:IPI00009899.4");
  row.unique.fromCellString("0");
  row.database.fromCellString("null");
  row.database_version.fromCellString("null");
  row.search_engine.fromCellString("[, , Percolator, ]");
  row.search_engine_score[0].fromCellString("9.55915773892318");

  e.first = "Percolator_score";
  s.fromCellString("0.157068");
  e.second = s;
  row.opt_.push_back(e);

  e.first = "Percolator_qvalue";
  s.fromCellString("0.00774619");
  e.second = s;
  row.opt_.push_back(e);

  e.first = "Percolator_PEP";
  s.fromCellString("0.0779777");
  e.second = s;
  row.opt_.push_back(e);

  e.first = "search_engine_sequence";
  s.fromCellString("IRRSSFS(Phospho)SK");
  e.second = s;
  row.opt_.push_back(e);

  e.first = "AScore_1";
  s.fromCellString("3.64384830671351");
  e.second = s;
  row.opt_.push_back(e);
  rows.push_back(row);

  mztab.setPSMSectionRows(rows);

  // Tests ///////////////////////////////
  vector<std::string> optional_columns = mztab.getPSMOptionalColumnNames();

  TEST_EQUAL(mztab.getPSMSectionRows().size(),2)
  TEST_EQUAL(optional_columns.size(),5)
}
END_SECTION

START_SECTION(static void addMetaInfoToOptionalColumns(const std::set<std::string>& keys, std::vector<MzTabOptionalColumnEntry>& opt, const std::string& id, const MetaInfoInterface& meta))
{
  // keys will have spaces replaced with underscore
  std::set<std::string> keys{ "FWHM", "with space", "ppm_errors" };
  // values should remain as they are
  MetaInfoInterface meta;
  meta.setMetaValue("FWHM", 34.5);
  DataValue dv(DoubleList{ 0.5, 1.4, -2.0, 0.1 });
  meta.setMetaValue("ppm_errors", dv);
  std::vector<MzTabOptionalColumnEntry> opt;
  MzTab::addMetaInfoToOptionalColumns(keys, opt, "global", meta);
  TEST_EQUAL(opt.size(), 3)
  TEST_EQUAL(opt[0].first, "opt_global_FWHM");
  TEST_EQUAL(opt[1].first, "opt_global_ppm_errors");
  TEST_EQUAL(opt[2].first, "opt_global_with_space");
  TEST_EQUAL(opt[0].second.toCellString(), "34.5");
  TEST_EQUAL(opt[1].second.toCellString(), "[0.5, 1.4, -2.0, 0.1]");
  TEST_EQUAL(opt[2].second.toCellString(), "null");
}
END_SECTION

START_SECTION([EXTRA] exportIdentificationsToMzTab terminates with export_all_psms on empty-hit PeptideIdentification)
{
  // regression: in IDMzTabStream::nextPSMRow the advance condition used
  // 'current_psm_idx_ == pid->getHits().size()-1', which underflows to SIZE_MAX for an
  // EMPTY hit list. With export_all_psms=true the caller then looped forever appending
  // rows (CI hang). The fix advances via 'current_psm_idx_ + 1 >= pid->getHits().size()'.
  // This test simply COMPLETING proves termination; before the fix it would hang.

  const std::string run_id = "run_0";

  // one protein run with a matching identifier and a primary MS run path
  ProteinIdentification prot;
  prot.setIdentifier(run_id);
  prot.setPrimaryMSRunPath(StringList{"file.mzML"});
  std::vector<ProteinIdentification> prot_ids{prot};

  // one spectrum-level PeptideIdentification with NO hits, matching the run identifier
  PeptideIdentification pep;
  pep.setIdentifier(run_id);
  pep.setRT(123.4);
  pep.setMZ(567.8);
  pep.setHits(std::vector<PeptideHit>{}); // empty hit list triggers the underflow pre-fix
  PeptideIdentificationList pep_ids;
  pep_ids.push_back(pep);

  MzTab mztab = MzTab::exportIdentificationsToMzTab(
      prot_ids,
      pep_ids,
      "test.idXML",
      /* first_run_inference_only */ false,
      /* export_empty_pep_ids */ false,
      /* export_all_psms */ true,
      /* title */ "regression test");

  // The key assertion is that we get here at all (finite, terminating): pre-fix the
  // caller looped forever appending rows. Post-fix the single empty-hit PeptideIdentification
  // yields at most one PSM row before the stream advances and finishes (the exact
  // empty-row count is an implementation detail; termination is what we regress on).
  TEST_TRUE(mztab.getPSMSectionRows().size() <= 1)
}
END_SECTION

START_SECTION([EXTRA] MzTabBoolean setNull / isNull polarity)
{
  // regression: setNull(true) must make the cell null, setNull(false) must make it not-null (the branches were inverted).
  MzTabBoolean b(true);
  TEST_EQUAL(b.isNull(), false)
  b.setNull(true);
  TEST_EQUAL(b.isNull(), true)
  TEST_EQUAL(b.toCellString(), "null")
  b.setNull(false);
  TEST_EQUAL(b.isNull(), false)
}
END_SECTION

START_SECTION([EXTRA] consensus-map assays use fraction-group/label grain and samples are study variables)
{
  ConsensusMap consensus_map;
  consensus_map.setExperimentType("label-free");

  // One sample owns two fraction groups. The first assay has two measured fractions.
  const vector<string> paths{"sample_f1.mzML", "sample_f2.mzML", "sample_replicate.mzML"};
  ConsensusMap::ColumnHeaders headers;
  for (Size i = 0; i < paths.size(); ++i)
  {
    ConsensusMap::ColumnHeader header;
    header.filename = paths[i];
    header.label = "label-free";
    header.setMetaValue("fraction_group", i < 2 ? 1 : 2);
    header.setMetaValue("fraction", i < 2 ? i + 1 : 1);
    header.setMetaValue("sample_name", "sample_A");
    headers[i] = header;
  }
  consensus_map.setColumnHeaders(headers);

  ConsensusFeature feature;
  feature.insert(FeatureHandle(0, Peak2D({10.0, 500.0}, 2.0), 1));
  feature.insert(FeatureHandle(1, Peak2D({11.0, 500.0}, 3.0), 2));
  feature.insert(FeatureHandle(2, Peak2D({12.0, 500.0}, 7.0), 3));
  consensus_map.push_back(feature);

  ProteinIdentification protein_id;
  protein_id.setIdentifier("run_0");
  protein_id.setPrimaryMSRunPath(paths);
  ProteinHit protein_hit;
  protein_hit.setAccession("P1");
  protein_id.setHits({protein_hit});

  ProteinIdentification::ProteinGroup group;
  group.accessions = {"P1"};
  group.getFloatDataArrays().resize(2);
  group.getFloatDataArrays()[0].setName("abundances");
  group.getFloatDataArrays()[0].push_back(30.0f);
  group.getFloatDataArrays()[1].setName("fraction_group_level_abundance");
  // Reverse the key order to verify matching through the parallel arrays, not by position.
  group.getFloatDataArrays()[1].push_back(20.0f);
  group.getFloatDataArrays()[1].push_back(10.0f);
  group.getIntegerDataArrays().resize(2);
  group.getIntegerDataArrays()[0].setName("fraction_group_level_fraction_group");
  group.getIntegerDataArrays()[0].push_back(2);
  group.getIntegerDataArrays()[0].push_back(1);
  group.getIntegerDataArrays()[1].setName("fraction_group_level_label");
  group.getIntegerDataArrays()[1].push_back(1);
  group.getIntegerDataArrays()[1].push_back(1);
  protein_id.insertIndistinguishableProteins(group);
  consensus_map.setProteinIdentifications({protein_id});

  const MzTab mz_tab = MzTab::exportConsensusMapToMzTab(
    consensus_map,
    "test.consensusXML",
    /* first_run_inference_only */ false,
    /* export_unidentified_features */ true,
    /* export_unassigned_ids */ false,
    /* export_subfeatures */ false);

  const MzTabMetaData& meta = mz_tab.getMetaData();
  TEST_EQUAL(meta.assay.size(), 2)
  if (meta.assay.size() == 2)
  {
    TEST_EQUAL(meta.assay.at(1).ms_run_ref.size(), 2)
    TEST_EQUAL(meta.assay.at(1).ms_run_ref[0], 1)
    TEST_EQUAL(meta.assay.at(1).ms_run_ref[1], 2)
    TEST_EQUAL(meta.assay.at(2).ms_run_ref.size(), 1)
    TEST_EQUAL(meta.assay.at(2).ms_run_ref[0], 3)
  }
  TEST_EQUAL(meta.study_variable.size(), 1)
  TEST_EQUAL(meta.study_variable.at(1).assay_refs.size(), 2)
  if (meta.study_variable.at(1).assay_refs.size() == 2)
  {
    TEST_EQUAL(meta.study_variable.at(1).assay_refs[0], 1)
    TEST_EQUAL(meta.study_variable.at(1).assay_refs[1], 2)
  }

  const MzTabProteinSectionRow& protein_row = mz_tab.getProteinSectionRows().back();
  TEST_EQUAL(protein_row.protein_abundance_assay.size(), 2)
  if (protein_row.protein_abundance_assay.size() == 2)
  {
    TEST_REAL_SIMILAR(protein_row.protein_abundance_assay.at(1).get(), 10.0)
    TEST_REAL_SIMILAR(protein_row.protein_abundance_assay.at(2).get(), 20.0)
  }
  TEST_TRUE(protein_row.protein_abundance_study_variable.at(1).isNull())

  const MzTabPeptideSectionRow& peptide_row = mz_tab.getPeptideSectionRows().front();
  TEST_EQUAL(peptide_row.peptide_abundance_assay.size(), 2)
  if (peptide_row.peptide_abundance_assay.size() == 2)
  {
    TEST_REAL_SIMILAR(peptide_row.peptide_abundance_assay.at(1).get(), 5.0)
    TEST_REAL_SIMILAR(peptide_row.peptide_abundance_assay.at(2).get(), 7.0)
  }
  TEST_TRUE(peptide_row.peptide_abundance_study_variable.at(1).isNull())
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
