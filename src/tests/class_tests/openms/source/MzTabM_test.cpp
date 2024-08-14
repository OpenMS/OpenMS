// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka$
// $Authors: Oliver Alka$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/MzTabM.h>
#include <OpenMS/FORMAT/OMSFile.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>
///////////////////////////

START_TEST(MzTabM, "$Id$")

using namespace OpenMS;
using namespace std;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MzTabM* ptr = nullptr;
MzTabM* null_ptr = nullptr;
START_SECTION(MzTabM())
{
  ptr = new MzTabM();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~MzTabM())
{
  delete ptr;
}
END_SECTION

START_SECTION(Fill data structure)
{
  MzTabM mztabm;

  // SML Small molecule section row
  MzTabMSmallMoleculeSectionRows sml_rows;
  MzTabMSmallMoleculeSectionRow sml_row;
  sml_row.sml_identifier.fromCellString(1);
  sml_row.smf_id_refs.fromCellString("1,2");
  sml_row.database_identifier.fromCellString("[HMDB:HMDB0001847]");
  sml_row.chemical_formula.fromCellString("[C17H20N4O2]");
  sml_row.smiles.fromCellString("[C1=CC=C(C=C1)CCNC(=O)CCNNC(=O)C2=CC=NC=C2]");
  sml_row.inchi.fromCellString("[InChI=1S/C17H20N4O2/c22-16(19-12-6-14-4-2-1-3-5-14)9-13-20-21-17(23)15-7-10-18-11-8-15/h1-5,7-8,10-11,20H,6,9,12-13H2,(H,19,22)(H,21,23)]");
  sml_row.chemical_name.fromCellString("[N-(2-phenylethyl)-3-[2-(pyridine-4-carbonyl)hydrazinyl]propanamide]");
  sml_row.uri.fromCellString("[http://www.hmdb.ca/metabolites/HMDB0001847]");
  vector<MzTabDouble> tnm = {MzTabDouble(312.17)};
  sml_row.theoretical_neutral_mass.set(tnm);
  sml_row.adducts.fromCellString("[[M+H]1+]");
  sml_row.reliability.set("3");
  sml_row.best_id_confidence_measure.fromCellString("[MS, MS:1000752, TOPP Software,]");
  sml_row.best_id_confidence_value.set(0.4);

  MzTabOptionalColumnEntry e;
  MzTabString s;
  e.first = "SIRIUS_TREE_score";
  s.fromCellString("-10.59083");
  e.second = s;
  sml_row.opt_.emplace_back(e);

  e.first = "SIRIUS_explained_intensity_score";
  s.fromCellString("96.67");
  e.second = s;
  sml_row.opt_.emplace_back(e);

  e.first = "SIRIUS_ISO_score";
  s.fromCellString("0.0649874");
  e.second = s;
  sml_row.opt_.emplace_back(e);

  sml_rows.emplace_back(sml_row);

  // SMF Small molecule feature section
  MzTabMSmallMoleculeFeatureSectionRows smf_rows;
  MzTabMSmallMoleculeFeatureSectionRow smf_row;
  smf_row.smf_identifier.fromCellString(1);
  smf_row.sme_id_refs.fromCellString("1");
  smf_row.sme_id_ref_ambiguity_code.fromCellString("null");
  smf_row.adduct.fromCellString("[M+H]1+");
  smf_row.isotopomer.setNull(true);
  smf_row.exp_mass_to_charge.set(313.1689);
  smf_row.charge.set(1);
  smf_row.retention_time.set(156.0); // is always in seconds
  smf_row.rt_start.set(152.2);
  smf_row.rt_end.set(163.4);
  smf_rows.emplace_back(smf_row);

  // SME Small molecule evidence section
  MzTabMSmallMoleculeEvidenceSectionRows sme_rows;
  MzTabMSmallMoleculeEvidenceSectionRow sme_row;
  sme_row.sme_identifier.set(1);
  sme_row.evidence_input_id.set("1234.5_156.0");
  sme_row.database_identifier.set("HMDB:HMDB0001847");
  sme_row.chemical_formula.set("C17H20N4O2");
  sme_row.smiles.set("C1=CC=C(C=C1)CCNC(=O)CCNNC(=O)C2=CC=NC=C2");
  sme_row.inchi.set("InChI=1S/C17H20N4O2/c22-16(19-12-6-14-4-2-1-3-5-14)9-13-20-21-17(23)15-7-10-18-11-8-15/h1-5,7-8,10-11,20H,6,9,12-13H2,(H,19,22)(H,21,23)");
  sme_row.chemical_name.set("N-(2-phenylethyl)-3-[2-(pyridine-4-carbonyl)hydrazinyl]propanamide");
  sme_row.uri.set("http://www.hmdb.ca/metabolites/HMDB0001847");
  sme_row.derivatized_form.isNull();
  sme_row.adduct.set("[M+H]1+");
  sme_row.exp_mass_to_charge.set(313.1689);
  sme_row.charge.set(1);
  sme_row.calc_mass_to_charge.set(313.1665);
  MzTabSpectraRef sp_ref;
  sp_ref.setMSFile(1);
  sp_ref.setSpecRef("index=5");
  sme_row.spectra_ref = sp_ref;
  sme_row.identification_method.fromCellString("[MS, MS:1000752, TOPP Software,]");
  sme_row.ms_level.fromCellString("[MS, MS:1000511, ms level, 1]");
  sme_row.id_confidence_measure[0] = MzTabDouble(123);
  sme_row.rank.set(1);

  e.first = "SIRIUS_TREE_score";
  s.fromCellString("-10.59083");
  e.second = s;
  sme_row.opt_.emplace_back(e);

  e.first = "SIRIUS_explained_intensity_score";
  s.fromCellString("96.67");
  e.second = s;
  sme_row.opt_.emplace_back(e);

  e.first = "SIRIUS_ISO_score";
  s.fromCellString("0.0649874");
  e.second = s;
  sme_row.opt_.emplace_back(e);

  sme_rows.emplace_back(sme_row);

  // Metadata for MzTab-M
  MzTabMMetaData mztabm_meta;
  mztabm_meta.mz_tab_id.set("local_identifier");
  mztabm_meta.title.set("SML_ROW_TEST");
  mztabm_meta.description.set("small_molecule_section_row_test");

  // sample proceessing
  MzTabParameterList sp;
  sp.fromCellString("[MS, MS:1000544, Conversion to mzML, ]|[MS, MS:1000035, Peak picking, ]|[MS, MS:1000594, Low intensity data point removal, ]");
  mztabm_meta.sample_processing[0] = sp;

  // instrument
  MzTabInstrumentMetaData meta_instrument;
  meta_instrument.name.fromCellString("[MS, MS:1000483, Thermo Fisher Scientific instrument model, LTQ Orbitrap Velos]");
  meta_instrument.source.fromCellString("[MS, MS:1000008, Ionization Type, ESI]");
  MzTabParameter ana;
  ana.fromCellString("[MS, MS:1000443, Mass Analyzer Type, Orbitrap]");
  meta_instrument.analyzer[0] = ana;
  meta_instrument.detector.fromCellString("[MS, MS:1000453, Detector, Dynode Detector]");
  mztabm_meta.instrument[0] = meta_instrument;

  // software
  MzTabSoftwareMetaData meta_software;
  MzTabParameter p_software;
  p_software.fromCellString("[MS, MS:1002205, ProteoWizard msconvert, ]");
  meta_software.software = p_software;
  meta_software.setting[0] = MzTabString("Peak Picking MS1");
  mztabm_meta.software[0] = meta_software;

  mztabm_meta.publication[0] = MzTabString("pubmed:21063943|doi:10.1007/978-1-60761-987-1_6");

  // contact
  MzTabContactMetaData meta_contact;
  meta_contact.name = MzTabString("Max MusterMann");
  meta_contact.affiliation = MzTabString("University of Musterhausen");
  meta_contact.email = MzTabString("MMM@please_do_not_try_to_write_an_email.com");

  mztabm_meta.contact[0] = meta_contact;
  mztabm_meta.uri[0] = MzTabString("https://www.ebi.ac.uk/metabolights/MTBLS");
  mztabm_meta.external_study_uri[0] = MzTabString("https://www.ebi.ac.uk/metabolights/MTBLS/files/i_Investigation.txt");
  mztabm_meta.quantification_method.fromCellString("[MS, MS:1001834, LC-MS label-free quantitation analysis, ]");

  // sample
  MzTabSampleMetaData meta_sample;
  meta_sample.description = MzTabString("Nice Sample");
  mztabm_meta.sample[0] = meta_sample;

  // ms-run
  MzTabMMSRunMetaData meta_msrun;
  meta_msrun.location = MzTabString("ftp://ftp.ebi.ac.uk/path/to/file");
  meta_msrun.instrument_ref = MzTabInteger(0); // only if different instruments are used.
  MzTabParameter p_format;
  p_format.fromCellString("[MS, MS:1000584, mzML file, ]");
  meta_msrun.format = p_format;
  MzTabParameter p_id_format;
  p_id_format.fromCellString("[MS, MS:1000584, mzML file, ]");
  meta_msrun.id_format = p_id_format;
  std::map<Size, MzTabParameter> pl_fragmentation_method;
  pl_fragmentation_method[0].fromCellString("[MS, MS:1000133, CID, ]");
  pl_fragmentation_method[1].fromCellString("[MS, MS:1000422, HCD, ]");
  meta_msrun.fragmentation_method = pl_fragmentation_method;
  std::map<Size, MzTabParameter> pl_scan_polarity;
  pl_scan_polarity[0].fromCellString("[MS, MS:1000130, positive scan, ]");
  pl_scan_polarity[1].fromCellString("[MS, MS:1000130, positive scan, ]");
  meta_msrun.scan_polarity = pl_scan_polarity;
  meta_msrun.hash = MzTabString("de9f2c7fd25e1b3afad3e85a0bd17d9b100db4b3");
  MzTabParameter p_hash_method;
  p_hash_method.fromCellString("[MS, MS:1000569, SHA-1, ]");
  meta_msrun.hash_method = p_hash_method;
  mztabm_meta.ms_run[0] = meta_msrun;

  // assay
  MzTabMAssayMetaData meta_assay;
  MzTabParameter p_custom;
  p_custom.fromCellString("[MS, , Assay operator, Blogs]");
  meta_assay.custom[0] = p_custom;
  meta_assay.external_uri = MzTabString("https://www.ebi.ac.uk/metabolights/MTBLS/files/i_Investigation.txt?STUDYASSAY=a_8pos.txt");
  meta_assay.sample_ref = MzTabInteger(1);
  meta_assay.ms_run_ref = MzTabInteger(1);
  mztabm_meta.assay[0] = meta_assay;

  // study variable
  MzTabMStudyVariableMetaData meta_study;
  std::vector<int> assay_refs{1};
  meta_study.assay_refs = assay_refs;
  MzTabParameter p_average_function;
  p_average_function.fromCellString("[MS, MS:1002883, median, ]");
  meta_study.average_function = p_average_function;
  MzTabParameter p_variation_function;
  p_variation_function.fromCellString("[MS, MS:1002885, standard error, ]"); // usually we will not average!
  meta_study.variation_function = p_variation_function;
  meta_study.description = MzTabString("control");
  MzTabParameterList pl_factors;
  pl_factors.fromCellString("[MS, MS:1000130, positive scan, ]");
  meta_study.factors = pl_factors;
  mztabm_meta.study_variable[0] = meta_study;

  // controlled vocabulary metadata
  MzTabCVMetaData meta_cv;
  meta_cv.label = MzTabString("MS");
  meta_cv.full_name = MzTabString("PSI-MS controlled vocabulary");
  meta_cv.version = MzTabString("4.1.155");
  meta_cv.url = MzTabString("share/OpenMS/CV/psi-ms.obo");
  mztabm_meta.cv[0] = meta_cv;

  // database
  MzTabMDatabaseMetaData meta_db;
  MzTabParameter p_db;
  p_db.fromCellString("[MIRIAM, MIR:00100079, HMDB, ]");
  meta_db.database = p_db;
  meta_db.prefix = MzTabString("HMDB");
  meta_db.version = MzTabString("4.0");
  meta_db.uri = MzTabString("null");
  mztabm_meta.database[0] = meta_db;

  MzTabParameter p_qunit;
  p_qunit.fromCellString("[MS, MS:1000042, peak intensity, ]");
  mztabm_meta.small_molecule_quantification_unit = p_qunit;
  MzTabParameter p_fqunit;
  p_fqunit.fromCellString("[MS, MS:1000042, peak intensity, ]");
  mztabm_meta.small_molecule_feature_quantification_unit = p_fqunit;
  MzTabParameter p_idre;
  p_idre.fromCellString("[MS, MS:1002955, hr-ms compound identification confidence level, ]");
  mztabm_meta.small_molecule_identification_reliability = p_idre;
  MzTabParameter p_confidence;
  p_confidence.fromCellString("[MS,MS:1002890,fragmentation score,]");
  mztabm_meta.id_confidence_measure[0] = p_confidence;

  // Fill mztab-m datastructure
  mztabm.setMetaData(mztabm_meta);
  mztabm.setMSmallMoleculeSectionRows(sml_rows);
  mztabm.setMSmallMoleculeFeatureSectionRows(smf_rows);
  mztabm.setMSmallMoleculeEvidenceSectionRows(sme_rows);

  // Tests ///////////////////////////////
  MzTabMSmallMoleculeSectionRow sml_test;
  sml_test = mztabm.getMSmallMoleculeSectionRows()[0];
  TEST_EQUAL(sml_test.smf_id_refs.toCellString(), "1,2")
  TEST_EQUAL(sml_test.adducts.toCellString(), "[[M+H]1+]")

  MzTabMSmallMoleculeFeatureSectionRow smf_test;
  smf_test = mztabm.getMSmallMoleculeFeatureSectionRows()[0];
  TEST_EQUAL(smf_test.exp_mass_to_charge.toCellString(), "313.168900000000008")
  TEST_EQUAL(smf_test.retention_time.toCellString(), "156.0")

  MzTabMSmallMoleculeEvidenceSectionRow sme_test;
  sme_test = mztabm.getMSmallMoleculeEvidenceSectionRows()[0];
  TEST_EQUAL(sme_test.database_identifier.toCellString(), "HMDB:HMDB0001847")
  TEST_EQUAL(sme_test.identification_method.toCellString(), "[MS, MS:1000752, TOPP Software, ]")

  MzTabMMetaData mtest;
  mtest = mztabm.getMetaData();
  TEST_EQUAL(mtest.mz_tab_version.toCellString(),"2.0.0-M") // set by constructor
  TEST_EQUAL(mtest.sample_processing[0].toCellString(), "[MS, MS:1000544, Conversion to mzML, ]|[MS, MS:1000035, Peak picking, ]|[MS, MS:1000594, Low intensity data point removal, ]")
  TEST_EQUAL(mtest.instrument[0].analyzer[0].toCellString(), "[MS, MS:1000443, Mass Analyzer Type, Orbitrap]")
  // meta_software.setting[0] = MzTabString("Peak Picking MS1");
  TEST_EQUAL(mtest.software[0].setting[0].toCellString(), "Peak Picking MS1")
  // meta_contact.affiliation = MzTabString("University of Musterhausen");
  TEST_EQUAL(mtest.contact[0].affiliation.toCellString(), "University of Musterhausen")
  // meta_sample.description = MzTabString("Nice Sample");
  TEST_EQUAL(mtest.sample[0].description.toCellString(), "Nice Sample")
  // p_format.fromCellString("[MS, MS:1000584, mzML file, ]");
  TEST_EQUAL(mtest.ms_run[0].format.toCellString(), "[MS, MS:1000584, mzML file, ]")
  // meta_study.description = MzTabString("control");
  TEST_EQUAL(mtest.study_variable[0].description.toCellString(), "control")
  // meta_db.prefix = MzTabString("HMDB");
  TEST_EQUAL(mtest.database[0].prefix.toCellString(), "HMDB")
  // p_qunit.fromCellString("[MS, MS:1000042, peak intensity, ]");
  TEST_EQUAL(mtest.small_molecule_quantification_unit.toCellString(), "[MS, MS:1000042, peak intensity, ]")

  vector<String> optional_sml_columns = mztabm.getMSmallMoleculeOptionalColumnNames();
  vector<String> optional_sme_columns = mztabm.getMSmallMoleculeEvidenceOptionalColumnNames();

  TEST_EQUAL(mztabm.getMSmallMoleculeSectionRows().size(),1)
  TEST_EQUAL(mztabm.getMSmallMoleculeFeatureSectionRows().size(), 1)
  TEST_EQUAL(mztabm.getMSmallMoleculeFeatureSectionRows().size(), 1)

  TEST_EQUAL(optional_sml_columns.size(), 3)
  TEST_EQUAL(optional_sme_columns.size(), 3)
}
END_SECTION

START_SECTION(MzTabM::exportFeatureMapToMzTabM(const FeatureMap& feature_map))
{
  FeatureMap feature_map;
  MzTabM mztabm;

  OMSFile().load(OPENMS_GET_TEST_DATA_PATH("MzTabMFile_input_1.oms"), feature_map);

  mztabm = mztabm.exportFeatureMapToMzTabM(feature_map);

  TEST_EQUAL(mztabm.getMSmallMoleculeSectionRows().size(), 83)
  TEST_EQUAL(mztabm.getMSmallMoleculeFeatureSectionRows().size(), 83)
  TEST_EQUAL(mztabm.getMSmallMoleculeEvidenceSectionRows().size(), 312)

  TEST_EQUAL(mztabm.getMSmallMoleculeOptionalColumnNames().size(), 0)
  TEST_EQUAL(mztabm.getMSmallMoleculeFeatureOptionalColumnNames().size(), 18)
  TEST_EQUAL(mztabm.getMSmallMoleculeEvidenceOptionalColumnNames().size(), 6)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
