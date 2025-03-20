// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/ANALYSIS/XLMS/OPXLHelper.h>

#include <OpenMS/ANALYSIS/XLMS/OPXLDataStructs.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGeneratorXLMS.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/ANALYSIS/XLMS/OPXLSpectrumProcessingAlgorithms.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CHEMISTRY/Tagger.h>
#include <QStringList>

using namespace OpenMS;

START_TEST(OPXLHelper, "$Id$")

// loading and building data structures required in several following tests
std::vector<FASTAFile::FASTAEntry> fasta_db;
FASTAFile file;
file.load(OPENMS_GET_TEST_DATA_PATH("FASTAFile_test.fasta"), fasta_db);


ProteaseDigestion digestor;
String enzyme_name = "Trypsin";
digestor.setEnzyme(enzyme_name);
digestor.setMissedCleavages(2);

Size min_peptide_length = 5;

QStringList q_str_list1;
QStringList q_str_list2;
q_str_list1 << "Carbamidomethyl (C)" << "Carbamidomethyl (T)";
q_str_list2 << "Oxidation (M)" << "Oxidation (Y)";
StringList fixedModNames = StringListUtils::fromQStringList(q_str_list1);
StringList varModNames = StringListUtils::fromQStringList(q_str_list2);
const ModifiedPeptideGenerator::MapToResidueType& fixed_modifications = ModifiedPeptideGenerator::getModifications(fixedModNames);
const ModifiedPeptideGenerator::MapToResidueType& variable_modifications = ModifiedPeptideGenerator::getModifications(varModNames);

QStringList q_str_list3;
QStringList q_str_list4;
q_str_list3 << "K" << "E";
q_str_list4 << "D" << "E" << "C-term";
StringList cross_link_residue1 = StringListUtils::fromQStringList(q_str_list3);
StringList cross_link_residue2 = StringListUtils::fromQStringList(q_str_list4);

Size max_variable_mods_per_peptide = 5;

START_SECTION(static std::vector<OPXLDataStructs::AASeqWithMass> digestDatabase(std::vector<FASTAFile::FASTAEntry> fasta_db, EnzymaticDigestion digestor, Size min_peptide_length, StringList cross_link_residue1, StringList cross_link_residue2, std::vector<const ResidueModification*> fixed_modifications, std::vector<const ResidueModification*> variable_modifications, Size max_variable_mods_per_peptide))

  std::vector<OPXLDataStructs::AASeqWithMass> peptides = OPXLHelper::digestDatabase(fasta_db, digestor, min_peptide_length, cross_link_residue1, cross_link_residue2, fixed_modifications, variable_modifications, max_variable_mods_per_peptide);

  TEST_EQUAL(peptides.size(), 886)
  TEST_EQUAL(peptides[5].peptide_mass > 5, true) // not an empty AASequence
  TEST_EQUAL(peptides[5].peptide_mass, peptides[5].peptide_seq.getMonoWeight())
  TEST_EQUAL(peptides[500].peptide_mass > 5, true) // not an empty AASequence
  TEST_EQUAL(peptides[500].peptide_mass, peptides[500].peptide_seq.getMonoWeight())
  TEST_EQUAL(peptides[668].position, OPXLDataStructs::INTERNAL)
  TEST_EQUAL(peptides[778].position, OPXLDataStructs::N_TERM)
END_SECTION

// building more data structures required in several following tests
std::vector<OPXLDataStructs::AASeqWithMass> peptides = OPXLHelper::digestDatabase(fasta_db, digestor, min_peptide_length, cross_link_residue1, cross_link_residue2, fixed_modifications, variable_modifications, max_variable_mods_per_peptide);

std::sort(peptides.begin(), peptides.end(), OPXLDataStructs::AASeqWithMassComparator());

double cross_link_mass = 150.0;
double precursor_mass_tolerance = 10;
bool precursor_mass_tolerance_unit_ppm = true;

std::vector<String> mono_masses;
mono_masses.push_back("50.0");
DoubleList cross_link_mass_mono_link = ListUtils::create<double>(mono_masses);

std::vector< double > spectrum_precursors;
double first_mass = peptides[700].peptide_mass + peptides[800].peptide_mass + cross_link_mass;

for (Size i = 0; i < 1000; i++)
{
  spectrum_precursors.push_back(first_mass + (i/4));
}
std::sort(spectrum_precursors.begin(), spectrum_precursors.end());

START_SECTION(static std::vector<OPXLDataStructs::XLPrecursor> enumerateCrossLinksAndMasses(const std::vector<OPXLDataStructs::AASeqWithMass>&  peptides, double cross_link_mass_light, const DoubleList& cross_link_mass_mono_link, const StringList& cross_link_residue1, const StringList& cross_link_residue2, std::vector< double >& spectrum_precursors, vector< int >& precursor_correction_positions, double precursor_mass_tolerance, bool precursor_mass_tolerance_unit_ppm))

  std::cout << std::endl;
  std::vector< int > spectrum_precursor_correction_positions;
  std::vector<OPXLDataStructs::XLPrecursor> precursors = OPXLHelper::enumerateCrossLinksAndMasses(peptides, cross_link_mass, cross_link_mass_mono_link, cross_link_residue1, cross_link_residue2, spectrum_precursors, spectrum_precursor_correction_positions, precursor_mass_tolerance, precursor_mass_tolerance_unit_ppm);
  // std::sort(precursors.begin(), precursors.end(), OPXLDataStructs::XLPrecursorComparator());

  TOLERANCE_ABSOLUTE(1e-3)
  TEST_EQUAL(precursors.size(), 9604)
  TEST_EQUAL(spectrum_precursor_correction_positions.size(), 9604)
  // sample about 1/15 of the data, since a lot of precursors are generated

  for (Size i = 0; i < precursors.size(); i += 2000)
  {
    if (precursors[i].beta_index > peptides.size())
    {
      // mono-link
      TEST_REAL_SIMILAR(peptides[precursors[i].alpha_index].peptide_mass + cross_link_mass_mono_link[0], precursors[i].precursor_mass)
    }
    else
    {
      // cross-link
      double computed_precursor = peptides[precursors[i].alpha_index].peptide_mass + peptides[precursors[i].beta_index].peptide_mass + cross_link_mass;
      TEST_REAL_SIMILAR(computed_precursor, precursors[i].precursor_mass)
    }
  }

END_SECTION

// building more data structures required in the following test
std::cout << std::endl;
std::vector< int > spectrum_precursor_correction_positions;
std::vector<OPXLDataStructs::XLPrecursor> precursors = OPXLHelper::enumerateCrossLinksAndMasses(peptides, cross_link_mass, cross_link_mass_mono_link, cross_link_residue1, cross_link_residue2, spectrum_precursors, spectrum_precursor_correction_positions, precursor_mass_tolerance, precursor_mass_tolerance_unit_ppm);
std::sort(precursors.begin(), precursors.end(), OPXLDataStructs::XLPrecursorComparator());

START_SECTION(static std::vector <OPXLDataStructs::ProteinProteinCrossLink> buildCandidates(const std::vector< OPXLDataStructs::XLPrecursor > & candidates, const std::vector< int > precursor_corrections, std::vector< int >& precursor_correction_positions, const std::vector<OPXLDataStructs::AASeqWithMass> & peptide_masses, const StringList & cross_link_residue1, const StringList & cross_link_residue2, double cross_link_mass, const DoubleList & cross_link_mass_mono_link, std::vector< double >& spectrum_precursor_vector, std::vector< double >& allowed_error_vector, String cross_link_name))
  double precursor_mass = 11814.50296;
  double allowed_error = 0.1;
  String cross_link_name = "MyLinker";

  std::vector< OPXLDataStructs::XLPrecursor > filtered_precursors;

  // determine MS2 precursors that match to the current peptide mass
  std::vector< OPXLDataStructs::XLPrecursor >::const_iterator low_it;
  std::vector< OPXLDataStructs::XLPrecursor >::const_iterator up_it;

  low_it = std::lower_bound(precursors.begin(), precursors.end(), precursor_mass - allowed_error, OPXLDataStructs::XLPrecursorComparator());
  up_it = std::upper_bound(precursors.begin(), precursors.end(), precursor_mass + allowed_error, OPXLDataStructs::XLPrecursorComparator());

  if (low_it != up_it) // no matching precursor in data
  {
    for (; low_it != up_it; ++low_it)
    {
      filtered_precursors.push_back(*low_it);
    }
  }
  TEST_EQUAL(precursors.size(), 9604)
  TEST_EQUAL(filtered_precursors.size(), 32)
  std::vector< int > precursor_corrections(59, 0);
  std::vector< int > precursor_correction_positions(59, 0);
  std::vector< double > spectrum_precursor_vector(1, 0.0);
  std::vector< double > allowed_error_vector(1, allowed_error);

  std::vector <OPXLDataStructs::ProteinProteinCrossLink> spectrum_candidates = OPXLHelper::buildCandidates(filtered_precursors, precursor_corrections, precursor_correction_positions, peptides, cross_link_residue1, cross_link_residue2, cross_link_mass, cross_link_mass_mono_link, spectrum_precursor_vector, allowed_error_vector, cross_link_name);

  TEST_EQUAL(spectrum_candidates.size(), 1152)
  TEST_EQUAL(spectrum_candidates[50].cross_linker_name, "MyLinker")
  for (Size i = 0; i < spectrum_candidates.size(); i += 200)
  {
    TEST_REAL_SIMILAR(spectrum_candidates[i].alpha->getMonoWeight() + spectrum_candidates[i].beta->getMonoWeight() + spectrum_candidates[i].cross_linker_mass, precursor_mass)
  }

END_SECTION

// prepare data for the next three tests
std::vector< PeptideIdentification > peptide_ids;
std::vector< ProteinIdentification > protein_ids;
IdXMLFile id_file;

// this is an old test file, that can not be easily reproduced anymore,
// since it represents an intermediate state of the data structures and is not written out
// in this form anymore
// But it is very useful to test the functions that change the old structure to the new one
id_file.load(OPENMS_GET_TEST_DATA_PATH("OPXLHelper_test.idXML"), protein_ids, peptide_ids);

for (auto& id : peptide_ids)  //OMS_CODING_TEST_EXCLUDE
{
  for (auto& hit : id.getHits())  //OMS_CODING_TEST_EXCLUDE
  {
    hit.removeMetaValue("XL_Protein_position_alpha");
    hit.removeMetaValue("XL_Protein_position_beta");
    hit.removeMetaValue("xl_target_decoy");
    hit.removeMetaValue("accessions_beta");
  }
}

START_SECTION(static void addProteinPositionMetaValues(std::vector< PeptideIdentification > & peptide_ids))

  // test that the MetaValues were removed
  for (const auto& id : peptide_ids)
  {
    for (const auto& hit : id.getHits())
    {
      TEST_EQUAL(hit.metaValueExists("XL_Protein_position_alpha"), false)
      TEST_EQUAL(hit.metaValueExists("XL_Protein_position_beta"), false)
      TEST_EQUAL(hit.metaValueExists("xl_target_decoy"), false)
      TEST_EQUAL(hit.metaValueExists("accessions_beta"), false)
      TEST_EQUAL(hit.metaValueExists(Constants::UserParam::OPENPEPXL_XL_POS1_PROT), false)
      TEST_EQUAL(hit.metaValueExists(Constants::UserParam::OPENPEPXL_XL_POS2_PROT), false)
      TEST_EQUAL(hit.metaValueExists(Constants::UserParam::OPENPEPXL_TARGET_DECOY_ALPHA), false)
      TEST_EQUAL(hit.metaValueExists(Constants::UserParam::OPENPEPXL_TARGET_DECOY_BETA), false)
      TEST_EQUAL(hit.metaValueExists(Constants::UserParam::OPENPEPXL_BETA_ACCESSIONS), false)
      TEST_EQUAL(hit.metaValueExists(Constants::UserParam::TARGET_DECOY), true)
    }
  }

  // add protein position MetaValues
  OPXLHelper::addProteinPositionMetaValues(peptide_ids);

  // check, that they were added to every PeptideHit
  for (const auto& id : peptide_ids)
  {
    const PeptideHit& hit = id.getHits()[0];
    TEST_EQUAL(hit.metaValueExists(Constants::UserParam::OPENPEPXL_XL_POS1_PROT), true)
    TEST_EQUAL(hit.metaValueExists(Constants::UserParam::OPENPEPXL_XL_POS2_PROT), true)
    TEST_EQUAL(hit.metaValueExists(Constants::UserParam::OPENPEPXL_TARGET_DECOY_ALPHA), false)
    TEST_EQUAL(hit.metaValueExists(Constants::UserParam::OPENPEPXL_TARGET_DECOY_BETA), false)
    TEST_EQUAL(hit.metaValueExists(Constants::UserParam::OPENPEPXL_BETA_ACCESSIONS), false)
  }

  // a few example values
  TEST_EQUAL(peptide_ids[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1_PROT), "1539")
  TEST_EQUAL(peptide_ids[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2_PROT), "182")
  TEST_EQUAL(peptide_ids[1].getHits()[1].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1_PROT), "1539")
  TEST_EQUAL(peptide_ids[1].getHits()[1].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2_PROT), "182")


END_SECTION

START_SECTION(static void addXLTargetDecoyMV(std::vector< PeptideIdentification > & peptide_ids))

  // add xl_target_decoy MetaValue
  OPXLHelper::addXLTargetDecoyMV(peptide_ids);
  // check, that they were added to every PeptideHit
  for (const auto& id : peptide_ids)
  {
      const PeptideHit& hit = id.getHits()[0];
      TEST_EQUAL(hit.metaValueExists(Constants::UserParam::OPENPEPXL_TARGET_DECOY_ALPHA), true)
      TEST_EQUAL(hit.metaValueExists(Constants::UserParam::OPENPEPXL_TARGET_DECOY_BETA), true)
  }

  // a few example values
  TEST_EQUAL(peptide_ids[0].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_TARGET_DECOY_ALPHA), "target")
  TEST_EQUAL(peptide_ids[0].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_TARGET_DECOY_BETA), "-")
  TEST_EQUAL(peptide_ids[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_TARGET_DECOY_ALPHA), "target")
  TEST_EQUAL(peptide_ids[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_TARGET_DECOY_BETA), "target")

END_SECTION

START_SECTION(static void addBetaAccessions(std::vector< PeptideIdentification > & peptide_ids))

  // add accessions_beta MV
  OPXLHelper::addBetaAccessions(peptide_ids);
  // check, that they were added to every PeptideHit
  for (const auto& id : peptide_ids)
  {
    for (const auto& hit : id.getHits())
    {
      TEST_EQUAL(hit.metaValueExists(Constants::UserParam::OPENPEPXL_BETA_ACCESSIONS), true)
    }
  }

  // a few example values
  TEST_EQUAL(peptide_ids[0].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_BETA_ACCESSIONS), "-")
  TEST_EQUAL(peptide_ids[1].getHits()[1].getMetaValue(Constants::UserParam::OPENPEPXL_BETA_ACCESSIONS), "Protein1")

END_SECTION

START_SECTION(static std::vector< PeptideIdentification > combineTopRanksFromPairs(std::vector< PeptideIdentification > & peptide_ids, Size number_top_hits))

  std::vector< PeptideIdentification > pep_ids = peptide_ids;
  // all hits are to separate spectra, so everything should be rank 1
  for (const auto& id : pep_ids)
  {
    for (const auto& hit : id.getHits())
    {
      TEST_EQUAL(hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_RANK), 1)
    }
  }

  // artificially assign one of the hits to the spectrum of another
  pep_ids[1].getHits()[0].setMetaValue("spectrum_index", pep_ids[0].getHits()[0].getMetaValue("spectrum_index"));
  pep_ids[1].getHits()[1].setMetaValue("spectrum_index", pep_ids[0].getHits()[0].getMetaValue("spectrum_index"));

  pep_ids = OPXLHelper::combineTopRanksFromPairs(pep_ids, 5);

  // there is one rank 2 now (in peptide_ids[2] now, because the order is not preserved)
  TEST_EQUAL(pep_ids[2].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_RANK), 2)

END_SECTION

START_SECTION(static void removeBetaPeptideHits(std::vector< PeptideIdentification > & peptide_ids))

  OPXLHelper::removeBetaPeptideHits(peptide_ids);

  TEST_EQUAL(peptide_ids.size(), 3)
  for (const auto& id : peptide_ids)
  {
    TEST_EQUAL(id.getHits().size(), 1)
  }

  // a few example values
  // mono-link
  TEST_EQUAL(peptide_ids[0].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1_PROT), "2078")
  TEST_EQUAL(peptide_ids[0].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2_PROT), "-")
  // cross-link
  TEST_EQUAL(peptide_ids[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1_PROT), "1539")
  TEST_EQUAL(peptide_ids[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2_PROT), "182")
  TEST_EQUAL(peptide_ids[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_BETA_PEPEV_PRE), "K")
  TEST_EQUAL(peptide_ids[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_BETA_PEPEV_END), "189")


END_SECTION

START_SECTION(static void buildFragmentAnnotations(std::vector<PeptideHit::PeakAnnotation> & frag_annotations, const std::vector< std::pair< Size, Size > > & matching, const PeakSpectrum & theoretical_spectrum, const PeakSpectrum & experiment_spectrum))
  TheoreticalSpectrumGeneratorXLMS specGen;
  Param param = specGen.getParameters();
  param.setValue("add_isotopes", "false");
  param.setValue("add_metainfo", "true");
  param.setValue("add_first_prefix_ion", "false");
  param.setValue("add_a_ions", "false");
  param.setValue("add_losses", "false");
  param.setValue("add_precursor_peaks", "false");
  param.setValue("add_k_linked_ions", "false");
  specGen.setParameters(param);

  PeakSpectrum theo_spec, exp_spec;
  // Theoretical Spec with metainfo
  AASequence peptedi = AASequence::fromString("PEPTEDI");
  specGen.getLinearIonSpectrum(theo_spec, peptedi, 4, true);

  param.setValue("add_metainfo", "false");
  specGen.setParameters(param);

  // Theoretical Spec without metainfo (Pseudo experimental spectrum)
  AASequence peptide = AASequence::fromString("PEPTIDE");
  specGen.getLinearIonSpectrum(exp_spec, peptide, 3, true);
  std::vector <std::pair <Size, Size> > alignment;

  DataArrays::FloatDataArray dummy_array;
  DataArrays::IntegerDataArray dummy_charge_array;
  OPXLSpectrumProcessingAlgorithms::getSpectrumAlignmentFastCharge(alignment, 50, true, theo_spec, exp_spec, dummy_charge_array, dummy_charge_array, dummy_array);

  std::vector<PeptideHit::PeakAnnotation> frag_annotations;

  // test, that additional annotations are added and do not replace existing ones
  PeptideHit::PeakAnnotation frag_anno;
  frag_anno.annotation = "TEST";
  frag_anno.charge = 50;
  frag_anno.mz = 1.0;
  frag_anno.intensity = 5.0;
  frag_annotations.push_back(frag_anno);

  OPXLHelper::buildFragmentAnnotations(frag_annotations, alignment, theo_spec, exp_spec);

  // number of annotations should be equal to number of aligned peaks (+ 1 for manual "TEST" annotation)
  TEST_EQUAL(frag_annotations.size(), alignment.size() + 1)
  TEST_EQUAL(frag_annotations[0].charge, 50)
  TEST_EQUAL(frag_annotations[0].mz, 1.0)
  TEST_EQUAL(frag_annotations[0].intensity, 5.0)
  TEST_EQUAL(frag_annotations[0].annotation, "TEST")

  TEST_EQUAL(frag_annotations[1].charge, 1)
  TEST_REAL_SIMILAR(frag_annotations[1].mz, 98.06004)
  TEST_EQUAL(frag_annotations[1].intensity, 1)
  TEST_EQUAL(frag_annotations[1].annotation, "[alpha|ci$b1]")

  TEST_EQUAL(frag_annotations[3].charge, 1)
  TEST_REAL_SIMILAR(frag_annotations[3].mz, 324.15539)
  TEST_EQUAL(frag_annotations[3].intensity, 1)
  TEST_EQUAL(frag_annotations[3].annotation, "[alpha|ci$b3]")

END_SECTION

START_SECTION(static std::vector <OPXLDataStructs::ProteinProteinCrossLink> OPXLHelper::collectPrecursorCandidates(IntList precursor_correction_steps, double precursor_mass, double precursor_mass_tolerance, bool precursor_mass_tolerance_unit_ppm, std::vector<OPXLDataStructs::AASeqWithMass> filtered_peptide_masses, double cross_link_mass, DoubleList cross_link_mass_mono_link, StringList cross_link_residue1, StringList cross_link_residue2, String cross_link_name, bool use_sequence_tags, std::vector<std::string>& tags))

  IntList precursor_correction_steps;
  precursor_correction_steps.push_back(2);
  precursor_correction_steps.push_back(1);

  double precursor_mass = 10668.85060;
  String cross_link_name = "MyLinker";
  precursor_mass_tolerance = 10;

  std::vector <OPXLDataStructs::ProteinProteinCrossLink> spectrum_candidates = OPXLHelper::collectPrecursorCandidates(precursor_correction_steps, precursor_mass, precursor_mass_tolerance, precursor_mass_tolerance_unit_ppm, peptides, cross_link_mass, cross_link_mass_mono_link, cross_link_residue1, cross_link_residue2, cross_link_name);

  TEST_EQUAL(spectrum_candidates.size(), 1050)
  TEST_EQUAL(spectrum_candidates[50].cross_linker_name, "MyLinker")
  for (Size i = 0; i < spectrum_candidates.size(); i += 100)
  {
    TEST_REAL_SIMILAR(spectrum_candidates[i].alpha->getMonoWeight() + spectrum_candidates[i].beta->getMonoWeight() + spectrum_candidates[i].cross_linker_mass, precursor_mass - 1 * Constants::C13C12_MASSDIFF_U)
  }

END_SECTION

START_SECTION(static double OPXLHelper::computePrecursorError(OPXLDataStructs::CrossLinkSpectrumMatch csm, double precursor_mz, int precursor_charge))

  OPXLDataStructs::ProteinProteinCrossLink ppcl;
  AASequence alpha = AASequence::fromString("TESTPEPTIDE");
  AASequence beta = AASequence::fromString("TESTTESTESTE");
  ppcl.alpha = &alpha;
  ppcl.beta = &beta;
  ppcl.cross_linker_mass = 150.0;

  OPXLDataStructs::CrossLinkSpectrumMatch csm;
  csm.cross_link = ppcl;
  csm.precursor_correction = 0;

  double precursor_charge = 3;
  double precursor_mz = (ppcl.alpha->getMonoWeight() + ppcl.beta->getMonoWeight() + ppcl.cross_linker_mass + precursor_charge * Constants::PROTON_MASS_U) / precursor_charge;

  double rel_error = OPXLHelper::computePrecursorError(csm, precursor_mz, precursor_charge);
  TEST_REAL_SIMILAR(rel_error, 0)

  precursor_mz += 0.05;
  rel_error = OPXLHelper::computePrecursorError(csm, precursor_mz, precursor_charge);
  TEST_REAL_SIMILAR(rel_error, 56.21777)

END_SECTION

START_SECTION(static void OPXLHelper::isoPeakMeans(OPXLDataStructs::CrossLinkSpectrumMatch& csm, DataArrays::IntegerDataArray& num_iso_peaks_array, std::vector< std::pair< Size, Size > >& matched_spec_linear_alpha, std::vector< std::pair< Size, Size > >& matched_spec_linear_beta, std::vector< std::pair< Size, Size > >& matched_spec_xlinks_alpha, std::vector< std::pair< Size, Size > >& matched_spec_xlinks_beta))

  DataArrays::IntegerDataArray iso_peaks;
  iso_peaks.push_back(3);
  iso_peaks.push_back(5);
  iso_peaks.push_back(2);
  iso_peaks.push_back(1);
  iso_peaks.push_back(1);
  iso_peaks.push_back(3);
  iso_peaks.push_back(1);
  iso_peaks.push_back(3);
  iso_peaks.push_back(2);

  std::vector< std::pair< Size, Size > > matched_spec_linear_alpha;
  matched_spec_linear_alpha.push_back(std::make_pair(1,1));
  matched_spec_linear_alpha.push_back(std::make_pair(2,2));
  matched_spec_linear_alpha.push_back(std::make_pair(4,3));
  matched_spec_linear_alpha.push_back(std::make_pair(6,4));
  matched_spec_linear_alpha.push_back(std::make_pair(7,5));
  std::vector< std::pair< Size, Size > > matched_spec_linear_beta;
  std::vector< std::pair< Size, Size > > matched_spec_xlinks_alpha;
  std::vector< std::pair< Size, Size > > matched_spec_xlinks_beta;
  matched_spec_xlinks_beta.push_back(std::make_pair(3,1));
  matched_spec_xlinks_beta.push_back(std::make_pair(5,2));
  matched_spec_xlinks_beta.push_back(std::make_pair(8,3));
  matched_spec_xlinks_beta.push_back(std::make_pair(0,4));

  OPXLDataStructs::CrossLinkSpectrumMatch csm;
  OPXLHelper::isoPeakMeans(csm, iso_peaks, matched_spec_linear_alpha, matched_spec_linear_beta, matched_spec_xlinks_alpha, matched_spec_xlinks_beta);

  TEST_REAL_SIMILAR(csm.num_iso_peaks_mean, 2.3333)
  TEST_REAL_SIMILAR(csm.num_iso_peaks_mean_linear_alpha, 2.4)
  TEST_REAL_SIMILAR(csm.num_iso_peaks_mean_linear_beta, 0)
  TEST_REAL_SIMILAR(csm.num_iso_peaks_mean_xlinks_alpha, 0)
  TEST_REAL_SIMILAR(csm.num_iso_peaks_mean_xlinks_beta, 2.25)
END_SECTION

START_SECTION(filterPrecursorsByTags(std::vector <OPXLDataStructs::XLPrecursor>& candidates, std::vector<std::string>& tags))

  std::cout << std::endl;
  std::vector< int > spectrum_precursor_correction_positions;
  std::vector<OPXLDataStructs::XLPrecursor> precursors = OPXLHelper::enumerateCrossLinksAndMasses(peptides, cross_link_mass, cross_link_mass_mono_link, cross_link_residue1, cross_link_residue2, spectrum_precursors, spectrum_precursor_correction_positions, precursor_mass_tolerance, precursor_mass_tolerance_unit_ppm);

  // set of tags
  std::vector<std::string> tags = {"DE", "PP", "FDA", "CIA", "FTC", "ESA", "ISRO", "NASA", "JAXA"};

  TEST_EQUAL(precursors.size(), 9604);

  // filter candidates
  OPXLHelper::filterPrecursorsByTags(precursors, spectrum_precursor_correction_positions, tags);
  TEST_EQUAL(precursors.size(), 4372);


  // // hasSubstring method runtime benchmark: search those 4372 candidates that do not contain the tags many times
  // // with 30000 iterations: 126.80 sec

  // std::cout << std::endl;
  // for (int i = 0; i < 30000; ++i)
  // {
  //   OPXLHelper::filterPrecursorsByTags(precursors, spectrum_precursor_correction_positions, tags);
  // }
  // TEST_EQUAL(precursors.size(), 4372);

  // // Aho-Corasick method runtime benchmark: search those 4372 candidates that do not contain the tags many times
  // // with 30000 iterations: Timeout after 1500.10 sec
  // // with 3000 iterations: 200.92 sec
  // for (int i = 0; i < 3000; ++i)
  // {
  //   OPXLHelper::filterPrecursorsByTagTrie(precursors, spectrum_precursor_correction_positions, tags);
  // }
  // TEST_EQUAL(precursors.size(), 4372);

END_SECTION

END_TEST
