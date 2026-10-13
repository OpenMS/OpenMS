// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
#define OPENMS_NASE_CLASS_TEST
#include "../NucleicAcidSearchEngine.cpp"
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include "RNAPrecursorCorrector.h"
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <algorithm>
#include <fstream>
#include <iterator>
#include <limits>

class TestNASE : public NucleicAcidSearchEngine
{
public:
  using NucleicAcidSearchEngine::getModifications_;
  using NucleicAcidSearchEngine::calculatePrecursorPercentage_;
  using NucleicAcidSearchEngine::flashMassUpperBound_;
  using NucleicAcidSearchEngine::updatePrecursorPercentage_;
  using NucleicAcidSearchEngine::filterBlacklistedIons_;
  using NucleicAcidSearchEngine::calculateFragmentCoverage_;
  void preprocess(MSExperiment& exp, bool remove, const std::string& method)
  {
    preprocessSpectra_(exp, 10.0, true, true, true, -1, -5, true,
      false, 100, 50, "slide", false, 1000, remove, 0.01, false, 2, method, 0.85, 0.25);
  }
  using NucleicAcidSearchEngine::resolveAmbiguousMods_;
  using NucleicAcidSearchEngine::HitsByScore;
  using NucleicAcidSearchEngine::AnnotatedHit;
  using NucleicAcidSearchEngine::PrecursorInfo;
  using NucleicAcidSearchEngine::retainTopHit_;
  using NucleicAcidSearchEngine::sortHitsByValue_;
  using NucleicAcidSearchEngine::generateLFQInput_;
  void enableAmbiguities() { resolve_ambiguous_mods_ = true; }
  bool ambiguitiesEnabled() const { return resolve_ambiguous_mods_; }
};

static Int runSearch(const std::vector<std::string>& arguments)
{
  std::vector<const char*> argv;
  argv.push_back("NucleicAcidSearchEngine");
  for (const auto& arg : arguments)
  {
    argv.push_back(arg.c_str());
  }
  TestNASE tool;
  return tool.main(static_cast<Int>(argv.size()), argv.data());
}

static Size hitCount(const std::string& path)
{
  std::vector<ProteinIdentification> proteins;
  PeptideIdentificationList peptides;
  IdXMLFile().load(path, proteins, peptides);
  Size count = 0;
  for (const auto& peptide : peptides)
  {
    count += peptide.getHits().size();
  }
  return count;
}

static MSSpectrum isotopeSpectrum(Int charge)
{
  MSSpectrum spectrum;
  spectrum.setMSLevel(2);
  spectrum.setNativeID("scan=1");
  spectrum.setRT(20);
  Precursor precursor;
  precursor.setMZ(500);
  precursor.setCharge(charge);
  spectrum.setPrecursors({precursor});
  for (Size isotope = 0; isotope < 3; ++isotope)
  {
    spectrum.push_back(Peak1D(500 + isotope * Constants::C13C12_MASSDIFF_U / 2, 100 - isotope * 20));
  }
  return spectrum;
}

namespace
{
  MSSpectrum rnaFragmentPattern(int charge, bool positive = false)
  {
    MSSpectrum s;
    s.setMSLevel(2);
    s.setType(SpectrumSettings::SpectrumType::CENTROID);
    s.setNativeID("scan=17");
    s.setRT(12.5);
    Precursor p;
    p.setMZ(1500.0);
    p.setCharge(positive ? 5 : -5);
    s.getPrecursors().push_back(p);
    auto formula = NASequence::fromString("ACGUACGU").getFormula();
    const double mass = formula.getMonoWeight();
    const double carrier = positive ? Constants::PROTON_MASS_U : -Constants::PROTON_MASS_U;
    auto isotopes = formula.getIsotopeDistribution(CoarseIsotopePatternGenerator(8));
    for (Size i = 0; i < isotopes.size(); ++i)
    {
      s.emplace_back((mass + i * Constants::C13C12_MASSDIFF_U) / charge + carrier,
                     isotopes[i].getIntensity() * 10000);
    }
    s.emplace_back(83.123, 9.0);
    s.emplace_back(1799.123, 11.0); // Unassigned signal must survive.
    s.sortByPosition();
    MSSpectrum::IntegerDataArray ids;
    ids.setName("original_index");
    for (Size i = 0; i < s.size(); ++i) ids.push_back(i);
    s.getIntegerDataArrays().push_back(ids);
    return s;
  }
}

namespace
{
  PeakGroup envelope(double mass = 2500, int charge = 2, bool positive = false, int first_isotope = 0)
  {
    PeakGroup group(charge, charge, positive);
    group.setMonoisotopicMass(mass);
    group.setChargeIsotopeCosine(charge, 0.95);
    group.setChargeSNR(charge, 10);
    const double carrier = positive ? Constants::PROTON_MASS_U : -Constants::PROTON_MASS_U;
    for (int i = first_isotope; i <= 8; ++i)
    {
      FLASHHelperClasses::LogMzPeak peak(Peak1D((mass + i * Constants::C13C12_MASSDIFF_U) / charge + carrier, 1000), positive);
      peak.abs_charge = charge;
      peak.isotopeIndex = i;
      group.push_back(peak);
    }
    return group;
  }

  MSSpectrum survey(bool positive = false, const std::string& sequence = "ACGUACGU")
  {
    MSSpectrum spectrum;
    spectrum.setMSLevel(1);
    spectrum.setNativeID("scan=1");
    spectrum.setRT(10);
    spectrum.setType(SpectrumSettings::SpectrumType::CENTROID);
    const auto formula = NASequence::fromString(sequence).getFormula();
    const auto isotopes = formula.getIsotopeDistribution(CoarseIsotopePatternGenerator(8));
    const double carrier = positive ? Constants::PROTON_MASS_U : -Constants::PROTON_MASS_U;
    for (Size i = 0; i < isotopes.size(); ++i)
    {
      spectrum.emplace_back((formula.getMonoWeight() + i * Constants::C13C12_MASSDIFF_U) / 2 + carrier,
                            isotopes[i].getIntensity() * 10000);
    }
    spectrum.sortByPosition();
    return spectrum;
  }

  MSSpectrum ms2(const MSSpectrum& ms1, int isotope = 1)
  {
    MSSpectrum spectrum;
    spectrum.setMSLevel(2);
    spectrum.setNativeID("scan=3");
    spectrum.setRT(12);
    spectrum.emplace_back(100.0, 500);
    Precursor precursor;
    precursor.setMZ(ms1[isotope].getMZ());
    precursor.setCharge(2);
    precursor.setIsolationWindowLowerOffset(0.5);
    precursor.setIsolationWindowUpperOffset(0.75);
    precursor.setMetaValue("spectrum_ref", "scan=1");
    spectrum.getPrecursors().push_back(precursor);
    return spectrum;
  }
}

START_TEST(NASEDigestion, "$Id$")

START_SECTION(([EXTRA] variable-only ambiguity setup preserves terminal hypotheses))
{
  TestNASE tool;
  tool.enableAmbiguities();
  tool.getModifications_({});
  TEST_TRUE(tool.ambiguitiesEnabled())
  TEST_EQUAL(tool.getModifications_({"mA?"}).size(), 2)
  IdentificationData data;
  TestNASE::HitsByScore hits;
  TestNASE::AnnotatedHit five, three;
  five.sequence = NASequence::fromString("pACGU");
  three.sequence = NASequence::fromString("ACGUp");
  five.oligo_ref = data.registerIdentifiedOligo(IdentificationData::IdentifiedOligo(five.sequence));
  three.oligo_ref = five.oligo_ref;
  hits.emplace(10.0, five);
  hits.emplace(10.0, three);
  tool.resolveAmbiguousMods_(hits);
  TEST_EQUAL(hits.size(), 2)
}
END_SECTION

START_SECTION(([EXTRA] searchable RNA candidates and cached digest provenance))
{
  std::string mzml, fasta, output, digest;
  NEW_TMP_FILE(mzml)
  NEW_TMP_FILE(fasta)
  NEW_TMP_FILE(output)
  NEW_TMP_FILE(digest)
  mzml += ".mzML";
  fasta += ".fasta";
  output += ".idXML";
  digest += ".oms";
  NucleicAcidSpectrumGenerator generator;
  Param parameters = generator.getParameters();
  parameters.setValue("add_first_prefix_ion", "true");
  parameters.setValue("add_a-B_ions", "true");
  parameters.setValue("add_metainfo", "true");
  generator.setParameters(parameters);
  const auto fixture = [&](const std::string& parent, const std::string& hit, Int charge)
  {
    std::ofstream(fasta) << ">TARGET\n" << parent << '\n';
    const auto sequence = NASequence::fromString(hit);
    MSSpectrum spectrum;
    generator.getSpectrum(spectrum, sequence, -1, -2);
    spectrum.setMSLevel(2);
    spectrum.setNativeID("scan=1");
    spectrum.setRT(20);
    Precursor precursor;
    precursor.setMZ(sequence.getMonoWeight(NASequence::Full, -2) / 2);
    precursor.setCharge(charge);
    spectrum.setPrecursors({precursor});
    MSExperiment exp;
    exp.addSpectrum(spectrum);
    FileHandler().storeExperiment(mzml, exp, {FileTypes::MZML});
  };
  const auto arguments = [&](const std::string& enzyme = "no cleavage")
  {
    return std::vector<std::string>{"-in", mzml, "-database", fasta, "-id_out", output, "-out", output + ".mzTab",
      "-oligo:enzyme", enzyme,
      "-threads", "1", "-oligo:min_size", "1", "-oligo:missed_cleavages", "0"};
  };
  fixture("ACGUAC", "ACGUAC", 2);
  auto args = arguments();
  args.insert(args.end(), {"-digest_out", digest});
  TEST_EQUAL(runSearch(args), 0)
  TEST_EQUAL(hitCount(output), 1)

  // Reuse the digest with different search tolerances while preserving its database.
  args = {"-in", mzml, "-digest", digest, "-id_out", output, "-out", output + ".mzTab",
    "-precursor:mass_tolerance", "20", "-fragment:mass_tolerance", "25", "-threads", "1",
    "-db_out", output + ".oms"};
  TEST_EQUAL(runSearch(args), 0)
  std::vector<ProteinIdentification> proteins;
  PeptideIdentificationList peptides;
  IdXMLFile().load(output, proteins, peptides);
  TEST_EQUAL(proteins.size(), 1)
  ABORT_IF(proteins.size() != 1)
  TEST_REAL_SIMILAR(proteins[0].getSearchParameters().precursor_mass_tolerance, 20)
  TEST_REAL_SIMILAR(proteins[0].getSearchParameters().fragment_mass_tolerance, 25)
  TEST_TRUE(proteins[0].getSearchParameters().db == fasta)
  IdentificationData result;
  OMSFile().load(output + ".oms", result);
  bool current_identification = false;
  for (const auto& match : result.getObservationMatches())
  {
    for (const auto& applied : match.steps_and_scores)
    {
      if (!applied.processing_step_opt ||
          !(*applied.processing_step_opt)->actions.count(DataProcessing::IDENTIFICATION)) continue;
      const auto mapping = result.getDBSearchSteps().find(*applied.processing_step_opt);
      TEST_TRUE(mapping != result.getDBSearchSteps().end())
      if (mapping == result.getDBSearchSteps().end()) continue;
      current_identification = true;
      TEST_REAL_SIMILAR(mapping->second->precursor_mass_tolerance, 20)
      TEST_REAL_SIMILAR(mapping->second->fragment_mass_tolerance, 25)
    }
  }
  TEST_TRUE(current_identification)

  for (const std::string credit : {"true", "false"})
  {
    auto scoring_args = arguments();
    scoring_args.insert(scoring_args.end(), {"-scoring:method", "mvh", "-scoring:single_peak_credit", credit});
    TEST_EQUAL(runSearch(scoring_args), 0)
    TEST_EQUAL(hitCount(output), 1)
    vector<ProteinIdentification> proteins;
    PeptideIdentificationList peptides;
    IdXMLFile().load(output, proteins, peptides);
    TEST_STRING_EQUAL(peptides[0].getScoreType(), "MVH score")
  }

  auto flash_args = arguments();
  flash_args.insert(flash_args.end(), {"-preprocessing:deisotoping", "flash"});
  TEST_EQUAL(runSearch(flash_args), 0)
  TEST_EQUAL(hitCount(output), 1)
  auto incompatible_args = arguments();
  incompatible_args.insert(incompatible_args.end(), {"-preprocessing:deisotoping", "none", "-decharge_ms2"});
  TEST_NOT_EQUAL(runSearch(incompatible_args), 0)
  incompatible_args = arguments();
  incompatible_args.insert(incompatible_args.end(), {"-preprocessing:deisotoping", "flash", "-fragment:mass_tolerance_unit", "Da"});
  TEST_NOT_EQUAL(runSearch(incompatible_args), 0)

  fixture("ACGUAC", "ACGUAC", -2);
  TEST_EQUAL(runSearch(arguments()), 0)
  TEST_EQUAL(hitCount(output), 1)

  MSExperiment raw;
  MzMLFile().load(mzml, raw);
  raw[0].push_back(Peak1D(raw[0].getPrecursors()[0].getMZ(), 10000.0));
  for (auto& array : raw[0].getStringDataArrays()) array.push_back("precursor");
  for (auto& array : raw[0].getIntegerDataArrays()) array.push_back(0);
  for (auto& array : raw[0].getFloatDataArrays()) array.push_back(0.0f);
  raw[0].sortByPosition();
  MzMLFile().store(mzml, raw);
  std::string oms_output;
  NEW_TMP_FILE(oms_output)
  oms_output += ".oms";
  auto report_args = arguments();
  report_args.insert(report_args.end(), {"-report:precursor_percentage", "-preprocessing:remove_precursor_peak", "-db_out", oms_output});
  TEST_EQUAL(runSearch(report_args), 0)
  vector<ProteinIdentification> report_proteins;
  PeptideIdentificationList report_peptides;
  IdXMLFile().load(output, report_proteins, report_peptides);
  TEST_TRUE(!report_peptides.empty())
  ABORT_IF(report_peptides.empty())
  TEST_TRUE(!report_peptides[0].getHits().empty())
  ABORT_IF(report_peptides[0].getHits().empty())
  TEST_REAL_SIMILAR(double(report_peptides[0].getMetaValue("precursor_percentage")), 100.0)
  TEST_FALSE(report_peptides[0].getHits()[0].metaValueExists("precursor_percentage"))
  IdentificationData report_data;
  OMSFile().load(oms_output, report_data);
  TEST_FALSE(report_data.getObservationMatches().empty())
  ABORT_IF(report_data.getObservationMatches().empty())
  TEST_REAL_SIMILAR(double(report_data.getObservations().begin()->getMetaValue("precursor_percentage")), 100.0)
  TEST_FALSE(report_data.getObservationMatches().begin()->metaValueExists("precursor_percentage"))
  TextFile mztab_output(output + ".mzTab");
  TEST_TRUE(std::any_of(mztab_output.begin(), mztab_output.end(), [](const std::string& line)
    { return line.find("opt_global_precursor_percentage") != std::string::npos; }))

  fixture("ACGUAC", "ACGUAC", 2);
  MSExperiment correction_input;
  MzMLFile().load(mzml, correction_input);
  auto survey_scan = survey(false, "ACGUAC");
  survey_scan.setNativeID("survey=1");
  survey_scan.setRT(correction_input[0].getRT() - 1.0);
  correction_input[0].getPrecursors()[0].setMZ(
    correction_input[0].getPrecursors()[0].getMZ() + Constants::C13C12_MASSDIFF_U / 2);
  correction_input[0].getPrecursors()[0].setMetaValue("spectrum_ref", "survey=1");
  correction_input.setSpectra({survey_scan, correction_input[0]});
  MzMLFile().store(mzml, correction_input);
  auto correction_args = arguments();
  correction_args.insert(correction_args.end(), {"-precursor:isotopes", "-1", "-precursor:max_charge", "-5"});
  TEST_EQUAL(runSearch(correction_args), 0)
  TEST_EQUAL(hitCount(output), 0)
  std::string audit_output;
  NEW_TMP_FILE(audit_output)
  audit_output += ".tsv";
  correction_args.insert(correction_args.end(), {"-precursor:correct_with_flash", "-out_precursor_corrections", audit_output});
  TEST_EQUAL(runSearch(correction_args), 0)
  TEST_EQUAL(hitCount(output), 1)
  vector<ProteinIdentification> corrected_proteins;
  PeptideIdentificationList corrected_peptides;
  IdXMLFile().load(output, corrected_proteins, corrected_peptides);
  TEST_FALSE(corrected_peptides.empty())
  ABORT_IF(corrected_peptides.empty())
  TEST_FALSE(corrected_peptides[0].getHits().empty())
  ABORT_IF(corrected_peptides[0].getHits().empty())
  TEST_STRING_EQUAL(corrected_peptides[0].getHits()[0].getMetaValue("flash_precursor_status").toString(), "corrected")
  TEST_EQUAL(Int(corrected_peptides[0].getHits()[0].getMetaValue("flash_precursor_isotope")), 1)
  TextFile correction_report(audit_output);
  TEST_EQUAL(std::distance(correction_report.begin(), correction_report.end()), 2)

  // mzML preserves earlier FLASH userParams. They must neither collapse this
  // run's configured isotope list when disabled nor leak stale provenance.
  const double selected_mz = correction_input[1].getPrecursors()[0].getMZ();
  correction_input[1].setMetaValue("flash_precursor_status", "corrected");
  correction_input[1].setMetaValue("flash_precursor_original_mz", 42.0);
  correction_input[1].setMetaValue("flash_precursor_mz", 42.0);
  correction_input[1].setMetaValue("flash_precursor_mono_mass", 42.0);
  correction_input[1].setMetaValue("flash_precursor_isotope", 77);
  correction_input[1].setMetaValue("flash_precursor_cosine", 0.99);
  correction_input[1].setMetaValue("flash_precursor_snr", 123.0);
  correction_input[1].setMetaValue("flash_precursor_mono_observed", "true");
  correction_input[1].setMetaValue("flash_precursor_old_counter", 999);
  correction_input[1].setMetaValue("sample_note", "retain");
  MzMLFile().store(mzml, correction_input);
  std::string stale_preprocessed, stale_oms;
  NEW_TMP_FILE(stale_preprocessed)
  NEW_TMP_FILE(stale_oms)
  stale_preprocessed += ".mzML";
  stale_oms += ".oms";
  auto stale_args = arguments();
  stale_args.insert(stale_args.end(), {"-precursor:isotopes", "1", "-precursor:max_charge", "-5",
    "-out_preprocessed", stale_preprocessed, "-db_out", stale_oms});
  TEST_EQUAL(runSearch(stale_args), 0)
  TEST_EQUAL(hitCount(output), 1)
  MSExperiment processed;
  MzMLFile().load(stale_preprocessed, processed);
  TEST_EQUAL(processed.size(), 1)
  ABORT_IF(processed.size() != 1)
  TEST_EQUAL(processed[0].getPrecursors()[0].getMZ(), selected_mz)
  TEST_EQUAL(processed[0].getMetaValue("sample_note").toString(), "retain")
  std::vector<std::string> metadata_keys;
  processed[0].getKeys(metadata_keys);
  for (const auto& key : metadata_keys) TEST_FALSE(key.starts_with("flash_precursor_"))
  IdentificationData stale_data;
  OMSFile().load(stale_oms, stale_data);
  TEST_FALSE(stale_data.getObservationMatches().empty())
  ABORT_IF(stale_data.getObservationMatches().empty())
  for (const auto& match : stale_data.getObservationMatches())
  {
    TEST_EQUAL(Int(match.getMetaValue("isotope_offset")), 1)
    metadata_keys.clear();
    match.getKeys(metadata_keys);
    for (const auto& key : metadata_keys) TEST_FALSE(key.starts_with("flash_precursor_"))
  }

  // Enabled correction derives a new decision from the actual selected ion,
  // without restoring the unrelated old original-m/z annotation.
  auto fresh_args = correction_args;
  fresh_args.insert(fresh_args.end(), {"-out_preprocessed", stale_preprocessed});
  TEST_EQUAL(runSearch(fresh_args), 0)
  TEST_EQUAL(hitCount(output), 1)
  MzMLFile().load(stale_preprocessed, processed);
  TEST_EQUAL(processed[0].getMetaValue("flash_precursor_status").toString(), "corrected")
  TEST_EQUAL(double(processed[0].getMetaValue("flash_precursor_original_mz")), selected_mz)
  TEST_EQUAL(Int(processed[0].getMetaValue("flash_precursor_isotope")), 1)
  TEST_FALSE(processed[0].metaValueExists("flash_precursor_old_counter"))
  TEST_EQUAL(processed[0].getMetaValue("sample_note").toString(), "retain")

  // Missing MS1 evidence makes the same annotated input unresolved. The
  // static list still finds the hit, and accepted-only annotations disappear.
  correction_input.setSpectra({correction_input[1]});
  MzMLFile().store(mzml, correction_input);
  stale_args.push_back("-precursor:correct_with_flash");
  TEST_EQUAL(runSearch(stale_args), 0)
  TEST_EQUAL(hitCount(output), 1)
  MzMLFile().load(stale_preprocessed, processed);
  TEST_EQUAL(processed[0].getMetaValue("flash_precursor_status").toString(), "no_ms1")
  TEST_EQUAL(processed[0].getPrecursors()[0].getMZ(), selected_mz)
  for (const auto& key : {"flash_precursor_mz", "flash_precursor_mono_mass", "flash_precursor_isotope",
       "flash_precursor_cosine", "flash_precursor_snr", "flash_precursor_mono_observed"})
  {
    TEST_FALSE(processed[0].metaValueExists(key))
  }
  stale_data.clear();
  OMSFile().load(stale_oms, stale_data);
  TEST_FALSE(stale_data.getObservationMatches().empty())
  ABORT_IF(stale_data.getObservationMatches().empty())
  for (const auto& match : stale_data.getObservationMatches())
  {
    TEST_EQUAL(Int(match.getMetaValue("isotope_offset")), 1)
    TEST_EQUAL(match.getMetaValue("flash_precursor_status").toString(), "no_ms1")
    for (const auto& key : {"flash_precursor_mz", "flash_precursor_mono_mass", "flash_precursor_isotope",
         "flash_precursor_cosine", "flash_precursor_snr", "flash_precursor_mono_observed"})
    {
      TEST_FALSE(match.metaValueExists(key))
    }
  }

  fixture("ACGUAC", "AC[m2,2G]UAC", 2);
  args = arguments("RNase_T1");
  args.insert(args.end(), {"-modifications:fixed", "m22G"});
  TEST_EQUAL(runSearch(args), 0)
  TEST_EQUAL(hitCount(output), 1)

  fixture("ACUACG", "ACUAC[m2,2G]", 2);
  args = arguments("RNase_T1");
  args.insert(args.end(), {"-modifications:variable", "m22G"});
  TEST_EQUAL(runSearch(args), 0)
  TEST_EQUAL(hitCount(output), 1)

  fixture("ACGUAC", "pACGUAC", 2);
  args = arguments();
  args.insert(args.end(), {"-modifications:variable", "5'-p", "-modifications:variable_max_per_oligo", "1"});
  TEST_EQUAL(runSearch(args), 0)
  TEST_EQUAL(hitCount(output), 1)
}
END_SECTION

START_SECTION(([EXTRA] imported digests use only their recorded enzyme for modification inference and budgets))
{
  std::string mzml, digest, output;
  NEW_TMP_FILE(mzml)
  NEW_TMP_FILE(digest)
  NEW_TMP_FILE(output)
  mzml += ".mzML";
  digest += ".oms";
  output += ".idXML";
  NucleicAcidSpectrumGenerator generator;
  Param parameters = generator.getParameters();
  parameters.setValue("add_first_prefix_ion", "true");
  parameters.setValue("add_a-B_ions", "true");
  parameters.setValue("add_metainfo", "true");
  generator.setParameters(parameters);
  const auto fixture = [&](const std::string& oligo_sequence, const std::string& hit,
                           const DigestionEnzyme* enzyme)
  {
    IdentificationData data;
    data.registerScoreType(IdentificationData::ScoreType("hyperscore", true));
    IdentificationData::DBSearchParam search_param;
    search_param.molecule_type = IdentificationData::MoleculeType::RNA;
    search_param.digestion_enzyme = enzyme;
    const auto search_ref = data.registerDBSearchParam(search_param);
    const auto software_ref = data.registerProcessingSoftware(
      IdentificationData::ProcessingSoftware("RNADigestor", "test"));
    const auto step_ref = data.registerProcessingStep(
      IdentificationData::ProcessingStep(software_ref), search_ref);
    data.setCurrentProcessingStep(step_ref);
    const auto parent_ref = data.registerParentSequence(
      IdentificationData::ParentSequence("TARGET", IdentificationData::MoleculeType::RNA, "ACGUAC"));
    IdentificationData::IdentifiedOligo oligo(NASequence::fromString(oligo_sequence));
    oligo.parent_matches[parent_ref].insert(IdentificationData::ParentMatch(0, oligo.sequence.size() - 1));
    data.registerIdentifiedOligo(oligo);
    OMSFile().store(digest, data);

    // Exercise the serialized metadata, including a genuinely absent enzyme.
    IdentificationData imported;
    OMSFile().load(digest, imported);
    const auto imported_enzyme = imported.getDBSearchParams().begin()->digestion_enzyme;
    TEST_TRUE((imported_enzyme == nullptr) == (enzyme == nullptr))
    if (enzyme && imported_enzyme)
    {
      TEST_EQUAL(imported_enzyme->getName(), enzyme->getName())
    }

    const auto sequence = NASequence::fromString(hit);
    MSSpectrum spectrum;
    generator.getSpectrum(spectrum, sequence, -1, -2);
    spectrum.setMSLevel(2);
    spectrum.setNativeID("scan=1");
    spectrum.setRT(20);
    Precursor precursor;
    precursor.setMZ(sequence.getMonoWeight(NASequence::Full, -2) / 2);
    precursor.setCharge(2);
    spectrum.setPrecursors({precursor});
    MSExperiment exp;
    exp.addSpectrum(spectrum);
    FileHandler().storeExperiment(mzml, exp, {FileTypes::MZML});
  };
  const auto arguments = [&](const std::string& cli_enzyme, const std::string& budget)
  {
    return std::vector<std::string>{"-in", mzml, "-digest", digest, "-id_out", output, "-out", output + ".mzTab",
      "-oligo:enzyme", cli_enzyme, "-threads", "1", "-oligo:min_size", "1",
      "-modifications:variable", "m22G", "5'-p", "-modifications:variable_max_per_oligo", budget};
  };
  const auto hit_labels = [&]()
  {
    std::vector<ProteinIdentification> proteins;
    PeptideIdentificationList peptides;
    IdXMLFile().load(output, proteins, peptides);
    StringList labels;
    for (const auto& peptide : peptides)
    {
      for (const auto& hit : peptide.getHits())
      {
        labels.push_back(hit.getMetaValue("label").toString());
      }
    }
    return labels;
  };

  // An unknown digest enzyme must not suppress a modification inferred from the CLI enzyme.
  fixture("ACGUAC", "AC[m2,2G]UAC", nullptr);
  TEST_EQUAL(runSearch(arguments("RNase_T1", "1")), 0)
  TEST_TRUE(hit_labels() == StringList{"AC[m2,2G]UAC"})

  // Existing modified residues do not consume the budget without enzyme metadata.
  fixture("AC[m2,2G]UAC", "pAC[m2,2G]UAC", nullptr);
  TEST_EQUAL(runSearch(arguments("RNase_T1", "1")), 0)
  TEST_TRUE(hit_labels() == StringList{"pAC[m2,2G]UAC"})

  // A recorded enzyme still governs the budget even when the CLI enzyme differs.
  fixture("AC[m2,2G]UAC", "pAC[m2,2G]UAC", RNaseDB::getInstance()->getEnzyme("RNase_T1"));
  TEST_EQUAL(runSearch(arguments("no cleavage", "1")), 0)
  TEST_EQUAL(hit_labels().size(), 0)
  TEST_EQUAL(runSearch(arguments("no cleavage", "2")), 0)
  TEST_TRUE(hit_labels() == StringList{"pAC[m2,2G]UAC"})

  fixture("ACGUAC", "AC[m2,2G]UAC", RNaseDB::getInstance()->getEnzyme("no cleavage"));
  TEST_EQUAL(runSearch(arguments("RNase_T1", "1")), 0)
  TEST_TRUE(hit_labels() == StringList{"AC[m2,2G]UAC"})
}
END_SECTION

START_SECTION(([EXTRA] complete cutoff ties are retained for every arrival order))
{
  TestNASE tool;
  IdentificationData data;
  auto oligo = data.registerIdentifiedOligo(IdentificationData::IdentifiedOligo(NASequence::fromString("ACGU")));
  TestNASE::PrecursorInfo precursor{0, -2, 0, nullopt};
  std::vector<std::string> sequences{"AAA", "AAC", "AAG", "AAU"};
  std::vector<double> scores{12.0, 10.0, 10.0, 9.0};
  std::vector<Size> order{0, 1, 2, 3};
  do
  {
    TestNASE::HitsByScore hits;
    for (Size index : order)
    {
      TestNASE::AnnotatedHit hit;
      hit.oligo_ref = oligo;
      hit.sequence = NASequence::fromString(sequences[index]);
      hit.precursor_error_ppm = 0;
      hit.precursor_ref = &precursor;
      tool.retainTopHit_(hits, scores[index], hit, 2);
    }
    tool.sortHitsByValue_(hits);
    TEST_EQUAL(hits.size(), 3)
    std::vector<std::string> retained;
    for (const auto& hit : hits)
    {
      retained.push_back(hit.second.sequence.toString());
    }
    TEST_TRUE(retained == std::vector<std::string>({"AAA", "AAC", "AAG"}))
  }
  while (std::next_permutation(order.begin(), order.end()));
}
END_SECTION

START_SECTION(([EXTRA] ambiguity variants merge within each precursor across tied-hit permutations))
{
  TestNASE tool;
  tool.enableAmbiguities();
  TEST_EQUAL(tool.getModifications_({"mA?"}).size(), 2)
  IdentificationData data;
  const auto oligo = data.registerIdentifiedOligo(IdentificationData::IdentifiedOligo(NASequence::fromString("AUCAGp")));
  const auto sodium = data.registerAdduct(AdductInfo("Na+", EmpiricalFormula("Na"), 1));
  TestNASE::PrecursorInfo first{0, -2, 0, nullopt};
  TestNASE::PrecursorInfo second{0, -3, 1, sodium};
  const std::vector<std::string> sequences{"AUC[Am]Gp", "AUC[Am]Gp", "AUC[mA]Gp", "AUC[mA]Gp", "AUC[m6A]Gp"};
  const std::vector<const TestNASE::PrecursorInfo*> precursors{&first, &second, &first, &second, &first};
  std::vector<Size> order{0, 1, 2, 3, 4};
  do
  {
    TestNASE::HitsByScore hits;
    for (const Size index : order)
    {
      TestNASE::AnnotatedHit hit;
      hit.oligo_ref = oligo;
      hit.sequence = NASequence::fromString(sequences[index]);
      hit.precursor_ref = precursors[index];
      hit.precursor_error_ppm = index;
      tool.retainTopHit_(hits, 10.0, hit, 1);
    }
    // Before the fix, value sorting interleaves the two precursor groups and
    // places the unrelated m6A hit between the variants of the first precursor.
    tool.sortHitsByValue_(hits);
    tool.resolveAmbiguousMods_(hits);
    tool.sortHitsByValue_(hits);
    TEST_EQUAL(hits.size(), 3)
    std::vector<std::string> retained;
    Size ambiguous_count = 0;
    for (const auto& hit : hits)
    {
      retained.push_back(hit.second.sequence.toString());
      TEST_EQUAL(hit.first, 10.0)
      TEST_TRUE(hit.second.oligo_ref == oligo)
      if (hit.second.sequence.toString() == "AUC[mA?]Gp")
      {
        ++ambiguous_count;
        TEST_TRUE(hit.second.precursor_ref == &first || hit.second.precursor_ref == &second)
        // Preserve the deterministically first specific variant's metadata.
        TEST_REAL_SIMILAR(hit.second.precursor_error_ppm, hit.second.precursor_ref == &first ? 0.0 : 1.0)
      }
      else
      {
        TEST_TRUE(hit.second.precursor_ref == &first)
        TEST_REAL_SIMILAR(hit.second.precursor_error_ppm, 4.0)
      }
    }
    TEST_EQUAL(ambiguous_count, 2)
    TEST_TRUE(retained == std::vector<std::string>({"AUC[m6A]Gp", "AUC[mA?]Gp", "AUC[mA?]Gp"}))
  }
  while (std::next_permutation(order.begin(), order.end()));
}
END_SECTION

START_SECTION(([EXTRA] ambiguity resolution preserves score source terminal and precursor identities))
{
  TestNASE tool;
  tool.enableAmbiguities();
  tool.getModifications_({"mA?"});
  IdentificationData data;
  const auto oligo = data.registerIdentifiedOligo(IdentificationData::IdentifiedOligo(NASequence::fromString("AUCAGp")));
  const auto other_oligo = data.registerIdentifiedOligo(IdentificationData::IdentifiedOligo(NASequence::fromString("AAUAGp")));
  const auto sodium = data.registerAdduct(AdductInfo("Na+", EmpiricalFormula("Na"), 1));
  TestNASE::PrecursorInfo precursor{0, -2, 0, nullopt};
  TestNASE::PrecursorInfo equal_values{0, -2, 0, nullopt};
  TestNASE::PrecursorInfo other_charge{0, -3, 0, nullopt};
  TestNASE::PrecursorInfo other_isotope{0, -2, 1, nullopt};
  TestNASE::PrecursorInfo other_adduct{0, -2, 0, sodium};
  TestNASE::AnnotatedHit first, second;
  first.oligo_ref = second.oligo_ref = oligo;
  first.precursor_ref = second.precursor_ref = &precursor;
  first.sequence = NASequence::fromString("AUC[Am]Gp");
  second.sequence = NASequence::fromString("AUC[mA]Gp");
  const auto check_distinct = [&](const TestNASE::AnnotatedHit& other, double other_score = 10.0)
  {
    TestNASE::HitsByScore hits;
    hits.emplace(10.0, first);
    hits.emplace(other_score, other);
    tool.sortHitsByValue_(hits);
    tool.resolveAmbiguousMods_(hits);
    TEST_EQUAL(hits.size(), 2)
    for (const auto& hit : hits)
    {
      TEST_TRUE(hit.second.sequence == first.sequence || hit.second.sequence == other.sequence)
    }
  };
  check_distinct(second, 9.0);
  auto different = second;
  different.oligo_ref = other_oligo;
  check_distinct(different);
  different = second;
  different.sequence = NASequence::fromString("pAUC[mA]Gp");
  check_distinct(different);
  different.sequence = NASequence::fromString("AUC[mA]G");
  check_distinct(different);
  different.sequence = NASequence::fromString("AUCC[mA]Gp");
  check_distinct(different);
  for (const auto* other_precursor : {&equal_values, &other_charge, &other_isotope, &other_adduct})
  {
    different = second;
    different.precursor_ref = other_precursor;
    check_distinct(different);
  }
}
END_SECTION

START_SECTION(([EXTRA] LFQ uses finite weights or the median of finite observations))
{
  TestNASE tool;
  IdentificationData data;
  auto input = data.registerInputFile(IdentificationData::InputFile("lfq.mzML"));
  Size index = 0;
  const auto add = [&](const std::string& sequence, double rt, double weight, Int charge = -2)
  {
    auto oligo = data.registerIdentifiedOligo(IdentificationData::IdentifiedOligo(NASequence::fromString(sequence)));
    IdentificationData::Observation obs("scan=" + std::to_string(++index), input, rt);
    obs.setMetaValue("precursor_intensity", weight);
    auto observation = data.registerObservation(obs);
    data.registerObservationMatch(IdentificationData::ObservationMatch(oligo, observation, charge));
  };
  const double nan = std::numeric_limits<double>::quiet_NaN();
  add("ACGU", 10, 0);
  add("ACGU", 20, -1);
  add("ACGU", 30, nan);
  add("ACGU", nan, 100);
  add("ACGU", nan, 100, -3);
  add("AGGU", 10, 1);
  add("AGGU", 30, 3);
  add("AGGU", 100, nan);
  add("AGGU", nan, 100);
  add("AUUU", nan, 1);
  std::string path;
  NEW_TMP_FILE(path)
  tool.generateLFQInput_(data, path);
  TextFile file(path);
  const std::vector<std::string> table(file.begin(), file.end());
  TEST_EQUAL(table.size(), 3)
  std::map<std::string, double> rts;
  for (Size row = 1; row < table.size(); ++row)
  {
    std::vector<std::string> fields;
    StringUtils::split(table[row], '\t', fields);
    rts[fields[0]] = StringUtils::toDouble(fields[4]);
    TEST_EQUAL(fields[3], "-2")
  }
  TEST_REAL_SIMILAR(rts.at("ACGU"), 20)
  TEST_REAL_SIMILAR(rts.at("AGGU"), 25)
  TEST_EQUAL(rts.count("AUUU"), 0)
}
END_SECTION

START_SECTION(([EXTRA] inferred and signed charges, array-aware removal before decharging))
{
  TestNASE tool;
  for (Int charge : {2, -2, 0})
  {
    MSExperiment exp;
    auto spectrum = isotopeSpectrum(charge);
    if (charge == 0)
    {
      Precursor second;
      second.setMZ(1000 + Constants::PROTON_MASS_U);
      spectrum.getPrecursors().push_back(second);
    }
    exp.addSpectrum(spectrum);
    tool.preprocess(exp, false, "legacy");
    TEST_EQUAL(exp[0].size(), 1)
    TEST_REAL_SIMILAR(exp[0][0].getMZ(), 1000 + Constants::PROTON_MASS_U)
  }
  for (const std::string method : {"none", "legacy"})
  {
    MSExperiment exp;
    auto spectrum = isotopeSpectrum(2);
    MSSpectrum::IntegerDataArray array;
    array.setName("peak_index");
    array.assign({0, 1, 2});
    spectrum.getIntegerDataArrays().push_back(array);
    exp.addSpectrum(spectrum);
    tool.preprocess(exp, true, method);
    TEST_EQUAL(exp[0].size(), 0)
    TEST_EQUAL(exp[0].getIntegerDataArrays()[0].size(), 0)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(([EXTRA] theoretical blacklist preserves arrays and coverage counts distinct internal positions))
{
  TestNASE tool;
  MSSpectrum spectrum;
  spectrum.push_back(Peak1D(100.0, 1.0));
  spectrum.push_back(Peak1D(200.0, 2.0));
  MSSpectrum::FloatDataArray floats;
  floats.assign({10.0, 20.0});
  spectrum.getFloatDataArrays().push_back(floats);
  MSSpectrum::StringDataArray names;
  names.assign({"c1", "y1"});
  spectrum.getStringDataArrays().push_back(names);
  MSSpectrum::IntegerDataArray charges;
  charges.assign({1, 2});
  spectrum.getIntegerDataArrays().push_back(charges);
  tool.filterBlacklistedIons_(spectrum, {100.0005}, 10.0);
  TEST_EQUAL(spectrum.size(), 1)
  TEST_EQUAL(spectrum.getFloatDataArrays()[0].size(), 1)
  TEST_REAL_SIMILAR(spectrum.getFloatDataArrays()[0][0], 20.0)
  TEST_STRING_EQUAL(spectrum.getStringDataArrays()[0][0], "y1")
  TEST_EQUAL(spectrum.getIntegerDataArrays()[0][0], 2)

  vector<PeptideHit::PeakAnnotation> annotations(4);
  annotations[0].annotation = "c1";
  annotations[1].annotation = "a1-B";
  annotations[2].annotation = "y1";
  annotations[3].annotation = "c5";
  TEST_REAL_SIMILAR(tool.calculateFragmentCoverage_(NASequence::fromString("ACGU"), annotations), 2.0 / 3.0)
  TEST_REAL_SIMILAR(tool.calculateFragmentCoverage_(NASequence::fromString("A"), {}), 1.0)
}
END_SECTION

START_SECTION(([EXTRA] residual precursor percentage uses raw intensity and strongest matching peak))
{
  MSSpectrum spectrum;
  Precursor precursor;
  precursor.setMZ(500.0);
  spectrum.setPrecursors({precursor});
  spectrum.push_back(Peak1D(100.0, 1000.0));
  spectrum.push_back(Peak1D(499.999, 100.0));
  spectrum.push_back(Peak1D(500.001, 250.0));
  spectrum.push_back(Peak1D(500.0 + Constants::C13C12_MASSDIFF_U / 2, 900.0));
  spectrum.sortByPosition(); // a narrow window does not include the isotope peak
  const auto percentage = TestNASE::calculatePrecursorPercentage_(spectrum, 10.0, true);
  TEST_TRUE(percentage.has_value())
  ABORT_IF(!percentage.has_value())
  TEST_REAL_SIMILAR(*percentage, 25.0)
  TEST_REAL_SIMILAR(*TestNASE::calculatePrecursorPercentage_(spectrum, 1.0, false), 90.0)
  TEST_REAL_SIMILAR(*TestNASE::calculatePrecursorPercentage_(spectrum, 0.0, false), 0.0)
  spectrum.setMetaValue("precursor_percentage", 999.0);
  spectrum.setMetaValue("unrelated_annotation", "retained");
  TestNASE::updatePrecursorPercentage_(spectrum, false, 10.0, true);
  TEST_FALSE(spectrum.metaValueExists("precursor_percentage"))
  TEST_EQUAL(spectrum.getMetaValue("unrelated_annotation"), "retained")
  TestNASE::updatePrecursorPercentage_(spectrum, true, 10.0, true);
  TEST_REAL_SIMILAR(double(spectrum.getMetaValue("precursor_percentage")), 25.0)
  spectrum.setPrecursors({precursor, precursor});
  TEST_FALSE(TestNASE::calculatePrecursorPercentage_(spectrum, 0.1, false).has_value())
  spectrum.setMetaValue("precursor_percentage", 999.0);
  TestNASE::updatePrecursorPercentage_(spectrum, true, 10.0, true);
  TEST_FALSE(spectrum.metaValueExists("precursor_percentage"))
  TEST_EQUAL(spectrum.getMetaValue("unrelated_annotation"), "retained")
  spectrum.setPrecursors({});
  TEST_FALSE(TestNASE::calculatePrecursorPercentage_(spectrum, 0.1, false).has_value())
}
END_SECTION

START_SECTION(([EXTRA] FLASH MS2: void filter(MSSpectrum& spectrum, int scan_number, bool decharge)))
{
  RNADeisotoper detector(10, 5, false, 8000);
  MSSpectrum empty;
  empty.setMSLevel(2);
  detector.filter(empty, 0, false);
  TEST_TRUE(empty.empty())

  auto original = rnaFragmentPattern(2);
  auto filtered = original;
  detector.filter(filtered, 17, false);
  TEST_EQUAL(filtered.getPrecursors()[0].getMZ(), original.getPrecursors()[0].getMZ())
  TEST_EQUAL(filtered.getPrecursors()[0].getCharge(), -5)
  TEST_EQUAL(filtered.getNativeID(), original.getNativeID())
  TEST_EQUAL(filtered.getRT(), original.getRT())
  TEST_TRUE(Int(filtered.getMetaValue("flash_isotope_peaks_removed")) > 0)
  TEST_EQUAL(filtered.getIntegerDataArrays()[0].size(), filtered.size())
  TEST_EQUAL(filtered.getIntegerDataArrays()[1].size(), filtered.size())
  for (Size i = 0; i < filtered.size(); ++i)
  {
    const auto index = filtered.getIntegerDataArrays()[0][i];
    TEST_REAL_SIMILAR(filtered[i].getMZ(), original[index].getMZ())
    TEST_REAL_SIMILAR(filtered[i].getIntensity(), original[index].getIntensity())
  }
  TEST_EQUAL(filtered.back().getMZ(), 1799.123)
  TEST_EQUAL(filtered.back().getIntensity(), 11.0)

  auto single = original;
  detector.filter(single, 17, true);
  const double mono_mass = NASequence::fromString("ACGUACGU").getFormula().getMonoWeight();
  bool found_single = false;
  for (const auto& peak : single)
    if (abs(peak.getMZ() - (mono_mass - Constants::PROTON_MASS_U)) < 1e-6) found_single = true;
  TEST_TRUE(found_single)
  TEST_EQUAL(single.getPrecursors()[0].getMZ(), 1500.0)

  auto unknown = original;
  unknown.getPrecursors()[0].setCharge(0);
  detector.filter(unknown, 17, false);
  TEST_TRUE(Int(unknown.getMetaValue("flash_isotope_peaks_removed")) > 0)

  // The real z=2 envelope must veto the z=1 harmonic explanation even when
  // the precursor permits only z=1 assignments. All original peaks survive.
  auto bounded = original;
  bounded.getPrecursors()[0].setCharge(-1);
  detector.filter(bounded, 17, false);
  TEST_EQUAL(bounded.size(), original.size())
  TEST_EQUAL(Int(bounded.getMetaValue("flash_isotope_peaks_removed")), 0)
  const auto& bounded_charges = bounded.getIntegerDataArrays().back();
  TEST_EQUAL(bounded_charges.getName(), "flash_charge")
  for (Size i = 0; i < bounded.size(); ++i)
  {
    TEST_EQUAL(bounded_charges[i], 0)
    TEST_EQUAL(bounded[i], original[i])
  }

  auto isolated = original;
  isolated.resize(1);
  isolated.getIntegerDataArrays()[0].resize(1);
  detector.filter(isolated, 17, false);
  TEST_EQUAL(isolated.size(), 1)
  TEST_EQUAL(isolated[0], original[0])

  // Missing M0 can be misindexed by the averagine model; never invent an observed peak.
  auto missing_mono = original;
  vector<Size> keep;
  for (Size i = 0; i < original.size(); ++i)
    if (i != 1) keep.push_back(i); // index 0 is the unassigned low-mass signal
  missing_mono.select(keep);
  const Size missing_size = missing_mono.size();
  detector.filter(missing_mono, 17, false);
  TEST_TRUE(missing_mono.size() <= missing_size)
  for (const auto& peak : missing_mono)
  {
    bool observed = false;
    for (const auto& old_peak : original)
      if (peak == old_peak) observed = true;
    TEST_TRUE(observed)
  }

  auto profile = original;
  profile.setType(SpectrumSettings::SpectrumType::PROFILE);
  detector.filter(profile, 17, false);
  TEST_EQUAL(profile.size(), filtered.size())
}
END_SECTION

START_SECTION(([EXTRA] FLASH MS2: RNADeisotoper(double tolerance_ppm, int max_abs_charge, bool positive, double max_mass, double min_cos = 0.85, double min_snr = 0.25)))
{
  RNADeisotoper positive(10, 5, true, 8000);
  auto s = rnaFragmentPattern(3, true);
  const Size size = s.size();
  positive.filter(s, 18, false);
  TEST_TRUE(s.size() < size)
  TEST_EQUAL(s.getPrecursors()[0].getCharge(), 5)

  // A supplied upper limit must not be silently raised to the default.
  MSSpectrum low_mass;
  low_mass.setMSLevel(2);
  low_mass.setType(SpectrumSettings::SpectrumType::CENTROID);
  const auto formula = NASequence::fromString("ACG").getFormula();
  const auto isotopes = formula.getIsotopeDistribution(CoarseIsotopePatternGenerator(8));
  for (Size i = 0; i < isotopes.size(); ++i)
  {
    low_mass.emplace_back((formula.getMonoWeight() + i * Constants::C13C12_MASSDIFF_U) / 2 + Constants::PROTON_MASS_U,
                          isotopes[i].getIntensity() * 10000);
  }
  RNADeisotoper limited(10, 5, true, 700);
  auto restricted = low_mass;
  limited.filter(restricted, 19, false);
  TEST_EQUAL(restricted.size(), low_mass.size())
  for (Size i = 0; i < low_mass.size(); ++i) TEST_EQUAL(restricted[i], low_mass[i])
  RNADeisotoper default_limit(10, 5, true);
  default_limit.filter(low_mass, 19, false);
  TEST_TRUE(low_mass.size() < restricted.size())

  TEST_EXCEPTION(Exception::IllegalArgument, RNADeisotoper(0, 2, true))
  TEST_EXCEPTION(Exception::IllegalArgument, RNADeisotoper(-10, 2, true))
  TEST_EXCEPTION(Exception::IllegalArgument, RNADeisotoper(std::numeric_limits<double>::quiet_NaN(), 2, true))
  TEST_EXCEPTION(Exception::IllegalArgument, RNADeisotoper(10, 0, true))
  TEST_EXCEPTION(Exception::IllegalArgument, RNADeisotoper(10, -1, true))
  for (double invalid : {std::numeric_limits<double>::quiet_NaN(),
                         std::numeric_limits<double>::infinity(),
                         -std::numeric_limits<double>::infinity()})
  {
    TEST_EXCEPTION(Exception::IllegalArgument, RNADeisotoper(10, 2, true, 1000, invalid, 0.25))
    TEST_EXCEPTION(Exception::IllegalArgument, RNADeisotoper(10, 2, true, 1000, 0.85, invalid))
  }
  // Finite out-of-range values retain the existing FLASH parameter contract.
  TEST_EXCEPTION(Exception::InvalidParameter, RNADeisotoper(10, 2, true, 1000, 0.49, 0.25))
  TEST_EXCEPTION(Exception::InvalidParameter, RNADeisotoper(10, 2, true, 1000, 1.01, 0.25))
  TEST_EXCEPTION(Exception::InvalidParameter, RNADeisotoper(10, 2, true, 1000, 0.85, 0.24))
  TEST_EXCEPTION(Exception::IllegalArgument, RNADeisotoper(10, 2, true, 49))
  TEST_EXCEPTION(Exception::IllegalArgument, RNADeisotoper(10, 2, true, std::numeric_limits<double>::quiet_NaN()))
  TEST_EXCEPTION(Exception::IllegalArgument, RNADeisotoper(10, 2, true, std::numeric_limits<double>::infinity()))
}
END_SECTION
START_SECTION(([EXTRA] double flashMassUpperBound_(double mz, Int charge, bool negative_mode) const))
{
  TestNASE tool;
  for (bool negative : {false, true})
  {
    const double carrier = negative ? -Constants::PROTON_MASS_U : Constants::PROTON_MASS_U;
    TEST_REAL_SIMILAR(tool.flashMassUpperBound_(1000.0, 3, negative), (1000.0 - carrier) * 3 + 100.0)
    TEST_EQUAL(tool.flashMassUpperBound_(1000.0, -3, negative), tool.flashMassUpperBound_(1000.0, 3, negative))
  }
}
END_SECTION
START_SECTION(([EXTRA] FLASH MS2: preprocessing retains raw-scale decisions with duplicate or empty native IDs))
{
  TestNASE tool;
  auto first = rnaFragmentPattern(2);
  auto second = first;
  for (auto& peak : second) peak.setIntensity(peak.getIntensity() * 100.0);
  RNADeisotoper detector(10, 5, false, 8000);
  auto expected_first = first, expected_second = second;
  const auto normalize = [](MSSpectrum& spectrum, const MSSpectrum& original)
  {
    double scale = 0.0;
    for (const auto& peak : original) scale = std::max(scale, double(peak.getIntensity()));
    for (auto& peak : spectrum) peak.setIntensity(peak.getIntensity() / scale);
  };
  detector.filter(expected_first, 0, true);
  detector.filter(expected_second, 1, true);
  TEST_TRUE(expected_first.size() < first.size())
  TEST_TRUE(expected_second.size() < second.size())
  normalize(expected_first, first);
  normalize(expected_second, second);
  for (const std::string& id : {std::string(), std::string("duplicate")})
  {
    first.setNativeID(id);
    second.setNativeID(id);
    MSExperiment experiment;
    experiment.setSpectra({first, second});
    tool.preprocess(experiment, false, "flash");
    TEST_EQUAL(experiment[0].size(), expected_first.size())
    TEST_EQUAL(experiment[1].size(), expected_second.size())
    ABORT_IF(experiment[0].size() != expected_first.size() || experiment[1].size() != expected_second.size())
    for (Size i = 0; i < expected_first.size(); ++i) TEST_EQUAL(experiment[0][i], expected_first[i])
    for (Size i = 0; i < expected_second.size(); ++i) TEST_EQUAL(experiment[1][i], expected_second[i])
    TEST_TRUE(experiment[0].getIntegerDataArrays() == expected_first.getIntegerDataArrays())
    TEST_TRUE(experiment[1].getIntegerDataArrays() == expected_second.getIntegerDataArrays())
  }
}
END_SECTION
START_SECTION(([EXTRA] FLASH MS2: shared monoisotopes and harmonic RNA fragments))
{
  MSSpectrum s;
  s.setMSLevel(2);
  s.setType(SpectrumSettings::SpectrumType::CENTROID);
  Precursor p;
  p.setMZ(1000.0);
  p.setCharge(-5);
  s.getPrecursors().push_back(p);
  vector<Peak1D> peaks;
  for (int z = 1; z <= 3; ++z)
  {
    const auto formula = NASequence::fromString(std::string(z, 'C')).getFormula(NASequence::CIon);
    const auto iso = formula.getIsotopeDistribution(CoarseIsotopePatternGenerator(8));
    for (Size i = 0; i < iso.size(); ++i)
      peaks.emplace_back((formula.getMonoWeight() + i * Constants::C13C12_MASSDIFF_U) / z - Constants::PROTON_MASS_U,
                         iso[i].getIntensity() * 10000);
  }
  sort(peaks.begin(), peaks.end(), [](const auto& a, const auto& b) { return a.getMZ() < b.getMZ(); });
  for (const auto& peak : peaks)
  {
    if (!s.empty() && abs(s.back().getMZ() - peak.getMZ()) < 1e-8)
      s.back().setIntensity(s.back().getIntensity() + peak.getIntensity());
    else s.push_back(peak);
  }
  const auto shared = s.front();
  s.emplace_back(83.123, 9.0);
  s.emplace_back(1799.123, 11.0);
  s.sortByPosition();
  const Size before = s.size();
  RNADeisotoper detector(10, 5, false, 6000);
  detector.filter(s, 19, false);
  int shared_count = 0;
  for (const auto& peak : s)
  {
    if (abs(peak.getMZ() - shared.getMZ()) < 1e-8)
    {
      ++shared_count;
      TEST_REAL_SIMILAR(peak.getIntensity(), shared.getIntensity())
    }
  }
  TEST_EQUAL(shared_count, 1)
  TEST_TRUE(s.size() <= before)
  TEST_EQUAL(s.front().getMZ(), 83.123)
  STATUS("Harmonic mixture: " << before << " -> " << s.size() << " peaks; detected groups " << s.getMetaValue("flash_detected_peak_groups"))
}
END_SECTION
START_SECTION(([EXTRA] FLASH MS2: profile samples belonging to monoisotopes are retained individually))
{
  auto sticks = rnaFragmentPattern(2);
  MSSpectrum profile = sticks;
  profile.clear(false);
  profile.getIntegerDataArrays().clear();
  profile.setType(SpectrumSettings::SpectrumType::PROFILE);
  // Sample each isotope at five points, with a symmetric intensity profile.
  const vector<double> offsets{-0.001, -0.0005, 0.0, 0.0005, 0.001};
  const vector<double> weights{0.2, 0.6, 1.0, 0.6, 0.2};
  for (Size i = 0; i < sticks.size(); ++i)
  {
    if (i == 0 || i + 1 == sticks.size()) profile.push_back(sticks[i]);
    else
    {
      for (Size j = 0; j < offsets.size(); ++j)
        profile.emplace_back(sticks[i].getMZ() + offsets[j], sticks[i].getIntensity() * weights[j]);
    }
  }
  profile.sortByPosition();
  MSSpectrum::IntegerDataArray ids;
  ids.setName("original_index");
  for (Size i = 0; i < profile.size(); ++i) ids.push_back(i);
  profile.getIntegerDataArrays().push_back(ids);
  const MSSpectrum original = profile;
  RNADeisotoper detector(10, 5, false, 8000);
  detector.filter(profile, 20, false);
  TEST_TRUE(profile.size() < original.size())
  TEST_TRUE(profile.getType() == SpectrumSettings::SpectrumType::PROFILE)
  TEST_EQUAL(profile.getNativeID(), original.getNativeID())
  TEST_EQUAL(profile.getPrecursors()[0].getCharge(), -5)
  Size mono_samples = 0;
  for (Size i = 0; i < profile.size(); ++i)
  {
    const int original_index = profile.getIntegerDataArrays()[0][i];
    TEST_EQUAL(profile[i], original[original_index])
    if (abs(profile[i].getMZ() - sticks[1].getMZ()) < 0.002)
    {
      ++mono_samples;
      TEST_EQUAL(profile.getIntegerDataArrays()[1][i], -2)
    }
  }
  TEST_EQUAL(mono_samples, 5)
  TEST_EQUAL(profile.front(), original.front())
  TEST_EQUAL(profile.back(), original.back())
  TEST_TRUE(Int(profile.getMetaValue("flash_isotope_peaks_removed")) > 0)
}
END_SECTION
START_SECTION(([EXTRA] FLASH MS1: RNAPrecursorCorrector(double tolerance_ppm, int min_abs_charge, int max_abs_charge, bool positive, double max_mass = 1000.0, double min_cos = 0.85, double min_snr = 0.25)))
{
  const auto parent = survey(true, "ACG");
  MSExperiment original;
  original.setSpectra({parent, ms2(parent)});
  RNAPrecursorCorrector limited(10, 1, 5, true, 700);
  auto restricted = original;
  const auto limited_rows = limited.correct(restricted);
  TEST_EQUAL(limited_rows.size(), 1)
  TEST_FALSE(limited_rows[0].accepted())
  TEST_EQUAL(restricted[1].getPrecursors()[0].getMZ(), original[1].getPrecursors()[0].getMZ())
  RNAPrecursorCorrector default_limit(10, 1, 5, true);
  auto unrestricted = original;
  const auto default_rows = default_limit.correct(unrestricted);
  TEST_EQUAL(default_rows.size(), 1)
  TEST_TRUE(default_rows[0].accepted())

  TEST_EXCEPTION(Exception::IllegalArgument, RNAPrecursorCorrector(0, 1, 2, true))
  TEST_EXCEPTION(Exception::IllegalArgument, RNAPrecursorCorrector(-10, 1, 2, true))
  TEST_EXCEPTION(Exception::IllegalArgument, RNAPrecursorCorrector(std::numeric_limits<double>::infinity(), 1, 2, true))
  TEST_EXCEPTION(Exception::IllegalArgument, RNAPrecursorCorrector(10, 0, 2, true))
  TEST_EXCEPTION(Exception::IllegalArgument, RNAPrecursorCorrector(10, 3, 2, true))
  for (double invalid : {std::numeric_limits<double>::quiet_NaN(),
                         std::numeric_limits<double>::infinity(),
                         -std::numeric_limits<double>::infinity()})
  {
    TEST_EXCEPTION(Exception::IllegalArgument, RNAPrecursorCorrector(10, 1, 2, true, 1000, invalid, 0.25))
    TEST_EXCEPTION(Exception::IllegalArgument, RNAPrecursorCorrector(10, 1, 2, true, 1000, 0.85, invalid))
  }
  TEST_EXCEPTION(Exception::InvalidParameter, RNAPrecursorCorrector(10, 1, 2, true, 1000, 0.49, 0.25))
  TEST_EXCEPTION(Exception::InvalidParameter, RNAPrecursorCorrector(10, 1, 2, true, 1000, 1.01, 0.25))
  TEST_EXCEPTION(Exception::InvalidParameter, RNAPrecursorCorrector(10, 1, 2, true, 1000, 0.85, 0.24))
  TEST_EXCEPTION(Exception::IllegalArgument, RNAPrecursorCorrector(10, 1, 2, true, 49))
  TEST_EXCEPTION(Exception::IllegalArgument, RNAPrecursorCorrector(10, 1, 2, true, std::numeric_limits<double>::quiet_NaN()))
  TEST_EXCEPTION(Exception::IllegalArgument, RNAPrecursorCorrector(10, 1, 2, true, std::numeric_limits<double>::infinity()))
}
END_SECTION

START_SECTION(([EXTRA] FLASH MS1: Result match(const DeconvolvedSpectrum& groups, const Precursor& precursor) const))
{
  RNAPrecursorCorrector corrector(10, 1, 5, false, 8000);
  DeconvolvedSpectrum groups;
  groups.push_back(envelope());
  Precursor precursor;
  precursor.setCharge(2);
  precursor.setMZ(groups[0][6].mz);
  auto result = corrector.match(groups, precursor);
  TEST_TRUE(result.accepted())
  TEST_EQUAL(result.status, "corrected")
  TEST_EQUAL(result.isotope, 6) // beyond the previous default 0..4 search hypotheses
  TEST_REAL_SIMILAR(result.corrected_mz, 1250 - Constants::PROTON_MASS_U)
  TEST_EQUAL(result.mono_observed, true)
  precursor.setCharge(-2);
  const auto signed_result = corrector.match(groups, precursor);
  TEST_TRUE(signed_result.accepted())
  TEST_EQUAL(signed_result.original_charge, -2)
  TEST_EQUAL(signed_result.charge, -2)
  TEST_EQUAL(signed_result.corrected_mz, result.corrected_mz)
  TEST_EQUAL(signed_result.mono_mass, result.mono_mass)
  precursor.setCharge(2);

  // The aggregate envelope mass must not replace the instrument's selected-ion precision.
  precursor.setMZ(groups[0][6].mz + 0.001);
  result = corrector.match(groups, precursor);
  TEST_TRUE(result.accepted())
  TEST_EQUAL(result.isotope, 6)
  TEST_EQUAL(result.corrected_mz, precursor.getMZ() - 6 * Constants::C13C12_MASSDIFF_U / 2)
  TEST_EQUAL(result.mono_mass, (result.corrected_mz + Constants::PROTON_MASS_U) * 2)

  precursor.setMZ(groups[0][0].mz);
  result = corrector.match(groups, precursor);
  TEST_EQUAL(result.status, "monoisotopic")
  TEST_EQUAL(result.isotope, 0)

  precursor.setCharge(3);
  TEST_EQUAL(corrector.match(groups, precursor).status, "no_envelope")
  precursor.setCharge(0);
  TEST_EQUAL(corrector.match(groups, precursor).charge, 2)
  precursor.setMZ(123.456);
  TEST_EQUAL(corrector.match(groups, precursor).status, "no_envelope")
  precursor.setMZ(0);
  TEST_EQUAL(corrector.match(groups, precursor).status, "invalid_precursor")

  groups.clear();
  groups.push_back(envelope(2500, 2, false, 1));
  precursor.setMZ(groups[0][5].mz);
  precursor.setCharge(2);
  result = corrector.match(groups, precursor);
  TEST_TRUE(result.accepted())
  TEST_EQUAL(result.isotope, 6)
  TEST_EQUAL(result.mono_observed, false)
  TEST_REAL_SIMILAR(result.corrected_mz, 1250 - Constants::PROTON_MASS_U)

  groups.push_back(groups[0]);
  TEST_TRUE(corrector.match(groups, precursor).accepted()) // duplicate detection, one mass
  groups.push_back(envelope(2500 + Constants::C13C12_MASSDIFF_U));
  result = corrector.match(groups, precursor);
  TEST_EQUAL(result.status, "ambiguous")
  TEST_REAL_SIMILAR(result.corrected_mz, precursor.getMZ())

  RNAPrecursorCorrector positive(10, 1, 5, true, 8000);
  groups.clear();
  groups.push_back(envelope(2500, 2, true));
  precursor.setMZ(groups[0][2].mz);
  result = positive.match(groups, precursor);
  TEST_REAL_SIMILAR(result.corrected_mz, 1250 + Constants::PROTON_MASS_U)
}
END_SECTION

START_SECTION(([EXTRA] FLASH MS1: std::vector<Result> correct(MSExperiment& experiment)))
{
  RNAPrecursorCorrector corrector(10, 1, 5, false, 8000);
  const auto parent = survey();
  auto later = parent;
  later.setNativeID("scan=2");
  later.setRT(11);
  later.clear(false); // an explicit reference must select the earlier, populated MS1
  auto fragmentation = ms2(parent);
  MSExperiment experiment;
  experiment.addSpectrum(parent);
  experiment.addSpectrum(later);
  experiment.addSpectrum(fragmentation);
  auto rows = corrector.correct(experiment);
  TEST_EQUAL(rows.size(), 1)
  TEST_TRUE(rows[0].accepted())
  TEST_EQUAL(rows[0].parent, "scan=1")
  TEST_EQUAL(rows[0].parent_source, "spectrum_ref")
  TEST_EQUAL(rows[0].isotope, 1)
  TEST_REAL_SIMILAR(experiment[2].getPrecursors()[0].getMZ(), parent[0].getMZ())
  TEST_EQUAL(experiment[2][0], fragmentation[0])
  TEST_EQUAL(experiment[2].getPrecursors()[0].getIsolationWindowLowerOffset(), 0.5)
  TEST_REAL_SIMILAR(double(experiment[2].getPrecursors()[0].getMetaValue("isolation window target m/z")), fragmentation.getPrecursors()[0].getMZ())
  TEST_EQUAL(experiment[0], parent)

  // A previously accepted result becomes unresolved on a later run with no
  // envelope. Keep the selected-ion m/z actually present in that input.
  auto unresolved = experiment;
  unresolved[0].clear(false);
  const double rerun_mz = unresolved[2].getPrecursors()[0].getMZ();
  const auto unresolved_rows = corrector.correct(unresolved);
  TEST_FALSE(unresolved_rows[0].accepted())
  TEST_EQUAL(unresolved_rows[0].status, "no_envelope")
  TEST_EQUAL(unresolved[2].getPrecursors()[0].getMZ(), rerun_mz)
  TEST_EQUAL(double(unresolved[2].getMetaValue("flash_precursor_original_mz")), rerun_mz)
  for (const auto& key : {"flash_precursor_mz", "flash_precursor_mono_mass", "flash_precursor_isotope",
       "flash_precursor_cosine", "flash_precursor_snr", "flash_precursor_mono_observed"})
  {
    TEST_FALSE(unresolved[2].metaValueExists(key))
  }

  // Round-trip selected-ion correction without moving the acquisition window.
  std::string mzml;
  NEW_TMP_FILE(mzml)
  MzMLFile().store(mzml, experiment);
  MSExperiment reread;
  MzMLFile().load(mzml, reread);
  TEST_REAL_SIMILAR(reread[2].getPrecursors()[0].getMZ(), parent[0].getMZ())
  TEST_REAL_SIMILAR(double(reread[2].getPrecursors()[0].getMetaValue("isolation window target m/z")), fragmentation.getPrecursors()[0].getMZ())

  // Signed negative charges survive correction, audit output, and mzML I/O.
  auto signed_fragmentation = fragmentation;
  signed_fragmentation.getPrecursors()[0].setCharge(-2);
  MSExperiment signed_experiment;
  signed_experiment.setSpectra({parent, signed_fragmentation});
  const auto signed_rows = corrector.correct(signed_experiment);
  TEST_TRUE(signed_rows[0].accepted())
  TEST_EQUAL(signed_rows[0].charge, -2)
  TEST_EQUAL(signed_experiment[1].getPrecursors()[0].getCharge(), -2)
  MzMLFile().store(mzml, signed_experiment);
  MzMLFile().load(mzml, reread);
  TEST_EQUAL(reread[1].getPrecursors()[0].getCharge(), -2)
  std::string signed_tsv;
  NEW_TMP_FILE(signed_tsv)
  RNAPrecursorCorrector::storeReport(signed_tsv, signed_rows);
  TextFile signed_report(signed_tsv);
  std::vector<std::string> signed_fields;
  StringUtils::split(*std::next(signed_report.begin()), '\t', signed_fields);
  TEST_EQUAL(signed_fields[6], "-2")
  TEST_EQUAL(signed_fields[7], "-2")

  // Some converted files incorrectly refer to the MS2 itself: document fallback.
  fragmentation.getPrecursors()[0].setMetaValue("spectrum_ref", "scan=3");
  experiment.setSpectra({parent, fragmentation});
  rows = corrector.correct(experiment);
  TEST_TRUE(rows[0].accepted())
  TEST_EQUAL(rows[0].parent_source, "invalid_reference_fallback")
  TEST_EQUAL(rows[0].parent, "scan=1")

  fragmentation.getPrecursors()[0].removeMetaValue("spectrum_ref");
  experiment.setSpectra({parent, fragmentation});
  TEST_EQUAL(corrector.correct(experiment)[0].parent_source, "preceding_ms1")

  experiment.setSpectra({fragmentation});
  rows = corrector.correct(experiment);
  TEST_EQUAL(rows[0].status, "no_ms1")
  TEST_REAL_SIMILAR(experiment[0].getPrecursors()[0].getMZ(), fragmentation.getPrecursors()[0].getMZ())

  fragmentation.getPrecursors().push_back(Precursor());
  experiment.setSpectra({parent, fragmentation});
  TEST_EQUAL(corrector.correct(experiment)[0].status, "invalid_precursor")

  std::string tsv;
  NEW_TMP_FILE(tsv)
  RNAPrecursorCorrector::storeReport(tsv, rows);
  TextFile report(tsv);
  TEST_EQUAL(std::distance(report.begin(), report.end()), 2)
}
END_SECTION
START_SECTION(([EXTRA] FLASH MS1: std::vector<Result> correct(MSExperiment& experiment) preserves unsorted parent scans))
{
  for (bool positive : {false, true})
  {
    auto sorted_parent = survey(positive);
    MSSpectrum::IntegerDataArray indices;
    indices.setName("original_index");
    MSSpectrum::FloatDataArray samples;
    samples.setName("sample_value");
    MSSpectrum::StringDataArray labels;
    labels.setName("sample_label");
    for (Size i = 0; i < sorted_parent.size(); ++i)
    {
      indices.push_back(i);
      samples.push_back(float(i) + 0.5f);
      labels.push_back("sample_" + std::to_string(i));
    }
    sorted_parent.getIntegerDataArrays().push_back(indices);
    sorted_parent.getFloatDataArrays().push_back(samples);
    sorted_parent.getStringDataArrays().push_back(labels);
    const auto fragmentation = ms2(sorted_parent);
    auto unsorted_parent = sorted_parent;
    // Leave the first/last m/z unchanged so only interior ordering differs.
    std::reverse(unsorted_parent.begin() + 1, unsorted_parent.end() - 1);
    for (auto& array : unsorted_parent.getIntegerDataArrays())
    {
      std::reverse(array.begin() + 1, array.end() - 1);
    }
    for (auto& array : unsorted_parent.getFloatDataArrays())
    {
      std::reverse(array.begin() + 1, array.end() - 1);
    }
    for (auto& array : unsorted_parent.getStringDataArrays())
    {
      std::reverse(array.begin() + 1, array.end() - 1);
    }
    TEST_TRUE(sorted_parent.isSorted())
    TEST_FALSE(unsorted_parent.isSorted())
    RNAPrecursorCorrector corrector(10, 1, 5, positive, 8000);
    MSExperiment sorted_experiment, unsorted_experiment;
    sorted_experiment.setSpectra({sorted_parent, fragmentation});
    unsorted_experiment.setSpectra({unsorted_parent, fragmentation});
    const auto expected_rows = corrector.correct(sorted_experiment);
    const auto actual_rows = corrector.correct(unsorted_experiment);
    TEST_EQUAL(expected_rows.size(), 1)
    TEST_EQUAL(actual_rows.size(), 1)
    ABORT_IF(expected_rows.size() != 1 || actual_rows.size() != 1)
    const auto& expected = expected_rows[0];
    const auto& actual = actual_rows[0];
    TEST_TRUE(expected.accepted())
    TEST_TRUE(actual.accepted())
    TEST_EQUAL(actual.status, expected.status)
    TEST_EQUAL(actual.spectrum, expected.spectrum)
    TEST_EQUAL(actual.parent, expected.parent)
    TEST_EQUAL(actual.parent_source, expected.parent_source)
    TEST_EQUAL(actual.isotope, expected.isotope)
    TEST_EQUAL(actual.charge, expected.charge)
    TEST_EQUAL(actual.original_mz, expected.original_mz)
    TEST_EQUAL(actual.original_charge, expected.original_charge)
    TEST_EQUAL(actual.corrected_mz, expected.corrected_mz)
    TEST_EQUAL(actual.mono_mass, expected.mono_mass)
    TEST_EQUAL(actual.cosine, expected.cosine)
    TEST_EQUAL(actual.snr, expected.snr)
    TEST_EQUAL(actual.selected_error_ppm, expected.selected_error_ppm)
    TEST_TRUE(actual.mono_observed == expected.mono_observed)
    TEST_EQUAL(actual.parent_peak_count, expected.parent_peak_count)
    TEST_EQUAL(actual.parent_query_count, expected.parent_query_count)
    TEST_EQUAL(sorted_experiment[0], sorted_parent)
    TEST_EQUAL(unsorted_experiment[0], unsorted_parent)
    TEST_EQUAL(unsorted_experiment[1], sorted_experiment[1])
  }
}
END_SECTION
START_SECTION(([EXTRA] directed FLASH MS1: Different precursors sharing an MS1 use one batch and retain the full scan))
{
  RNAPrecursorCorrector corrector(10, 1, 5, false, 8000);
  const auto first = survey();
  auto second = first;
  for (auto& peak : second) peak.setMZ(peak.getMZ() + 30.0);
  auto parent = first;
  parent.insert(parent.end(), second.begin(), second.end());
  parent.sortByPosition();
  auto a = ms2(first, 1);
  auto b = ms2(second, 3);
  b.setNativeID("scan=4");
  auto repeat = a;
  repeat.setNativeID("scan=5");
  auto unknown = a;
  unknown.setNativeID("scan=6");
  unknown.getPrecursors()[0].setCharge(0);
  MSExperiment experiment;
  experiment.setSpectra({parent, a, b, repeat, unknown});
  const auto rows = corrector.correct(experiment);
  TEST_EQUAL(rows.size(), 4)
  for (Size i = 0; i < rows.size(); ++i)
  {
    TEST_TRUE(rows[i].accepted())
    TEST_EQUAL(rows[i].parent_peak_count, parent.size())
    TEST_EQUAL(rows[i].charge, 2)
    TEST_EQUAL(rows[i].parent_query_count, 3)
  }
  TEST_EQUAL(rows[0].isotope, 1)
  TEST_EQUAL(rows[1].isotope, 3)
  TEST_REAL_SIMILAR(rows[0].corrected_mz, first[0].getMZ())
  TEST_REAL_SIMILAR(rows[1].corrected_mz, second[0].getMZ())
  TEST_REAL_SIMILAR(rows[2].corrected_mz, rows[0].corrected_mz)
  TEST_EQUAL(rows[3].original_charge, 0)
  TEST_EQUAL(experiment[0], parent)
}
END_SECTION
START_SECTION(([EXTRA] directed FLASH MS1: std::vector<Result> correct(MSExperiment& experiment) retains overlapping co-batch corrections))
{
  MSExperiment fixture;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("RNAPrecursorCorrector_competition.mzML"), fixture);
  RNAPrecursorCorrector corrector(10, 1, 20, false, 8800.181396484375, 0.85, 0.25);
  for (Size i = 0; i < 2; ++i)
  {
    const auto& parent = fixture[2 * i];
    const auto& child = fixture[2 * i + 1];
    MSExperiment alone;
    alone.setSpectra({parent, child});
    const auto expected = corrector.correct(alone)[0];
    auto overlapping = child;
    overlapping.setNativeID("overlapping");
    auto& query = overlapping.getPrecursors()[0];
    query.setMZ(query.getMZ() + Constants::C13C12_MASSDIFF_U / std::max(1, std::abs(query.getCharge())));
    auto unknown = child;
    unknown.setNativeID("unknown_charge");
    unknown.getPrecursors()[0].setCharge(0);
    for (bool reverse : {false, true})
    {
      MSExperiment batch;
      batch.setSpectra(reverse ? std::vector<MSSpectrum>{parent, child, unknown, overlapping} :
                                std::vector<MSSpectrum>{parent, child, overlapping, unknown});
      const auto rows = corrector.correct(batch);
      TEST_EQUAL(rows.size(), 3)
      ABORT_IF(rows.size() != 3)
      const auto& actual = rows[0];
      TEST_EQUAL(actual.status, expected.status)
      TEST_EQUAL(actual.isotope, expected.isotope)
      TEST_EQUAL(actual.charge, expected.charge)
      TEST_EQUAL(actual.original_mz, expected.original_mz)
      TEST_EQUAL(actual.original_charge, expected.original_charge)
      TEST_EQUAL(actual.corrected_mz, expected.corrected_mz)
      TEST_EQUAL(actual.mono_mass, expected.mono_mass)
      TEST_EQUAL(actual.cosine, expected.cosine)
      TEST_EQUAL(actual.snr, expected.snr)
      TEST_EQUAL(actual.selected_error_ppm, expected.selected_error_ppm)
      TEST_TRUE(actual.mono_observed == expected.mono_observed)
      TEST_EQUAL(batch[0], parent)
    }
  }
}
END_SECTION
START_SECTION(([EXTRA] directed FLASH MS1: Full-MS1 evidence at distant charge states))
{
  for (bool positive : {false, true})
  {
    auto parent = survey(positive);
    const double carrier = positive ? Constants::PROTON_MASS_U : -Constants::PROTON_MASS_U;
    const auto charge_two = parent;
    for (const auto& peak : charge_two)
    {
      parent.emplace_back((peak.getMZ() - carrier) * 2 / 3 + carrier, peak.getIntensity() * 0.8);
    }
    parent.sortByPosition();
    SpectralDeconvolution full;
    Param parameters = full.getParameters();
    parameters.setValue("tol", DoubleList{10, 10});
    parameters.setValue("min_charge", positive ? 1 : -5);
    parameters.setValue("max_charge", positive ? 5 : -1);
    parameters.setValue("max_mass", 8000.0);
    parameters.setValue("min_cos", DoubleList{0.85, 0.85});
    parameters.setValue("min_snr", DoubleList{0.25, 0.25});
    full.setParameters(parameters);
    full.calculateAveragine(true);
    for (int charge : {0, 2})
    {
      auto fragmentation = ms2(charge_two, 1);
      fragmentation.getPrecursors()[0].setCharge(charge);
      RNAPrecursorCorrector corrector(10, 1, 5, positive, 8000);
      full.performSpectrumDeconvolution(parent, 0, PeakGroup());
      const auto expected = corrector.match(full.getDeconvolvedSpectrum(), fragmentation.getPrecursors()[0]);
      MSExperiment experiment;
      experiment.setSpectra({parent, fragmentation});
      const auto result = corrector.correct(experiment)[0];
      TEST_TRUE(expected.accepted())
      TEST_EQUAL(result.status, expected.status)
      TEST_EQUAL(result.isotope, expected.isotope)
      TEST_EQUAL(result.charge, expected.charge)
      TEST_EQUAL(result.corrected_mz, expected.corrected_mz)
      TEST_EQUAL(result.original_mz, expected.original_mz)
      TEST_EQUAL(result.original_charge, expected.original_charge)
      TEST_EQUAL(result.mono_mass, expected.mono_mass)
      TEST_EQUAL(result.cosine, expected.cosine)
      TEST_EQUAL(result.snr, expected.snr)
      TEST_EQUAL(result.selected_error_ppm, expected.selected_error_ppm)
      TEST_EQUAL(result.mono_observed, expected.mono_observed)
      TEST_EQUAL(experiment[0], parent)
    }
  }
}
END_SECTION
START_SECTION(([EXTRA] directed FLASH MS1: Refined overlapping envelopes retain the full-scan competition result))
{
  // Two reduced Kaiser MS1 cases: an overlapping high-SNR envelope that fails
  // the final quality cutoff, and a true envelope whose initial apex shifted.
  MSExperiment experiment;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("RNAPrecursorCorrector_competition.mzML"), experiment);
  RNAPrecursorCorrector corrector(10, 1, 20, false, 8800.181396484375, 0.85, 0.25);
  const auto original = experiment;
  const auto rows = corrector.correct(experiment);
  TEST_EQUAL(rows.size(), 2)
  ABORT_IF(rows.size() != 2)
  SpectralDeconvolution full;
  Param parameters = full.getParameters();
  parameters.setValue("tol", DoubleList{10, 10});
  parameters.setValue("min_charge", -20);
  parameters.setValue("max_charge", -1);
  parameters.setValue("max_mass", 8800.181396484375);
  parameters.setValue("min_cos", DoubleList{0.85, 0.85});
  parameters.setValue("min_snr", DoubleList{0.25, 0.25});
  full.setParameters(parameters);
  full.calculateAveragine(true);
  for (Size i = 0; i < rows.size(); ++i)
  {
    full.performSpectrumDeconvolution(original[2 * i], int(2 * i), PeakGroup());
    const auto expected = corrector.match(full.getDeconvolvedSpectrum(), original[2 * i + 1].getPrecursors()[0]);
    TEST_EQUAL(rows[i].status, expected.status)
    TEST_EQUAL(rows[i].isotope, expected.isotope)
    TEST_EQUAL(rows[i].charge, expected.charge)
    TEST_EQUAL(rows[i].original_mz, expected.original_mz)
    TEST_EQUAL(rows[i].original_charge, expected.original_charge)
    TEST_EQUAL(rows[i].corrected_mz, expected.corrected_mz)
    TEST_EQUAL(rows[i].mono_mass, expected.mono_mass)
    TEST_EQUAL(rows[i].cosine, expected.cosine)
    TEST_EQUAL(rows[i].snr, expected.snr)
    TEST_EQUAL(rows[i].selected_error_ppm, expected.selected_error_ppm)
    TEST_EQUAL(rows[i].mono_observed, expected.mono_observed)
  }
  TEST_EQUAL(rows[0].status, "no_envelope")
  TEST_EQUAL(rows[1].status, "monoisotopic")
  TEST_EQUAL(rows[1].isotope, 0)
  TEST_EQUAL(rows[1].charge, 4)
  TEST_EQUAL(rows[1].corrected_mz, rows[1].original_mz)
  TEST_EQUAL(experiment[0], original[0])
  TEST_EQUAL(experiment[2], original[2])
}
END_SECTION
END_TEST
