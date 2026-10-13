// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
#define OPENMS_NASE_CLASS_TEST
#include "../NucleicAcidSearchEngine.cpp"
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <fstream>
#include <limits>

class TestNASE : public NucleicAcidSearchEngine
{
public:
  using NucleicAcidSearchEngine::getModifications_;
  using NucleicAcidSearchEngine::filterBlacklistedIons_;
  using NucleicAcidSearchEngine::calculateFragmentCoverage_;
  void preprocess(MSExperiment& exp, bool remove, const std::string& method)
  {
    preprocessSpectra_(exp, 10.0, true, true, true, -1, -5, true,
      false, 100, 50, "slide", false, 1000, remove, 0.01, false, 2, method);
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

START_TEST(NASEDigestion, "$Id$")

START_SECTION(([EXTRA] variable-only ambiguity setup preserves terminal hypotheses))
{
  TestNASE tool;
  tool.enableAmbiguities();
  tool.getModifications_({});
  TEST_TRUE(tool.ambiguitiesEnabled())
  TEST_EQUAL(tool.getModifications_({"mA?"}).size(), 2)
  TestNASE::HitsByScore hits;
  TestNASE::AnnotatedHit five, three;
  five.sequence = NASequence::fromString("pACGU");
  three.sequence = NASequence::fromString("ACGUp");
  RNaseDigestion::DigestedOligo oligo;
  oligo.sequence = five.sequence;
  five.oligo_ref = &oligo;
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
  digest += ".idparquet";
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
    "-db_out", output + ".db.idparquet"};
  TEST_EQUAL(runSearch(args), 0)
  std::vector<ProteinIdentification> proteins;
  PeptideIdentificationList peptides;
  IdXMLFile().load(output, proteins, peptides);
  TEST_EQUAL(proteins.size(), 1)
  ABORT_IF(proteins.size() != 1)
  TEST_REAL_SIMILAR(proteins[0].getSearchParameters().precursor_mass_tolerance, 20)
  TEST_REAL_SIMILAR(proteins[0].getSearchParameters().fragment_mass_tolerance, 25)
  TEST_TRUE(proteins[0].getSearchParameters().db == fasta)
  // The search run records the tolerances of this search, and the database and digestion of the digest
  // (it was made with -oligo:min_size 1; this search uses the default).
  IdentificationData result;
  FileHandler().loadIdentifications(output + ".db.idparquet", result, {FileTypes::IDPARQUET});
  TEST_EQUAL(result.getRuns().size(), 1)
  ABORT_IF(result.getRuns().size() != 1)
  const auto& search_run = result.getRuns().front();
  TEST_REAL_SIMILAR(search_run.getSettings().search.precursor_mass_tolerance, 20)
  TEST_REAL_SIMILAR(search_run.getSettings().search.fragment_mass_tolerance, 25)
  TEST_EQUAL(static_cast<Int>(search_run.getSettings().search.getMetaValue("rna_min_length")), 1)
  TEST_EQUAL(search_run.getDatabases().size(), 1)
  ABORT_IF(search_run.getDatabases().empty())
  TEST_EQUAL(search_run.getDatabases()[0].path, fasta)

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

  fixture("ACGUAC", "ACGUAC", -2);
  TEST_EQUAL(runSearch(arguments()), 0)
  TEST_EQUAL(hitCount(output), 1)

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
  digest += ".idparquet";
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
    // a digest as NucleicAcidSearchEngine -digest_out stores it: a catalog run with the candidate in its database sequence
    IdentificationData data;
    auto& run = data.addRun("digestion", IdentificationData::MoleculeKind::OLIGONUCLEOTIDE);
    IdentificationData::RunSettings settings;
    settings.software = "RNADigestor";
    settings.software_version = "test";
    settings.setMetaValue("identification:catalog", "true");
    if (enzyme) settings.search.setMetaValue("rna_enzyme", enzyme->getName());
    run.setSettings(settings);
    IdentificationData::Database database;
    database.path = "test.fasta";
    const auto database_id = run.addDatabase(database);
    IdentificationData::DatabaseSequence parent {database_id, "TARGET", IdentificationData::TargetDecoy::TARGET};
    parent.sequence = "ACGUAC";
    run.setDatabaseSequences(std::vector {parent});
    IdentificationData::SourceFile source;
    source.path = database.path;
    IdentificationData::Observation observation;
    observation.data_id = "digest=1";
    const auto query = run.addIdentification(run.addSource(source), observation);
    IdentificationData::MatchData oligo;
    oligo.encoding = IdentificationData::Encoding::NA_SEQUENCE;
    oligo.representation = oligo_sequence;
    oligo.target_decoy = IdentificationData::TargetDecoy::TARGET;
    oligo.sequence_evidence.push_back({database_id, "TARGET", 0, static_cast<UInt32>(NASequence::fromString(oligo_sequence).size() - 1), '[', ']'});
    run.addMatch(query, oligo);
    FileHandler().storeIdentifications(digest, data, {FileTypes::IDPARQUET});

    // Exercise the serialized metadata, including a genuinely absent enzyme.
    IdentificationData imported;
    FileHandler().loadIdentifications(digest, imported, {FileTypes::IDPARQUET});
    const auto& search = imported.getRuns().front().getSettings().search;
    TEST_EQUAL(search.metaValueExists("rna_enzyme"), enzyme != nullptr)
    if (enzyme && search.metaValueExists("rna_enzyme"))
    {
      TEST_EQUAL(search.getMetaValue("rna_enzyme").toString(), enzyme->getName())
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
  RNaseDigestion::DigestedOligo digested;
  digested.sequence = NASequence::fromString("ACGU");
  const auto* oligo = &digested;
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
  RNaseDigestion::DigestedOligo digested;
  digested.sequence = NASequence::fromString("AUCAGp");
  const auto* oligo = &digested;
  const AdductInfo sodium("Na+", EmpiricalFormula("Na"), 1);
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
  RNaseDigestion::DigestedOligo digested, other_digested;
  digested.sequence = NASequence::fromString("AUCAGp");
  other_digested.sequence = NASequence::fromString("AAUAGp");
  const auto* oligo = &digested;
  const auto* other_oligo = &other_digested;
  const AdductInfo sodium("Na+", EmpiricalFormula("Na"), 1);
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
  auto& run = data.addRun("lfq", IdentificationData::MoleculeKind::OLIGONUCLEOTIDE);
  IdentificationData::ScoreDefinition score;
  score.name = "hyperscore";
  run.setPrimaryScore(run.addScore(score));
  IdentificationData::SourceFile input;
  input.path = "lfq.mzML";
  const auto source = run.addSource(input);
  Size index = 0;
  const auto add = [&](const std::string& sequence, double rt, double weight, Int charge = -2)
  {
    IdentificationData::Observation observation;
    observation.data_id = "scan=" + std::to_string(++index);
    // observation coordinates are finite: a retention time that is not known is missing
    if (std::isfinite(rt)) observation.rt = rt;
    observation.setMetaValue("precursor_intensity", weight);
    const auto query = run.addIdentification(source, observation);
    IdentificationData::MatchData match;
    match.encoding = IdentificationData::Encoding::NA_SEQUENCE;
    match.representation = sequence;
    match.charge = charge;
    run.addMatch(query, match, {1.0});
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

END_TEST
