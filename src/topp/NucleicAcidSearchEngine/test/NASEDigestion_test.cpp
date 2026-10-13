// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
#define OPENMS_NASE_CLASS_TEST
#include "../NucleicAcidSearchEngine.cpp"
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <fstream>

class TestNASE : public NucleicAcidSearchEngine
{
public:
  using NucleicAcidSearchEngine::getModifications_;
  using NucleicAcidSearchEngine::reshuffleDecoyIfCollision_;
  using NucleicAcidSearchEngine::registerDigestedOligo_;
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

START_SECTION(([EXTRA] collision retries are bounded and preserve chain termini))
{
  TestNASE tool;
  const NASequence original = NASequence::fromString("pACGUp");
  const set<NASequence> targets{original};
  bool resolved = false;
  Size attempts = 0;
  const auto first = tool.reshuffleDecoyIfCollision_(original, targets, 5, resolved, attempts);
  TEST_TRUE(resolved)
  TEST_TRUE(attempts > 0 && attempts <= 5)
  TEST_EQUAL(targets.count(first), 0)
  TEST_EQUAL(first.getFivePrimeMod(), original.getFivePrimeMod())
  TEST_EQUAL(first.getThreePrimeMod(), original.getThreePrimeMod())
  const auto second = tool.reshuffleDecoyIfCollision_(original, targets, 5, resolved, attempts);
  TEST_EQUAL(first, second)
  tool.reshuffleDecoyIfCollision_(original, targets, 0, resolved, attempts);
  TEST_FALSE(resolved)
  TEST_EQUAL(attempts, 0)
  const auto unique = tool.reshuffleDecoyIfCollision_(first, targets, 0, resolved, attempts);
  TEST_TRUE(resolved)
  TEST_EQUAL(attempts, 0)
  TEST_EQUAL(unique, first)
  const auto homopolymer = NASequence::fromString("AAAA");
  const auto unresolved = tool.reshuffleDecoyIfCollision_(homopolymer, {homopolymer}, 5, resolved, attempts);
  TEST_FALSE(resolved)
  TEST_EQUAL(attempts, 5)
  TEST_EQUAL(unresolved, homopolymer)
}
END_SECTION

START_SECTION(([EXTRA] opt-in repaired decoys refer to actual durable synthetic parents))
{
  TestNASE tool;
  IdentificationData data;
  IdentificationData::ParentSequence target("TARGET", IdentificationData::MoleculeType::RNA, "pACGUp");
  IdentificationData::ParentSequence decoy("DECOY", IdentificationData::MoleculeType::RNA, "pACGUp");
  decoy.is_decoy = true;
  auto target_ref = data.registerParentSequence(target);
  auto decoy_ref = data.registerParentSequence(decoy);
  const auto parent = NASequence::fromString(target.sequence);
  set<NASequence> targets;
  Size collisions = 0, resolved = 0, unresolved = 0, exhausted = 0;
  tool.registerDigestedOligo_(data, target_ref, parent, {0, 4}, parent, targets, 5, collisions, resolved, unresolved, exhausted);
  tool.registerDigestedOligo_(data, decoy_ref, parent, {0, 4}, parent, targets, 5, collisions, resolved, unresolved, exhausted);
  TEST_EQUAL(collisions, 1)
  TEST_EQUAL(resolved, 1)
  TEST_EQUAL(unresolved, 0)
  TEST_EQUAL(data.getIdentifiedOligos().size(), 2)
  for (const auto& oligo : data.getIdentifiedOligos())
  {
    TEST_EQUAL(oligo.parent_matches.size(), 1)
    ABORT_IF(oligo.parent_matches.size() != 1)
    const auto& mapping = *oligo.parent_matches.begin();
    if (mapping.first->is_decoy)
    {
      TEST_NOT_EQUAL(oligo.sequence, parent)
      TEST_EQUAL(mapping.second.size(), 1)
      ABORT_IF(mapping.second.size() != 1)
      const auto& match = *mapping.second.begin();
      TEST_EQUAL(match.start_pos, 0)
      TEST_EQUAL(match.end_pos, 3)
      TEST_EQUAL(mapping.first->getMetaValue("reshuffled_decoy").toString(), "true")
      TEST_EQUAL(mapping.first->getMetaValue("source_parent_accession").toString(), "DECOY")
      TEST_EQUAL(NASequence::fromString(mapping.first->sequence).getSubsequence(match.start_pos, 4), oligo.sequence)
    }
  }
  tool.registerDigestedOligo_(data, decoy_ref, parent, {0, 4}, parent, targets, 0, collisions, resolved, unresolved, exhausted);
  TEST_EQUAL(unresolved, 1)
  TEST_EQUAL(data.getIdentifiedOligos().size(), 2)
  TEST_EQUAL(data.getParentSequences().size(), 3)
  TEST_EQUAL(decoy_ref->sequence, "pACGUp")
  bool shared = false;
  for (const auto& oligo : data.getIdentifiedOligos())
  {
    if (oligo.sequence == parent)
    {
      shared = true;
      TEST_EQUAL(oligo.parent_matches.size(), 2)
      TEST_FALSE(oligo.allParentsAreDecoys())
    }
  }
  TEST_TRUE(shared)
  std::string persisted;
  NEW_TMP_FILE(persisted)
  OMSFile().store(persisted, data);
  IdentificationData loaded;
  OMSFile().load(persisted, loaded);
  bool synthetic = false;
  for (const auto& source : loaded.getParentSequences())
  {
    if (source.metaValueExists("reshuffled_decoy"))
    {
      synthetic = true;
      TEST_EQUAL(source.getMetaValue("source_parent_accession").toString(), "DECOY")
      TEST_TRUE(source.accession.find("__NASE_reshuffled_") != std::string::npos)
    }
  }
  TEST_TRUE(synthetic)
  bool success = false;
  Size successful_budget = 0;
  tool.reshuffleDecoyIfCollision_(parent, targets, 5, success, successful_budget);
  TEST_TRUE(success)
  collisions = resolved = unresolved = exhausted = 0;
  tool.registerDigestedOligo_(data, decoy_ref, parent, {0, 4}, parent, targets, successful_budget,
    collisions, resolved, unresolved, exhausted);
  TEST_EQUAL(resolved, 1)
  TEST_EQUAL(exhausted, 0)

  IdentificationData contextual_data;
  const auto contextual_parent = NASequence::fromString("pGAACGUCp");
  IdentificationData::ParentSequence contextual_source("source_DECOY", IdentificationData::MoleculeType::RNA,
    contextual_parent.toString());
  contextual_source.is_decoy = true;
  auto contextual_ref = contextual_data.registerParentSequence(contextual_source);
  const auto contextual_product = contextual_parent.getSubsequence(1, 4);
  set<NASequence> contextual_targets{contextual_product};
  collisions = resolved = unresolved = exhausted = 0;
  tool.registerDigestedOligo_(contextual_data, contextual_ref, contextual_parent, {1, 4}, contextual_product,
    contextual_targets, 20, collisions, resolved, unresolved, exhausted, "RNase_T1");
  TEST_EQUAL(resolved, 1)
  TEST_EQUAL(unresolved, 0)
  TEST_EQUAL(contextual_ref->sequence, contextual_parent.toString())
  RNaseDigestion t1;
  t1.setEnzyme("RNase_T1");
  for (const auto& oligo : contextual_data.getIdentifiedOligos())
  {
    TEST_EQUAL(oligo.sequence[3]->getCode(), "G")
    const auto& mapping = *oligo.parent_matches.begin();
    const auto actual_parent = NASequence::fromString(mapping.first->sequence);
    TEST_TRUE(t1.getFragmentPositions(actual_parent) == t1.getFragmentPositions(contextual_parent))
    TEST_EQUAL(actual_parent.getSubsequence(1, 4), oligo.sequence)
    TEST_TRUE(mapping.first->accession.ends_with("source_DECOY"))
    const auto& match = *mapping.second.begin();
    TEST_EQUAL(match.start_pos, 1)
    TEST_EQUAL(match.end_pos, 4)
    TEST_EQUAL(match.left_neighbor, "G")
    TEST_EQUAL(match.right_neighbor, "U")
  }
  TEST_EXCEPTION(Exception::IllegalArgument, tool.registerDigestedOligo_(contextual_data, contextual_ref,
    contextual_parent, {6, 4}, contextual_product, contextual_targets, 0, collisions, resolved, unresolved, exhausted))
}
END_SECTION

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

START_SECTION(([EXTRA] fixed and variable modifications produce searchable RNA candidates))
{
  std::string mzml, fasta, output, cached;
  NEW_TMP_FILE(mzml)
  NEW_TMP_FILE(fasta)
  NEW_TMP_FILE(output)
  NEW_TMP_FILE(cached)
  cached += ".oms";
  mzml += ".mzML";
  fasta += ".fasta";
  output += ".idXML";
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
  auto producer_args = arguments();
  producer_args.insert(producer_args.end(), {"-digest_out", cached});
  TEST_EQUAL(runSearch(producer_args), 0)
  IdentificationData cached_data;
  OMSFile().load(cached, cached_data);
  TEST_FALSE(cached_data.getDBSearchParams().empty())
  ABORT_IF(cached_data.getDBSearchParams().empty())
  TEST_EQUAL(static_cast<int>(cached_data.getDBSearchParams().begin()->getMetaValue("digestion:max_decoy_reshuffle_attempts")), 0)
  const std::vector<std::string> cache_args{"-in", mzml, "-digest", cached, "-id_out", output, "-threads", "1"};
  TEST_EQUAL(runSearch(cache_args), 0)
  auto mismatch_args = cache_args;
  mismatch_args.insert(mismatch_args.end(), {"-fdr:max_decoy_reshuffle_attempts", "1"});
  TEST_EQUAL(runSearch(mismatch_args), 6)
  auto matching_args = cache_args;
  matching_args.insert(matching_args.end(), {"-fdr:max_decoy_reshuffle_attempts", "0"});
  TEST_EQUAL(runSearch(matching_args), 0)
  TEST_EQUAL(hitCount(output), 1)

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
  auto args = arguments("RNase_T1");
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
