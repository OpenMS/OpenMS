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

START_SECTION(([EXTRA] collision-free decoys retain annotated source coordinates))
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
      TEST_EQUAL(match.getMetaValue("reshuffled_decoy").toString(), "true")
    }
  }
  tool.registerDigestedOligo_(data, decoy_ref, parent, {0, 4}, parent, targets, 0, collisions, resolved, unresolved, exhausted);
  TEST_EQUAL(unresolved, 1)
  TEST_EQUAL(data.getIdentifiedOligos().size(), 2)
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
  std::string mzml, fasta, output;
  NEW_TMP_FILE(mzml)
  NEW_TMP_FILE(fasta)
  NEW_TMP_FILE(output)
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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
