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

END_TEST
