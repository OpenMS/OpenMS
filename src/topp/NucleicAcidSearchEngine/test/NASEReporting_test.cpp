// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
#define OPENMS_NASE_CLASS_TEST
#include "../NucleicAcidSearchEngine.cpp"
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <fstream>
#include <limits>

class TestNASEReporting : public NucleicAcidSearchEngine
{
public:
  using NucleicAcidSearchEngine::HitsByScore;
  using NucleicAcidSearchEngine::AnnotatedHit;
  using NucleicAcidSearchEngine::PrecursorInfo;
  using NucleicAcidSearchEngine::retainTopHit_;
  using NucleicAcidSearchEngine::sortHitsByValue_;
  using NucleicAcidSearchEngine::generateLFQInput_;
};

static Int runSearch(const std::vector<std::string>& arguments)
{
  std::vector<const char*> argv{"NucleicAcidSearchEngine"};
  for (const auto& arg : arguments)
  {
    argv.push_back(arg.c_str());
  }
  TestNASEReporting tool;
  return tool.main(static_cast<Int>(argv.size()), argv.data());
}

static std::vector<std::tuple<std::string, double, Int>> readHits(const std::string& path)
{
  std::vector<ProteinIdentification> proteins;
  PeptideIdentificationList peptides;
  IdXMLFile().load(path, proteins, peptides);
  std::vector<std::tuple<std::string, double, Int>> result;
  for (const auto& peptide : peptides)
  {
    for (const auto& hit : peptide.getHits())
    {
      result.emplace_back(hit.getSequence().toString(), hit.getScore(), hit.getCharge());
    }
  }
  std::sort(result.begin(), result.end());
  return result;
}

START_TEST(NASEReporting, "$Id$")

START_SECTION(([EXTRA] complete cutoff ties are retained for every arrival order))
{
  TestNASEReporting tool;
  IdentificationData data;
  auto oligo = data.registerIdentifiedOligo(IdentificationData::IdentifiedOligo(NASequence::fromString("ACGU")));
  TestNASEReporting::PrecursorInfo precursor{0, -2, 0, nullopt};
  std::vector<std::string> sequences{"AAA", "AAC", "AAG", "AAU"};
  std::vector<double> scores{12.0, 10.0, 10.0, 9.0};
  std::vector<Size> order{0, 1, 2, 3};
  do
  {
    TestNASEReporting::HitsByScore hits;
    for (Size index : order)
    {
      TestNASEReporting::AnnotatedHit hit;
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
    TEST_EQUAL(retained == std::vector<std::string>({"AAA", "AAC", "AAG"}), true)
  }
  while (std::next_permutation(order.begin(), order.end()));
}
END_SECTION

START_SECTION(([EXTRA] LFQ uses finite weights or the median of finite observations))
{
  TestNASEReporting tool;
  IdentificationData data;
  auto input = data.registerInputFile(IdentificationData::InputFile("lfq.mzML"));
  Size index = 0;
  const auto add = [&](const std::string& sequence, double rt, double weight)
  {
    auto oligo = data.registerIdentifiedOligo(IdentificationData::IdentifiedOligo(NASequence::fromString(sequence)));
    IdentificationData::Observation obs("scan=" + std::to_string(++index), input, rt);
    obs.setMetaValue("precursor_intensity", weight);
    auto observation = data.registerObservation(obs);
    data.registerObservationMatch(IdentificationData::ObservationMatch(oligo, observation, -2));
  };
  const double nan = std::numeric_limits<double>::quiet_NaN();
  add("ACGU", 10, 0);
  add("ACGU", 20, -1);
  add("ACGU", 30, nan);
  add("ACGU", nan, 100);
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
  }
  TEST_REAL_SIMILAR(rts.at("ACGU"), 20)
  TEST_REAL_SIMILAR(rts.at("AGGU"), 25)
  TEST_EQUAL(rts.count("AUUU"), 0)
}
END_SECTION

START_SECTION(([EXTRA] full search and cached digest report current parameters))
{
  std::string mzml, fasta, output, digest, lfq;
  NEW_TMP_FILE(mzml)
  NEW_TMP_FILE(fasta)
  NEW_TMP_FILE(output)
  NEW_TMP_FILE(digest)
  NEW_TMP_FILE(lfq)
  mzml += ".mzML";
  fasta += ".fasta";
  output += ".idXML";
  digest += ".oms";
  lfq += ".tsv";
  std::ofstream(fasta) << ">TARGET\nACGUAC\n";
  NucleicAcidSpectrumGenerator generator;
  Param parameters = generator.getParameters();
  parameters.setValue("add_first_prefix_ion", "true");
  parameters.setValue("add_a-B_ions", "true");
  parameters.setValue("add_metainfo", "true");
  generator.setParameters(parameters);
  const auto sequence = NASequence::fromString("ACGUAC");
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
  std::vector<std::string> args{"-in", mzml, "-database", fasta, "-id_out", output, "-out", output + ".mzTab",
    "-oligo:enzyme", "no cleavage", "-threads", "1", "-digest_out", digest, "-lfq_out", lfq};
  TEST_EQUAL(runSearch(args), 0)
  const auto serial = readHits(output);
  TEST_EQUAL(serial.size(), 1)
  args[11] = "32";
  TEST_EQUAL(runSearch(args), 0)
  TEST_EQUAL(serial == readHits(output), true)
  args = {"-in", mzml, "-digest", digest, "-id_out", output, "-out", output + ".mzTab",
    "-precursor:mass_tolerance", "20", "-fragment:mass_tolerance", "25", "-threads", "1", "-db_out", output + ".oms"};
  TEST_EQUAL(runSearch(args), 0)
  std::vector<ProteinIdentification> proteins;
  PeptideIdentificationList peptides;
  IdXMLFile().load(output, proteins, peptides);
  TEST_REAL_SIMILAR(proteins[0].getSearchParameters().precursor_mass_tolerance, 20)
  TEST_REAL_SIMILAR(proteins[0].getSearchParameters().fragment_mass_tolerance, 25)
  TEST_EQUAL(proteins[0].getSearchParameters().db == fasta, true)
  IdentificationData result;
  OMSFile().load(output + ".oms", result);
  bool current_identification = false;
  for (const auto& mapping : result.getDBSearchSteps())
  {
    if (mapping.first->actions.count(DataProcessing::IDENTIFICATION))
    {
      current_identification = true;
      TEST_REAL_SIMILAR(mapping.second->precursor_mass_tolerance, 20)
      TEST_REAL_SIMILAR(mapping.second->fragment_mass_tolerance, 25)
    }
  }
  TEST_TRUE(current_identification)
}
END_SECTION

END_TEST
