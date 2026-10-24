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
#include <OpenMS/FORMAT/PercolatorInfile.h>
///////////////////////////

#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <vector>

using namespace OpenMS;
using namespace std;

namespace
{
  // PSMs for buildFeatureMatrix: 300 spectra of two files with up to 3 hits; hits without evidence or target/decoy
  // type, decoys, an input CalcMass, both isotope error meta values, an extra feature as a number, as a string and
  // missing, and an input meta value under a PIN column name ('score')
  PeptideIdentificationList makeMatrixInput()
  {
    const std::vector<std::string> peptides{"SAMPLER", "PEPTIDEK", "ACDEFGHIK", "M(Oxidation)LVEK", ".(Acetyl)SAMPLEK"};
    PeptideIdentificationList pids;
    for (int i = 0; i < 300; ++i)
    {
      PeptideIdentification pid;
      pid.setMZ(400.0 + 0.37 * i);
      pid.setRT(10.0 + i);
      pid.setSpectrumReference("scan=" + std::to_string(1000 + i % 150));
      pid.setMetaValue("file_origin", i < 150 ? "a.mzML" : "b.mzML");
      for (int k = 0; k < 1 + i % 3; ++k)
      {
        const int n = i * 3 + k;
        PeptideHit hit;
        hit.setSequence(AASequence::fromString(peptides[n % peptides.size()]));
        hit.setCharge(2 + n % 3);
        hit.setScore(0.1 * n - 3.0);
        if (n % 17 != 5) hit.setTargetDecoyType(n % 4 == 1 ? PeptideHit::TargetDecoyType::DECOY : PeptideHit::TargetDecoyType::TARGET);
        if (n % 23 != 7)
        {
          PeptideEvidence ev;
          ev.setProteinAccession(n % 4 == 1 ? "DECOY_P" : "P" + std::to_string(n % 9));
          ev.setAABefore(n % 5 == 0 ? '[' : 'K');
          ev.setAAAfter(n % 7 == 0 ? ']' : 'A');
          hit.addPeptideEvidence(ev);
        }
        if (n % 11 == 3) hit.setMetaValue("CalcMass", 777.25 + n);
        if (n % 13 == 2) hit.setMetaValue("IsotopeError", "1");
        if (n % 6 == 4) hit.setMetaValue(Constants::UserParam::ISOTOPE_ERROR, -1);
        if (n % 19 != 8) hit.setMetaValue("engine_x", n % 2 ? DataValue(0.5 * n) : DataValue(std::to_string(n) + ".25"));
        if (n % 10 == 9) hit.setMetaValue("score", 42.0);
        pid.insertHit(hit);
      }
      pids.push_back(pid);
    }
    return pids;
  }
}

START_TEST(PercolatorInfile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PercolatorInfile* ptr = nullptr;
PercolatorInfile* null_pointer = nullptr;

START_SECTION(PercolatorInfile())
{
  ptr = new PercolatorInfile();
  TEST_NOT_EQUAL(ptr, null_pointer);
}
END_SECTION

START_SECTION(~PercolatorInfile())
{
  delete ptr;
}
END_SECTION

START_SECTION(PeptideIdentificationList PercolatorInfile::load(const std::string& pin_file, bool higher_score_better, const std::string& score_name, std::string decoy_prefix))
{
  StringList filenames;
  // test loading of pin file with automatic update of target/decoy annotation based on decoy prefix in protein accessions

  // test some extra scores
  StringList extra_scores = {"ln(delta_next)", "ln(delta_best)", "matched_peaks"};

  auto pids = PercolatorInfile::load(OPENMS_GET_TEST_DATA_PATH("sage.pin"), 
    true, 
    "ln(hyperscore)", 
    extra_scores, 
    filenames, 
    "DECOY_");
  TEST_EQUAL(pids.size(), 9)
  TEST_EQUAL(filenames.size(), 2)
  TEST_EQUAL(pids[0].getSpectrumReference(), "30381")
  TEST_EQUAL(pids[6].getSpectrumReference(), "spectrum=2041")
  TEST_EQUAL(pids[7].getHits()[0].getMetaValue("target_decoy"),"decoy") // 8th entry is annotated as target in pin file but only maps to decoy proteins with prefix "DECOY_" -> set to decoy
}
END_SECTION

START_SECTION((static StringList getStandardFeatureSet(int min_charge, int max_charge)))
{
  // The standard Percolator feature set is the three mandatory header columns
  // (SpecId, Label, ScanNr), then the mass/length features, then one "charge<c>"
  // column per charge state in [min,max], then the enzyme / mass-delta features.
  // Callers append any extra_features and finally Peptide, Proteins before writing
  // the .pin header that Percolator expects (#9195). This validates that contract.

  // exact ordered feature set for charges 2..4: three mandatory columns, then
  // mass/length features, then one "charge<c>" column per charge state, then the
  // enzyme / mass-delta features.
  StringList fs = PercolatorInfile::getStandardFeatureSet(2, 4);
  StringList expected = {"SpecId", "Label", "ScanNr", "ExpMass", "CalcMass", "mass", "peplen",
                         "charge2", "charge3", "charge4",
                         "enzN", "enzC", "enzInt", "dm", "absdm"};
  TEST_EQUAL(fs.size(), expected.size())
  for (Size i = 0; i < expected.size(); ++i)
  {
    TEST_STRING_EQUAL(fs[i], expected[i])
  }

  // the three mandatory leading columns
  TEST_STRING_EQUAL(fs[0], "SpecId")
  TEST_STRING_EQUAL(fs[1], "Label")
  TEST_STRING_EQUAL(fs[2], "ScanNr")

  // a single charge state yields exactly one "charge3" column and no other charge column
  StringList single = PercolatorInfile::getStandardFeatureSet(3, 3);
  StringList expected_single = {"SpecId", "Label", "ScanNr", "ExpMass", "CalcMass", "mass", "peplen",
                                "charge3", "enzN", "enzC", "enzInt", "dm", "absdm"};
  TEST_EQUAL(single.size(), expected_single.size())
  for (Size i = 0; i < expected_single.size(); ++i)
  {
    TEST_STRING_EQUAL(single[i], expected_single[i])
  }

  // the assembled .pin header (features + Peptide + Proteins) begins with SpecId
  // and ends with Peptide, Proteins -- the exact shape Percolator parses (#9195).
  StringList header = fs;
  header.push_back("Peptide");
  header.push_back("Proteins");
  TEST_STRING_EQUAL(header.front(), "SpecId")
  TEST_STRING_EQUAL(header[header.size() - 2], "Peptide")
  TEST_STRING_EQUAL(header.back(), "Proteins")
}
END_SECTION

START_SECTION((static void store(const std::string& pin_file, const PeptideIdentificationList& peptide_ids, const StringList& feature_set, const std::string& enz, int min_charge, int max_charge)))
{
  // Validates the -out_pin output that PercolatorAdapter/ProSE/OpenNuXL emit:
  // the header line must be the exact column contract Percolator parses
  // (SpecId Label ScanNr ... Peptide Proteins) and store() must compute and
  // write the standard per-PSM features (stampPinFeaturesOnHits). Build one
  // fully-annotated target PSM so it is not skipped and a data row is written.
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("SAMPLER"));
  hit.setCharge(2);
  hit.setScore(1.0);
  hit.setTargetDecoyType(PeptideHit::TargetDecoyType::TARGET);
  PeptideEvidence ev;
  ev.setProteinAccession("PROT1");
  ev.setAABefore('K'); // tryptic flanks so enzN/enzC are true
  ev.setAAAfter('S');
  hit.setPeptideEvidences(std::vector<PeptideEvidence>{ev});

  PeptideIdentification pid;
  pid.setMZ(500.25);
  pid.setRT(123.4);
  pid.setSpectrumReference("scan=529");
  pid.setHits(std::vector<PeptideHit>{hit});

  PeptideIdentificationList pids;
  pids.push_back(pid);

  // the real .pin column contract: standard features + Peptide + Proteins
  StringList feature_set = PercolatorInfile::getStandardFeatureSet(2, 3);
  feature_set.push_back("Peptide");
  feature_set.push_back("Proteins");

  std::string pin_file;
  NEW_TMP_FILE(pin_file);
  PercolatorInfile::store(pin_file, pids, feature_set, "trypsin", 2, 3);

  // read the written .pin back
  std::ifstream is(pin_file.c_str());
  std::vector<std::string> lines;
  std::string line;
  while (std::getline(is, line)) lines.push_back(line);

  // header + exactly one data row (the PSM must not be skipped)
  TEST_EQUAL(lines.size(), 2)
  ABORT_IF(lines.empty())

  // --- header: column names and order Percolator expects ---
  StringList header;
  StringUtils::split(lines[0], '\t', header);
  TEST_EQUAL(header.size(), feature_set.size())
  TEST_STRING_EQUAL(header[0], "SpecId")
  TEST_STRING_EQUAL(header[1], "Label")
  TEST_STRING_EQUAL(header[2], "ScanNr")
  TEST_STRING_EQUAL(header[header.size() - 2], "Peptide")
  TEST_STRING_EQUAL(header.back(), "Proteins")

  auto has = [&header](const std::string& n) {
    return std::find(header.begin(), header.end(), n) != header.end();
  };
  // mass + enzyme feature columns are present
  TEST_EQUAL(has("ExpMass"), true)
  TEST_EQUAL(has("CalcMass"), true)
  TEST_EQUAL(has("mass"), true)
  TEST_EQUAL(has("peplen"), true)
  TEST_EQUAL(has("charge2"), true)
  TEST_EQUAL(has("charge3"), true)
  TEST_EQUAL(has("enzN"), true)
  TEST_EQUAL(has("enzC"), true)
  TEST_EQUAL(has("enzInt"), true)
  TEST_EQUAL(has("dm"), true)
  TEST_EQUAL(has("absdm"), true)

  // --- data row: same column count + the deterministic computed feature values ---
  ABORT_IF(lines.size() < 2)
  StringList row;
  StringUtils::split(lines[1], '\t', row);
  TEST_EQUAL(row.size(), header.size())

  auto col = [&header](const std::string& n) -> Size {
    return (Size)(std::find(header.begin(), header.end(), n) - header.begin());
  };
  TEST_STRING_EQUAL(row[col("SpecId")], "scan=529")
  TEST_STRING_EQUAL(row[col("ScanNr")], "529")
  TEST_STRING_EQUAL(row[col("Label")], "1")          // target -> 1
  TEST_STRING_EQUAL(row[col("peplen")], "7")         // SAMPLER
  TEST_STRING_EQUAL(row[col("charge2")], "1")        // charge == 2 (one-hot)
  TEST_STRING_EQUAL(row[col("charge3")], "0")
  TEST_STRING_EQUAL(row[col("enzN")], "1")           // tryptic N-terminus
  TEST_STRING_EQUAL(row[col("enzC")], "1")           // tryptic C-terminus
  TEST_STRING_EQUAL(row[col("Peptide")], "K.SAMPLER.S")
  TEST_STRING_EQUAL(row[col("Proteins")], "PROT1")
}
END_SECTION

START_SECTION(([EXTRA] stampPinFeaturesOnHits: the isotope error corrects dm/absdm only; ExpMass is the precursor m/z of the spectrum))
{
  // One spectrum whose precursor was selected at the first 13C isotope peak of SAMPLER (2+), with two PSMs: SAMPLER
  // at isotope_error +1 (observed minus theoretical) and SAMPLEK at isotope_error 0. Percolator identifies a spectrum
  // by ScanNr and ExpMass, so both rows keep the observed precursor m/z there and in 'mass'.
  const AASequence sampler = AASequence::fromString("SAMPLER");
  const AASequence samplek = AASequence::fromString("SAMPLEK");
  const double observed_mz = sampler.getMZ(2) + Constants::C13C12_MASSDIFF_U / 2.0;
  auto make_hit = [](const AASequence& sequence, int isotope_error, double score)
  {
    PeptideHit hit;
    hit.setSequence(sequence);
    hit.setCharge(2);
    hit.setScore(score);
    hit.setTargetDecoyType(PeptideHit::TargetDecoyType::TARGET);
    PeptideEvidence ev;
    ev.setProteinAccession("PROT1");
    ev.setAABefore('K');
    ev.setAAAfter('S');
    hit.setPeptideEvidences(std::vector<PeptideEvidence>{ev});
    hit.setMetaValue(Constants::UserParam::ISOTOPE_ERROR, isotope_error);
    return hit;
  };
  PeptideIdentification pid;
  pid.setMZ(observed_mz);
  pid.setRT(10.0);
  pid.setSpectrumReference("scan=7");
  pid.setHits(std::vector<PeptideHit>{make_hit(sampler, 1, 2.0), make_hit(samplek, 0, 1.0)});
  PeptideIdentificationList pids;
  pids.push_back(pid);

  TEST_EQUAL(PercolatorInfile::stampPinFeaturesOnHits(pids, "trypsin", 2, 3).size(), 0)
  const std::vector<PeptideHit>& hits = pids[0].getHits();
  ABORT_IF(hits.size() != 2)
  TOLERANCE_ABSOLUTE(1e-9)
  for (const PeptideHit& hit : hits)
  {
    TEST_REAL_SIMILAR(static_cast<double>(hit.getMetaValue("ExpMass")), observed_mz)
    TEST_REAL_SIMILAR(static_cast<double>(hit.getMetaValue("mass")), observed_mz)
  }
  // SAMPLER: the observed precursor lies one 13C spacing above it; the correction removes that
  TEST_REAL_SIMILAR(static_cast<double>(hits[0].getMetaValue("dm")), 0.0)
  TEST_REAL_SIMILAR(static_cast<double>(hits[0].getMetaValue("absdm")), 0.0)
  TEST_REAL_SIMILAR(static_cast<double>(hits[0].getMetaValue("deltamass")), 0.0)
  // SAMPLEK: no isotope error, the plain difference
  TEST_REAL_SIMILAR(static_cast<double>(hits[1].getMetaValue("dm")), observed_mz - samplek.getMZ(2))
  TEST_REAL_SIMILAR(static_cast<double>(hits[1].getMetaValue("absdm")), std::abs(observed_mz - samplek.getMZ(2)))
}
END_SECTION

START_SECTION(([EXTRA] stampPinFeaturesOnHits: an error stops at the failing hit, as in a serial loop))
{
  // 300 spectra (stamped in parallel); the hit of spectrum 150 carries a CalcMass that is not a number. Every hit has
  // a 'score' meta value before, which stamping overwrites.
  PeptideIdentificationList pids;
  for (int i = 0; i < 300; ++i)
  {
    PeptideHit hit;
    hit.setSequence(AASequence::fromString("SAMPLER"));
    hit.setCharge(2);
    hit.setScore(1.0 + i);
    hit.setTargetDecoyType(PeptideHit::TargetDecoyType::TARGET);
    PeptideEvidence ev;
    ev.setProteinAccession("PROT1");
    ev.setAABefore('K');
    ev.setAAAfter('S');
    hit.setPeptideEvidences(std::vector<PeptideEvidence>{ev});
    hit.setMetaValue("score", 42.0);
    if (i == 150) hit.setMetaValue("CalcMass", "bad");
    PeptideIdentification pid;
    pid.setMZ(400.0);
    pid.setRT(10.0 + i);
    pid.setSpectrumReference("scan=" + std::to_string(i + 1));
    pid.setHits(std::vector<PeptideHit>{hit});
    pids.push_back(pid);
  }
  PercolatorInfile::PinFeatureMetaValueMap added;
  TEST_EXCEPTION(Exception::InvalidValue, PercolatorInfile::stampPinFeaturesOnHits(pids, "trypsin", 2, 3, added))
  // stamped and recorded up to the failing hit, including the keys it added before the error
  TEST_EQUAL(added.size(), 151)
  TEST_EQUAL(added.count({149, 0}), 1)
  TEST_EQUAL(static_cast<double>(pids[149].getHits()[0].getMetaValue("score")), 150.0)
  ABORT_IF(added.count({150, 0}) != 1)
  TEST_EQUAL(added.at({150, 0}).size(), 3)  // SpecId, ScanNr, Label
  TEST_EQUAL(pids[150].getHits()[0].metaValueExists("SpecId"), true)
  // the spectra after it are as before
  bool untouched = true;
  for (Size i = 151; i < pids.size(); ++i)
  {
    const PeptideHit& hit = pids[i].getHits()[0];
    untouched = untouched && !hit.metaValueExists("SpecId") && !hit.metaValueExists("CalcMass")
                && static_cast<double>(hit.getMetaValue("score")) == 42.0;
  }
  TEST_EQUAL(untouched, true)
}
END_SECTION

START_SECTION(([EXTRA] isotope correction preserves spectrum identity in stamped features and exported PIN rows))
{
  const AASequence sequence = AASequence::fromString(".(Acetyl)SAM(Oxidation)PLER");
  for (const std::string& key : {std::string(Constants::UserParam::ISOTOPE_ERROR), std::string("IsotopeError")})
  {
    for (const int charge : {2, 3})
    {
      const double observed_mz = sequence.getMZ(charge) + 0.02;
      PeptideIdentification pid;
      pid.setMZ(observed_mz);
      pid.setRT(10.0);
      pid.setSpectrumReference("scan=7");
      for (const int isotope_error : {-1, 0, 1, 2})
      {
        PeptideHit hit;
        hit.setSequence(sequence);
        hit.setCharge(charge);
        hit.setScore(1.0);
        hit.setTargetDecoyType(isotope_error == 0 ? PeptideHit::TargetDecoyType::DECOY : PeptideHit::TargetDecoyType::TARGET);
        PeptideEvidence evidence;
        evidence.setProteinAccession(isotope_error == 0 ? "DECOY_PROT1" : "PROT1");
        evidence.setAABefore('K');
        evidence.setAAAfter('S');
        hit.addPeptideEvidence(evidence);
        // Zero without a key also covers inputs without isotope annotation.
        if (isotope_error != 0) hit.setMetaValue(key, isotope_error);
        pid.insertHit(hit);
      }
      PeptideIdentificationList pids;
      pids.push_back(pid);
      std::string pin_file;
      NEW_TMP_FILE(pin_file);
      StringList features = PercolatorInfile::getStandardFeatureSet(2, 3);
      features.insert(features.end(), {"deltamass", "Peptide", "Proteins"});
      PercolatorInfile::store(pin_file, pids, features, "trypsin", 2, 3);
      // store() must leave its input untouched.
      TEST_FALSE(pids[0].getHits()[0].metaValueExists("ExpMass"))
      TEST_EQUAL(PercolatorInfile::stampPinFeaturesOnHits(pids, "trypsin", 2, 3).size(), 0)
      TextFile txt(pin_file);
      const StringList lines(txt.begin(), txt.end());
      ABORT_IF(lines.size() != 5)
      auto col = [&features](const std::string& name) -> Size
      {
        return std::find(features.begin(), features.end(), name) - features.begin();
      };
      const std::vector<int> isotope_errors = {-1, 0, 1, 2};
      for (Size i = 0; i < isotope_errors.size(); ++i)
      {
        const auto& hit = pids[0].getHits()[i];
        const double delta_mass = observed_mz - isotope_errors[i] * Constants::C13C12_MASSDIFF_U / charge - sequence.getMZ(charge);
        TEST_REAL_SIMILAR(static_cast<double>(hit.getMetaValue("ExpMass")), observed_mz)
        TEST_REAL_SIMILAR(static_cast<double>(hit.getMetaValue("mass")), observed_mz)
        TEST_REAL_SIMILAR(static_cast<double>(hit.getMetaValue("CalcMass")), sequence.getMZ(charge))
        TEST_REAL_SIMILAR(static_cast<double>(hit.getMetaValue("dm")), delta_mass)
        TEST_REAL_SIMILAR(static_cast<double>(hit.getMetaValue("absdm")), std::abs(delta_mass))
        TEST_REAL_SIMILAR(static_cast<double>(hit.getMetaValue("deltamass")), delta_mass)
        StringList row;
        StringUtils::split(lines[i + 1], '\t', row);
        ABORT_IF(row.size() != features.size())
        TEST_STRING_EQUAL(row[col("ScanNr")], "7")
        TEST_REAL_SIMILAR(StringUtils::toDouble(row[col("ExpMass")]), observed_mz)
        TEST_REAL_SIMILAR(StringUtils::toDouble(row[col("mass")]), observed_mz)
        TEST_REAL_SIMILAR(StringUtils::toDouble(row[col("dm")]), delta_mass)
        TEST_REAL_SIMILAR(StringUtils::toDouble(row[col("absdm")]), std::abs(delta_mass))
        TEST_REAL_SIMILAR(StringUtils::toDouble(row[col("deltamass")]), delta_mass)
      }
    }
  }
}
END_SECTION

START_SECTION((static std::string getFileIdentifier(const PeptideIdentification& pid)))
{
  PeptideIdentification pid;
  TEST_STRING_EQUAL(PercolatorInfile::getFileIdentifier(pid), "")
  pid.setMetaValue("file_origin", "a.idXML");
  TEST_STRING_EQUAL(PercolatorInfile::getFileIdentifier(pid), "a.idXML")
  pid.setMetaValue("id_merge_index", 1);
  TEST_STRING_EQUAL(PercolatorInfile::getFileIdentifier(pid), "a.idXML|1")
  // an index never runs into the file origin
  PeptideIdentification other;
  other.setMetaValue("file_origin", "a.idXML1");
  TEST_STRING_EQUAL(PercolatorInfile::getFileIdentifier(other), "a.idXML1")
  other.setMetaValue("file_origin", "run1");
  other.setMetaValue("id_merge_index", 2);
  pid.setMetaValue("file_origin", "run");
  pid.setMetaValue("id_merge_index", 12);
  TEST_NOT_EQUAL(PercolatorInfile::getFileIdentifier(pid), PercolatorInfile::getFileIdentifier(other))
  PeptideIdentification merged;
  merged.setMetaValue("id_merge_index", 0);
  TEST_STRING_EQUAL(PercolatorInfile::getFileIdentifier(merged), "|0")
}
END_SECTION

START_SECTION(([EXTRA] store: PSMs of several spectrum files get a FileName column among the optional columns))
{
  // Percolator identifies a spectrum by spectrum file, ScanNr and ExpMass. Without a FileName column all PSMs
  // count as one file, so equal scan numbers and precursor m/z of two files would be one spectrum.
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("SAMPLER"));
  hit.setCharge(2);
  hit.setScore(1.0);
  hit.setTargetDecoyType(PeptideHit::TargetDecoyType::TARGET);
  PeptideEvidence ev;
  ev.setProteinAccession("PROT1");
  ev.setAABefore('K');
  ev.setAAAfter('S');
  hit.setPeptideEvidences(std::vector<PeptideEvidence>{ev});
  PeptideIdentification pid;
  pid.setMZ(500.25);
  pid.setRT(123.4);
  pid.setSpectrumReference("scan=529");
  pid.setHits(std::vector<PeptideHit>{hit});

  StringList feature_set = PercolatorInfile::getStandardFeatureSet(2, 3);
  feature_set.push_back("Peptide");
  feature_set.push_back("Proteins");
  auto store = [&](const PeptideIdentificationList& pids) {
    std::string pin_file;
    NEW_TMP_FILE(pin_file);
    PercolatorInfile::store(pin_file, pids, feature_set, "trypsin", 2, 3);
    std::ifstream is(pin_file.c_str());
    std::vector<StringList> lines;
    std::string line;
    while (std::getline(is, line))
    {
      lines.emplace_back();
      StringUtils::split(line, '\t', lines.back());
    }
    return lines;
  };

  // one file: no FileName column (the .pin file is unchanged)
  PeptideIdentificationList one_file;
  pid.setMetaValue("file_origin", "a.mzML");
  one_file.push_back(pid);
  one_file.push_back(pid);
  auto lines = store(one_file);
  ABORT_IF(lines.size() != 3)
  TEST_EQUAL(lines[0] == feature_set, true)

  // two files: FileName after SpecId, Label, ScanNr, ExpMass and CalcMass, before the first feature
  PeptideIdentificationList two_files = one_file;
  two_files[1].setMetaValue("file_origin", "b.mzML");
  lines = store(two_files);
  ABORT_IF(lines.size() != 3)
  StringList expected = feature_set;
  expected.insert(expected.begin() + 5, "FileName");
  TEST_EQUAL(lines[0] == expected, true)
  TEST_STRING_EQUAL(lines[0][4], "CalcMass")
  TEST_STRING_EQUAL(lines[0][6], "mass")
  TEST_EQUAL(lines[1].size(), expected.size())
  TEST_EQUAL(lines[2].size(), expected.size())
  TEST_STRING_EQUAL(lines[1][5], "a.mzML")
  TEST_STRING_EQUAL(lines[2][5], "b.mzML")
  // the same scan number, told apart by the file
  TEST_STRING_EQUAL(lines[1][2], "529")
  TEST_STRING_EQUAL(lines[2][2], "529")
  TEST_STRING_EQUAL(lines[1][0], "a.mzMLscan=529")
  TEST_STRING_EQUAL(lines[2][0], "b.mzMLscan=529")

  // a merged file: 'id_merge_index' tells the files apart as well
  PeptideIdentificationList merged = one_file;
  merged[0].setMetaValue("id_merge_index", 0);
  merged[1].setMetaValue("id_merge_index", 1);
  lines = store(merged);
  ABORT_IF(lines.size() != 3)
  TEST_EQUAL(lines[0] == expected, true)
  TEST_STRING_EQUAL(lines[1][5], "a.mzML|0")
  TEST_STRING_EQUAL(lines[2][5], "a.mzML|1")
}
END_SECTION

START_SECTION((static double getFeatureValue(const DataValue& value, const std::string& feature)))
{
  // numeric meta values count as they are
  TEST_REAL_SIMILAR(PercolatorInfile::getFeatureValue(DataValue(2.5), "f"), 2.5)
  TEST_REAL_SIMILAR(PercolatorInfile::getFeatureValue(DataValue(-7), "f"), -7.0)

  // Search engine scores read from text are strings (e.g. SageAdapter stores 'SAGE:ln(-poisson)'
  // like this). They count as the numbers they spell, as for the executable parsing the .pin.
  TEST_REAL_SIMILAR(PercolatorInfile::getFeatureValue(DataValue("1.527394838114684"), "SAGE:ln(-poisson)"), 1.527394838114684)
  TEST_REAL_SIMILAR(PercolatorInfile::getFeatureValue(DataValue(" -3e-2 "), "f"), -0.03)
  TEST_REAL_SIMILAR(PercolatorInfile::getFeatureValue(DataValue("48"), "SAGE:scored_candidates"), 48.0)

  // text with all the digits of a double reads back as that double
  const double value = 0.1 + 0.2;
  std::ostringstream value_text;
  value_text << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
  TEST_EQUAL(PercolatorInfile::getFeatureValue(DataValue(value_text.str()), "f"), value)

  // no numeric value: an error, not an unrelated number
  TEST_EXCEPTION(Exception::InvalidValue, PercolatorInfile::getFeatureValue(DataValue("abc"), "f"))
  TEST_EXCEPTION(Exception::InvalidValue, PercolatorInfile::getFeatureValue(DataValue("1.5 abc"), "f"))
  TEST_EXCEPTION(Exception::InvalidValue, PercolatorInfile::getFeatureValue(DataValue(""), "f"))
  TEST_EXCEPTION(Exception::InvalidValue, PercolatorInfile::getFeatureValue(DataValue(), "f"))
  TEST_EXCEPTION(Exception::InvalidValue, PercolatorInfile::getFeatureValue(DataValue(DoubleList{1.0, 2.0}), "f"))

  // not finite, as text or as a number: the executable stops at such a feature as well
  TEST_EXCEPTION(Exception::InvalidValue, PercolatorInfile::getFeatureValue(DataValue("nan"), "f"))
  TEST_EXCEPTION(Exception::InvalidValue, PercolatorInfile::getFeatureValue(DataValue("inf"), "f"))
  TEST_EXCEPTION(Exception::InvalidValue, PercolatorInfile::getFeatureValue(DataValue("-inf"), "f"))
  TEST_EXCEPTION(Exception::InvalidValue, PercolatorInfile::getFeatureValue(DataValue(std::numeric_limits<double>::quiet_NaN()), "f"))
  TEST_EXCEPTION(Exception::InvalidValue, PercolatorInfile::getFeatureValue(DataValue(std::numeric_limits<double>::infinity()), "f"))
  TEST_EXCEPTION(Exception::InvalidValue, PercolatorInfile::getFeatureValue(DataValue(-std::numeric_limits<double>::infinity()), "f"))
}
END_SECTION

START_SECTION((static Size buildFeatureMatrix(const PeptideIdentificationList& peptide_ids, const StringList& features, const std::string& enz, int min_charge, int max_charge, RescoreInput& input, std::vector<std::pair<size_t, size_t>>& rows)))
{
  const PeptideIdentificationList pids = makeMatrixInput();
  StringList features = PercolatorInfile::getStandardFeatureSet(2, 3);
  features.erase(std::remove_if(features.begin(), features.end(), [](const std::string& f)
    { return f == "SpecId" || f == "Label" || f == "ScanNr" || f == "ExpMass"; }), features.end());
  features.insert(features.end(), {"deltamass", "retentiontime", "score", "engine_x"});

  // reference: the stamped meta values, read back as in-process rescoring formerly did
  PeptideIdentificationList stamped = pids;
  const auto skipped = PercolatorInfile::stampPinFeaturesOnHits(stamped, "trypsin", 2, 3);
  std::vector<std::vector<double>> ref_rows;
  std::vector<std::pair<size_t, size_t>> ref_locations;
  std::vector<int> ref_scans, ref_files;
  std::vector<double> ref_exp, ref_calc;
  std::vector<bool> ref_decoy;
  std::map<std::string, int> files;
  for (size_t i = 0; i < stamped.size(); ++i)
  {
    for (size_t j = 0; j < stamped[i].getHits().size(); ++j)
    {
      if (skipped.contains({i, j})) continue;
      const PeptideHit& hit = stamped[i].getHits()[j];
      std::vector<double> row;
      bool ok = true;
      for (const std::string& f : features)
      {
        if (!hit.metaValueExists(f)) { ok = false; break; }
        row.push_back(PercolatorInfile::getFeatureValue(hit.getMetaValue(f), f));
      }
      if (!ok) continue;
      ref_rows.push_back(row);
      ref_locations.emplace_back(i, j);
      ref_decoy.push_back(static_cast<int>(hit.getMetaValue("Label")) == -1);
      ref_scans.push_back(static_cast<int>(hit.getMetaValue("ScanNr")));
      ref_exp.push_back(static_cast<double>(hit.getMetaValue("ExpMass")));
      ref_calc.push_back(PercolatorInfile::getFeatureValue(hit.getMetaValue("CalcMass"), "CalcMass"));
      ref_files.push_back(files.try_emplace(PercolatorInfile::getFileIdentifier(stamped[i]), static_cast<int>(files.size())).first->second);
    }
  }
  TEST_EQUAL(skipped.empty(), false)
  TEST_EQUAL(ref_rows.size() > 300, true)

  RescoreInput input;
  std::vector<std::pair<size_t, size_t>> rows;
  TEST_EQUAL(PercolatorInfile::buildFeatureMatrix(pids, features, "trypsin", 2, 3, input, rows), skipped.size())
  TEST_EQUAL(rows == ref_locations, true)
  TEST_EQUAL(input.features.empty(), true)
  TEST_EQUAL(input.feature_names == features, true)
  ABORT_IF(input.rowCount() != ref_rows.size() || input.featureCount() != features.size())
  bool same_values = true; // bit for bit
  for (size_t r = 0; r < ref_rows.size(); ++r)
  {
    same_values = same_values && std::equal(ref_rows[r].begin(), ref_rows[r].end(), input.row(r));
  }
  TEST_EQUAL(same_values, true)
  TEST_EQUAL(input.is_decoy == ref_decoy, true)
  TEST_EQUAL(input.scan_numbers == ref_scans, true)
  TEST_EQUAL(input.spec_file_numbers == ref_files, true)
  TEST_EQUAL(input.exp_masses == ref_exp, true)
  TEST_EQUAL(input.calc_masses == ref_calc, true)
  // no meta value was added to the hits
  TEST_EQUAL(pids == makeMatrixInput(), true)

  // errors as by stamping, then by reading the rows: an unknown target/decoy value, a feature that is not a number
  PeptideIdentificationList bad_td = pids;
  bad_td[200].getHits()[0].setMetaValue("target_decoy", "neither");
  TEST_EXCEPTION(Exception::InvalidValue, PercolatorInfile::buildFeatureMatrix(bad_td, features, "trypsin", 2, 3, input, rows))
  PeptideIdentificationList bad_feature = pids;
  bad_feature[100].getHits()[0].setMetaValue("engine_x", "abc");
  TEST_EXCEPTION(Exception::InvalidValue, PercolatorInfile::buildFeatureMatrix(bad_feature, features, "trypsin", 2, 3, input, rows))
  // a string column is not a feature
  TEST_EXCEPTION(Exception::InvalidValue, PercolatorInfile::buildFeatureMatrix(pids, StringList{"mass", "Peptide"}, "trypsin", 2, 3, input, rows))
}
END_SECTION

START_SECTION((static void refreshPinMetaValues(PeptideIdentificationList& peptide_ids, const std::string& enz, int min_charge, int max_charge)))
{
  // as stamping: an input meta value under a PIN column name gets its PIN value; no other meta value is added
  const PeptideIdentificationList pids = makeMatrixInput();
  PeptideIdentificationList refreshed = pids;
  PercolatorInfile::refreshPinMetaValues(refreshed, "trypsin", 2, 3);
  PeptideIdentificationList stamped = pids;
  PercolatorInfile::PinFeatureMetaValueMap added;
  PercolatorInfile::stampPinFeaturesOnHits(stamped, "trypsin", 2, 3, added);
  for (const auto& [location, keys] : added)
  {
    for (const std::string& key : keys) stamped[location.first].getHits()[location.second].removeMetaValue(key);
  }
  TEST_EQUAL(refreshed == stamped, true)
  TEST_EQUAL(refreshed == pids, false)
  TEST_EQUAL(static_cast<double>(refreshed[3].getHits()[0].getMetaValue("score")), refreshed[3].getHits()[0].getScore())
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
