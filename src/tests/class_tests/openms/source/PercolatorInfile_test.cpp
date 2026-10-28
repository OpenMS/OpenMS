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
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <tuple>
#include <vector>

using namespace OpenMS;
using namespace std;

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

START_SECTION(([EXTRA] stampPinFeaturesOnHits: a warning names the spectra without a scan number))
{
  // The scan number pattern comes from the first spectrum reference ("scan=5"). The second one ("index=3", e.g. from
  // another input) does not match it, so its PSM gets ScanNr -1, and Percolator puts all such PSMs into one
  // cross-validation fold. That is reported, not silent (#10719).
  auto make_pid = [](const std::string& reference, const std::string& sequence)
  {
    PeptideHit hit;
    hit.setSequence(AASequence::fromString(sequence));
    hit.setCharge(2);
    hit.setScore(1.0);
    hit.setTargetDecoyType(PeptideHit::TargetDecoyType::TARGET);
    PeptideEvidence ev;
    ev.setProteinAccession("PROT1");
    ev.setAABefore('K');
    ev.setAAAfter('S');
    hit.setPeptideEvidences(std::vector<PeptideEvidence>{ev});
    PeptideIdentification pid;
    pid.setMZ(AASequence::fromString(sequence).getMZ(2));
    pid.setRT(10.0);
    pid.setSpectrumReference(reference);
    pid.setHits(std::vector<PeptideHit>{hit});
    return pid;
  };
  PeptideIdentificationList pids;
  pids.push_back(make_pid("scan=5", "SAMPLER"));
  pids.push_back(make_pid("index=3", "SAMPLEK"));
  // a spectrum whose only PSM is skipped (no target/decoy information) gets no ScanNr, so the warning leaves it out
  pids.push_back(make_pid("index=4", "SAMPLEG"));
  std::vector<PeptideHit> unknown_type = pids[2].getHits();
  unknown_type[0].setTargetDecoyType(PeptideHit::TargetDecoyType::UNKNOWN);
  pids[2].setHits(unknown_type);

  ostringstream warnings;
  getGlobalLogWarn().insert(warnings);
  TEST_EQUAL(PercolatorInfile::stampPinFeaturesOnHits(pids, "trypsin", 2, 3).size(), 1)
  getGlobalLogWarn().remove(warnings);
  TEST_EQUAL(static_cast<int>(pids[0].getHits()[0].getMetaValue("ScanNr")), 5)
  TEST_EQUAL(static_cast<int>(pids[1].getHits()[0].getMetaValue("ScanNr")), -1)
  TEST_EQUAL(pids[2].getHits()[0].metaValueExists("ScanNr"), false)
  TEST_TRUE(warnings.str().find("no scan number found in 1 of 2 spectrum references (e.g. 'index=3')") != std::string::npos)

  // all scan numbers found: no warning
  PeptideIdentificationList scans;
  scans.push_back(make_pid("scan=5", "SAMPLER"));
  scans.push_back(make_pid("scan=6", "SAMPLEK"));
  ostringstream none;
  getGlobalLogWarn().insert(none);
  PercolatorInfile::stampPinFeaturesOnHits(scans, "trypsin", 2, 3);
  getGlobalLogWarn().remove(none);
  TEST_EQUAL(none.str().find("no scan number") == std::string::npos, true)
}
END_SECTION

START_SECTION(([EXTRA] stampPinFeaturesOnHits: the enzyme features with ProteaseDB and Percolator enzyme names))
{
  // The enzyme features of a tryptic peptide with a missed cleavage, preceded by M and followed by S: the C-terminus
  // and one internal site are tryptic
  auto stamp_hit = [](const std::string& enz, char aa_before, const std::string& peptide, char aa_after, Int start)
  {
    PeptideHit hit;
    hit.setSequence(AASequence::fromString(peptide));
    hit.setCharge(2);
    hit.setScore(1.0);
    hit.setTargetDecoyType(PeptideHit::TargetDecoyType::TARGET);
    PeptideEvidence ev;
    ev.setProteinAccession("PROT1");
    ev.setStart(start);
    ev.setAABefore(aa_before);
    ev.setAAAfter(aa_after);
    hit.setPeptideEvidences(std::vector<PeptideEvidence>{ev});
    PeptideIdentification pid;
    pid.setMZ(hit.getSequence().getMZ(2));
    pid.setRT(10.0);
    pid.setSpectrumReference("scan=1");
    pid.setHits(std::vector<PeptideHit>{hit});
    PeptideIdentificationList pids;
    pids.push_back(pid);
    PercolatorInfile::stampPinFeaturesOnHits(pids, enz, 2, 3);
    return pids[0].getHits()[0];
  };
  auto stamp = [&stamp_hit](const std::string& enz, char aa_before = 'M', const std::string& peptide = "PEPKTIDER", char aa_after = 'S',
                            Int start = PeptideEvidence::UNKNOWN_POSITION)
  {
    const PeptideHit stamped = stamp_hit(enz, aa_before, peptide, aa_after, start);
    return std::make_tuple(static_cast<int>(stamped.getMetaValue("enzN")), static_cast<int>(stamped.getMetaValue("enzC")),
                           static_cast<int>(stamped.getMetaValue("enzInt")));
  };
  // ProSE -out_pin passed 'Trypsin/P', which used to count every site (1, 1, 8)
  for (const std::string enzyme : {"Trypsin/P", "trypsinp", "Trypsin", "trypsin"})
  {
    TEST_TRUE(stamp(enzyme) == std::make_tuple(0, 1, 1))
  }
  TEST_TRUE(stamp("Lys-C") == std::make_tuple(0, 0, 1)) // K|T inside; M|P and R|S are not Lys-C sites
  TEST_TRUE(stamp("lys-c") == std::make_tuple(0, 0, 1))
  // ProteaseDB enzymes without a Percolator name use the rule of ProteaseDB (they used to count every site: 1, 1, 8)
  TEST_TRUE(stamp("Lys-C/P") == std::make_tuple(0, 0, 1))
  TEST_TRUE(stamp("Arg-C/P") == std::make_tuple(0, 1, 0)) // R|S
  TEST_TRUE(stamp("glutamyl endopeptidase") == std::make_tuple(0, 0, 3)) // E|P, D|E and E|R inside
  TEST_TRUE(stamp("V8-E") == std::make_tuple(0, 0, 1)) // E|R (E|P is not a site)
  TEST_TRUE(stamp("unspecific cleavage") == std::make_tuple(1, 1, 8))
  TEST_TRUE(stamp("no cleavage") == std::make_tuple(0, 0, 0))
  // a rule that looks at two residues (after HP, KP, RP, unless before P): the residue before the flanking P is unknown,
  // so P|A counts; K-P|E inside the peptide does, with the flanking AA as context
  TEST_TRUE(stamp("proline endopeptidase", 'P', "AKPEPK") == std::make_tuple(1, 0, 1))
  TEST_TRUE(stamp("proline endopeptidase", 'K', "PEPK") == std::make_tuple(0, 0, 1)) // K|P is no site, K-P|E is
  // an enzyme the features do not know counts every site
  TEST_TRUE(stamp("no_enzyme") == std::make_tuple(1, 1, 8))
  TEST_TRUE(stamp("not an enzyme") == std::make_tuple(1, 1, 8))
  // a flank of '-' is a terminus, by either rule
  TEST_TRUE(stamp("Trypsin/P", '-') == std::make_tuple(1, 1, 1))
  TEST_TRUE(stamp("Lys-C/P", '-') == std::make_tuple(1, 0, 1))
  // so is a protein terminus as OpenMS marks it ('[', ']'), which used to count as no cleavage site
  TEST_TRUE(stamp("Trypsin/P", PeptideEvidence::N_TERMINAL_AA) == std::make_tuple(1, 1, 1))
  TEST_TRUE(stamp("Lys-C/P", PeptideEvidence::N_TERMINAL_AA) == std::make_tuple(1, 0, 1))
  TEST_TRUE(stamp("Lys-C", 'M', "PEPKTIDER", PeptideEvidence::C_TERMINAL_AA) == std::make_tuple(0, 1, 1))
  TEST_TRUE(stamp("Lys-C/P", 'M', "PEPKTIDER", PeptideEvidence::C_TERMINAL_AA) == std::make_tuple(0, 1, 1))
  // the peptide after the initial Met (initial-Met clipping) is protein N-terminal as well, by either rule: a reversed
  // decoy protein rarely starts with M, so its N-terminal peptides keep '[' and must count the same (#10726)
  TEST_TRUE(stamp("Trypsin/P", 'M', "PEPKTIDER", 'S', 1) == std::make_tuple(1, 1, 1))
  TEST_TRUE(stamp("Trypsin/P", PeptideEvidence::N_TERMINAL_AA, "PEPKTIDER", 'S', 0) == std::make_tuple(1, 1, 1))
  TEST_TRUE(stamp("Lys-C/P", 'M', "PEPKTIDER", 'S', 1) == std::make_tuple(1, 0, 1))
  TEST_TRUE(stamp("Lys-C/P", PeptideEvidence::N_TERMINAL_AA, "PEPKTIDER", 'S', 0) == std::make_tuple(1, 0, 1))
  // an M further inside the protein is no terminus
  TEST_TRUE(stamp("Trypsin/P", 'M', "PEPKTIDER", 'S', 2) == std::make_tuple(0, 1, 1))
  TEST_TRUE(stamp("Lys-C/P", 'M', "PEPKTIDER", 'S', 2) == std::make_tuple(0, 0, 1))
  // the Peptide column keeps the flanking M
  TEST_STRING_EQUAL(stamp_hit("Trypsin/P", 'M', "PEPKTIDER", 'S', 1).getMetaValue("Peptide").toString(), "M.PEPKTIDER.S")
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

START_SECTION(([EXTRA] stampPinFeaturesOnHits: ScanNr tells the spectra of a file apart))
{
  // Percolator takes the PSMs of a file with the same ScanNr for one spectrum: train_best_positive keeps one positive
  // per ScanNr, and the cross-validation splits by it. Thermo keeps its scan numbers. Bruker TDF DDA spectra get their
  // precursor; Waters, SCIEX and other Bruker TDF native IDs, which number spectra only within a part of the run, a
  // running number per file.
  auto scan_numbers = [](const std::vector<std::pair<std::string, std::string>>& spectra) // (file, native ID)
  {
    PeptideIdentificationList pids;
    for (const auto& [file, native_id] : spectra)
    {
      PeptideHit hit;
      hit.setSequence(AASequence::fromString("PEPTIDEK"));
      hit.setCharge(2);
      hit.setScore(1.0);
      hit.setTargetDecoyType(PeptideHit::TargetDecoyType::TARGET);
      PeptideEvidence ev;
      ev.setProteinAccession("PROT1");
      ev.setAABefore('K');
      ev.setAAAfter('A');
      hit.setPeptideEvidences(std::vector<PeptideEvidence>{ev});
      PeptideIdentification pid;
      pid.setMZ(465.7);
      pid.setRT(10.0);
      pid.setSpectrumReference(native_id);
      pid.setMetaValue("file_origin", file);
      pid.setHits(std::vector<PeptideHit>{hit});
      pids.push_back(pid);
    }
    PercolatorInfile::stampPinFeaturesOnHits(pids, "trypsin", 2, 3);
    std::string numbers;
    for (const PeptideIdentification& pid : pids)
    {
      if (!numbers.empty()) numbers += " ";
      numbers += pid.getHits()[0].getMetaValue("ScanNr").toString();
    }
    return numbers;
  };

  // Thermo: the scan numbers
  TEST_STRING_EQUAL(scan_numbers({{"a", "controllerType=0 controllerNumber=1 scan=1001"},
                                  {"a", "controllerType=0 controllerNumber=1 scan=1002"}}), "1001 1002")
  // Bruker TDF DDA: the precursors, not the first ion mobility scan of the precursor (475 in all three frames)
  TEST_STRING_EQUAL(scan_numbers({{"a", "frame=120 scan=475 precursor=2"},
                                  {"a", "frame=388 scan=475 precursor=61"},
                                  {"a", "frame=1203 scan=475 precursor=230"}}), "2 61 230")
  // SCIEX: the experiment (the last number) recurs in every cycle. A running number per file, the same for the
  // PSMs of one spectrum.
  TEST_STRING_EQUAL(scan_numbers({{"a", "sample=1 period=1 cycle=100 experiment=4"},
                                  {"a", "sample=1 period=1 cycle=101 experiment=4"},
                                  {"a", "sample=1 period=1 cycle=100 experiment=4"},
                                  {"b", "sample=1 period=1 cycle=100 experiment=4"}}), "1 2 1 1")
  // Waters: scans are numbered within a function
  TEST_STRING_EQUAL(scan_numbers({{"a", "function=2 process=0 scan=345"},
                                  {"a", "function=3 process=0 scan=345"}}), "1 2")
  // Bruker TDF without a precursor (diaPASEF): scans are numbered within a frame
  TEST_STRING_EQUAL(scan_numbers({{"a", "frame=120 windowGroup=3 scan=475"},
                                  {"a", "frame=121 windowGroup=3 scan=475"}}), "1 2")
  // "frame=F" alone (Bruker TDF MS1, diaTracer pseudo-DDA): the frame
  TEST_STRING_EQUAL(scan_numbers({{"a", "frame=10"}, {"a", "frame=11"}}), "10 11")
  // files from different instruments: each file numbers its spectra after the native ID format of its first one
  TEST_STRING_EQUAL(scan_numbers({{"a", "frame=120 scan=475 precursor=2"},
                                  {"b", "controllerType=0 controllerNumber=1 scan=1001"},
                                  {"c", "sample=1 period=1 cycle=100 experiment=4"},
                                  {"a", "frame=388 scan=475 precursor=61"}}), "2 1001 1 61")
  TEST_STRING_EQUAL(scan_numbers({{"b", "controllerType=0 controllerNumber=1 scan=1001"},
                                  {"a", "frame=120 scan=475 precursor=2"}}), "1001 2")
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

START_SECTION(([EXTRA] load: a single-file PIN with no FileName column loads without crashing))
{
  // CPP-161: load() unconditionally read map_filename_to_idx.at(raw_file_name) to set id_merge_index,
  // but that map is only populated when a FileName column exists in the header. A single-file PIN,
  // the exact shape PercolatorInfile::store() itself writes (see the "one file: no FileName column"
  // section above), has no such column, so raw_file_name stays "UNKNOWN" and .at() threw
  // std::out_of_range on the very first PSM -- OpenMS could not read back its own PIN output.
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
  pid.setMetaValue("file_origin", "a.mzML");

  PeptideIdentificationList one_file;
  one_file.push_back(pid);

  StringList feature_set = PercolatorInfile::getStandardFeatureSet(2, 3);
  feature_set.push_back("retentiontime"); // required by load(), not by store()'s own internal checks
  feature_set.push_back("Peptide");
  feature_set.push_back("Proteins");

  std::string pin_file;
  NEW_TMP_FILE(pin_file);
  PercolatorInfile::store(pin_file, one_file, feature_set, "trypsin", 2, 3);

  StringList filenames;
  StringList extra_scores;
  PeptideIdentificationList loaded = PercolatorInfile::load(pin_file, true, "dm", extra_scores, filenames, "DECOY_");
  TEST_EQUAL(loaded.size(), 1)
  TEST_EQUAL(filenames.empty(), true)
  TEST_EQUAL(loaded[0].metaValueExists("id_merge_index"), false)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
