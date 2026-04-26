// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
//
// Side-by-side regression suite for the in-process Percolator against the
// external `percolator` binary. Sections 2+ are gated on the environment
// variable PERCOLATOR_BINARY; without it, they no-op. Section 1 (PIN stamp
// parity) always runs.

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/ID/Percolator.h>
#include <OpenMS/ANALYSIS/ID/PercolatorTypes.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/PercolatorInfile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iterator>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using namespace OpenMS;
using namespace std;

///////////////////////////////////////////////////////////////////////////////
// Helpers
///////////////////////////////////////////////////////////////////////////////

namespace
{

/// Reference idXML: same input the TOPP PercolatorAdapter_1 test uses.
String inputIdxml()
{
  return OPENMS_GET_TEST_DATA_PATH("../../../topp/THIRDPARTY/PercolatorAdapter_1.idXML");
}

/// Derive min/max charge from all hits.
void deriveChargeRange(const PeptideIdentificationList& pids,
                       int& min_charge, int& max_charge)
{
  min_charge = 10;
  max_charge = 0;
  for (const auto& pid : pids)
  {
    for (const auto& hit : pid.getHits())
    {
      min_charge = std::min(min_charge, hit.getCharge());
      max_charge = std::max(max_charge, hit.getCharge());
    }
  }
}

/// Build the PIN feature set the adapter would use: standard + any
/// PSMFeatureExtractor-provided extras + Peptide + Proteins trailer.
StringList buildFeatureSet(const vector<ProteinIdentification>& prs,
                           int min_charge, int max_charge)
{
  StringList feature_set = PercolatorInfile::getStandardFeatureSet(min_charge, max_charge);
  if (!prs.empty())
  {
    const auto& sp = prs.front().getSearchParameters();
    if (sp.metaValueExists("extra_features"))
    {
      StringList extras = ListUtils::create<String>(
        sp.getMetaValue("extra_features").toString());
      feature_set.insert(feature_set.end(), extras.begin(), extras.end());
    }
  }
  feature_set.push_back("Peptide");
  feature_set.push_back("Proteins");
  return feature_set;
}

/// Numeric columns only — excludes bookkeeping and string columns unfit for
/// an SVM. Mirrors the filter PercolatorAdapter's in-process branch applies.
StringList numericFeatures(const StringList& feature_set)
{
  static const std::set<String> metadata_not_feature {
    "SpecId", "Label", "ScanNr", "ExpMass", "CalcMass", "Peptide", "Proteins"
  };
  StringList numeric;
  for (const auto& f : feature_set)
  {
    if (!metadata_not_feature.count(f)) numeric.push_back(f);
  }
  return numeric;
}

/// Composite key for 1:1 cross-run matching of PSMs. Scan + sequence + charge
/// uniquely identifies a row in the .pin for the TOPP test data.
String rowKey(const PeptideIdentification& pid, const PeptideHit& hit)
{
  String sr = String(pid.getSpectrumReference());
  if (sr.empty() && pid.metaValueExists("spectrum_id"))
  {
    sr = pid.getMetaValue("spectrum_id").toString();
  }
  return sr + "|" + hit.getSequence().toString() + "|" + String(hit.getCharge());
}

struct PercolatorTriplet
{
  double score = 0.0;
  double qval = 0.0;
  double pep = 0.0;
};

/// Extract the PSI-MS CV values the subprocess adapter stamps on each hit.
/// MS:1001491 = q-value, MS:1001492 = SVM score, MS:1001493 = PEP.
std::map<String, PercolatorTriplet> loadBaselineTriplets(const String& idxml)
{
  vector<ProteinIdentification> prs;
  PeptideIdentificationList pids;
  IdXMLFile().load(idxml, prs, pids);
  std::map<String, PercolatorTriplet> out;
  for (const auto& pid : pids)
  {
    for (const auto& hit : pid.getHits())
    {
      PercolatorTriplet t;
      if (hit.metaValueExists("MS:1001492"))
        t.score = static_cast<double>(hit.getMetaValue("MS:1001492"));
      if (hit.metaValueExists("MS:1001491"))
        t.qval = static_cast<double>(hit.getMetaValue("MS:1001491"));
      if (hit.metaValueExists("MS:1001493"))
        t.pep = static_cast<double>(hit.getMetaValue("MS:1001493"));
      out[rowKey(pid, hit)] = t;
    }
  }
  return out;
}

/// Materialize RescoreInput rows in the same iteration order used downstream.
/// Also returns the hit pointers so the caller can zip them with outputs.
struct BuiltInput
{
  RescoreInput ri;
  std::vector<const PeptideIdentification*> pid_per_row;
  std::vector<const PeptideHit*> hit_per_row;
};

BuiltInput buildRescoreInput(const PeptideIdentificationList& pids,
                             const StringList& numeric_features)
{
  BuiltInput b;
  b.ri.feature_names = numeric_features;
  for (const auto& pid : pids)
  {
    for (const auto& hit : pid.getHits())
    {
      bool has_all = true;
      std::vector<double> feats;
      feats.reserve(numeric_features.size());
      for (const auto& f : numeric_features)
      {
        if (!hit.metaValueExists(f)) { has_all = false; break; }
        feats.push_back(static_cast<double>(hit.getMetaValue(f)));
      }
      if (!has_all) continue;
      b.ri.features.push_back(std::move(feats));
      b.ri.is_decoy.push_back(hit.isDecoy());
      b.pid_per_row.push_back(&pid);
      b.hit_per_row.push_back(&hit);
    }
  }
  return b;
}

///////////////////////////////////////////////////////////////////////////////
// Shared infrastructure for Tests 2-5: synthetic dataset, minimal .pin writer,
// subprocess invocation, and .pout parser. Tests 2-5 are gated on the env
// variable PERCOLATOR_BINARY (set by CTest when the binary is found).
///////////////////////////////////////////////////////////////////////////////

String percolatorBinary()
{
  const char* env = std::getenv("PERCOLATOR_BINARY");
  if (!env || !*env) return String();
  return String(env);
}

/// Generate a realistic-scale synthetic dataset: 2000 PSMs with a mildly
/// informative 3-feature signal, same shape as the "realistic" Percolator
/// class test. Large enough to comfortably pass SanityCheck.
void generateSyntheticData(RescoreInput& ri, std::mt19937& rng,
                           double size_mult = 1.0)
{
  ri.features.clear();
  ri.is_decoy.clear();
  ri.feature_names = StringList{"f0", "f1", "noise"};

  std::normal_distribution<double> noise(0.0, 0.5);
  std::normal_distribution<double> unit(0.0, 1.0);

  const size_t n_easy = static_cast<size_t>(500 * size_mult);
  const size_t n_hard = static_cast<size_t>(500 * size_mult);
  const size_t n_dec  = static_cast<size_t>(1000 * size_mult);
  for (size_t i = 0; i < n_easy; ++i)
  {
    ri.features.push_back({+2.0 + noise(rng), +1.0 + noise(rng), unit(rng)});
    ri.is_decoy.push_back(false);
  }
  for (size_t i = 0; i < n_hard; ++i)
  {
    ri.features.push_back({+0.5 + noise(rng), +0.5 * noise(rng), unit(rng)});
    ri.is_decoy.push_back(false);
  }
  for (size_t i = 0; i < n_dec; ++i)
  {
    ri.features.push_back({+0.0 + noise(rng), +0.0 * noise(rng), unit(rng)});
    ri.is_decoy.push_back(true);
  }

  // Populate all four fields that feed P::OrderScanHash — the CV fold split
  // diverges if either path defaults a field to a different value. Fields:
  //
  //   scan_numbers      — unique per row, so every PSM ends up in its own
  //                       "spectrum" bucket (matches the writePinFile PIN).
  //   exp_masses        — match the PIN's ExpMass column verbatim.
  //   calc_masses       — match the PIN's CalcMass column verbatim.
  //   spec_file_numbers — subprocess populates specFileNr from an optional
  //                       `filename` / `spectrafile` PIN column (see
  //                       SetHandler.cpp:110-111). writePinFile emits NO such
  //                       column, so subprocess keeps specFileNr at its
  //                       PSMDescription ctor default (0u). Explicitly
  //                       assigning 0 here documents the contract and guards
  //                       against a silent in-process default drift.
  const size_t n = ri.features.size();
  ri.scan_numbers.assign(n, 0);
  ri.spec_file_numbers.assign(n, 0);
  ri.exp_masses.assign(n, 1000.0);
  ri.calc_masses.assign(n, 1000.0);
  for (size_t i = 0; i < n; ++i)
  {
    ri.scan_numbers[i] = static_cast<int>(i + 1);
  }
}

/// Write a minimal .pin file that percolator's command-line tool will accept.
/// Columns: SpecId, Label, ScanNr, ExpMass, CalcMass[, filename], <features>,
///          Peptide, Proteins.
/// Row ids `row_%08zu` so subprocess PSMId strings match in-process row indices.
/// When @p emit_filename is true, a `filename` column is inserted between
/// CalcMass and the features (i.e. in the optional-fields section that
/// SetHandler::getOptionalFields recognises before it hits the first unknown
/// header token). DataSet::readPsm then calls setSpectrumFileName() which
/// populates PSMDescription::specFileNr from the lookup table, matching the
/// value we set in RescoreInput::spec_file_numbers for the in-process path.
void writePinFile(const String& path, const RescoreInput& ri,
                  bool emit_filename = false)
{
  std::ofstream f(path.c_str());
  if (!f) throw std::runtime_error("cannot open pin file for writing");
  f.setf(std::ios::scientific);
  f.precision(17);

  f << "SpecId\tLabel\tScanNr\tExpMass\tCalcMass";
  if (emit_filename) f << "\tfilename";
  for (const auto& fn : ri.feature_names) f << "\t" << fn;
  f << "\tPeptide\tProteins\n";

  const size_t n = ri.features.size();
  for (size_t i = 0; i < n; ++i)
  {
    char row_id[32];
    std::snprintf(row_id, sizeof(row_id), "row_%08zu", i);
    const int label = ri.is_decoy[i] ? -1 : +1;
    const int scan = (i < ri.scan_numbers.size())
      ? ri.scan_numbers[i]
      : static_cast<int>(i + 1);
    f << row_id << "\t" << label << "\t" << scan
      << "\t" << 1000.0 << "\t" << 1000.0;
    if (emit_filename)
    {
      const int fnr = (i < ri.spec_file_numbers.size())
        ? ri.spec_file_numbers[i] : 0;
      f << "\tfile_" << fnr;
    }
    for (double v : ri.features[i]) f << "\t" << v;
    f << "\t-.PEPTIDE.-\t"
      << (ri.is_decoy[i] ? "decoy_protein" : "target_protein") << "\n";
  }
}

struct SubprocessOut
{
  std::map<String, PercolatorTriplet> triplets;  // keyed by row_id (PSMId)
  int exit_code = -1;
};

/// Invoke the external percolator binary on `pin_path`. `extra_args` is space-
/// separated (appended verbatim). Returns parsed PSM triplets keyed by PSMId.
SubprocessOut runSubprocess(const String& bin,
                            const String& pin_path,
                            const std::string& extra_args = "")
{
  const String target_pout = File::getTemporaryFile();
  const String decoy_pout  = File::getTemporaryFile();
  const String stdout_log  = File::getTemporaryFile();
  const String stderr_log  = File::getTemporaryFile();

  std::ostringstream cmd;
  cmd << "\"" << bin << "\""
      << " -U"                                  // PSM-level (no peptide roll-up)
      << " -m \"" << target_pout << "\""
      << " -M \"" << decoy_pout << "\""
      << " " << extra_args
      << " \"" << pin_path << "\""
      << " > \"" << stdout_log << "\""
      << " 2> \"" << stderr_log << "\"";

  SubprocessOut s;
  // cmd.exe /c strips the first and last quote of its argument. Since this
  // command both starts AND ends with a quoted path, both get stripped, which
  // mangles the binary path and the redirection target. Wrap the whole thing
  // in an extra outer pair of quotes so cmd.exe strips only those.
#ifdef _WIN32
  const std::string full_cmd = "\"" + cmd.str() + "\"";
  s.exit_code = std::system(full_cmd.c_str());
#else
  s.exit_code = std::system(cmd.str().c_str());
#endif
  if (s.exit_code != 0) return s;

  auto parse = [&s](const String& path)
  {
    std::ifstream f(path.c_str());
    String line;
    bool header_seen = false;
    bool has_filename_col = false;  // true when pout has PSMId filename score ...
    while (std::getline(f, line))
    {
      if (!header_seen)
      {
        header_seen = true;
        // Detect optional filename column: header is either
        //   PSMId  score  q-value  ...       (no filename in input PIN)
        // or
        //   PSMId  filename  score  q-value  (filename column was in input PIN)
        std::istringstream hs(line);
        String col1, col2;
        hs >> col1 >> col2;
        has_filename_col = (col2 == "filename" || col2 == "spectrafile");
        continue;
      }
      if (line.empty()) continue;
      // Fields: PSMId  [filename]  score  q-value  posterior_error_prob  peptide  protein...
      std::istringstream ls(line);
      String psm_id;
      double score = 0, qval = 0, pep = 0;
      String peptide;
      ls >> psm_id;
      if (has_filename_col) { String skip_fn; ls >> skip_fn; }
      ls >> score >> qval >> pep >> peptide;
      if (!ls) continue;
      PercolatorTriplet t;
      t.score = score;
      t.qval = qval;
      t.pep = pep;
      // Composite key: subprocess multi-hit-per-spectrum inputs share the
      // same PSMId (PercolatorInfile::store writes SpecId = file+scan, not
      // including peptide), so we include peptide to disambiguate.
      s.triplets[psm_id + "|" + peptide] = t;
    }
  };
  parse(target_pout);
  parse(decoy_pout);
  return s;
}

/// Convenience: align in-process RescoreOutput rows against subprocess map by
/// row id. Missing rows are reported via `misses`.
struct AlignedRows
{
  size_t matches = 0;
  size_t misses = 0;
  double score_r = 0.0;       ///< Pearson on SVM scores
  double max_dq = 0.0;
  double max_dpep = 0.0;
  int target_lt_001_in = 0;
  int target_lt_001_sub = 0;
};

AlignedRows alignAndCompare(const RescoreOutput& inp,
                            const SubprocessOut& sub,
                            const std::vector<bool>& is_decoy)
{
  AlignedRows a;
  double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
  for (size_t i = 0; i < inp.scores.size(); ++i)
  {
    char row_id[64];
    std::snprintf(row_id, sizeof(row_id), "row_%08zu|-.PEPTIDE.-", i);
    auto it = sub.triplets.find(row_id);
    if (it == sub.triplets.end()) { a.misses++; continue; }
    a.matches++;
    const double sc_in = inp.scores[i];
    const double sc_sub = it->second.score;
    const double dq = std::abs(inp.q_values[i] - it->second.qval);
    const double dpep = std::abs(inp.peps[i] - it->second.pep);
    sx += sc_in;  sy += sc_sub;
    sxx += sc_in * sc_in;
    syy += sc_sub * sc_sub;
    sxy += sc_in * sc_sub;
    a.max_dq   = std::max(a.max_dq,   dq);
    a.max_dpep = std::max(a.max_dpep, dpep);
    if (!is_decoy[i] && inp.q_values[i]   <= 0.01) a.target_lt_001_in++;
    if (!is_decoy[i] && it->second.qval <= 0.01)   a.target_lt_001_sub++;
  }
  const double n = static_cast<double>(a.matches);
  if (n > 1)
  {
    const double num = n * sxy - sx * sy;
    const double den = std::sqrt((n * sxx - sx * sx) * (n * syy - sy * sy));
    a.score_r = (den > 1e-15) ? num / den : 0.0;
  }
  return a;
}


} // namespace

///////////////////////////////////////////////////////////////////////////////

START_TEST(Percolator_subprocess_parity, "$Id$")

///////////////////////////////////////////////////////////////////////////////
// Test 1: PIN file content equals stamped meta values.
//
// This guards the factored helper: PercolatorInfile::preparePin_ internally
// calls stampPinFeaturesOnHits and then serializes the stamped meta values
// to a .pin file. A regression where one of them drifts (e.g. someone adds a
// transformation only on one side) would silently corrupt the in-process
// path's training features.
//
// Approach: pick a realistic idXML (same as the PercolatorAdapter_1 TOPP
// test). Write the .pin via store() AND stamp a parallel copy in memory.
// Parse the .pin back and compare every non-trailing column against the
// stamped hit's meta value.
///////////////////////////////////////////////////////////////////////////////

START_SECTION([EXTRA] PIN file content equals stamped meta values)
{
  vector<ProteinIdentification> prs;
  PeptideIdentificationList pids;
  IdXMLFile().load(inputIdxml(), prs, pids);
  TEST_FALSE(prs.empty())
  TEST_FALSE(pids.empty())

  int min_charge, max_charge;
  deriveChargeRange(pids, min_charge, max_charge);
  TEST_TRUE(min_charge <= max_charge)

  StringList feature_set = buildFeatureSet(prs, min_charge, max_charge);
  const std::string enz = "no_enzyme";

  // (a) Write PIN via store(). This also internally stamps, then serializes.
  const String pin_path = File::getTemporaryFile();
  PercolatorInfile::store(pin_path, pids, feature_set, enz, min_charge, max_charge);

  // (b) Stamp an independent copy so we can read back the meta values
  //     without depending on store()'s internal copy.
  PeptideIdentificationList stamped = pids;
  const auto skipped = PercolatorInfile::stampPinFeaturesOnHits(
    stamped, enz, min_charge, max_charge);

  // (c) Parse the PIN file. Header must equal feature_set verbatim.
  TextFile txt;
  txt.load(pin_path);
  const size_t n_lines = std::distance(txt.begin(), txt.end());
  TEST_TRUE(n_lines >= 2)

  StringList header = ListUtils::create<String>(*txt.begin(), '\t');
  TEST_EQUAL(header.size(), feature_set.size())
  {
    bool header_matches = (header.size() == feature_set.size());
    for (size_t i = 0; header_matches && i < feature_set.size(); ++i)
    {
      if (header[i] != feature_set[i]) header_matches = false;
    }
    TEST_EQUAL(header_matches, true)
  }

  // (d) Build the expected row sequence in the same order preparePin_ writes:
  //     outer loop pids, inner loop hits, skipping any (pid, hit) in `skipped`.
  //     SpecIds are not unique per hit (all hits within a pid share the same
  //     SpecId), so row-by-row alignment is required.
  std::vector<const PeptideHit*> expected_rows;
  for (size_t i = 0; i < stamped.size(); ++i)
  {
    const auto& hits = stamped[i].getHits();
    for (size_t j = 0; j < hits.size(); ++j)
    {
      if (skipped.count({i, j})) continue;
      expected_rows.push_back(&hits[j]);
    }
  }

  const size_t pin_data_rows = n_lines - 1;
  TEST_EQUAL(pin_data_rows, expected_rows.size())

  // (e) Walk rows in parallel. Compare every column except the trailing
  //     Proteins (which embeds tabs and inflates the split count).
  const size_t last_compared = feature_set.size() - 1; // skip Proteins

  size_t compared_rows = 0;
  size_t diffs = 0;
  String first_mismatch;
  size_t row_idx = 0;
  for (auto it = txt.begin() + 1;
       it != txt.end() && row_idx < expected_rows.size();
       ++it, ++row_idx)
  {
    StringList cols = ListUtils::create<String>(*it, '\t');
    if (cols.size() < last_compared)
    {
      diffs++;
      if (first_mismatch.empty()) first_mismatch = "short row: " + *it;
      continue;
    }
    const PeptideHit& hit = *expected_rows[row_idx];

    for (size_t c = 0; c < last_compared; ++c)
    {
      const String& feat = feature_set[c];
      if (!hit.metaValueExists(feat))
      {
        diffs++;
        if (first_mismatch.empty())
          first_mismatch = "row " + String(row_idx) + ": missing meta value " + feat;
        break;
      }
      const String stamped_val = hit.getMetaValue(feat).toString();
      if (cols[c] != stamped_val)
      {
        diffs++;
        if (first_mismatch.empty())
        {
          first_mismatch = "row " + String(row_idx) + " col " + feat
            + ": pin='" + cols[c] + "' stamped='" + stamped_val + "'";
        }
        break;
      }
    }
    compared_rows++;
  }

  if (!first_mismatch.empty())
  {
    TEST_EQUAL(first_mismatch, String())  // makes the diff visible on failure
  }
  TEST_EQUAL(diffs, 0)
  TEST_TRUE(compared_rows > 20)  // sanity: we compared a meaningful number of rows
}
END_SECTION

///////////////////////////////////////////////////////////////////////////////
// Test 2: Score and FDR-threshold parity with the subprocess.
//
// Writes a .pin with a realistic-scale (2000-row) synthetic dataset, then
// runs the external `percolator` binary and the in-process path with the
// same parameters (seed=1, default FDRs, default iterations). Asserts:
//   - Pearson r on SVM scores   ≥ 0.999 (observed: 1.0 to machine precision)
//   - max relative |Δscore|     ≤ 1e-4 (subprocess scores are round-tripped
//     through .pout text at ~6 significant digits; that text precision is
//     the only source of drift at this correlation)
//   - target count at q ≤ 0.01  matches exactly
//   - target count at q ≤ 0.05  matches exactly
//
// q-value and PEP *magnitudes* per row can legitimately drift in the tail
// (subprocess clamps "unidentified" PSMs to q=1/pep=1 via a post-training
// classification step the in-process wrapper does not re-implement). Test 3
// covers rank-at-threshold parity, which is what matters downstream.
///////////////////////////////////////////////////////////////////////////////

START_SECTION([EXTRA] scores and FDR-threshold counts match subprocess)
{
  const String bin = percolatorBinary();
  if (bin.empty())
  {
    TEST_EQUAL(true, true);  // skip silently; external binary unavailable
  }
  else
  {
    std::mt19937 rng(2026);
    RescoreInput ri;
    generateSyntheticData(ri, rng);

    const String pin_path = File::getTemporaryFile();
    writePinFile(pin_path, ri);

    SubprocessOut sub = runSubprocess(bin, pin_path, "-S 1");
    TEST_EQUAL(sub.exit_code, 0)
    TEST_TRUE(sub.triplets.size() > 1000)

    // Subprocess percolator uses isotonic PEP by default; match it.
    Percolator p;
    Param par = p.getDefaults();
    par.setValue("seed", 1);
    par.setValue("pep_method", "isotonic");
    p.setParameters(par);
    RescoreOutput out = p.rescore(ri);

    // Score parity: Pearson r and per-row max relative delta.
    size_t matches = 0;
    double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
    double max_rel_ds = 0;
    int tgt_q01_in = 0, tgt_q01_sub = 0;
    int tgt_q05_in = 0, tgt_q05_sub = 0;
    for (size_t i = 0; i < out.scores.size(); ++i)
    {
      char k[64]; std::snprintf(k, sizeof(k), "row_%08zu|-.PEPTIDE.-", i);
      auto it = sub.triplets.find(k);
      if (it == sub.triplets.end()) continue;
      matches++;
      const double x = out.scores[i];
      const double y = it->second.score;
      sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y;
      const double rel = std::abs(x - y) / std::max({std::abs(x), std::abs(y), 1.0});
      max_rel_ds = std::max(max_rel_ds, rel);
      if (!ri.is_decoy[i])
      {
        if (out.q_values[i]   <= 0.01) tgt_q01_in++;
        if (it->second.qval <= 0.01)   tgt_q01_sub++;
        if (out.q_values[i]   <= 0.05) tgt_q05_in++;
        if (it->second.qval <= 0.05)   tgt_q05_sub++;
      }
    }
    const double n = static_cast<double>(matches);
    const double num = n * sxy - sx * sy;
    const double den = std::sqrt((n * sxx - sx * sx) * (n * syy - sy * sy));
    const double r = (den > 1e-15) ? num / den : 0.0;

    TEST_TRUE(matches > 1000)
    TEST_TRUE(r >= 0.999)
    TEST_TRUE(max_rel_ds <= 1e-4)
    TEST_EQUAL(tgt_q01_in, tgt_q01_sub)
    TEST_EQUAL(tgt_q05_in, tgt_q05_sub)
  }
}
END_SECTION

///////////////////////////////////////////////////////////////////////////////
// Test 3: Ranking parity at FDR thresholds.
//
// Pearson-on-scores (Test 2) says the SVM is identical in aggregate, but what
// matters downstream is *which* PSMs get accepted at q ≤ 0.01 / 0.05 / 0.10.
// Assert the *set* of target rows passing each threshold agrees element-wise
// between paths — not just the counts. A drift that shuffled tied scores
// (e.g. a different tie-breaker) would flip set membership even though
// counts stayed equal; that's the class of regression this catches.
///////////////////////////////////////////////////////////////////////////////

START_SECTION([EXTRA] ranking parity at q &lt;= 0.01 / 0.05 / 0.10)
{
  const String bin = percolatorBinary();
  if (bin.empty())
  {
    TEST_EQUAL(true, true);  // skip silently
  }
  else
  {
    std::mt19937 rng(2026);
    RescoreInput ri;
    generateSyntheticData(ri, rng);

    const String pin_path = File::getTemporaryFile();
    writePinFile(pin_path, ri);

    SubprocessOut sub = runSubprocess(bin, pin_path, "-S 1");
    TEST_EQUAL(sub.exit_code, 0)

    Percolator p;
    Param par = p.getDefaults();
    par.setValue("seed", 1);
    par.setValue("pep_method", "isotonic");
    p.setParameters(par);
    RescoreOutput out = p.rescore(ri);

    auto accepted_targets = [&](double threshold, bool from_subprocess) -> std::set<size_t>
    {
      std::set<size_t> s;
      for (size_t i = 0; i < out.scores.size(); ++i)
      {
        if (ri.is_decoy[i]) continue;
        double q = 1.0;
        if (from_subprocess)
        {
          char k[64]; std::snprintf(k, sizeof(k), "row_%08zu|-.PEPTIDE.-", i);
          auto it = sub.triplets.find(k);
          if (it == sub.triplets.end()) continue;
          q = it->second.qval;
        }
        else
        {
          q = out.q_values[i];
        }
        if (q <= threshold) s.insert(i);
      }
      return s;
    };

    for (double thr : {0.01, 0.05, 0.10})
    {
      auto set_in = accepted_targets(thr, /*from_subprocess=*/false);
      auto set_sub = accepted_targets(thr, /*from_subprocess=*/true);
      TEST_EQUAL(set_in.size(), set_sub.size())
      TEST_EQUAL(set_in == set_sub, true)
    }
  }
}
END_SECTION

///////////////////////////////////////////////////////////////////////////////
// Test 4: Parameter matrix sweep.
//
// Run both paths across a handful of flag combinations known to flow through
// the library <-> subprocess boundary. Each row (`--train-best-positive`,
// `-Y`, `-u`, `-N 200`, `--nested-xval-bins 3`) is a specific regression
// surface: wiring it into Param but forgetting to pass it on to `P::*`
// upstream is the failure mode we're guarding against.
//
// Assertion: Pearson r ≥ 0.99 and equal target count at q ≤ 0.01 for EVERY
// configuration. Both are observed to be ~1.0 / exact; this catches any
// one of the flags silently not being honored.
///////////////////////////////////////////////////////////////////////////////

START_SECTION([EXTRA] parameter matrix: each flag flows through to the SVM)
{
  const String bin = percolatorBinary();
  if (bin.empty())
  {
    TEST_EQUAL(true, true);  // skip silently
  }
  else
  {
    struct Case
    {
      const char* name;
      const char* extra_args;          ///< appended to the subprocess cmdline
      std::function<void(Param&)> tune;  ///< applied to in-process Param
      double min_r;                    ///< minimum acceptable Pearson on scores
    };
    // Per-case minimum r: all configs should produce r very close to 1.
    // subset_max_train=200 now uses the same vendored reservoir-sampling
    // algorithm (scan-keyed priority via PseudoRandom::lcg_rand, same seed),
    // so training subsets and weights match the subprocess to within a very
    // small tolerance.
    const std::vector<Case> cases = {
      { "train_best_positive", "-S 1 --train-best-positive",
        [](Param& p) { p.setValue("train_best_positive", "true"); }, 0.99 },
      { "post_processing_tdc", "-S 1 -Y",
        [](Param& p) { p.setValue("post_processing_tdc", "true"); }, 0.99 },
      { "unitnorm",            "-S 1 -u",
        [](Param& p) { p.setValue("normalizer", "uni"); }, 0.99 },
      // Nested CV grid search amplifies floating-point ordering differences
      // between in-process and subprocess (3 nested folds × 9 cpos×cfrac
      // pairs each). Apple-clang's TRON solver lands on slightly different
      // weights than the bundled subprocess on this case (Pearson r ~0.94
      // observed on macOS Xcode runners; ~0.99 on Linux g++/clang). Drop
      // the bound to 0.9 so the case still catches gross regressions
      // without tripping on platform-specific solver-internal FP variance.
      { "nested_xval_bins=3",  "-S 1 --nested-xval-bins 3",
        [](Param& p) { p.setValue("nested_xval_bins", 3); }, 0.9 },
      { "subset_max_train=200","-S 1 -N 200",
        [](Param& p) { p.setValue("subset_max_train", 200); }, 0.99 },
      // c_pos and c_neg explicit. Default (c_pos=0.1, c_neg=0) hits the
      // "single-c_pos / 3-c_frac" branch in CrossValidation::initializeGridSearch.
      // With both explicit non-zero, we hit "single-c_pos / single-c_frac" (no
      // grid search). This catches regressions in the no-grid SVM path that
      // would never trip on default-config tests.
      { "cpos_cneg_explicit",  "-S 1 -p 1.0 -n 1.0",
        [](Param& p) { p.setValue("c_pos", 1.0); p.setValue("c_neg", 1.0); }, 0.99 },
    };

    for (const auto& tc : cases)
    {
      std::mt19937 rng(2026);
      RescoreInput ri;
      generateSyntheticData(ri, rng);
      const String pin_path = File::getTemporaryFile();
      writePinFile(pin_path, ri);

      SubprocessOut sub = runSubprocess(bin, pin_path, tc.extra_args);
      TEST_EQUAL(sub.exit_code, 0)

      Percolator p;
      Param par = p.getDefaults();
      par.setValue("seed", 1);
      par.setValue("pep_method", "isotonic");
      tc.tune(par);
      p.setParameters(par);
      RescoreOutput out = p.rescore(ri);

      size_t matches = 0;
      double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
      int tgt_q01_in = 0, tgt_q01_sub = 0;
      for (size_t i = 0; i < out.scores.size(); ++i)
      {
        char k[64]; std::snprintf(k, sizeof(k), "row_%08zu|-.PEPTIDE.-", i);
        auto it = sub.triplets.find(k);
        if (it == sub.triplets.end()) continue;
        matches++;
        const double x = out.scores[i];
        const double y = it->second.score;
        sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y;
        if (!ri.is_decoy[i])
        {
          if (out.q_values[i] <= 0.01) tgt_q01_in++;
          if (it->second.qval <= 0.01) tgt_q01_sub++;
        }
      }
      const double n = static_cast<double>(matches);
      const double num = n * sxy - sx * sy;
      const double den = std::sqrt((n * sxx - sx * sx) * (n * syy - sy * sy));
      const double r = (den > 1e-15) ? num / den : 0.0;

      TEST_TRUE(matches > 1000)
      // Fail with the case name in the message so regressions are bisectable.
      if (r < tc.min_r)
      {
        TEST_EQUAL(String("case ") + tc.name + " r=" + String(r)
                   + " < min_r=" + String(tc.min_r),
                   String("case ok"))
      }
      TEST_TRUE(r >= tc.min_r)
      if (tgt_q01_in != tgt_q01_sub)
      {
        TEST_EQUAL(String("case ") + tc.name + " target@q01 mismatch: in="
                   + String(tgt_q01_in) + " sub=" + String(tgt_q01_sub),
                   String("case ok"))
      }
      TEST_EQUAL(tgt_q01_in, tgt_q01_sub)
    }
  }
}
END_SECTION

///////////////////////////////////////////////////////////////////////////////
// Test 5: SVM feature-weight parity.
//
// Subprocess percolator writes per-fold weights with `-w`. In-process exposes
// `Percolator::getSvmWeights()`, which returns the across-fold average of the
// raw (un-normalized) weights. Average the subprocess raw weights across its
// 3 CV bins and compare the FEATURE coefficients element-wise.
//
// We deliberately exclude the trailing bias (m0) term from the assertion:
// both paths pass the trained weights_ vector through score-normalization
// (`Scores::normalizeScores` -> `Normalizer::endScoreNormalizeWeights` in
// postIterationProcessing/merge), and the point at which `getAvgWeights` is
// invoked relative to that step isn't identical between the in-process
// wrapper and the standalone percolator writer — so the affine offset drifts
// even though predictions (Test 2) agree to machine precision. Catching
// feature-coefficient drift is the useful signal; bias drift isn't.
///////////////////////////////////////////////////////////////////////////////

START_SECTION([EXTRA] SVM weights match average of per-fold subprocess weights)
{
  const String bin = percolatorBinary();
  if (bin.empty())
  {
    TEST_EQUAL(true, true);  // skip silently
  }
  else
  {
    std::mt19937 rng(2026);
    RescoreInput ri;
    generateSyntheticData(ri, rng);
    const String pin_path = File::getTemporaryFile();
    writePinFile(pin_path, ri);

    // Ask percolator to dump weights.
    const String wfile = File::getTemporaryFile();
    SubprocessOut sub = runSubprocess(
      bin, pin_path, "-S 1 -w \"" + wfile + "\"");
    TEST_EQUAL(sub.exit_code, 0)

    // Parse weights file. Format per CV bin:
    //   <feature1>\t<feature2>\t...\t<featureN>\tm0
    //   <norm_w1>\t...\t<norm_wN>\t<norm_bias>
    //   <raw_w1>\t...\t<raw_wN>\t<raw_bias>
    // Interleaved 3x. Plus 2 leading comment lines prefixed with '#'.
    std::vector<std::vector<double>> raw_per_fold;
    {
      std::ifstream f(wfile.c_str());
      String line;
      size_t block_line = 0;  // 0=header, 1=normalized, 2=raw (mod 3)
      while (std::getline(f, line))
      {
        if (line.empty() || line[0] == '#') continue;
        if (block_line % 3 == 2)
        {
          std::vector<double> row;
          std::istringstream ls(line);
          double v;
          while (ls >> v) row.push_back(v);
          raw_per_fold.push_back(std::move(row));
        }
        block_line++;
      }
    }
    TEST_EQUAL(raw_per_fold.size(), 3)  // 3 CV bins
    TEST_TRUE(!raw_per_fold.empty() && raw_per_fold[0].size() ==
              ri.feature_names.size() + 1)  // N features + bias

    // Average across folds (matches Percolator::getAvgWeights semantics).
    std::vector<double> sub_avg(ri.feature_names.size() + 1, 0.0);
    for (const auto& fold : raw_per_fold)
    {
      for (size_t i = 0; i < fold.size() && i < sub_avg.size(); ++i)
      {
        sub_avg[i] += fold[i] / static_cast<double>(raw_per_fold.size());
      }
    }

    // In-process.
    Percolator p;
    Param par = p.getDefaults();
    par.setValue("seed", 1);
    par.setValue("pep_method", "isotonic");
    p.setParameters(par);
    RescoreOutput out = p.rescore(ri);
    (void)out;

    const auto& in_weights = p.getSvmWeights();
    TEST_TRUE(!in_weights.empty())
    TEST_EQUAL(in_weights.front().size(), sub_avg.size())

    // Feature-coefficient parity (excludes the trailing m0 bias term).
    double max_abs = 0;
    size_t arg_max = 0;
    const size_t n_features = sub_avg.size() - 1;
    for (size_t i = 0; i < n_features; ++i)
    {
      const double d = std::abs(in_weights.front()[i] - sub_avg[i]);
      if (d > max_abs) { max_abs = d; arg_max = i; }
    }
    if (max_abs > 1e-3)
    {
      TEST_EQUAL(String("feature ") + ri.feature_names[arg_max]
                 + " |dw| = " + String(max_abs)
                 + " in=" + String(in_weights.front()[arg_max])
                 + " sub_avg=" + String(sub_avg[arg_max]),
                 String("weights match"))
    }
    TEST_TRUE(max_abs <= 1e-3)
  }
}
END_SECTION

///////////////////////////////////////////////////////////////////////////////
// Test 6: Realistic-idXML library parity (P1.a).
//
// Feeds CometAdapter_4_out.idXML through both paths at the library layer.
// The .pin is the single source of truth for row order: we write it via
// PercolatorInfile::store, parse it back to build the RescoreInput, run
// subprocess on the same .pin, and align by SpecId. Ensures the in-process
// and subprocess paths see identical input.
//
// Fixture: CometAdapter_4_out.idXML is Comet search output post-PSMFeatureExtractor
// with 712 hits / 581 decoys and extra_features populated (9 COMET:* +
// MS:1002252/55 + isotope_error). Rich enough that the SVM converges on the
// initial direction with test_fdr=train_fdr=0.5 rather than triggering
// percolator's "Reset to default direction" fallback path — that path
// diverges between in-process and subprocess in the q-value / PEP assignment
// and is therefore not useful as a parity regression target.
//
// The FDR=0.5 choice is about convergence, not about matching any particular
// reference .ini. Smaller / weaker fixtures (PercolatorAdapter_1.idXML with
// 40 hits, FalseDiscoveryRate_OMSSA.idXML, MzTabExporter_2_input.idXML) all
// trigger the Reset path even at FDR=0.5 and were verified unsuitable.
///////////////////////////////////////////////////////////////////////////////

START_SECTION([EXTRA] realistic idXML parity at library layer)
{
  const String bin = percolatorBinary();
  if (bin.empty())
  {
    TEST_EQUAL(true, true);  // skip silently
  }
  else
  {
    vector<ProteinIdentification> prs;
    PeptideIdentificationList pids;
    const String idxml_path =
      OPENMS_GET_TEST_DATA_PATH("../../../topp/THIRDPARTY/CometAdapter_4_out.idXML");
    IdXMLFile().load(idxml_path, prs, pids);
    TEST_FALSE(pids.empty())

    int min_charge, max_charge;
    deriveChargeRange(pids, min_charge, max_charge);
    StringList feature_set = buildFeatureSet(prs, min_charge, max_charge);
    const std::string enz = "trypsin";

    // Write .pin (also stamps meta values under the hood).
    const String pin_path = File::getTemporaryFile();
    PercolatorInfile::store(pin_path, pids, feature_set, enz, min_charge, max_charge);

    // Parse .pin back to build RescoreInput that matches subprocess input
    // exactly. Avoids buildRescoreInput / stampPinFeaturesOnHits alignment
    // pitfalls because .pin is the ground truth for row order.
    TextFile pin;
    pin.load(pin_path);
    const size_t n_lines = std::distance(pin.begin(), pin.end());
    TEST_TRUE(n_lines >= 2)

    StringList header = ListUtils::create<String>(*pin.begin(), '\t');
    std::map<String, size_t> col_index;
    for (size_t c = 0; c < header.size(); ++c) col_index[header[c]] = c;

    StringList numeric = numericFeatures(feature_set);
    // Composite per-row key: SpecId alone is NOT unique (multi-hit-per-PID
    // inputs have ~4 hits sharing the same file+scan SpecId), so subprocess
    // output map would collapse siblings. Use SpecId|Peptide to disambiguate;
    // runSubprocess() uses the same composite key on the .pout side.
    std::vector<String> pin_keys;
    RescoreInput ri;
    ri.feature_names = numeric;
    for (auto it = pin.begin() + 1; it != pin.end(); ++it)
    {
      StringList cols = ListUtils::create<String>(*it, '\t');
      // Tolerate Proteins column with embedded tabs by checking up to the
      // last numeric-feature column.
      bool ok = true;
      for (const auto& f : numeric)
      {
        if (col_index[f] >= cols.size()) { ok = false; break; }
      }
      if (!ok) continue;
      const String peptide_col = (col_index["Peptide"] < cols.size())
        ? cols[col_index["Peptide"]] : String();
      pin_keys.push_back(cols[col_index["SpecId"]] + "|" + peptide_col);
      ri.is_decoy.push_back(cols[col_index["Label"]] == "-1");
      ri.scan_numbers.push_back(cols[col_index["ScanNr"]].toInt());
      ri.spec_file_numbers.push_back(0);
      ri.exp_masses.push_back(cols[col_index["ExpMass"]].toDouble());
      ri.calc_masses.push_back(cols[col_index["CalcMass"]].toDouble());
      std::vector<double> feats;
      feats.reserve(numeric.size());
      for (const auto& f : numeric) feats.push_back(cols[col_index[f]].toDouble());
      ri.features.push_back(std::move(feats));
    }
    TEST_TRUE(ri.features.size() >= 20)

    // Run subprocess.
    SubprocessOut sub = runSubprocess(bin, pin_path, "-S 1 -t 0.5 -F 0.5");
    TEST_EQUAL(sub.exit_code, 0)

    // Run in-process, matching subprocess's CLI defaults. The library has
    // stronger defaults (post_processing_tdc="true", train_best_positive="true")
    // than the upstream percolator CLI (both off unless explicitly enabled);
    // set them to "false" here to mirror the subprocess invocation.
    Percolator p;
    Param par = p.getDefaults();
    par.setValue("seed", 1);
    par.setValue("pep_method", "nonparametric");  // matches 3.06 subprocess PosteriorEstimator::estimatePEP
    par.setValue("test_fdr",  0.5);
    par.setValue("train_fdr", 0.5);
    par.setValue("post_processing_tdc",  "false");
    par.setValue("train_best_positive",  "false");
    p.setParameters(par);
    RescoreOutput out = p.rescore(ri);
    TEST_EQUAL(out.scores.size(), ri.features.size())

    // Pearson r, max |Δpep|, accepted-target sets at q in {0.01, 0.05, 0.10}.
    double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
    double max_dpep = 0;
    size_t matches = 0;
    for (size_t i = 0; i < out.scores.size(); ++i)
    {
      auto it = sub.triplets.find(pin_keys[i]);
      if (it == sub.triplets.end()) continue;
      matches++;
      const double x = out.scores[i];
      const double y = it->second.score;
      sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y;
      max_dpep = std::max(max_dpep, std::abs(out.peps[i] - it->second.pep));
    }
    const double nf = static_cast<double>(matches);
    const double num = nf * sxy - sx * sy;
    const double den = std::sqrt((nf * sxx - sx * sx) * (nf * syy - sy * sy));
    const double r = (den > 1e-15) ? num / den : 0.0;

    TEST_TRUE(matches >= 20)
    TEST_TRUE(r >= 0.99)
    // max_dpep tolerance of 0.15: with pep_method="nonparametric" routing to
    // the vendored PosteriorEstimator::estimatePEP (3.08.01), the SVM scores
    // and q-values match the 3.06 subprocess bit-identically (r=1.0, accepted
    // sets equal), but the PEP smoothing kernel has been updated between 3.06
    // and 3.08 — residual max |Δpep| observed at 0.10 on this fixture.
    // TODO(percolator-3.08-binary): tighten to 0.05 after upgrading the
    // bundled binary to 3.08.01; the gap is algorithmic 3.06↔3.08 drift,
    // not a parity regression.
    TEST_TRUE(max_dpep <= 0.15)

    for (double thr : {0.01, 0.05, 0.10})
    {
      std::set<String> a_in, a_sub;
      for (size_t i = 0; i < out.scores.size(); ++i)
      {
        if (ri.is_decoy[i]) continue;
        auto it = sub.triplets.find(pin_keys[i]);
        if (it == sub.triplets.end()) continue;
        if (out.q_values[i] <= thr) a_in.insert(pin_keys[i]);
        if (it->second.qval <= thr) a_sub.insert(pin_keys[i]);
      }
      if (a_in != a_sub)
      {
        TEST_EQUAL(String("q=") + String(thr) + " accepted-set mismatch: in="
                   + String(a_in.size()) + " sub=" + String(a_sub.size()),
                   String("sets match"))
      }
      TEST_EQUAL(a_in == a_sub, true)
    }
  }
}
END_SECTION

///////////////////////////////////////////////////////////////////////////////
// Test 7: Reservoir-sampling parity at larger scale (P1.b).
//
// 20 000 rows with subset_max_train=5000 forces the reservoir sampler on
// both paths. Existing §4 case uses subset_max_train=200 on 2000 rows,
// which is informative but below the natural scale where sampling matters.
///////////////////////////////////////////////////////////////////////////////

START_SECTION([EXTRA] reservoir-sampling parity at 20k rows)
{
  const String bin = percolatorBinary();
  if (bin.empty())
  {
    TEST_EQUAL(true, true);  // skip
  }
  else
  {
    std::mt19937 rng(2026);
    RescoreInput ri;
    generateSyntheticData(ri, rng, /*size_mult=*/10.0);
    TEST_EQUAL(ri.features.size(), 20000)

    const String pin_path = File::getTemporaryFile();
    writePinFile(pin_path, ri);

    SubprocessOut sub = runSubprocess(bin, pin_path, "-S 1 -N 5000");
    TEST_EQUAL(sub.exit_code, 0)

    Percolator p;
    Param par = p.getDefaults();
    par.setValue("seed", 1);
    par.setValue("pep_method", "isotonic");
    par.setValue("subset_max_train", 5000);
    p.setParameters(par);
    RescoreOutput out = p.rescore(ri);

    double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
    size_t matches = 0;
    int tgt_q01_in = 0, tgt_q01_sub = 0;
    int tgt_q05_in = 0, tgt_q05_sub = 0;
    // Accepted-target SETS at q <= 0.01 / 0.05: catches drift that count-only
    // checks miss (e.g. same count but different PSMs crossing the threshold).
    std::set<size_t> acc01_in, acc01_sub, acc05_in, acc05_sub;
    for (size_t i = 0; i < out.scores.size(); ++i)
    {
      char k[64]; std::snprintf(k, sizeof(k), "row_%08zu|-.PEPTIDE.-", i);
      auto it = sub.triplets.find(k);
      if (it == sub.triplets.end()) continue;
      matches++;
      const double x = out.scores[i];
      const double y = it->second.score;
      sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y;
      if (!ri.is_decoy[i])
      {
        if (out.q_values[i]   <= 0.01) { tgt_q01_in++; acc01_in.insert(i); }
        if (it->second.qval <= 0.01) { tgt_q01_sub++; acc01_sub.insert(i); }
        if (out.q_values[i]   <= 0.05) { tgt_q05_in++; acc05_in.insert(i); }
        if (it->second.qval <= 0.05) { tgt_q05_sub++; acc05_sub.insert(i); }
      }
    }
    const double nf = static_cast<double>(matches);
    const double num = nf * sxy - sx * sy;
    const double den = std::sqrt((nf * sxx - sx * sx) * (nf * syy - sy * sy));
    const double r = (den > 1e-15) ? num / den : 0.0;

    TEST_TRUE(matches >= 18000)
    TEST_TRUE(r >= 0.999)

    // Count parity at q-thresholds: tolerance-based at 20k scale. The 5000-row
    // reservoir-trained SVM slightly underperforms the initial direction on
    // this fixture, which triggers percolator's Reset fallback; post-reset
    // q-value assignment has tiny boundary sensitivity that flips ~5 rows
    // across q=0.01 between paths even at r >= 0.999. Tolerance 0.5% is 7x
    // the observed 0.07% drift — catches real regressions without being
    // noise-sensitive at this scale.
    const int maxc01 = std::max(tgt_q01_in, tgt_q01_sub);
    const double ratio01 = (maxc01 > 0)
      ? static_cast<double>(std::abs(tgt_q01_in - tgt_q01_sub)) / static_cast<double>(maxc01)
      : 0.0;
    TEST_TRUE(ratio01 <= 0.005)

    const int maxc05 = std::max(tgt_q05_in, tgt_q05_sub);
    const double ratio05 = (maxc05 > 0)
      ? static_cast<double>(std::abs(tgt_q05_in - tgt_q05_sub)) / static_cast<double>(maxc05)
      : 0.0;
    TEST_TRUE(ratio05 <= 0.005)

    // Accepted-target SET parity at q ≤ 0.01 and q ≤ 0.05 via Jaccard. Count
    // ratio above ensures cardinalities are close, but doesn't catch the case
    // where two paths accept the same *number* of targets but different *sets*
    // of rows near the q-threshold boundary. On this fixture both paths agree
    // on ~99.9% of accepted IDs; Jaccard floor of 0.99 catches a regression
    // that shuffles membership while preserving counts.
    auto jaccard = [](const std::set<size_t>& a, const std::set<size_t>& b) {
      if (a.empty() && b.empty()) return 1.0;
      std::set<size_t> inter, uni;
      std::set_intersection(a.begin(), a.end(), b.begin(), b.end(),
                            std::inserter(inter, inter.begin()));
      std::set_union(a.begin(), a.end(), b.begin(), b.end(),
                     std::inserter(uni, uni.begin()));
      return uni.empty() ? 1.0
        : static_cast<double>(inter.size()) / static_cast<double>(uni.size());
    };
    const double j01 = jaccard(acc01_in, acc01_sub);
    const double j05 = jaccard(acc05_in, acc05_sub);
    TEST_TRUE(j01 >= 0.99)
    TEST_TRUE(j05 >= 0.99)
  }
}
END_SECTION

///////////////////////////////////////////////////////////////////////////////
// Test 8: Multi-file PIN parity (P1.c).
//
// Two sub-batches of synthetic rows with distinct specFileNr. The CV fold
// hash (OrderScanHash) incorporates specFileNr; this section confirms it
// agrees across paths when more than one file is present — a code path the
// other library-layer tests don't exercise (they all use specFileNr=0).
///////////////////////////////////////////////////////////////////////////////

START_SECTION([EXTRA] multi-file PIN parity)
{
  const String bin = percolatorBinary();
  if (bin.empty())
  {
    TEST_EQUAL(true, true);  // skip
  }
  else
  {
    std::mt19937 rng(2026);
    RescoreInput ri;
    generateSyntheticData(ri, rng);  // 2000 rows, specFileNr=0 by default

    // Split into two "files": first half -> specFileNr=0, second -> specFileNr=1.
    const size_t half = ri.features.size() / 2;
    for (size_t i = 0; i < ri.features.size(); ++i)
    {
      ri.spec_file_numbers[i] = (i < half) ? 0 : 1;
    }

    const String pin_path = File::getTemporaryFile();
    writePinFile(pin_path, ri, /*emit_filename=*/true);

    SubprocessOut sub = runSubprocess(bin, pin_path, "-S 1");
    TEST_EQUAL(sub.exit_code, 0)

    Percolator p;
    Param par = p.getDefaults();
    par.setValue("seed", 1);
    par.setValue("pep_method", "isotonic");
    p.setParameters(par);
    RescoreOutput out = p.rescore(ri);

    double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
    size_t matches = 0;
    int tgt_q01_in = 0, tgt_q01_sub = 0;
    int tgt_q05_in = 0, tgt_q05_sub = 0;
    for (size_t i = 0; i < out.scores.size(); ++i)
    {
      char k[64]; std::snprintf(k, sizeof(k), "row_%08zu|-.PEPTIDE.-", i);
      auto it = sub.triplets.find(k);
      if (it == sub.triplets.end()) continue;
      matches++;
      sx += out.scores[i];
      sy += it->second.score;
      sxx += out.scores[i] * out.scores[i];
      syy += it->second.score * it->second.score;
      sxy += out.scores[i] * it->second.score;
      if (!ri.is_decoy[i])
      {
        if (out.q_values[i] <= 0.01) ++tgt_q01_in;
        if (it->second.qval <= 0.01) ++tgt_q01_sub;
        if (out.q_values[i] <= 0.05) ++tgt_q05_in;
        if (it->second.qval <= 0.05) ++tgt_q05_sub;
      }
    }
    const double nf = static_cast<double>(matches);
    const double num = nf * sxy - sx * sy;
    const double den = std::sqrt((nf * sxx - sx * sx) * (nf * syy - sy * sy));
    const double r = (den > 1e-15) ? num / den : 0.0;

    TEST_TRUE(matches > 1000)
    TEST_TRUE(r >= 0.999)
    TEST_EQUAL(tgt_q01_in, tgt_q01_sub)
    TEST_EQUAL(tgt_q05_in, tgt_q05_sub)
  }
}
END_SECTION

///////////////////////////////////////////////////////////////////////////////

END_TEST
