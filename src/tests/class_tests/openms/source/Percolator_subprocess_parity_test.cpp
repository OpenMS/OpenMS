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
    "SpecId", "Label", "ScanNr", "ExpMass", "Peptide", "Proteins"
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
void generateSyntheticData(RescoreInput& ri, std::mt19937& rng)
{
  ri.features.clear();
  ri.is_decoy.clear();
  ri.feature_names = StringList{"f0", "f1", "noise"};

  std::normal_distribution<double> noise(0.0, 0.5);
  std::normal_distribution<double> unit(0.0, 1.0);

  const size_t n_easy = 500, n_hard = 500, n_dec = 1000;
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

  // scan_numbers power the CV fold split via Percolator's scan-hash. The hash
  // also includes exp_mass/calc_mass (see OrderScanHash), so these must match
  // the values written in the .pin or folds will diverge between paths.
  const size_t n = ri.features.size();
  ri.scan_numbers.assign(n, 0);
  ri.exp_masses.assign(n, 1000.0);
  ri.calc_masses.assign(n, 1000.0);
  for (size_t i = 0; i < n; ++i)
  {
    ri.scan_numbers[i] = static_cast<int>(i + 1);
  }
}

/// Write a minimal .pin file that percolator's command-line tool will accept.
/// Columns: SpecId, Label, ScanNr, ExpMass, CalcMass, <features>, Peptide, Proteins.
/// Row ids `row_%08zu` so subprocess PSMId strings match in-process row indices.
void writePinFile(const String& path, const RescoreInput& ri)
{
  std::ofstream f(path.c_str());
  if (!f) throw std::runtime_error("cannot open pin file for writing");
  f.setf(std::ios::scientific);
  f.precision(17);

  f << "SpecId\tLabel\tScanNr\tExpMass\tCalcMass";
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
  s.exit_code = std::system(cmd.str().c_str());
  if (s.exit_code != 0) return s;

  auto parse = [&s](const String& path)
  {
    std::ifstream f(path.c_str());
    String line;
    bool header_seen = false;
    while (std::getline(f, line))
    {
      if (!header_seen) { header_seen = true; continue; }
      if (line.empty()) continue;
      // Fields: PSMId  score  q-value  posterior_error_prob  peptide  protein...
      std::istringstream ls(line);
      String psm_id;
      double score = 0, qval = 0, pep = 0;
      ls >> psm_id >> score >> qval >> pep;
      if (!ls) continue;
      PercolatorTriplet t;
      t.score = score;
      t.qval = qval;
      t.pep = pep;
      s.triplets[psm_id] = t;
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
    char row_id[32];
    std::snprintf(row_id, sizeof(row_id), "row_%08zu", i);
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
      char k[32]; std::snprintf(k, sizeof(k), "row_%08zu", i);
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
          char k[32]; std::snprintf(k, sizeof(k), "row_%08zu", i);
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
    // Per-case minimum r: most configs produce r=1 to machine precision.
    // subset_max_train=200 does random sub-sampling that uses independent RNG
    // streams between paths, so score magnitudes drift even at fixed seed —
    // the target count at q<=0.01 is still the right invariant to check there.
    const std::vector<Case> cases = {
      { "train_best_positive", "-S 1 --train-best-positive",
        [](Param& p) { p.setValue("train_best_positive", "true"); }, 0.99 },
      { "post_processing_tdc", "-S 1 -Y",
        [](Param& p) { p.setValue("post_processing_tdc", "true"); }, 0.99 },
      { "unitnorm",            "-S 1 -u",
        [](Param& p) { p.setValue("normalizer", "uni"); }, 0.99 },
      { "nested_xval_bins=3",  "-S 1 --nested-xval-bins 3",
        [](Param& p) { p.setValue("nested_xval_bins", 3); }, 0.99 },
      { "subset_max_train=200","-S 1 -N 200",
        [](Param& p) { p.setValue("subset_max_train", 200); }, 0.60 },
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
        char k[32]; std::snprintf(k, sizeof(k), "row_%08zu", i);
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

END_TEST
