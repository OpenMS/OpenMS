// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
//
// Side-by-side parity test for the PercolatorAdapter TOPP tool: runs it twice
// on the same idXML with -use_subprocess true/false, aligns output hits by
// (spectrum_reference, sequence, charge), and compares the stamped Percolator
// outputs (score, q-value, PEP).
//
// Gated on two env vars set by CMake:
//   PERCOLATOR_BINARY         -> forwarded as -percolator_executable
//   PERCOLATOR_ADAPTER_BINARY -> absolute path to bin/PercolatorAdapter
// Skips silently when either is unset.

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace OpenMS;
using namespace std;

namespace
{

/// Composite key for 1:1 cross-run matching of adapter output hits.
String rowKey(const PeptideIdentification& pid, const PeptideHit& hit)
{
  String sr = String(pid.getSpectrumReference());
  if (sr.empty() && pid.metaValueExists("spectrum_id"))
  {
    sr = pid.getMetaValue("spectrum_id").toString();
  }
  return sr + "|" + hit.getSequence().toString() + "|" + String(hit.getCharge());
}

struct OutputTriplet
{
  double score = 0.0;
  double qval  = 0.0;
  double pep   = 0.0;
  bool   is_decoy = false;
  String key_for_score = "";  ///< meta key actually found on the hit
};

/// Extract the Percolator triplet from a hit, probing both meta-key schemes.
OutputTriplet extractTriplet(const PeptideIdentification& pid, const PeptideHit& hit)
{
  OutputTriplet t;
  t.is_decoy = hit.isDecoy();
  if (hit.metaValueExists("MS:1001492"))
  {
    t.score = static_cast<double>(hit.getMetaValue("MS:1001492"));
    t.qval  = hit.metaValueExists("MS:1001491")
      ? static_cast<double>(hit.getMetaValue("MS:1001491")) : 0.0;
    t.pep   = hit.metaValueExists("MS:1001493")
      ? static_cast<double>(hit.getMetaValue("MS:1001493")) : 0.0;
    t.key_for_score = "MS:1001492";
  }
  else if (hit.metaValueExists("percolator_score"))
  {
    t.score = static_cast<double>(hit.getMetaValue("percolator_score"));
    t.qval  = hit.metaValueExists("percolator_q_value")
      ? static_cast<double>(hit.getMetaValue("percolator_q_value")) : 0.0;
    t.pep   = hit.metaValueExists("percolator_pep")
      ? static_cast<double>(hit.getMetaValue("percolator_pep")) : 0.0;
    t.key_for_score = "percolator_score";
  }
  (void)pid;
  return t;
}

/// Invoke PercolatorAdapter on @p in_idxml, writing @p out_idxml.
/// Passes extra_args verbatim (e.g. "-testFDR 0.5 -trainFDR 0.5") so the
/// caller can match the FDR settings used by the TOPP_PercolatorAdapter_1 test.
int runAdapter(const String& adapter_bin,
               const String& percolator_bin,
               bool use_subprocess,
               const String& in_idxml,
               const String& out_idxml,
               const std::string& extra_args = "")
{
  const String stdout_log = File::getTemporaryFile();
  const String stderr_log = File::getTemporaryFile();

  std::ostringstream cmd;
  cmd << "\"" << adapter_bin << "\""
      << " -test"
      << " -in \""           << in_idxml     << "\""
      << " -out \""          << out_idxml    << "\""
      << " -out_type idXML"
      << " -percolator_executable \"" << percolator_bin << "\""
      << " -use_subprocess " << (use_subprocess ? "true" : "false");
  if (!extra_args.empty()) cmd << " " << extra_args;
  cmd << " > \"" << stdout_log << "\""
      << " 2> \"" << stderr_log << "\"";
  // cmd.exe /c strips the first and last quote of its argument. Since this
  // command both starts AND ends with a quoted path, both get stripped, which
  // mangles the binary path and the redirection target. Wrap the whole thing
  // in an extra outer pair of quotes so cmd.exe strips only those.
#ifdef _WIN32
  const std::string full_cmd = "\"" + cmd.str() + "\"";
  return std::system(full_cmd.c_str());
#else
  return std::system(cmd.str().c_str());
#endif
}

/// Build a map[rowKey -> triplet] from an adapter output idXML.
std::map<String, OutputTriplet> loadTripletsByRowKey(const String& idxml)
{
  vector<ProteinIdentification> prs;
  PeptideIdentificationList pids;
  IdXMLFile().load(idxml, prs, pids);
  std::map<String, OutputTriplet> out;
  for (const auto& pid : pids)
  {
    for (const auto& hit : pid.getHits())
    {
      out[rowKey(pid, hit)] = extractTriplet(pid, hit);
    }
  }
  return out;
}

} // namespace

///////////////////////////////////////////////////////////////////////////////

START_TEST(PercolatorAdapter_parity, "$Id$")

START_SECTION([EXTRA] adapter parity: -use_subprocess true vs false on same idXML)
{
  const char* perc = std::getenv("PERCOLATOR_BINARY");
  const char* adap = std::getenv("PERCOLATOR_ADAPTER_BINARY");
  if (!perc || !*perc || !adap || !*adap)
  {
    TEST_EQUAL(true, true);  // skip silently
  }
  else
  {
    const String percolator_bin = String(perc);
    const String adapter_bin    = String(adap);
    const String in_idxml =
      OPENMS_GET_TEST_DATA_PATH("../../../topp/THIRDPARTY/CometAdapter_4_out.idXML");

    const String out_sub = File::getTemporaryFile() + ".idxml";
    const String out_inp = File::getTemporaryFile() + ".idxml";

    TEST_EQUAL(runAdapter(adapter_bin, percolator_bin, /*use_subprocess=*/true,
                          in_idxml, out_sub,
                          "-testFDR 0.5 -trainFDR 0.5"), 0)
    TEST_EQUAL(runAdapter(adapter_bin, percolator_bin, /*use_subprocess=*/false,
                          in_idxml, out_inp,
                          "-testFDR 0.5 -trainFDR 0.5"), 0)

    auto sub = loadTripletsByRowKey(out_sub);
    auto inp = loadTripletsByRowKey(out_inp);

    TEST_TRUE(!sub.empty())
    TEST_TRUE(!inp.empty())

    // Prelude: both paths must stamp the same meta-key scheme. If one stamps
    // MS:1001492 and the other stamps percolator_score, comparison would
    // silently read zeros and pass spuriously.
    std::set<String> sub_keys, inp_keys;
    for (const auto& kv : sub) if (!kv.second.key_for_score.empty())
      sub_keys.insert(kv.second.key_for_score);
    for (const auto& kv : inp) if (!kv.second.key_for_score.empty())
      inp_keys.insert(kv.second.key_for_score);
    if (sub_keys != inp_keys)
    {
      std::ostringstream msg;
      msg << "meta-key mismatch: sub={ ";
      for (const auto& k : sub_keys) msg << k << " ";
      msg << "} vs inp={ ";
      for (const auto& k : inp_keys) msg << k << " ";
      msg << "}";
      TEST_EQUAL(String(msg.str()), String("meta keys match"))
    }
    TEST_EQUAL(sub_keys == inp_keys, true)

    // Same surviving hit set (symmetric difference = empty).
    std::set<String> only_sub, only_inp;
    for (const auto& kv : sub) if (!inp.count(kv.first)) only_sub.insert(kv.first);
    for (const auto& kv : inp) if (!sub.count(kv.first)) only_inp.insert(kv.first);
    TEST_EQUAL(only_sub.size(), 0)
    TEST_EQUAL(only_inp.size(), 0)

    // Pearson r, max |Δpep|, accepted-target sets at q in {0.01, 0.05, 0.10}.
    double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
    double max_dpep = 0;
    size_t matches = 0;
    struct Mismatch { String key; double s_sub, s_inp, q_sub, q_inp, p_sub, p_inp; };
    std::vector<Mismatch> mismatches;

    for (const auto& kv : inp)
    {
      auto it = sub.find(kv.first);
      if (it == sub.end()) continue;
      matches++;
      const double x = kv.second.score;
      const double y = it->second.score;
      sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y;
      const double dp = std::abs(kv.second.pep - it->second.pep);
      max_dpep = std::max(max_dpep, dp);
      if (std::abs(x - y) > 0.5 * std::max({std::abs(x), std::abs(y), 1.0})
          && mismatches.size() < 5)
      {
        mismatches.push_back({kv.first, y, x,
                              it->second.qval, kv.second.qval,
                              it->second.pep,  kv.second.pep});
      }
    }
    const double nf = static_cast<double>(matches);
    const double num = nf * sxy - sx * sy;
    const double den = std::sqrt((nf * sxx - sx * sx) * (nf * syy - sy * sy));
    const double r = (den > 1e-15) ? num / den : 0.0;

    TEST_TRUE(matches >= 20)
    if (r < 0.99)
    {
      std::ostringstream msg;
      msg << "r=" << r << " (expected >= 0.99); first mismatches: ";
      for (const auto& m : mismatches)
      {
        msg << "[" << m.key << " s_sub=" << m.s_sub << " s_inp=" << m.s_inp
            << " q_sub=" << m.q_sub << " q_inp=" << m.q_inp
            << " p_sub=" << m.p_sub << " p_inp=" << m.p_inp << "] ";
      }
      TEST_EQUAL(String(msg.str()), String("r ok"))
    }
    TEST_TRUE(r >= 0.99)  // TODO(percolator-3.08): tighten to 0.999 after upgrade
    // max_dpep tolerance 0.25 at adapter layer: SVM scores match (r=1.0),
    // accepted-target sets match, but the PEP smoothing kernel has been
    // updated between the bundled 3.06 binary and our vendored 3.08.01
    // PosteriorEstimator::estimatePEP — residual max |Δpep| ≈ 0.10 at the
    // library layer (see Percolator_subprocess_parity_test §6) and ≈ 0.23 at
    // the adapter layer (where adapter-side post-processing amplifies it).
    // Algorithmic 3.06↔3.08 drift, not a parity regression.
    // TODO(percolator-3.08-binary): tighten to 0.05 after upgrading the
    // bundled binary.
    TEST_TRUE(max_dpep <= 0.25)

    // Accepted-target set equality at q in {0.01, 0.05, 0.10}.
    for (double thr : {0.01, 0.05, 0.10})
    {
      std::set<String> acc_sub, acc_inp;
      for (const auto& kv : inp)
      {
        if (kv.second.is_decoy) continue;
        auto it = sub.find(kv.first);
        if (it == sub.end()) continue;
        if (kv.second.qval <= thr)  acc_inp.insert(kv.first);
        if (it->second.qval <= thr) acc_sub.insert(kv.first);
      }
      if (acc_sub != acc_inp)
      {
        TEST_EQUAL(String("q=") + String(thr) + " accepted-set mismatch: sub="
                   + String(acc_sub.size()) + " inp=" + String(acc_inp.size()),
                   String("sets match"))
      }
      TEST_EQUAL(acc_sub == acc_inp, true)
    }
  }
}
END_SECTION

END_TEST
