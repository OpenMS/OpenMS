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
#include <OpenMS/ANALYSIS/MAPMATCHING/PipEchoAlgorithm.h>
///////////////////////////

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>

using namespace OpenMS;
using namespace std;

// --------------------------------------------------------------------------
// A "realistic" PIP-ECHO match-between-runs (MBR) test.
//
// The hand-built PipEchoAlgorithm_test exercises only the two degenerate
// regimes (the fdr=1.0 escape hatch and the 0-decoy conservative drop). It
// never drives the scientific core: the SVM-based PEP estimation and the
// target/decoy transfer-FDR filtering, which only run once there are
// >= 100 acceptor transfers and >= 20 decoy transfers.
//
// Here we synthesise a ground-truth LC-MS experiment (~10k features, 3
// pre-aligned runs) and run group() on it. Because we *know* which features
// are true peptide signals and which are noise, we can verify the property
// that matters and that no other test checks: the realised transfer FDR
// respects the requested cutoff, while a useful fraction of the true
// transfers is recovered.
//
// SVM (LIBSVM) output is not bit-reproducible across platforms, so we assert
// on structure and statistics (counts, monotonicity, an FDR bound, a
// sensitivity floor, and same-seed determinism) rather than on a reference
// file or exact q-values.
//
// The full-scale run drives a LIBSVM grid search over thousands of acceptors
// across ten rounds, so it is slow (seconds to minutes). It is therefore
// gated behind the OPENMS_RUN_SLOW_TESTS environment variable and skips
// silently when unset, so it does not bloat the default CI run. Enable with:
//   OPENMS_RUN_SLOW_TESTS=1 ./PipEchoAlgorithm_realistic_test
// --------------------------------------------------------------------------

namespace
{
  // Ground truth attached to every generated feature, keyed by its (globally
  // unique) feature id. A ConsensusFeature's FeatureHandle preserves the
  // unique id but NOT meta values, so the unique id is the only reliable way
  // to map an output handle back to the truth.
  struct Truth
  {
    std::string seq;  // peptide sequence ("" for pure noise features)
    Int charge {0};
    bool is_donor {false}; // identified (carries a PeptideIdentification)
    bool noise {false};    // random feature with no underlying peptide
  };

  struct Dataset
  {
    std::vector<FeatureMap> maps;
    std::map<UInt64, Truth> truth;
    // Bookkeeping for diagnostics / sensitivity denominator.
    std::size_t n_donor_features {};
    std::size_t n_true_acceptors {};   // present-but-unidentified true peptides
    std::size_t n_noise_features {};
    std::size_t n_recoverable {};      // true acceptors whose peptide has a
                                       // donor in a *different* run
  };

  // A master peptide: a sequence, a fixed charge, a true RT and a base
  // (log2) intensity. Its theoretical m/z is derived from sequence+charge
  // exactly as the algorithm does, so identified donors have ~zero mass
  // error (plus the small ppm jitter we add per observation).
  struct MasterPeptide
  {
    std::string seq;
    Int charge;
    double rt;
    double theo_mz;
    double log2_intensity;
  };

  Feature makeFeatureNoId(double rt, double mz, Int charge, double intensity, UInt64 uid)
  {
    Feature f;
    f.setRT(rt);
    f.setMZ(mz);
    f.setCharge(charge);
    f.setIntensity(intensity);
    f.setUniqueId(uid);
    return f;
  }

  Feature makeDonorFeature(double rt, double mz, Int charge, double intensity, UInt64 uid,
                           const AASequence& seq)
  {
    Feature f = makeFeatureNoId(rt, mz, charge, intensity, uid);
    PeptideHit hit;
    hit.setSequence(seq);
    hit.setCharge(charge);
    hit.setTargetDecoyType(PeptideHit::TargetDecoyType::TARGET);
    PeptideIdentification pid;
    pid.setHits({hit});
    f.setPeptideIdentifications({pid});
    return f;
  }

  // Build ~target_total features across `n_runs` pre-aligned runs.
  Dataset buildSyntheticDataset(uint32_t seed,
                                std::size_t n_peptides,
                                std::size_t n_runs,
                                std::size_t target_total,
                                double p_present,
                                double p_ident)
  {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    // A moderately dense RT span and a tight m/z band (mimicking co-eluting
    // tryptic peptides) so the grid neighbourhood actually contains the
    // mass-similar partners that MBR decoy generation needs (see below).
    std::uniform_real_distribution<double> u_rt(300.0, 2300.0);
    std::uniform_real_distribution<double> u_noise_mz(380.0, 700.0);
    std::uniform_int_distribution<int> u_len(9, 20);
    std::normal_distribution<double> rt_jitter(0.0, 1.5);      // aligned runs
    std::normal_distribution<double> ppm_jitter(0.0, 2.0);     // m/z, ppm
    std::normal_distribution<double> log2_int(20.0, 1.8);      // ~1e6 median
    std::normal_distribution<double> int_mult(0.0, 0.3);       // per-obs (log2)

    const std::string aas = "ACDEFGHIKLMNPQRSTVWY";
    // Runs are pre-aligned: only a tiny systematic shift remains.
    std::vector<double> run_shift(n_runs, 0.0);
    for (std::size_t r = 0; r < n_runs; ++r) { run_shift[r] = (double(r) - double(n_runs) / 2.0) * 2.0; }

    auto draw_charge = [&]() -> Int {
      double u = u01(rng);
      if (u < 0.60) return 2;
      if (u < 0.90) return 3;
      return 4;
    };

    // 1) Master peptide list (distinct sequences).
    std::vector<MasterPeptide> peptides;
    peptides.reserve(n_peptides);
    std::set<std::string> seen_seq;
    while (peptides.size() < n_peptides)
    {
      int len = u_len(rng);
      std::string s;
      s.reserve(len);
      for (int i = 0; i < len - 1; ++i) { s.push_back(aas[std::uniform_int_distribution<size_t>(0, aas.size() - 1)(rng)]); }
      s.push_back(u01(rng) < 0.5 ? 'K' : 'R'); // tryptic-ish C-terminus
      if (!seen_seq.insert(s).second) continue; // ensure uniqueness

      // Charge follows mass (~1 charge per 500 Da, clamped to 2..4) as for real
      // tryptic peptides. This keeps m/z in a tight ~[380,750] band AND yields a
      // realistic 3+/4+ majority -- decoy pairs need |dmass| in [5,11] Da within
      // the grid's +-2 Da m/z reach, i.e. |dm/z| <= 2 Da, which is only reachable
      // at charge >= 3 (a 2+ pair 5-11 Da apart is >= 2.5 Da in m/z).
      AASequence aa = AASequence::fromString(s);
      double mass = aa.getMonoWeight();
      Int charge = static_cast<Int>(std::lround(mass / 500.0));
      charge = std::min<Int>(4, std::max<Int>(2, charge));

      MasterPeptide mp;
      mp.seq = s;
      mp.charge = charge;
      mp.rt = u_rt(rng);
      mp.theo_mz = aa.getMZ(charge);
      mp.log2_intensity = log2_int(rng);
      peptides.push_back(mp);
    }

    Dataset ds;
    ds.maps.resize(n_runs);
    for (std::size_t r = 0; r < n_runs; ++r)
    {
      ds.maps[r].setPrimaryMSRunPath({"run" + std::to_string(r) + ".mzML"});
    }

    UInt64 uid = 1;

    // donor_runs[seq] = set of runs in which this peptide is identified.
    std::map<std::string, std::set<std::size_t>> donor_runs;
    // For the recoverable denominator we remember each true acceptor's
    // (uid, seq, run) so we can resolve "has a donor elsewhere" afterwards.
    struct AcceptorRec { UInt64 uid; std::string seq; std::size_t run; };
    std::vector<AcceptorRec> true_acceptor_recs;

    // 2) Per run, per peptide: maybe present, maybe identified.
    for (const auto& mp : peptides)
    {
      const AASequence aaseq = AASequence::fromString(mp.seq);
      for (std::size_t r = 0; r < n_runs; ++r)
      {
        if (u01(rng) >= p_present) continue; // not detected in this run

        double rt = mp.rt + run_shift[r] + rt_jitter(rng);
        double mz = mp.theo_mz * (1.0 + ppm_jitter(rng) * 1e-6);
        double intensity = std::pow(2.0, mp.log2_intensity + int_mult(rng));

        const bool identified = u01(rng) < p_ident;
        Truth t;
        t.seq = mp.seq;
        t.charge = mp.charge;
        t.is_donor = identified;
        t.noise = false;

        if (identified)
        {
          ds.maps[r].push_back(makeDonorFeature(rt, mz, mp.charge, intensity, uid, aaseq));
          donor_runs[mp.seq].insert(r);
          ++ds.n_donor_features;
        }
        else
        {
          ds.maps[r].push_back(makeFeatureNoId(rt, mz, mp.charge, intensity, uid));
          true_acceptor_recs.push_back({uid, mp.seq, r});
          ++ds.n_true_acceptors;
        }
        ds.truth[uid] = t;
        ++uid;
      }
    }

    // 3) Noise: unidentified features with no underlying peptide. They raise
    //    acceptor density and are the population false transfers land on.
    std::size_t current = ds.n_donor_features + ds.n_true_acceptors;
    std::size_t n_noise = (target_total > current) ? (target_total - current) : 0;
    std::uniform_int_distribution<size_t> pick_pep(0, peptides.size() - 1);
    for (std::size_t i = 0; i < n_noise; ++i)
    {
      std::size_t r = std::uniform_int_distribution<size_t>(0, n_runs - 1)(rng);
      double rt = u_rt(rng);
      double mz;
      Int charge;
      // A realistic fraction of background ions interfere at a real peptide's
      // m/z but at an unrelated RT. These are the main fuel for randomized-RT
      // decoys and for FALSE transfers, so the transfer FDR has something to
      // estimate. They carry no true peptide, so any transfer onto them is wrong.
      if (u01(rng) < 0.6)
      {
        const MasterPeptide& mp = peptides[pick_pep(rng)];
        mz = mp.theo_mz * (1.0 + ppm_jitter(rng) * 1e-6);
        charge = mp.charge;
      }
      else
      {
        mz = u_noise_mz(rng);
        charge = draw_charge();
      }
      double intensity = std::pow(2.0, log2_int(rng) + int_mult(rng));
      ds.maps[r].push_back(makeFeatureNoId(rt, mz, charge, intensity, uid));
      Truth t;
      t.seq = "";        // no true peptide here -> any transfer is a false transfer
      t.charge = charge;
      t.is_donor = false;
      t.noise = true;
      ds.truth[uid] = t;
      ++uid;
      ++ds.n_noise_features;
    }

    // 4) Recoverable denominator: true acceptors whose peptide is a donor in
    //    a different run (an MBR transfer is at least possible).
    for (const auto& rec : true_acceptor_recs)
    {
      auto it = donor_runs.find(rec.seq);
      if (it == donor_runs.end()) continue;
      for (std::size_t dr : it->second)
      {
        if (dr != rec.run) { ++ds.n_recoverable; break; }
      }
    }

    return ds;
  }

  // Result of classifying the transfers in a consensus map against truth.
  struct TransferStats
  {
    std::size_t n_consensus {};
    std::size_t n_groups_with_transfer {};
    std::size_t n_transfers {};        // total transferred acceptor handles
    std::size_t n_correct {};          // transferred id matches the acceptor's true peptide
    std::size_t n_false {};            // wrong peptide or noise
    std::size_t n_meta_listed {};      // sum of mbr_transfer_map_indices sizes
    std::set<double> distinct_q;       // distinct reported q-values
    bool all_q_within_cutoff {true};
    double max_q {0.0};
    std::size_t n_donor_groups {};     // consensus features carrying >=1 donor
  };

  TransferStats analyzeTransfers(const ConsensusMap& cm, const std::map<UInt64, Truth>& truth,
                                 double cutoff)
  {
    TransferStats s;
    s.n_consensus = cm.size();

    for (const ConsensusFeature& cf : cm)
    {
      // Find the donor identity of this group (every group seeded from a
      // donor, so at least one handle is a donor).
      bool have_donor = false;
      std::string donor_seq;
      Int donor_charge = 0;
      for (const FeatureHandle& h : cf)
      {
        auto it = truth.find(h.getUniqueId());
        if (it != truth.end() && it->second.is_donor)
        {
          have_donor = true;
          donor_seq = it->second.seq;
          donor_charge = it->second.charge;
          break;
        }
      }
      if (have_donor) ++s.n_donor_groups;

      const bool has_meta = cf.metaValueExists("mbr_transfer_map_indices");
      if (has_meta)
      {
        ++s.n_groups_with_transfer;
        DoubleList qvals = cf.getMetaValue("mbr_transfer_qvalues");
        IntList midx = cf.getMetaValue("mbr_transfer_map_indices");
        s.n_meta_listed += midx.size();
        for (double q : qvals)
        {
          s.distinct_q.insert(q);
          s.max_q = std::max(s.max_q, q);
          if (q > cutoff + 1e-9) s.all_q_within_cutoff = false;
        }
      }

      // A handle is a transferred acceptor iff its uid is a non-donor
      // (acceptors only ever enter the consensus via the transfer path).
      for (const FeatureHandle& h : cf)
      {
        auto it = truth.find(h.getUniqueId());
        if (it == truth.end()) continue;            // shouldn't happen
        if (it->second.is_donor) continue;          // a donor handle
        ++s.n_transfers;
        const bool correct = have_donor && !it->second.noise
                             && it->second.seq == donor_seq
                             && it->second.charge == donor_charge;
        if (correct) ++s.n_correct; else ++s.n_false;
      }
    }
    return s;
  }
}

START_TEST(PipEchoAlgorithm_realistic, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION([EXTRA] realistic MBR controls the transfer FDR and recovers true transfers)
{
  const char* run_slow = std::getenv("OPENMS_RUN_SLOW_TESTS");
  if (!run_slow || !*run_slow)
  {
    std::cout << "\n[PipEchoAlgorithm_realistic] SKIPPED (slow). "
                 "Set OPENMS_RUN_SLOW_TESTS=1 to run the full ~10k-feature MBR test.\n";
    TEST_EQUAL(true, true) // skip silently; keeps default CI fast
  }
  else
  {
    // ---- Generation parameters (tuned so the SVM/FDR path is driven) ----
    const uint32_t GEN_SEED = 20260624u;
    const std::size_t N_PEPTIDES = 2750;
    const std::size_t N_RUNS = 3;
    const std::size_t TARGET_TOTAL = 10000;
    const double P_PRESENT = 0.85;
    const double P_IDENT = 0.5;

    Dataset ds = buildSyntheticDataset(GEN_SEED, N_PEPTIDES, N_RUNS, TARGET_TOTAL, P_PRESENT, P_IDENT);

    std::size_t total_features = 0;
    for (const auto& m : ds.maps) total_features += m.size();
    std::cout << "\n[gen] runs=" << ds.maps.size()
              << " total_features=" << total_features
              << " donors=" << ds.n_donor_features
              << " true_acceptors=" << ds.n_true_acceptors
              << " noise=" << ds.n_noise_features
              << " recoverable=" << ds.n_recoverable << "\n";
    for (std::size_t r = 0; r < ds.maps.size(); ++r)
    {
      std::cout << "[gen] run" << r << " features=" << ds.maps[r].size() << "\n";
    }

    // ---- Run group() at several FDR cutoffs ----
    // 0.001 is unresolvable for any realistic decoy count (it would need >= 1000
    // decoys), so it exercises the conservative drop-all gate at scale. 0.05 and
    // 0.10 are resolvable, so they exercise the SVM/PEP path and real FDR control.
    const std::vector<double> cutoffs = {0.001, 0.05, 0.10};
    std::vector<TransferStats> stats;
    std::vector<ConsensusMap> outputs;
    for (double fdr : cutoffs)
    {
      PipEchoAlgorithm algo;
      Param p = algo.getParameters();
      p.setValue("fdr", fdr);
      p.setValue("random_seed", 0); // deterministic decoy selection
      algo.setParameters(p);

      ConsensusMap cm;
      algo.group(ds.maps, cm);

      TransferStats s = analyzeTransfers(cm, ds.truth, fdr);
      double realized = s.n_transfers > 0 ? double(s.n_false) / double(s.n_transfers) : 0.0;
      std::cout << "[fdr=" << fdr << "] consensus=" << s.n_consensus
                << " donor_groups=" << s.n_donor_groups
                << " groups_with_transfer=" << s.n_groups_with_transfer
                << " transfers=" << s.n_transfers
                << " correct=" << s.n_correct
                << " false=" << s.n_false
                << " realized_FDR=" << realized
                << " distinct_q=" << s.distinct_q.size()
                << " max_q=" << s.max_q
                << " meta_listed=" << s.n_meta_listed << "\n";
      stats.push_back(s);
      outputs.push_back(std::move(cm));
    }

    const TransferStats& drop = stats[0]; // fdr = 0.001 (unresolvable)
    const TransferStats& s05 = stats[1];  // fdr = 0.05
    const TransferStats& s10 = stats[2];  // fdr = 0.10

    // ---- Structural invariants (hold at every cutoff) ----
    for (std::size_t i = 0; i < cutoffs.size(); ++i)
    {
      const TransferStats& s = stats[i];
      // Output is non-empty and every group is anchored by a donor.
      TEST_EQUAL(s.n_consensus > 0, true)
      TEST_EQUAL(s.n_donor_groups, s.n_consensus)
      // The transfer-provenance meta agrees with the actual transferred handles.
      TEST_EQUAL(s.n_meta_listed, s.n_transfers)
      // The FDR filter must never report a q-value above its own cutoff.
      TEST_EQUAL(s.all_q_within_cutoff, true)
    }

    // ---- Conservative gate at scale ----
    // An unresolvable cutoff drops ALL transfers but keeps every direct
    // identification (donor). This is the realistic-data regime that exposed the
    // decoy-generation bug: too few decoys -> FDR not estimable -> drop.
    TEST_EQUAL(drop.n_transfers, 0)
    TEST_EQUAL(drop.n_groups_with_transfer, 0)

    // ---- The SVM/FDR path actually ran (not the toy fallbacks) ----
    // Transfers surviving at fdr=0.05 prove the conservative estimability gate
    // passed (>= min_decoys decoys, decoys*fdr >= 1), i.e. real FDR control ran.
    TEST_EQUAL(s05.n_transfers > 0, true)
    // Non-degenerate q-values: the toy test collapses to a single q=0.5; real
    // FDR estimation over genuine decoys produces a spread of values.
    TEST_EQUAL(s05.distinct_q.size() > 5, true)

    // ---- FDR control (the headline) ----
    // Target/decoy FDR control is a finite-sample estimate, and the MBR-decoy
    // null only approximates the true false-transfer distribution (here made
    // deliberately adversarial with interfering ions), so we allow generous
    // slack but still require the realised FDR to be bounded -- a broken filter
    // lets the false fraction blow well past the cutoff.
    for (std::size_t i = 0; i < cutoffs.size(); ++i)
    {
      const TransferStats& s = stats[i];
      double realized = s.n_transfers > 0 ? double(s.n_false) / double(s.n_transfers) : 0.0;
      double bound = std::max(3.0 * cutoffs[i], cutoffs[i] + 0.05);
      TEST_EQUAL(realized <= bound, true)
    }

    // ---- Sensitivity: MBR must do useful work, not just be conservative ----
    // (A "drop everything" implementation would trivially have FDR 0.)
    TEST_EQUAL(s05.n_correct > 0, true)
    // Recover a substantial fraction of the transfers that are even possible.
    TEST_EQUAL(double(s05.n_correct) >= 0.3 * double(ds.n_recoverable), true)

    // ---- Monotonicity: a looser cutoff keeps at least as many transfers ----
    TEST_EQUAL(s05.n_transfers >= drop.n_transfers, true) // 0.05 >= 0.001
    TEST_EQUAL(s10.n_transfers >= s05.n_transfers, true)  // 0.10 >= 0.05

    // ---- Determinism: same seeds -> identical consensus structure ----
    // Re-run the fdr=0.05 case and compare against the stored output.
    {
      PipEchoAlgorithm a2;
      Param p = a2.getParameters();
      p.setValue("fdr", 0.05);
      p.setValue("random_seed", 0);
      a2.setParameters(p);
      ConsensusMap cm2;
      a2.group(ds.maps, cm2);

      const ConsensusMap& cm1 = outputs[1];
      TEST_EQUAL(cm1.size(), cm2.size())
      bool identical = cm1.size() == cm2.size();
      for (Size i = 0; identical && i < cm1.size(); ++i)
      {
        identical = (cm1[i].getRT() == cm2[i].getRT())
                    && (cm1[i].getMZ() == cm2[i].getMZ())
                    && (cm1[i].size() == cm2[i].size());
      }
      TEST_EQUAL(identical, true)
    }

    // ---- Training-point cap: FDR control + sensitivity survive subsampling ----
    // This dataset produces thousands of acceptors per fold, so a small
    // max_training_points forces the stratified-subsample path in round().
    // Because prediction and q-values still use EVERY acceptor, the realised
    // FDR must stay bounded and a useful fraction of true transfers must still
    // be recovered -- the cap may only change the fitted boundary, not the set
    // of scored transfers.
    {
      PipEchoAlgorithm algo;
      Param p = algo.getParameters();
      p.setValue("fdr", 0.05);
      p.setValue("random_seed", 0);
      p.setValue("max_training_points", 500); // well below the per-fold training size
      algo.setParameters(p);

      ConsensusMap cm;
      algo.group(ds.maps, cm);
      TransferStats s = analyzeTransfers(cm, ds.truth, 0.05);
      double realized = s.n_transfers > 0 ? double(s.n_false) / double(s.n_transfers) : 0.0;

      // FDR control is not broken by the cap (same generous bound as above).
      TEST_TRUE(realized <= std::max(3.0 * 0.05, 0.05 + 0.05))
      // The cap must not collapse the method to "drop everything".
      TEST_TRUE(s.n_transfers > 0)
      TEST_TRUE(s.n_correct > 0)
      // The q-value filter still respects its own cutoff.
      TEST_TRUE(s.all_q_within_cutoff)

      // Deterministic subsample (uniform stride, no RNG) -> reproducible.
      PipEchoAlgorithm a2;
      Param p2 = a2.getParameters();
      p2.setValue("fdr", 0.05);
      p2.setValue("random_seed", 0);
      p2.setValue("max_training_points", 500);
      a2.setParameters(p2);
      ConsensusMap cm2;
      a2.group(ds.maps, cm2);
      bool same = (cm.size() == cm2.size());
      for (Size i = 0; same && i < cm.size(); ++i)
      {
        same = (cm[i].getRT() == cm2[i].getRT())
               && (cm[i].getMZ() == cm2[i].getMZ())
               && (cm[i].size() == cm2[i].size());
      }
      TEST_TRUE(same)
    }
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
