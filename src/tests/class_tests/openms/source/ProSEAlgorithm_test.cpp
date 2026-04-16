// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/LogStream.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/ProSEAlgorithm.h>
///////////////////////////

#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/ANALYSIS/ID/OpenSearchModificationAnalysis.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>

#include <random>
#include <set>

using namespace OpenMS;
using namespace std;

// Test subclass exposing internal state for white-box assertions on
// asymmetric bounds, calibration-pass results, and the mod-match tolerance helper.
//
// The `friend class ProSEAlgorithm_test;` declaration in the
// production header makes the `using` re-exposure below legal.
//
// Note: `fragment_index_` is NOT a member of ProSEAlgorithm — it
// lives on the SearchContext and is a local reference during search(). Tests that
// need to observe FragmentIndex state must do so via prepareContext() + the
// context-taking search() overload before/after the restore hook fires.
class ProSEAlgorithm_test : public ProSEAlgorithm
{
public:
  using ProSEAlgorithm::precursor_mass_tolerance_lower_;
  using ProSEAlgorithm::precursor_mass_tolerance_upper_;
  using ProSEAlgorithm::precursor_mass_tolerance_unit_;
  using ProSEAlgorithm::computeModMatchTolerance_;
  using ProSEAlgorithm::last_calibration_result_;
  using ProSEAlgorithm::last_mod_match_tolerance_used_;
  using ProSEAlgorithm::CalibrationResult_;
};

// --- Shared calibration fixture -------------------------------------------------
//
// Tests 7 and 8 share: a small protein database, a list of synthetic MS2 spectra
// with per-spectrum scattered precursor errors whose median is ~+7 ppm, and
// identical algorithm parameters. The free functions below build the fixture to
// keep the tests' top-level assertions focused.
//
// Design:
//   1. Digest the test protein into tryptic peptides, keep those >= 8 residues
//      so the TSG + preprocess pipeline produces enough fragment peaks to score.
//   2. For each kept peptide, generate a synthetic MS2 spectrum and apply a
//      per-spectrum ppm-level precursor error from a scattered distribution
//      centered on +7 ppm.
//   3. User window [20, 30] ppm is wide enough that:
//        - the candidate look-up in computeMassWindow_() finds the theoretical
//          peptide (max applied error < 30)
//        - the wrong-match filter at ProSEAlgorithm.cpp:1584 passes every hit
//      and the calibration result is a genuine tightening relative to user
//      bounds rather than a clamp-to-user artifact.
//   4. calibration:min_psms = 5 bypasses the top-50% score crop (ProSEAlgorithm.cpp:1571)
//      for our small fixture (only ~12 hits) so the full scattered distribution
//      reaches the median/MAD estimator.
//
// Error distribution rationale:
//   The 12-element error vector below has median = 7.0 and a spread such that
//   precursor_spread = median(|e-7|) + 3 * MAD(|e-7|) > 7, keeping
//   extreme_bias = false. See comments in ProSEAlgorithm.cpp runCalibrationPass_.
static vector<FASTAFile::FASTAEntry> calibration_fasta_db_()
{
  return {
    {"P01", "Test", "MSDEREKVLGFHQRMPNASTICYWDLKEGFVRTHQPSANLDIKCMYKWTE"
                    "RHASGDFLKPIVEQNCTMYRGWSADELKHPFNQGTICMSYREWDAVLKPH"
                    "GITNSEYRQWDLKAPMFHCVSITGNREYWDKLMPAHFQCSTVINEYRWDLK"
                    "APMHSCFTGQNVIREYWDKLMSPAHCFQNTSGIVREYWDKLHMPASCFQGN"},
  };
}

static PeakMap build_calibration_spectra_(const vector<double>& ppm_shifts)
{
  // Digest the test protein into tryptic peptides >= 8 residues.
  ProteaseDigestion digester;
  digester.setEnzyme("Trypsin");
  digester.setMissedCleavages(1);

  vector<AASequence> peptides;
  for (const auto& entry : calibration_fasta_db_())
  {
    AASequence protein = AASequence::fromString(entry.sequence);
    digester.digest(protein, peptides, 8, 40);
  }

  TheoreticalSpectrumGenerator tsg;
  Param tsg_param = tsg.getParameters();
  tsg_param.setValue("add_first_prefix_ion", "true");
  tsg_param.setValue("add_metainfo", "true");
  tsg.setParameters(tsg_param);

  PeakMap spectra;
  double rt = 100.0;
  Size emitted = 0;
  for (const auto& pep : peptides)
  {
    if (emitted >= ppm_shifts.size()) break;
    if (pep.size() < 8) continue;

    int charge = 2;
    MSSpectrum spec;
    // Two fragment charges (b+y at z=1 and z=2) — matches the working pattern
    // used by the Synthetic modification discovery test. Short peptides with
    // only z=1 ions get stripped to nothing by the downstream preprocessing
    // pipeline (Deisotoper + WindowMower + NLargest).
    tsg.getSpectrum(spec, pep, 1, std::min<int>(charge - 1, 2));
    spec.sortByPosition();
    if (spec.size() < 10) continue;

    spec.setMSLevel(2);
    spec.setRT(rt);
    rt += 1.0;

    Precursor prec;
    double mz = pep.getMZ(charge);
    prec.setMZ(mz * (1.0 + ppm_shifts[emitted] * 1e-6));
    prec.setCharge(charge);
    spec.setPrecursors({prec});
    spec.setNativeID("spectrum=" + String(spectra.size()));

    spectra.addSpectrum(std::move(spec));
    ++emitted;
  }
  return spectra;
}

static void configure_calibration_params_(ProSEAlgorithm& algo,
                                          double lower_ppm,
                                          double upper_ppm,
                                          Size min_psms = 3)
{
  Param p = algo.getParameters();
  p.setValue("precursor:mass_tolerance_lower", lower_ppm);
  p.setValue("precursor:mass_tolerance_upper", upper_ppm);
  p.setValue("precursor:mass_tolerance_unit", "ppm");
  p.setValue("fragment:mass_tolerance", 20.0);
  p.setValue("fragment:mass_tolerance_unit", "ppm");
  p.setValue("fragment:min_ion_index", 0);
  p.setValue("calibration:enabled", "true");
  p.setValue("calibration:subset_ratio", 1.0);
  // min_psms is chosen per-test. When min_psms > cal_hits/2, the top-50% score
  // crop at ProSEAlgorithm.cpp:1571 is skipped and every collected error reaches the
  // estimator. For our small fixture that's what we want.
  p.setValue("calibration:min_psms", static_cast<Int>(min_psms));
  p.setValue("decoys", "false");
  p.setValue("peptide:min_size", 7);
  p.setValue("peptide:max_size", 40);
  p.setValue("peptide:missed_cleavages", 1);
  algo.setParameters(p);
}

START_TEST(ProSEAlgorithm, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ProSEAlgorithm* ptr = nullptr;
ProSEAlgorithm* null_ptr = nullptr;

START_SECTION(ProSEAlgorithm())
{
  ptr = new ProSEAlgorithm();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~ProSEAlgorithm())
{
  delete ptr;
}
END_SECTION

START_SECTION(([EXTRA] Synthetic modification discovery - open search))
{
  // =========================================================================
  // Strategy: Generate fragments from unmodified+Carbamidomethyl sequences
  //   (ensuring perfect fragment matching), but shift precursor m/z by the
  //   modification mass (creating the delta mass for open search discovery).
  //   This mimics real open search data where fragment ions match the
  //   unmodified backbone and the precursor reveals the mass shift.
  // =========================================================================

  // Modifications to test: {name, mass_shift_Da}
  struct ModDef
  {
    string name;
    double mass;
  };

  vector<ModDef> test_mods = {
    {"Oxidation (M)",            15.9949},
    {"Phospho (S)",              79.9663},
    {"Phospho (T)",              79.9663},
    {"Acetyl (K)",               42.0106},
    {"Deamidated (N)",            0.9840},
    {"Methyl (R)",               14.0157},
    {"Dimethyl (K)",             28.0314},
    {"Carbamyl (K)",             43.0058},
    {"Dioxidation (M)",          31.9898},
    {"HexNAc (S)",              203.0794},
    {"Acetyl (Protein N-term)",  42.0106},
    {"Formyl (Protein N-term)",  27.9949},
    {"Unknown (artificial)",    123.4560},  // artificial mass not matching any known modification
  };

  // =========================================================================
  // Create synthetic protein database
  // =========================================================================
  vector<FASTAFile::FASTAEntry> fasta_db = {
    {"P01", "Protein01",
     "MSDEREKVLGFHQRMPNASTICYWDLKEGFVRTHQPSANLDIKCMYKWTE"
     "RHASGDFLKPIVEQNCTMYRGWSADELKHPFNQGTICMSYREWDAVLKPH"
     "GITNSEYRQWDLKAPMFHCVSITGNREYWDKLMPAHFQCSTVINEYRWDLK"
     "APMHSCFTGQNVIREYWDKLMSPAHCFQNTSGIVREYWDKLHMPASCFQGN"},
    {"P02", "Protein02",
     "MKAILNHVGSTFREDWQCPYLKMISGDTFNHRVAWQECPLKYMTGISNHFR"
     "DVEWAQCPLKTMIYGSNHFRDVEWAQCPKLIMTGSYNHFRDVEWAQCKPLIM"
     "TGSYHFNRDVEWAQCPKLITMGSYHFNRDVEWAQCPKLITMGSYHNFRDVEW"
     "AQCPKLITMGSYHFNRDVEWAQCPKLITMGSYHFNRDVEWAQCPKLITMGSY"},
    {"P03", "Protein03",
     "MGHYIKLTPNRESWDVAFQCKMHGYILKTPNRESWDVAFQCKMHGLIYTKP"
     "NRESWDVAFQCKHMGIYLTKPNRESWDVAFQCKMHGIYLKTNPRESWDVAFQ"
     "CKMHGYLIKTPNRESWDVAFQCKMHGIYLTKPNRESWDVAFQCKMHGIYLTK"
     "PNRESWDVAFQCKMHGIYLTKPNRESWDVAFQCKMHGIYLTKPNRESWDVAF"},
    {"P04", "Protein04",
     "MSVDNKTHFRGECAWYPILQMSDKTNHFRGEVAWCYQPILKMSDETKNHFRG"
     "VAWCEQYPILKMSDTKENHFRGVAWCEYQPILKMSDETKHNFRGVACWEYQPI"
     "LKMSDETKNHFRGVAWCEYQPILKMSDETKNHFRGVAWCEYQPILKMSDETNK"
     "HFRGVAWCEYQPILKMSDETKNHFRGVAWCEYQPILKMSDETKNHFRGVAWCE"},
    {"P05", "Protein05",
     "MAKLFGYNRSTECWDIPHQVMKALYFGNRSTWECDIPHQVKMALGFYNRSTWE"
     "CDIPHQVKMALFGYNRSTEWCDIPHQVKMAFGLYNRSTWECDIPHQVKMALFY"
     "GNRSTEWCDIPHQVKMALFGYNRSTEWCDIPHQVKMALFGYNRSTEWCDIPHQ"
     "VKMALFGYNRSTEWCDIPHQVKMALFGYNRSTEWCDIPHQVKMALFGYNRSTE"},
    {"P06", "Protein06",
     "MTGYLSKFHERNDICWPAQVMTGLYSKHFERNDICWPAQVMTGLYSKFHRNDE"
     "ICWPAQVKMTGYLSKFHERNDICWPAQVKMTGLYSKHFERNDICWPAQVKMTG"
     "LYKSFHERNDICWPAQVKMTGLYSKHFERNDICWPAQVKMTGLYSKFHERNDIC"
     "WPAQVKMTGLYSKFHERNDICWPAQVKMTGLYSKFHERNDICWPAQVKMTGLY"},
    {"P07", "Protein07",
     "MDIKHWNRSYPLCTGEFAQVKMDIHKWNRSYPLCTGFAEQVKMDIKHWNRSYP"
     "LCTGFEAQVKMDIHKWNRSYPLCTGEFAQVKMDIKHWNRSYPLCTGFEAQVKM"
     "DIHKWNRSYPLCTGEFAQVKMDIHKWNRSYPLCTGFEAQVKMDIHKWNRSYPL"
     "CTGFEAQVKMDIHKWNRSYPLCTGFEAQVKMDIHKWNRSYPLCTGFEAQVKMD"},
    {"P08", "Protein08",
     "MEYKFADLHGSNTCRWQPIVKMEYFKADLHGSNTCRWQPVIKMEYFKADLHGS"
     "NTCRWQPIVKMEYFKADLHGSNTRWCQPIVKMEYFDKALGHSNTRWCQPIVKM"
     "EYFKADLGHSNTCRWQPIVKMEYFKADLHGSNTCRWQPIVKMEYFKADLHGSN"
     "TCRWQPIVKMEYFKADLHGSNTCRWQPIVKMEYFKADLHGSNTCRWQPIVKME"},
    {"P09", "Protein09",
     "MQHWVDESYRFTNGPILCKAMQHWVDESYRTFNGPILCKAMQHWVEDYSRTFNG"
     "PILCKAMQHWVEDYSRFTNGPILCKAMQHWVEDYSRFTNGPILCKAMQHWVEDY"
     "SRFTNGPILCKAMQHWVEDYSRFTNGPILCKAMQHWVEDYSRFTNGPILCKAMQ"
     "HWVEDYSRFTNGPILCKAMQHWVEDYSRFTNGPILCKAMQHWVEDYSRFTNGPI"},
    {"P10", "Protein10",
     "MTEFLNQGDKSYCRHWPIVAMTEFNLQGDKSYCRHWPIVAMTEFLNQGDKSYCR"
     "HWPIVAMTEFLNQGDKSYCHRRWPIVAMTEFLNQGDKSYCRHWPIVAMTEFLNQ"
     "GDKSYCRHWPIVAMTEFLNQGDKSYCRHWPIVAMTEFLNQGDKSYCRHWPIVAM"
     "TEFLNQGDKSYCRHWPIVAMTEFLNQGDKSYCRHWPIVAMTEFLNQGDKSYCR"},
  };

  // =========================================================================
  // Digest proteins
  // =========================================================================
  ProteaseDigestion digester;
  digester.setEnzyme("Trypsin");
  digester.setMissedCleavages(1);

  // Digest and apply Carbamidomethyl(C) - this is what the search engine will do
  ModifiedPeptideGenerator::MapToResidueType fixed_mods =
    ModifiedPeptideGenerator::getModifications({"Carbamidomethyl (C)"});

  vector<AASequence> all_peptides;
  for (const auto& entry : fasta_db)
  {
    AASequence protein = AASequence::fromString(entry.sequence);
    vector<AASequence> peptides;
    digester.digest(protein, peptides, 7, 40);
    for (auto& pep : peptides)
    {
      ModifiedPeptideGenerator::applyFixedModifications(fixed_mods, pep);
      all_peptides.push_back(std::move(pep));
    }
  }

  OPENMS_LOG_INFO << "[TEST] Total tryptic peptides (with Carbamidomethyl): " << all_peptides.size() << std::endl;
  TEST_TRUE(all_peptides.size() > 50)

  // =========================================================================
  // Generate spectra: unmodified fragments + shifted precursor
  // =========================================================================
  mt19937 rng(42);

  TheoreticalSpectrumGenerator tsg;
  Param tsg_param = tsg.getParameters();
  tsg_param.setValue("add_first_prefix_ion", "true");
  tsg_param.setValue("add_metainfo", "true");
  tsg.setParameters(tsg_param);

  PeakMap spectra;
  double rt = 100.0;
  Size total_spectra = 0;
  map<string, Size> mod_spectrum_count;

  const Size target_per_mod = 400;

  for (const auto& mod_def : test_mods)
  {
    // Shuffle peptides and take first target_per_mod
    vector<size_t> indices(all_peptides.size());
    iota(indices.begin(), indices.end(), 0);
    shuffle(indices.begin(), indices.end(), rng);

    Size created = 0;
    for (size_t idx : indices)
    {
      if (created >= target_per_mod) break;
      const AASequence& pep = all_peptides[idx];
      if (pep.size() < 8) continue; // need enough fragments

      int charge = 2 + (int)(rng() % 3); // charge 2-4

      // Generate theoretical spectrum from the unmodified+Carbamidomethyl peptide
      MSSpectrum spec;
      tsg.getSpectrum(spec, pep, 1, min(charge - 1, 2));
      spec.sortByPosition();
      if (spec.size() < 10) continue; // need enough fragment peaks

      spec.setMSLevel(2);
      spec.setRT(rt);
      rt += 0.1;

      // Set precursor: true m/z of unmodified peptide + modification mass shift
      double unmod_mz = pep.getMZ(charge);
      double shifted_mz = unmod_mz + mod_def.mass / (double)charge;

      Precursor prec;
      prec.setMZ(shifted_mz);
      prec.setCharge(charge);
      spec.setPrecursors({prec});
      spec.setNativeID("spectrum=" + String(spectra.size()));

      spectra.addSpectrum(std::move(spec));
      created++;
      total_spectra++;
    }

    mod_spectrum_count[mod_def.name] = created;
    OPENMS_LOG_INFO << "[TEST] " << mod_def.name << ": " << created << " spectra" << std::endl;
  }

  OPENMS_LOG_INFO << "[TEST] Total spectra: " << total_spectra << std::endl;
  TEST_TRUE(total_spectra > 2000)

  // =========================================================================
  // Configure and run open search
  // =========================================================================
  ProSEAlgorithm algo;
  Param p = algo.getParameters();
  p.setValue("precursor:mass_tolerance_lower", 500.0);
  p.setValue("precursor:mass_tolerance_upper", 500.0);
  p.setValue("precursor:mass_tolerance_unit", "Da");
  p.setValue("fragment:mass_tolerance", 20.0);
  p.setValue("fragment:mass_tolerance_unit", "ppm");
  p.setValue("modifications:fixed", vector<string>{"Carbamidomethyl (C)"});
  p.setValue("modifications:variable", vector<string>{});
  p.setValue("decoys", "false");
  p.setValue("peptide:min_size", 7);
  p.setValue("peptide:max_size", 40);
  p.setValue("peptide:missed_cleavages", 1);
  p.setValue("report:top_hits", 1);
  algo.setParameters(p);

  auto result = algo.searchWithModificationAnalysis(spectra, fasta_db, "");

  // =========================================================================
  // Verify results
  // =========================================================================
  TEST_EQUAL(result.exit_code == ProSEAlgorithm::ExitCodes::EXECUTION_OK, true)
  TEST_EQUAL(result.is_open_search, true)

  OPENMS_LOG_INFO << "[TEST] Total PSMs: " << result.peptide_ids.size() << std::endl;
  TEST_TRUE(result.peptide_ids.size() > 500)

  // Log PTM analysis results
  for (const auto& ptm : result.modification_analysis.ptm_stats.entries)
  {
    OPENMS_LOG_INFO << "[TEST] PTM: " << ptm.name
                    << " (count=" << ptm.count
                    << ", theo_mass=" << ptm.theoretical_mass
                    << ", obs_mass=" << ptm.observed_mass << ")" << std::endl;
  }

  // Log delta mass entries (the core of modification discovery)
  OPENMS_LOG_INFO << "[TEST] Delta mass entries: "
                  << result.modification_analysis.delta_mass_stats.entries.size() << std::endl;
  for (const auto& dm : result.modification_analysis.delta_mass_stats.entries)
  {
    if (dm.count >= 10)
    {
      OPENMS_LOG_INFO << "[TEST] DeltaMass=" << dm.delta_mass
                      << " count=" << dm.count
                      << " mapped=" << dm.mapped_modification << std::endl;
    }
  }

  // ===================================================================
  // Verify delta mass bins: the heart of modification discovery.
  // The delta mass histogram should contain bins at each modification mass
  // we injected, with counts close to the number of spectra generated.
  // ===================================================================
  const auto& dm_entries = result.modification_analysis.delta_mass_stats.entries;
  TEST_TRUE(!dm_entries.empty())
  TEST_TRUE(result.modification_analysis.delta_mass_stats.total_psms > 0)

  // Build a lookup: delta_mass -> count
  map<double, Size> dm_counts; // key = rounded delta mass
  for (const auto& dm : dm_entries)
  {
    dm_counts[dm.delta_mass] = dm.count;
  }

  // Expected delta masses and minimum expected counts
  vector<pair<string, double>> must_find = {
    {"Oxidation (M)",             15.9949},
    {"Phospho (S/T/Y)",           79.9663},
    {"Acetyl (K/N-term)",         42.0106},
    {"Deamidated (N/Q)",           0.9840},
    {"Methyl (R)",                14.0157},
    {"Dimethyl (K)",              28.0314},
    {"Carbamyl (K)",              43.0058},
    {"Dioxidation (M)",           31.9898},
    {"HexNAc (S)",              203.0794},
    {"Formyl (Protein N-term)",   27.9949},
  };

  for (const auto& [label, expected_mass] : must_find)
  {
    // Find a delta mass bin within 0.05 Da of the expected mass
    bool found = false;
    Size count = 0;
    for (const auto& [dm, cnt] : dm_counts)
    {
      if (fabs(dm - expected_mass) < 0.05)
      {
        found = true;
        count = cnt;
        break;
      }
    }
    if (!found)
    {
      OPENMS_LOG_INFO << "[TEST] MISSING delta mass bin for " << label
                      << " (expected ~" << expected_mass << " Da)" << std::endl;
    }
    else
    {
      OPENMS_LOG_INFO << "[TEST] Found delta mass for " << label
                      << ": count=" << count << std::endl;
      // Each modification should have significant counts (we generated ~378 per mod)
      TEST_TRUE(count >= 100)
    }
    TEST_EQUAL(found, true)
  }

  // ===================================================================
  // Verify that the artificial unknown modification (123.456 Da) is
  // present in the delta mass histogram but NOT mapped to a known mod.
  // ===================================================================
  {
    bool found_unknown_bin = false;
    for (const auto& dm : dm_entries)
    {
      if (fabs(dm.delta_mass - 123.456) < 0.05)
      {
        found_unknown_bin = true;
        OPENMS_LOG_INFO << "[TEST] Unknown delta mass bin at " << dm.delta_mass
                        << ": count=" << dm.count
                        << " mapped='" << dm.mapped_modification << "'"
                        << " is_known=" << dm.is_known_modification << std::endl;
        TEST_TRUE(dm.count >= 50)
        TEST_EQUAL(dm.is_known_modification, false)
        break;
      }
    }
    TEST_EQUAL(found_unknown_bin, true)
  }
}
END_SECTION

START_SECTION(([EXTRA] FDR-filtered modification discovery))
{
  // =========================================================================
  // Same synthetic data as the open search test, but with FDR filtering
  // before modification analysis. Demonstrates the workflow:
  //   1. Search with decoys
  //   2. Compute q-values via FalseDiscoveryRate
  //   3. Filter by FDR threshold
  //   4. Remove decoy hits
  //   5. Run modification analysis on filtered results
  // =========================================================================

  struct ModDef
  {
    string name;
    double mass;
  };

  vector<ModDef> test_mods = {
    {"Oxidation (M)",            15.9949},
    {"Phospho (S)",              79.9663},
    {"Phospho (T)",              79.9663},
    {"Acetyl (K)",               42.0106},
    {"Deamidated (N)",            0.9840},
    {"Methyl (R)",               14.0157},
    {"Dimethyl (K)",             28.0314},
    {"Carbamyl (K)",             43.0058},
    {"Dioxidation (M)",          31.9898},
    {"HexNAc (S)",              203.0794},
    {"Acetyl (Protein N-term)",  42.0106},
    {"Formyl (Protein N-term)",  27.9949},
    {"Unknown (artificial)",    123.4560},
  };

  // Reuse the same protein database
  vector<FASTAFile::FASTAEntry> fasta_db = {
    {"P01", "Protein01",
     "MSDEREKVLGFHQRMPNASTICYWDLKEGFVRTHQPSANLDIKCMYKWTE"
     "RHASGDFLKPIVEQNCTMYRGWSADELKHPFNQGTICMSYREWDAVLKPH"
     "GITNSEYRQWDLKAPMFHCVSITGNREYWDKLMPAHFQCSTVINEYRWDLK"
     "APMHSCFTGQNVIREYWDKLMSPAHCFQNTSGIVREYWDKLHMPASCFQGN"},
    {"P02", "Protein02",
     "MKAILNHVGSTFREDWQCPYLKMISGDTFNHRVAWQECPLKYMTGISNHFR"
     "DVEWAQCPLKTMIYGSNHFRDVEWAQCPKLIMTGSYNHFRDVEWAQCKPLIM"
     "TGSYHFNRDVEWAQCPKLITMGSYHFNRDVEWAQCPKLITMGSYHNFRDVEW"
     "AQCPKLITMGSYHFNRDVEWAQCPKLITMGSYHFNRDVEWAQCPKLITMGSY"},
    {"P03", "Protein03",
     "MGHYIKLTPNRESWDVAFQCKMHGYILKTPNRESWDVAFQCKMHGLIYTKP"
     "NRESWDVAFQCKHMGIYLTKPNRESWDVAFQCKMHGIYLKTNPRESWDVAFQ"
     "CKMHGYLIKTPNRESWDVAFQCKMHGIYLTKPNRESWDVAFQCKMHGIYLTK"
     "PNRESWDVAFQCKMHGIYLTKPNRESWDVAFQCKMHGIYLTKPNRESWDVAF"},
    {"P04", "Protein04",
     "MSVDNKTHFRGECAWYPILQMSDKTNHFRGEVAWCYQPILKMSDETKNHFRG"
     "VAWCEQYPILKMSDTKENHFRGVAWCEYQPILKMSDETKHNFRGVACWEYQPI"
     "LKMSDETKNHFRGVAWCEYQPILKMSDETKNHFRGVAWCEYQPILKMSDETNK"
     "HFRGVAWCEYQPILKMSDETKNHFRGVAWCEYQPILKMSDETKNHFRGVAWCE"},
    {"P05", "Protein05",
     "MAKLFGYNRSTECWDIPHQVMKALYFGNRSTWECDIPHQVKMALGFYNRSTWE"
     "CDIPHQVKMALFGYNRSTEWCDIPHQVKMAFGLYNRSTWECDIPHQVKMALFY"
     "GNRSTEWCDIPHQVKMALFGYNRSTEWCDIPHQVKMALFGYNRSTEWCDIPHQ"
     "VKMALFGYNRSTEWCDIPHQVKMALFGYNRSTEWCDIPHQVKMALFGYNRSTE"},
    {"P06", "Protein06",
     "MTGYLSKFHERNDICWPAQVMTGLYSKHFERNDICWPAQVMTGLYSKFHRNDE"
     "ICWPAQVKMTGYLSKFHERNDICWPAQVKMTGLYSKHFERNDICWPAQVKMTG"
     "LYKSFHERNDICWPAQVKMTGLYSKHFERNDICWPAQVKMTGLYSKFHERNDIC"
     "WPAQVKMTGLYSKFHERNDICWPAQVKMTGLYSKFHERNDICWPAQVKMTGLY"},
    {"P07", "Protein07",
     "MDIKHWNRSYPLCTGEFAQVKMDIHKWNRSYPLCTGFAEQVKMDIKHWNRSYP"
     "LCTGFEAQVKMDIHKWNRSYPLCTGEFAQVKMDIKHWNRSYPLCTGFEAQVKM"
     "DIHKWNRSYPLCTGEFAQVKMDIHKWNRSYPLCTGFEAQVKMDIHKWNRSYPL"
     "CTGFEAQVKMDIHKWNRSYPLCTGFEAQVKMDIHKWNRSYPLCTGFEAQVKMD"},
    {"P08", "Protein08",
     "MEYKFADLHGSNTCRWQPIVKMEYFKADLHGSNTCRWQPVIKMEYFKADLHGS"
     "NTCRWQPIVKMEYFKADLHGSNTRWCQPIVKMEYFDKALGHSNTRWCQPIVKM"
     "EYFKADLGHSNTCRWQPIVKMEYFKADLHGSNTCRWQPIVKMEYFKADLHGSN"
     "TCRWQPIVKMEYFKADLHGSNTCRWQPIVKMEYFKADLHGSNTCRWQPIVKME"},
    {"P09", "Protein09",
     "MQHWVDESYRFTNGPILCKAMQHWVDESYRTFNGPILCKAMQHWVEDYSRTFNG"
     "PILCKAMQHWVEDYSRFTNGPILCKAMQHWVEDYSRFTNGPILCKAMQHWVEDY"
     "SRFTNGPILCKAMQHWVEDYSRFTNGPILCKAMQHWVEDYSRFTNGPILCKAMQ"
     "HWVEDYSRFTNGPILCKAMQHWVEDYSRFTNGPILCKAMQHWVEDYSRFTNGPI"},
    {"P10", "Protein10",
     "MTEFLNQGDKSYCRHWPIVAMTEFNLQGDKSYCRHWPIVAMTEFLNQGDKSYCR"
     "HWPIVAMTEFLNQGDKSYCHRRWPIVAMTEFLNQGDKSYCRHWPIVAMTEFLNQ"
     "GDKSYCRHWPIVAMTEFLNQGDKSYCRHWPIVAMTEFLNQGDKSYCRHWPIVAM"
     "TEFLNQGDKSYCRHWPIVAMTEFLNQGDKSYCRHWPIVAMTEFLNQGDKSYCR"},
  };

  // Digest proteins
  ProteaseDigestion digester;
  digester.setEnzyme("Trypsin");
  digester.setMissedCleavages(1);

  ModifiedPeptideGenerator::MapToResidueType fixed_mods =
    ModifiedPeptideGenerator::getModifications({"Carbamidomethyl (C)"});

  vector<AASequence> all_peptides;
  for (const auto& entry : fasta_db)
  {
    AASequence protein = AASequence::fromString(entry.sequence);
    vector<AASequence> peptides;
    digester.digest(protein, peptides, 7, 40);
    for (auto& pep : peptides)
    {
      ModifiedPeptideGenerator::applyFixedModifications(fixed_mods, pep);
      all_peptides.push_back(std::move(pep));
    }
  }
  TEST_TRUE(all_peptides.size() > 50)

  // Generate spectra with shifted precursors
  mt19937 rng(42);
  TheoreticalSpectrumGenerator tsg;
  Param tsg_param = tsg.getParameters();
  tsg_param.setValue("add_first_prefix_ion", "true");
  tsg_param.setValue("add_metainfo", "true");
  tsg.setParameters(tsg_param);

  PeakMap spectra;
  double rt = 100.0;
  const Size target_per_mod = 400;

  for (const auto& mod_def : test_mods)
  {
    vector<size_t> indices(all_peptides.size());
    iota(indices.begin(), indices.end(), 0);
    shuffle(indices.begin(), indices.end(), rng);

    Size created = 0;
    for (size_t idx : indices)
    {
      if (created >= target_per_mod) break;
      const AASequence& pep = all_peptides[idx];
      if (pep.size() < 8) continue;

      int charge = 2 + (int)(rng() % 3);
      MSSpectrum spec;
      tsg.getSpectrum(spec, pep, 1, min(charge - 1, 2));
      spec.sortByPosition();
      if (spec.size() < 10) continue;

      spec.setMSLevel(2);
      spec.setRT(rt);
      rt += 0.1;

      double unmod_mz = pep.getMZ(charge);
      double shifted_mz = unmod_mz + mod_def.mass / (double)charge;

      Precursor prec;
      prec.setMZ(shifted_mz);
      prec.setCharge(charge);
      spec.setPrecursors({prec});
      spec.setNativeID("spectrum=" + String(spectra.size()));

      spectra.addSpectrum(std::move(spec));
      created++;
    }
  }
  TEST_TRUE(spectra.size() > 2000)

  // =========================================================================
  // Step 1: Search with decoys enabled
  // =========================================================================
  ProSEAlgorithm algo;
  Param p = algo.getParameters();
  p.setValue("precursor:mass_tolerance_lower", 500.0);
  p.setValue("precursor:mass_tolerance_upper", 500.0);
  p.setValue("precursor:mass_tolerance_unit", "Da");
  p.setValue("fragment:mass_tolerance", 20.0);
  p.setValue("fragment:mass_tolerance_unit", "ppm");
  p.setValue("modifications:fixed", vector<string>{"Carbamidomethyl (C)"});
  p.setValue("modifications:variable", vector<string>{});
  p.setValue("decoys", "true");  // Enable decoys for FDR
  p.setValue("peptide:min_size", 7);
  p.setValue("peptide:max_size", 40);
  p.setValue("peptide:missed_cleavages", 1);
  p.setValue("report:top_hits", 1);
  algo.setParameters(p);

  vector<ProteinIdentification> prot_ids;
  PeptideIdentificationList pep_ids;
  auto ec = algo.search(spectra, fasta_db, prot_ids, pep_ids);

  TEST_EQUAL(ec == ProSEAlgorithm::ExitCodes::EXECUTION_OK, true)
  Size total_before = pep_ids.size();
  OPENMS_LOG_INFO << "[TEST FDR] Total PSMs before filtering: " << total_before << std::endl;
  TEST_TRUE(total_before > 500)

  // =========================================================================
  // Step 2: Compute q-values via target-decoy FDR
  // =========================================================================
  FalseDiscoveryRate fdr_calculator;
  fdr_calculator.apply(pep_ids);

  // =========================================================================
  // Step 3: Filter at 5% FDR, remove decoys, clean up
  // =========================================================================
  IDFilter::filterHitsByScore(pep_ids, 0.05);
  IDFilter::removeDecoyHits(pep_ids);
  IDFilter::removeEmptyIdentifications(pep_ids);

  Size total_after = pep_ids.size();
  OPENMS_LOG_INFO << "[TEST FDR] PSMs after FDR filtering: " << total_after << std::endl;
  TEST_TRUE(total_after > 0)
  TEST_TRUE(total_after < total_before)

  // =========================================================================
  // Step 4: Run modification analysis on filtered results
  // =========================================================================
  OpenSearchModificationAnalysis mod_analyzer;
  auto mod_result = mod_analyzer.analyzeModificationsWithStatistics(
    pep_ids, 500.0, false, false, "");

  const auto& dm_entries = mod_result.delta_mass_stats.entries;
  TEST_TRUE(!dm_entries.empty())

  OPENMS_LOG_INFO << "[TEST FDR] Delta mass entries after filtering: "
                  << dm_entries.size() << std::endl;
  for (const auto& dm : dm_entries)
  {
    if (dm.count >= 10)
    {
      OPENMS_LOG_INFO << "[TEST FDR] DeltaMass=" << dm.delta_mass
                      << " count=" << dm.count
                      << " mapped=" << dm.mapped_modification << std::endl;
    }
  }

  // Verify that major modifications survive FDR filtering
  vector<pair<string, double>> must_find = {
    {"Oxidation (M)",     15.9949},
    {"Phospho (S/T/Y)",   79.9663},
    {"Acetyl (K/N-term)", 42.0106},
    {"Deamidated (N/Q)",   0.9840},
    {"Dimethyl (K)",      28.0314},
    {"Dioxidation (M)",   31.9898},
    {"HexNAc (S)",       203.0794},
  };

  Size found_count = 0;
  for (const auto& [label, expected_mass] : must_find)
  {
    for (const auto& dm : dm_entries)
    {
      if (fabs(dm.delta_mass - expected_mass) < 0.05 && dm.count >= 10)
      {
        OPENMS_LOG_INFO << "[TEST FDR] Found " << label
                        << ": count=" << dm.count << std::endl;
        found_count++;
        break;
      }
    }
  }
  // At least 4 out of 7 major modifications should survive FDR filtering
  OPENMS_LOG_INFO << "[TEST FDR] Found " << found_count << "/"
                  << must_find.size() << " major modifications" << std::endl;
  TEST_TRUE(found_count >= 4)
}
END_SECTION

START_SECTION(([EXTRA] Closed search baseline))
{
  vector<FASTAFile::FASTAEntry> fasta_db = {
    {"P01", "Test", "MSDEREKVLGFHQRMPNASTICYWDLKEGFVRTHQPSANLDIKCMYKWTE"
                    "RHASGDFLKPIVEQNCTMYRGWSADELKHPFNQGTICMSYREWDAVLKPH"},
  };

  TheoreticalSpectrumGenerator tsg;
  Param tsg_param = tsg.getParameters();
  tsg_param.setValue("add_first_prefix_ion", "true");
  tsg_param.setValue("add_metainfo", "true");
  tsg.setParameters(tsg_param);

  // Generate spectra for known peptides with fixed Carbamidomethyl(C) + variable Oxidation(M)
  vector<string> test_seqs = {
    "VLGFHQR",
    "M(Oxidation)PNASTIC(Carbamidomethyl)YWDLK",
    "EGFVRTHQPSANLDK",
    "PIVEQNC(Carbamidomethyl)TM(Oxidation)YR",
  };

  PeakMap spectra;
  double rt = 100.0;
  for (const auto& seq_str : test_seqs)
  {
    AASequence seq = AASequence::fromString(seq_str);
    int charge = 2;
    MSSpectrum spec;
    tsg.getSpectrum(spec, seq, 1, 1);
    spec.sortByPosition();
    spec.setMSLevel(2);
    spec.setRT(rt);
    rt += 1.0;

    Precursor prec;
    prec.setMZ(seq.getMZ(charge));
    prec.setCharge(charge);
    spec.setPrecursors({prec});
    spec.setNativeID("spectrum=" + String(spectra.size()));

    spectra.addSpectrum(std::move(spec));
  }

  ProSEAlgorithm algo;
  Param p = algo.getParameters();
  p.setValue("precursor:mass_tolerance_lower", 10.0);
  p.setValue("precursor:mass_tolerance_upper", 10.0);
  p.setValue("precursor:mass_tolerance_unit", "ppm");
  p.setValue("fragment:mass_tolerance", 20.0);
  p.setValue("fragment:mass_tolerance_unit", "ppm");
  p.setValue("modifications:fixed", vector<string>{"Carbamidomethyl (C)"});
  p.setValue("modifications:variable", vector<string>{"Oxidation (M)"});
  p.setValue("decoys", "false");
  p.setValue("peptide:min_size", 7);
  p.setValue("peptide:max_size", 40);
  p.setValue("peptide:missed_cleavages", 1);
  algo.setParameters(p);

  vector<ProteinIdentification> prot_ids;
  PeptideIdentificationList pep_ids;
  auto ec = algo.search(spectra, fasta_db, prot_ids, pep_ids);

  TEST_EQUAL(ec == ProSEAlgorithm::ExitCodes::EXECUTION_OK, true)
  TEST_TRUE(pep_ids.size() > 0)
  TEST_EQUAL(prot_ids.size(), 1)
  TEST_EQUAL(prot_ids[0].getSearchEngine(), "ProSE")
}
END_SECTION

START_SECTION(([EXTRA] Ion mobility annotation))
{
  // Create a small protein database
  vector<FASTAFile::FASTAEntry> fasta_db = {
    {"P01", "Test", "MSDEREKVLGFHQRMPNASTICYWDLKEGFVRTHQPSANLDIKCMYKWTE"
                    "RHASGDFLKPIVEQNCTMYRGWSADELKHPFNQGTICMSYREWDAVLKPH"},
  };

  TheoreticalSpectrumGenerator tsg;
  Param tsg_param = tsg.getParameters();
  tsg_param.setValue("add_first_prefix_ion", "true");
  tsg_param.setValue("add_metainfo", "true");
  tsg.setParameters(tsg_param);

  // Generate spectra with drift times set (simulating DDA-PASEF)
  vector<string> test_seqs = {"VLGFHQR", "EGFVRTHQPSANLDIK"};
  vector<double> test_ims = {0.85, 1.12}; // 1/K0 values

  PeakMap spectra;
  double rt = 100.0;
  for (Size i = 0; i < test_seqs.size(); ++i)
  {
    AASequence seq = AASequence::fromString(test_seqs[i]);
    int charge = 2;
    MSSpectrum spec;
    tsg.getSpectrum(spec, seq, 1, 1);
    spec.sortByPosition();
    spec.setMSLevel(2);
    spec.setRT(rt);
    rt += 1.0;

    // Set ion mobility (simulates BrukerTimsFile DDA-PASEF loading)
    spec.setDriftTime(test_ims[i]);
    spec.setDriftTimeUnit(DriftTimeUnit::VSSC);

    Precursor prec;
    prec.setMZ(seq.getMZ(charge));
    prec.setCharge(charge);
    spec.setPrecursors({prec});
    spec.setNativeID("spectrum=" + String(spectra.size()));

    spectra.addSpectrum(std::move(spec));
  }

  ProSEAlgorithm algo;
  Param p = algo.getParameters();
  p.setValue("precursor:mass_tolerance_lower", 10.0);
  p.setValue("precursor:mass_tolerance_upper", 10.0);
  p.setValue("precursor:mass_tolerance_unit", "ppm");
  p.setValue("fragment:mass_tolerance", 20.0);
  p.setValue("fragment:mass_tolerance_unit", "ppm");
  p.setValue("modifications:fixed", vector<string>{});
  p.setValue("modifications:variable", vector<string>{});
  p.setValue("decoys", "false");
  p.setValue("peptide:min_size", 7);
  p.setValue("peptide:max_size", 40);
  p.setValue("peptide:missed_cleavages", 1);
  algo.setParameters(p);

  vector<ProteinIdentification> prot_ids;
  PeptideIdentificationList pep_ids;
  auto ec = algo.search(spectra, fasta_db, prot_ids, pep_ids);

  TEST_EQUAL(ec == ProSEAlgorithm::ExitCodes::EXECUTION_OK, true)
  TEST_EQUAL(prot_ids.size(), 1)
  TEST_EQUAL(pep_ids.size(), 2) // one PSM per input spectrum

  // Verify IM annotation on every PeptideIdentification — both values must appear
  bool has_085 = false, has_112 = false;
  for (const auto& pid : pep_ids)
  {
    TEST_EQUAL(pid.metaValueExists(Constants::UserParam::IM), true)
    double im_val = pid.getMetaValue(Constants::UserParam::IM);
    TEST_TRUE(im_val > 0.0)
    if (fabs(im_val - 0.85) < 1e-6) has_085 = true;
    if (fabs(im_val - 1.12) < 1e-6) has_112 = true;
  }
  TEST_TRUE(has_085) // first spectrum's IM (0.85) found
  TEST_TRUE(has_112) // second spectrum's IM (1.12) found

  // Verify IM unit on ProteinIdentification
  TEST_EQUAL(prot_ids[0].metaValueExists(Constants::UserParam::IM), true)
  TEST_STRING_EQUAL(String(prot_ids[0].getMetaValue(Constants::UserParam::IM)), "1/K0")
}
END_SECTION

START_SECTION(([EXTRA] Edge cases - empty inputs))
{
  // Empty spectra
  {
    PeakMap empty_spectra;
    vector<FASTAFile::FASTAEntry> fasta_db = {{"P01", "Test", "MSDEREKVLGFHQR"}};

    ProSEAlgorithm algo;
    Param p = algo.getParameters();
    p.setValue("decoys", "false");
    algo.setParameters(p);

    vector<ProteinIdentification> prot_ids;
    PeptideIdentificationList pep_ids;
    auto ec = algo.search(empty_spectra, fasta_db, prot_ids, pep_ids);
    TEST_EQUAL(ec == ProSEAlgorithm::ExitCodes::EXECUTION_OK, true)
    TEST_EQUAL(pep_ids.size(), 0)
  }

  // Empty FASTA database: FragmentIndex does not handle empty databases gracefully,
  // so we skip this edge case (it would crash in FragmentIndex::build).
  // This is an existing limitation, not specific to the in-memory overload.
}
END_SECTION

START_SECTION((ExitCodes search(const String &, const String &, std::vector<ProteinIdentification> &, PeptideIdentificationList &) const))
{
  NOT_TESTABLE // tested via TOPP tool
}
END_SECTION

START_SECTION((SearchResult searchWithModificationAnalysis(const String &, const String &, const String &) const))
{
  NOT_TESTABLE // tested via TOPP tool
}
END_SECTION

START_SECTION(([EXTRA] prepareContext + context-based search produces same IDs as single-shot search))
{
  // Build a tiny synthetic dataset where we know the search returns hits.
  // Then verify that:
  //   1. search(spectra, fasta_db, ...) (single-shot, builds index internally)
  //   2. prepareContext(fasta_db) + search(spectra, ctx, ...) (context reuse)
  // produce identical PSM counts (same internal pipeline).

  vector<FASTAFile::FASTAEntry> fasta_db = {
    {"P01", "Test01",
     "MSDEREKVLGFHQRMPNASTICYWDLKEGFVRTHQPSANLDIKCMYKWTE"
     "RHASGDFLKPIVEQNCTMYRGWSADELKHPFNQGTICMSYREWDAVLKPH"},
    {"P02", "Test02",
     "MKAILNHVGSTFREDWQCPYLKMISGDTFNHRVAWQECPLKYMTGISNHFR"
     "DVEWAQCPLKTMIYGSNHFRDVEWAQCPKLIMTGSYNHFRDVEWAQCKPLIM"},
  };

  // Generate spectra from a few tryptic peptides (closed search, perfect matches).
  ProteaseDigestion digester;
  digester.setEnzyme("Trypsin");
  digester.setMissedCleavages(1);

  ModifiedPeptideGenerator::MapToResidueType fixed_mods =
    ModifiedPeptideGenerator::getModifications({"Carbamidomethyl (C)"});

  TheoreticalSpectrumGenerator tsg;
  Param tsg_param = tsg.getParameters();
  tsg_param.setValue("add_first_prefix_ion", "true");
  tsg_param.setValue("add_metainfo", "true");
  tsg.setParameters(tsg_param);

  PeakMap spectra;
  double rt = 100.0;
  for (const auto& entry : fasta_db)
  {
    AASequence protein = AASequence::fromString(entry.sequence);
    vector<AASequence> peptides;
    digester.digest(protein, peptides, 7, 40);
    for (auto& pep : peptides)
    {
      ModifiedPeptideGenerator::applyFixedModifications(fixed_mods, pep);
      if (pep.size() < 8) continue;

      MSSpectrum spec;
      tsg.getSpectrum(spec, pep, 1, 1);
      spec.sortByPosition();
      if (spec.size() < 10) continue;

      spec.setMSLevel(2);
      spec.setRT(rt);
      rt += 0.1;

      Precursor prec;
      prec.setMZ(pep.getMZ(2));
      prec.setCharge(2);
      spec.setPrecursors({prec});
      spec.setNativeID("spectrum=" + String(spectra.size()));

      spectra.addSpectrum(std::move(spec));
    }
  }
  TEST_TRUE(spectra.size() > 5)

  ProSEAlgorithm algo;
  Param p = algo.getParameters();
  p.setValue("precursor:mass_tolerance_lower", 20.0);
  p.setValue("precursor:mass_tolerance_upper", 20.0);
  p.setValue("precursor:mass_tolerance_unit", "ppm");
  p.setValue("fragment:mass_tolerance", 20.0);
  p.setValue("fragment:mass_tolerance_unit", "ppm");
  p.setValue("modifications:fixed", vector<string>{"Carbamidomethyl (C)"});
  p.setValue("modifications:variable", vector<string>{});
  p.setValue("decoys", "false");
  algo.setParameters(p);

  // Path A: single-shot search (builds + tears down the index internally).
  PeakMap spectra_a = spectra; // search() preprocesses in place
  vector<ProteinIdentification> prot_a;
  PeptideIdentificationList pep_a;
  auto ec_a = algo.search(spectra_a, fasta_db, prot_a, pep_a);
  TEST_EQUAL(ec_a == ProSEAlgorithm::ExitCodes::EXECUTION_OK, true)
  TEST_TRUE(pep_a.size() > 0)

  // Path B: prepareContext + context-based search.
  PeakMap spectra_b = spectra;
  ProSEAlgorithm::SearchContext ctx = algo.prepareContext(fasta_db);
  TEST_EQUAL(ctx.fragment_index.isBuild(), true)
  vector<ProteinIdentification> prot_b;
  PeptideIdentificationList pep_b;
  auto ec_b = algo.search(spectra_b, ctx, prot_b, pep_b);
  TEST_EQUAL(ec_b == ProSEAlgorithm::ExitCodes::EXECUTION_OK, true)

  // Both paths must yield the same number of PSMs (the search engine itself is
  // deterministic when decoys are disabled).
  TEST_EQUAL(pep_a.size(), pep_b.size())

  // Compare top-hit sequences spectrum-by-spectrum: they should match exactly.
  TEST_EQUAL(prot_a.size(), prot_b.size())
  for (Size i = 0; i < pep_a.size(); ++i)
  {
    TEST_EQUAL(pep_a[i].getHits().empty(), pep_b[i].getHits().empty())
    if (!pep_a[i].getHits().empty() && !pep_b[i].getHits().empty())
    {
      TEST_STRING_EQUAL(pep_a[i].getHits()[0].getSequence().toString(),
                        pep_b[i].getHits()[0].getSequence().toString())
    }
  }

  // Reusing the same context for a second search must also work.
  PeakMap spectra_c = spectra;
  vector<ProteinIdentification> prot_c;
  PeptideIdentificationList pep_c;
  auto ec_c = algo.search(spectra_c, ctx, prot_c, pep_c);
  TEST_EQUAL(ec_c == ProSEAlgorithm::ExitCodes::EXECUTION_OK, true)
  TEST_EQUAL(pep_c.size(), pep_b.size())
}
END_SECTION

START_SECTION((MultiFileSearchResult searchWithModificationAnalysis(const std::vector<String>&, const std::vector<FASTAFile::FASTAEntry>&, const std::vector<String>&, const String&) const))
{
  // Verify the multi-file in-memory FASTA overload validates input list lengths.
  ProSEAlgorithm algo;
  Param p = algo.getParameters();
  p.setValue("decoys", "false");
  algo.setParameters(p);

  vector<FASTAFile::FASTAEntry> fasta_db = {{"P01", "Test", "MSDEREKVLGFHQRMPNASTICYWDLK"}};
  vector<String> in_files = {"a.mzML", "b.mzML"};
  vector<String> mismatched_base_names = {"a"}; // wrong size

  TEST_EXCEPTION(Exception::InvalidParameter,
                 algo.searchWithModificationAnalysis(in_files, fasta_db, mismatched_base_names, ""))

  // Empty input file list returns INPUT_FILE_EMPTY (no exception).
  auto empty_res = algo.searchWithModificationAnalysis(std::vector<String>{}, fasta_db, std::vector<String>{}, "");
  TEST_EQUAL(empty_res.per_file.empty(), true)
  TEST_EQUAL(empty_res.aggregate.exit_code == ProSEAlgorithm::ExitCodes::INPUT_FILE_EMPTY, true)
}
END_SECTION

START_SECTION((MultiFileSearchResult searchWithModificationAnalysis(const std::vector<String>&, const String&, const std::vector<String>&, const String&) const))
{
  NOT_TESTABLE // tested via TOPP tool (multi-file integration test)
}
END_SECTION

START_SECTION(([EXTRA] PSM annotations - matched ion counts, longest run, fragment annotations))
{
  // Create a small FASTA database with one protein containing a single tryptic peptide
  vector<FASTAFile::FASTAEntry> fasta_db = {
    {"P01", "TestProtein", "PEPTIDEK"}
  };

  // Generate a synthetic MS2 spectrum from a known tryptic peptide
  AASequence peptide = AASequence::fromString("PEPTIDEK");
  TheoreticalSpectrumGenerator tsg;
  Param tsg_param(tsg.getParameters());
  tsg_param.setValue("add_metainfo", "true");
  tsg_param.setValue("add_first_prefix_ion", "true");
  tsg.setParameters(tsg_param);

  PeakSpectrum theo;
  tsg.getSpectrum(theo, peptide, 1, 1);

  // Build a PeakMap with one MS2 spectrum
  PeakMap exp;
  MSSpectrum ms2;
  ms2.setMSLevel(2);
  ms2.setRT(100.0);
  Precursor prec;
  prec.setMZ(peptide.getMZ(2));  // charge 2
  prec.setCharge(2);
  ms2.setPrecursors({prec});

  // Copy theoretical peaks to experimental (perfect match)
  for (const auto& p : theo)
  {
    ms2.emplace_back(p.getMZ(), p.getIntensity());
  }
  ms2.sortByPosition();
  exp.addSpectrum(std::move(ms2));

  // Configure search engine
  ProSEAlgorithm algo;
  Param p = algo.getParameters();
  p.setValue("precursor:mass_tolerance_lower", 10.0);
  p.setValue("precursor:mass_tolerance_upper", 10.0);
  p.setValue("precursor:mass_tolerance_unit", "ppm");
  p.setValue("fragment:mass_tolerance", 10.0);
  p.setValue("fragment:mass_tolerance_unit", "ppm");
  p.setValue("enzyme", "Trypsin");
  p.setValue("peptide:missed_cleavages", 0);
  p.setValue("peptide:min_size", 5);
  p.setValue("peptide:max_size", 40);
  p.setValue("annotate:PSM", vector<string>{"ALL"});
  p.setValue("modifications:fixed", vector<string>{});
  p.setValue("modifications:variable", vector<string>{});
  algo.setParameters(p);

  vector<ProteinIdentification> prot_ids;
  PeptideIdentificationList pep_ids;
  algo.search(exp, fasta_db, prot_ids, pep_ids);

  // Should have found our peptide
  TEST_EQUAL(pep_ids.size(), 1)
  TEST_EQUAL(pep_ids[0].getHits().size() >= 1, true)

  const PeptideHit& hit = pep_ids[0].getHits()[0];
  TEST_EQUAL(hit.getSequence(), peptide)

  // Verify matched ion count annotations exist and are positive
  TEST_EQUAL(hit.metaValueExists(Constants::UserParam::NUM_MATCHED_PEAKS), true)
  TEST_EQUAL(hit.metaValueExists(Constants::UserParam::MATCHED_B_IONS), true)
  TEST_EQUAL(hit.metaValueExists(Constants::UserParam::MATCHED_Y_IONS), true)
  TEST_EQUAL(hit.metaValueExists(Constants::UserParam::LONGEST_PEPTIDE_ION_SEQUENCE), true)

  int num_matched = hit.getMetaValue(Constants::UserParam::NUM_MATCHED_PEAKS);
  int b_ions = hit.getMetaValue(Constants::UserParam::MATCHED_B_IONS);
  int y_ions = hit.getMetaValue(Constants::UserParam::MATCHED_Y_IONS);
  int longest_run = hit.getMetaValue(Constants::UserParam::LONGEST_PEPTIDE_ION_SEQUENCE);

  TEST_EQUAL(num_matched > 0, true)
  TEST_EQUAL(num_matched, b_ions + y_ions)
  TEST_EQUAL(b_ions > 0, true)
  TEST_EQUAL(y_ions > 0, true)

  // Perfect match: longest run should be substantial (peptide length - 1 for one series)
  TEST_EQUAL(longest_run >= 3, true)

  // Verify fragment annotations
  const auto& annotations = hit.getPeakAnnotations();
  TEST_EQUAL(annotations.empty(), false)
  // Each annotation should have mz > 0, non-empty name, and charge >= 1
  for (const auto& ann : annotations)
  {
    TEST_EQUAL(ann.mz > 0.0, true)
    TEST_EQUAL(ann.annotation.empty(), false)
    TEST_EQUAL(ann.charge >= 1, true)
  }
}
END_SECTION

START_SECTION(([EXTRA] calibration preserves asymmetric bias - normal case))
{
  // User sets an asymmetric [20, 30] ppm window (skewed toward a known positive
  // bias). The scattered +7 ppm distribution keeps |shift| < spread so the
  // calibration writeback path runs and the algo-level members get rewritten
  // to the calibrated (cal_lower, cal_upper) window.
  //
  // PLAN DEVIATION: the plan asked for [20, 5] ppm + a uniform +7 ppm shift.
  //   - [lower=20, upper=5] = [-20, +5], so the wrong-match filter at
  //     ProSEAlgorithm.cpp:1584 rejects every +7 ppm hit.
  //   - A uniform shift gives residual spread ~ 1e-6, so |shift| >> spread
  //     and extreme_bias triggers (same pathology as test 9).
  //   - The plan's SimpleSearchEngine_1.mzML fixture does not exist in
  //     src/tests/class_tests/openms/data/. The existing tests in this file
  //     are all synthetic, so we follow that idiom.
  // We preserve the plan's intent (asymmetric user window + non-extreme
  // calibration result) by widening to [20, 30] and scattering the shift.
  //
  // Shift distribution (12 values): median = 7.0, spread ~ 8 ppm,
  //   |shift| (= 7) < spread -> extreme_bias = false.
  const vector<double> ppm_shifts = {
    0.0, 2.0, 4.0, 5.0, 6.0, 7.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0
  };
  PeakMap spectra = build_calibration_spectra_(ppm_shifts);
  auto fasta_db = calibration_fasta_db_();

  ProSEAlgorithm_test algo;
  configure_calibration_params_(algo, /*lower_ppm*/ 20.0, /*upper_ppm*/ 30.0,
                                /*min_psms*/ 3);

  vector<ProteinIdentification> prot_ids;
  PeptideIdentificationList pep_ids;
  auto ec = algo.search(spectra, fasta_db, prot_ids, pep_ids);
  TEST_EQUAL(ec == ProSEAlgorithm::ExitCodes::EXECUTION_OK, true)

  const auto& cal = algo.last_calibration_result_;
  TEST_EQUAL(cal.success, true)
  TEST_EQUAL(cal.extreme_bias, false)
  // Fixture's ppm_shifts are all positive, so the calibration direction must
  // come out positive too. Spread is strictly positive by construction.
  TEST_EQUAL(cal.precursor_shift > 0.0, true)
  TEST_EQUAL(cal.precursor_spread > 0.0, true)
  // |shift| < spread is the precondition for the writeback block (extreme_bias
  // already asserted false above, but state it as a positive numerical check).
  TEST_EQUAL(std::abs(cal.precursor_shift) < cal.precursor_spread, true)

  // Load-bearing check — guards the sign convention in runCalibrationPass_'s
  // writeback block. Under (lower, upper) where signed error e lies in
  // [-upper, +lower], the calibrated window [shift - spread, shift + spread]
  // maps to:
  //   cal_lower = spread + shift   (upper endpoint of the signed window)
  //   cal_upper = spread - shift   (|lower endpoint| of the signed window)
  // A regression of the swap would flip both of these.
  TEST_REAL_SIMILAR(cal.cal_lower, cal.precursor_spread + cal.precursor_shift)
  TEST_REAL_SIMILAR(cal.cal_upper, cal.precursor_spread - cal.precursor_shift)
  // Positive bias => cal_lower > cal_upper. A regression of the swap would
  // produce cal_upper > cal_lower — the assertion below catches that even
  // if the functional identities above are tautologically consistent with a
  // mislabeled spread/shift pair.
  TEST_EQUAL(cal.cal_lower > cal.cal_upper, true)
  // Both tightened from user-configured (20, 30); std::min cap inactive, so
  // the functional identities above are unconstrained.
  TEST_EQUAL(cal.cal_lower < 20.0, true)
  TEST_EQUAL(cal.cal_upper < 30.0, true)
  // Post-search, the tolerance members have been RESTORED to the user-configured
  // values to avoid per-file state leaks in the multi-file wrapper (which reuses a
  // single ProSEAlgorithm instance across files). The calibrated
  // values are observable via last_calibration_result_, which is checked above.
  TEST_REAL_SIMILAR(algo.precursor_mass_tolerance_lower_, 20.0)
  TEST_REAL_SIMILAR(algo.precursor_mass_tolerance_upper_, 30.0)
}
END_SECTION

START_SECTION(([EXTRA] OpenSearchModificationAnalysis received post-calibration tolerance))
{
  // Double-bookkeeping regression guard: OpenSearchModificationAnalysis must be
  // called with the CALIBRATED mod-match tolerance, not the pre-calibration
  // user-configured one.
  //
  // We observe this via the `last_mod_match_tolerance_used_` hook, which captures
  // what computeModMatchTolerance_() returned at the moment the OSMA call fired.
  // Post-search, the tolerance members are restored to user values (see previous
  // test), so calling computeModMatchTolerance_() directly after search() would
  // return the user-configured value — the opposite of what we want to check.
  const vector<double> ppm_shifts = {
    0.0, 2.0, 4.0, 5.0, 6.0, 7.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0
  };
  PeakMap spectra = build_calibration_spectra_(ppm_shifts);
  auto fasta_db = calibration_fasta_db_();

  ProSEAlgorithm_test algo;
  configure_calibration_params_(algo, /*lower_ppm*/ 20.0, /*upper_ppm*/ 30.0,
                                /*min_psms*/ 3);

  vector<ProteinIdentification> prot_ids;
  PeptideIdentificationList pep_ids;
  auto ec = algo.search(spectra, fasta_db, prot_ids, pep_ids);
  TEST_EQUAL(ec == ProSEAlgorithm::ExitCodes::EXECUTION_OK, true)

  const auto& cal = algo.last_calibration_result_;
  TEST_EQUAL(cal.success, true)
  TEST_EQUAL(cal.extreme_bias, false)

  // OSMA must have received the calibrated min(cal_lower, cal_upper) — NOT the
  // user-configured min(20, 30) = 20.
  const double expected = std::min(cal.cal_lower, cal.cal_upper);
  TEST_REAL_SIMILAR(algo.last_mod_match_tolerance_used_, expected)
  TEST_NOT_EQUAL(algo.last_mod_match_tolerance_used_, 20.0)
}
END_SECTION

START_SECTION(([EXTRA] calibration extreme-bias path preserves user bounds))
{
  // Uniform +50 ppm shift → residual median/MAD collapse to ~0 → spread = 1e-6
  // (the floor in runCalibrationPass_). |shift| = 50 >> spread, so extreme_bias
  // triggers and the writeback block is skipped: algo members stay at the
  // user-configured values.
  //
  // User window [100, 100] ppm is wide enough that (a) the candidate look-up
  // finds the theoretical peptide (the +50 ppm error is within the [-100, +100]
  // window), and (b) the wrong-match filter passes every hit.
  const vector<double> ppm_shifts(12, 50.0); // uniform
  PeakMap spectra = build_calibration_spectra_(ppm_shifts);
  auto fasta_db = calibration_fasta_db_();

  ProSEAlgorithm_test algo;
  configure_calibration_params_(algo, /*lower_ppm*/ 100.0, /*upper_ppm*/ 100.0,
                                /*min_psms*/ 3);

  vector<ProteinIdentification> prot_ids;
  PeptideIdentificationList pep_ids;
  auto ec = algo.search(spectra, fasta_db, prot_ids, pep_ids);
  TEST_EQUAL(ec == ProSEAlgorithm::ExitCodes::EXECUTION_OK, true)

  const auto& cal = algo.last_calibration_result_;
  TEST_EQUAL(cal.success, true)
  TEST_EQUAL(cal.extreme_bias, true)
  // User bounds unchanged — no writeback happened.
  TEST_REAL_SIMILAR(algo.precursor_mass_tolerance_lower_, 100.0)
  TEST_REAL_SIMILAR(algo.precursor_mass_tolerance_upper_, 100.0)
}
END_SECTION

START_SECTION(([EXTRA] computeModMatchTolerance_ returns min(lower, upper)))
{
  // Pure unit test — no search, no calibration. Pins the min() reduction rule
  // so a future change to max() or midpoint is caught.
  ProSEAlgorithm_test algo;
  Param p = algo.getParameters();
  p.setValue("precursor:mass_tolerance_unit", "ppm");

  p.setValue("precursor:mass_tolerance_lower", 5.0);
  p.setValue("precursor:mass_tolerance_upper", 50.0);
  algo.setParameters(p);
  TEST_REAL_SIMILAR(algo.computeModMatchTolerance_(), 5.0)

  p.setValue("precursor:mass_tolerance_lower", 50.0);
  p.setValue("precursor:mass_tolerance_upper", 5.0);
  algo.setParameters(p);
  TEST_REAL_SIMILAR(algo.computeModMatchTolerance_(), 5.0)

  // Da unit
  p.setValue("precursor:mass_tolerance_unit", "Da");
  p.setValue("precursor:mass_tolerance_lower", 0.5);
  p.setValue("precursor:mass_tolerance_upper", 2.0);
  algo.setParameters(p);
  TEST_REAL_SIMILAR(algo.computeModMatchTolerance_(), 0.5)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
