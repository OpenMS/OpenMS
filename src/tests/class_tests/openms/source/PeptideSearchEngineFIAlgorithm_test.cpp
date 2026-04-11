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
#include <OpenMS/ANALYSIS/ID/PeptideSearchEngineFIAlgorithm.h>
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
// The `friend class PeptideSearchEngineFIAlgorithm_test;` declaration in the
// production header makes the `using` re-exposure below legal.
//
// Note: `fragment_index_` is NOT a member of PeptideSearchEngineFIAlgorithm — it
// lives on the SearchContext and is a local reference during search(). Tests that
// need to observe FragmentIndex state must do so via prepareContext() + the
// context-taking search() overload before/after the restore hook fires.
class PeptideSearchEngineFIAlgorithm_test : public PeptideSearchEngineFIAlgorithm
{
public:
  using PeptideSearchEngineFIAlgorithm::precursor_mass_tolerance_lower_;
  using PeptideSearchEngineFIAlgorithm::precursor_mass_tolerance_upper_;
  using PeptideSearchEngineFIAlgorithm::precursor_mass_tolerance_unit_;
  using PeptideSearchEngineFIAlgorithm::computeModMatchTolerance_;
  using PeptideSearchEngineFIAlgorithm::last_calibration_result_;
  using PeptideSearchEngineFIAlgorithm::CalibrationResult_;
};

START_TEST(PeptideSearchEngineFIAlgorithm, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeptideSearchEngineFIAlgorithm* ptr = nullptr;
PeptideSearchEngineFIAlgorithm* null_ptr = nullptr;

START_SECTION(PeptideSearchEngineFIAlgorithm())
{
  ptr = new PeptideSearchEngineFIAlgorithm();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~PeptideSearchEngineFIAlgorithm())
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
  PeptideSearchEngineFIAlgorithm algo;
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
  TEST_EQUAL(result.exit_code == PeptideSearchEngineFIAlgorithm::ExitCodes::EXECUTION_OK, true)
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
  PeptideSearchEngineFIAlgorithm algo;
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

  TEST_EQUAL(ec == PeptideSearchEngineFIAlgorithm::ExitCodes::EXECUTION_OK, true)
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

  PeptideSearchEngineFIAlgorithm algo;
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

  TEST_EQUAL(ec == PeptideSearchEngineFIAlgorithm::ExitCodes::EXECUTION_OK, true)
  TEST_TRUE(pep_ids.size() > 0)
  TEST_EQUAL(prot_ids.size(), 1)
  TEST_EQUAL(prot_ids[0].getSearchEngine(), "PeptideDataBaseSearchFI")
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

  PeptideSearchEngineFIAlgorithm algo;
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

  TEST_EQUAL(ec == PeptideSearchEngineFIAlgorithm::ExitCodes::EXECUTION_OK, true)
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

    PeptideSearchEngineFIAlgorithm algo;
    Param p = algo.getParameters();
    p.setValue("decoys", "false");
    algo.setParameters(p);

    vector<ProteinIdentification> prot_ids;
    PeptideIdentificationList pep_ids;
    auto ec = algo.search(empty_spectra, fasta_db, prot_ids, pep_ids);
    TEST_EQUAL(ec == PeptideSearchEngineFIAlgorithm::ExitCodes::EXECUTION_OK, true)
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

  PeptideSearchEngineFIAlgorithm algo;
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
  TEST_EQUAL(ec_a == PeptideSearchEngineFIAlgorithm::ExitCodes::EXECUTION_OK, true)
  TEST_TRUE(pep_a.size() > 0)

  // Path B: prepareContext + context-based search.
  PeakMap spectra_b = spectra;
  PeptideSearchEngineFIAlgorithm::SearchContext ctx = algo.prepareContext(fasta_db);
  TEST_EQUAL(ctx.fragment_index.isBuild(), true)
  vector<ProteinIdentification> prot_b;
  PeptideIdentificationList pep_b;
  auto ec_b = algo.search(spectra_b, ctx, prot_b, pep_b);
  TEST_EQUAL(ec_b == PeptideSearchEngineFIAlgorithm::ExitCodes::EXECUTION_OK, true)

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
  TEST_EQUAL(ec_c == PeptideSearchEngineFIAlgorithm::ExitCodes::EXECUTION_OK, true)
  TEST_EQUAL(pep_c.size(), pep_b.size())
}
END_SECTION

START_SECTION((MultiFileSearchResult searchWithModificationAnalysis(const std::vector<String>&, const std::vector<FASTAFile::FASTAEntry>&, const std::vector<String>&, const String&) const))
{
  // Verify the multi-file in-memory FASTA overload validates input list lengths.
  PeptideSearchEngineFIAlgorithm algo;
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
  TEST_EQUAL(empty_res.aggregate.exit_code == PeptideSearchEngineFIAlgorithm::ExitCodes::INPUT_FILE_EMPTY, true)
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
  PeptideSearchEngineFIAlgorithm algo;
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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
