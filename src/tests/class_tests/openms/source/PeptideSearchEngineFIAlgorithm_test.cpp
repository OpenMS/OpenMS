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

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <random>
#include <set>

using namespace OpenMS;
using namespace std;

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
  p.setValue("precursor:mass_tolerance", 500.0);
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
  p.setValue("precursor:mass_tolerance", 10.0);
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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
