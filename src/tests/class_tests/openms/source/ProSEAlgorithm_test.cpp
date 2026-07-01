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
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>

#include <algorithm>
#include <numeric>
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
  using ProSEAlgorithm::preprocessSpectra_;
  using ProSEAlgorithm::resolveDecoyStrategy_;
  using ProSEAlgorithm::DecoyStrategy_;
  using ProSEAlgorithm::buildDecoyAugmentedDB_;
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
    spec.setNativeID("spectrum=" + StringUtils::toStr(spectra.size()));

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
  p.setValue("decoys", "ignore");
  p.setValue("peptide:min_size", 7);
  p.setValue("peptide:max_size", 40);
  p.setValue("peptide:missed_cleavages", 1);
  algo.setParameters(p);
}

// ---------------------------------------------------------------------------
// Shared synthetic search problem for the protein-FDR contract tests below.
// 10 proteins + many modified-precursor spectra under a wide precursor window, so
// ProSEAlgorithm with decoys=true reliably produces BOTH target and decoy protein
// hits — the prerequisite for exercising picked-protein FDR.
// ---------------------------------------------------------------------------
void buildSyntheticProteinFDRData(std::vector<FASTAFile::FASTAEntry>& fasta_db, PeakMap& spectra)
{
  fasta_db = {
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

  ProteaseDigestion digester;
  digester.setEnzyme("Trypsin");
  digester.setMissedCleavages(1);
  ModifiedPeptideGenerator::MapToResidueType fixed_mods =
    ModifiedPeptideGenerator::getModifications({"Carbamidomethyl (C)"});

  std::vector<AASequence> all_peptides;
  for (const auto& entry : fasta_db)
  {
    AASequence protein = AASequence::fromString(entry.sequence);
    std::vector<AASequence> peptides;
    digester.digest(protein, peptides, 7, 40);
    for (auto& pep : peptides)
    {
      ModifiedPeptideGenerator::applyFixedModifications(fixed_mods, pep);
      all_peptides.push_back(std::move(pep));
    }
  }

  const std::vector<double> shift_masses = {15.9949, 79.9663, 42.0106, 0.9840, 28.0314, 31.9898, 203.0794};
  std::mt19937 rng(42);
  TheoreticalSpectrumGenerator tsg;
  Param tsg_param = tsg.getParameters();
  tsg_param.setValue("add_first_prefix_ion", "true");
  tsg_param.setValue("add_metainfo", "true");
  tsg.setParameters(tsg_param);

  double rt = 100.0;
  const Size target_per_shift = 300;
  for (double shift : shift_masses)
  {
    std::vector<size_t> indices(all_peptides.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), rng);
    Size created = 0;
    for (size_t idx : indices)
    {
      if (created >= target_per_shift) break;
      const AASequence& pep = all_peptides[idx];
      if (pep.size() < 8) continue;
      int charge = 2 + (int)(rng() % 3);
      MSSpectrum spec;
      tsg.getSpectrum(spec, pep, 1, std::min(charge - 1, 2));
      spec.sortByPosition();
      if (spec.size() < 10) continue;
      spec.setMSLevel(2);
      spec.setRT(rt);
      rt += 0.1;
      double shifted_mz = pep.getMZ(charge) + shift / (double)charge;
      Precursor prec;
      prec.setMZ(shifted_mz);
      prec.setCharge(charge);
      spec.setPrecursors({prec});
      spec.setNativeID("spectrum=" + StringUtils::toStr(spectra.size()));
      spectra.addSpectrum(std::move(spec));
      created++;
    }
  }
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

START_SECTION(([EXTRA] default mass tolerances))
{
  ProSEAlgorithm algo;
  Param p = algo.getParameters();
  TEST_REAL_SIMILAR((double)p.getValue("precursor:mass_tolerance_lower"), 10.0)
  TEST_REAL_SIMILAR((double)p.getValue("precursor:mass_tolerance_upper"), 10.0)
  TEST_STRING_EQUAL(p.getValue("precursor:mass_tolerance_unit").toString(), "ppm")
  TEST_REAL_SIMILAR((double)p.getValue("fragment:mass_tolerance"), 20.0)
  TEST_STRING_EQUAL(p.getValue("fragment:mass_tolerance_unit").toString(), "ppm")
}
END_SECTION

START_SECTION(([EXTRA] resolveDecoyStrategy_ / buildDecoyAugmentedDB_: auto/generate/ignore))
{
  // Target+decoy database (50% decoys, conventional DECOY_ prefix).
  const std::vector<FASTAFile::FASTAEntry> td_db = {
    FASTAFile::FASTAEntry("sp|P1|A", "", "PEPTIDEKAAR"),
    FASTAFile::FASTAEntry("sp|P2|B", "", "SAMPLERPEPTIDEK"),
    FASTAFile::FASTAEntry("DECOY_sp|P1|A", "", "RAAKEDITPEP"),
    FASTAFile::FASTAEntry("DECOY_sp|P2|B", "", "KEDITPEPRELPMAS") };
  // Target-only database.
  const std::vector<FASTAFile::FASTAEntry> t_db = {
    FASTAFile::FASTAEntry("sp|P1|A", "", "PEPTIDEKAAR"),
    FASTAFile::FASTAEntry("sp|P2|B", "", "SAMPLERPEPTIDEK") };

  auto count_prefix = [](const std::vector<FASTAFile::FASTAEntry>& db, const std::string& pre)
  {
    Size n = 0;
    for (const auto& e : db) if (e.identifier.rfind(pre, 0) == 0) ++n;
    return n;
  };

  // --- auto: reuse existing decoys (detected), do not generate -------------
  {
    ProSEAlgorithm_test algo;
    Param p = algo.getParameters();
    p.setValue("decoys", "auto");
    algo.setParameters(p);
    ProSEAlgorithm_test::DecoyStrategy_ s = algo.resolveDecoyStrategy_(td_db);
    TEST_EQUAL(s.generate, false)
    TEST_EQUAL(s.strip_existing, false)
    TEST_EQUAL(s.have_decoys, true)
    TEST_STRING_EQUAL(s.decoy_string, "DECOY_")
    TEST_EQUAL(s.is_prefix, true)
    // DB is searched unchanged.
    std::vector<FASTAFile::FASTAEntry> built = algo.buildDecoyAugmentedDB_(td_db, s);
    TEST_EQUAL(built.size(), 4)
    TEST_EQUAL(count_prefix(built, "DECOY_"), 2)
  }

  // --- auto: no decoys present -> generate them ---------------------------
  {
    ProSEAlgorithm_test algo;
    Param p = algo.getParameters();
    p.setValue("decoys", "auto");
    algo.setParameters(p);
    ProSEAlgorithm_test::DecoyStrategy_ s = algo.resolveDecoyStrategy_(t_db);
    TEST_EQUAL(s.generate, true)
    TEST_EQUAL(s.strip_existing, false)
    TEST_EQUAL(s.have_decoys, true)
    TEST_STRING_EQUAL(s.decoy_string, "DECOY_")
    std::vector<FASTAFile::FASTAEntry> built = algo.buildDecoyAugmentedDB_(t_db, s);
    TEST_EQUAL(built.size(), 4)             // 2 targets + 2 generated decoys
    TEST_EQUAL(count_prefix(built, "DECOY_"), 2)
  }

  // --- ignore: strip existing decoys, search targets only -----------------
  {
    ProSEAlgorithm_test algo;
    Param p = algo.getParameters();
    p.setValue("decoys", "ignore");
    algo.setParameters(p);
    ProSEAlgorithm_test::DecoyStrategy_ s = algo.resolveDecoyStrategy_(td_db);
    TEST_EQUAL(s.generate, false)
    TEST_EQUAL(s.strip_existing, true)
    TEST_EQUAL(s.have_decoys, false)
    std::vector<FASTAFile::FASTAEntry> built = algo.buildDecoyAugmentedDB_(td_db, s);
    TEST_EQUAL(built.size(), 2)             // decoys removed
    TEST_EQUAL(count_prefix(built, "DECOY_"), 0)
  }

  // --- generate: strip pre-existing decoys, then regenerate from targets ---
  {
    ProSEAlgorithm_test algo;
    Param p = algo.getParameters();
    p.setValue("decoys", "generate");
    algo.setParameters(p);
    ProSEAlgorithm_test::DecoyStrategy_ s = algo.resolveDecoyStrategy_(td_db);
    TEST_EQUAL(s.generate, true)
    TEST_EQUAL(s.strip_existing, true)
    TEST_EQUAL(s.have_decoys, true)
    std::vector<FASTAFile::FASTAEntry> built = algo.buildDecoyAugmentedDB_(td_db, s);
    TEST_EQUAL(built.size(), 4)             // 2 targets + 2 freshly generated
    TEST_EQUAL(count_prefix(built, "DECOY_"), 2)
  }

  // --- custom marker outside the common vocabulary: literal fall-back -----
  {
    ProSEAlgorithm_test algo;
    Param p = algo.getParameters();
    p.setValue("decoys", "auto");
    p.setValue("decoy_prefix", "BOGUS_");
    algo.setParameters(p);
    const std::vector<FASTAFile::FASTAEntry> custom_db = {
      FASTAFile::FASTAEntry("sp|P1|A", "", "PEPTIDEKAAR"),
      FASTAFile::FASTAEntry("BOGUS_sp|P1|A", "", "RAAKEDITPEP") };
    ProSEAlgorithm_test::DecoyStrategy_ s = algo.resolveDecoyStrategy_(custom_db);
    TEST_EQUAL(s.generate, false)           // existing decoys recognised via fall-back
    TEST_EQUAL(s.have_decoys, true)
    TEST_STRING_EQUAL(s.decoy_string, "BOGUS_")
    TEST_EQUAL(s.is_prefix, true)
  }

  // --- auto: reuse decoys detected by a SUFFIX marker (prefix/suffix aware) -----
  // Headline #9634 feature: decoys can be marked as a suffix (e.g. from DecoyDatabase
  // with -decoy_string_position suffix). DecoyHelper detects it; ProSE must reuse them
  // (not double-generate) and thread is_prefix=false through the whole FDR chain.
  {
    ProSEAlgorithm_test algo;
    Param p = algo.getParameters();
    p.setValue("decoys", "auto");
    algo.setParameters(p);
    const std::vector<FASTAFile::FASTAEntry> suffix_db = {
      FASTAFile::FASTAEntry("sp|P1|A", "", "PEPTIDEKAAR"),
      FASTAFile::FASTAEntry("sp|P2|B", "", "SAMPLERPEPTIDEK"),
      FASTAFile::FASTAEntry("sp|P1|A_decoy", "", "RAAKEDITPEP"),
      FASTAFile::FASTAEntry("sp|P2|B_decoy", "", "KEDITPEPRELPMAS") };
    ProSEAlgorithm_test::DecoyStrategy_ s = algo.resolveDecoyStrategy_(suffix_db);
    TEST_EQUAL(s.generate, false)           // reuse existing, do not generate
    TEST_EQUAL(s.strip_existing, false)
    TEST_EQUAL(s.have_decoys, true)
    TEST_EQUAL(s.is_prefix, false)          // detected as a SUFFIX marker
    // searched unchanged (no double-generation -> no *_decoy_decoy entries).
    std::vector<FASTAFile::FASTAEntry> built = algo.buildDecoyAugmentedDB_(suffix_db, s);
    TEST_EQUAL(built.size(), 4)
  }
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
      spec.setNativeID("spectrum=" + StringUtils::toStr(spectra.size()));

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
  p.setValue("decoys", "ignore");
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
      spec.setNativeID("spectrum=" + StringUtils::toStr(spectra.size()));

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
  p.setValue("decoys", "auto");  // Enable decoys for FDR
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
    spec.setNativeID("spectrum=" + StringUtils::toStr(spectra.size()));

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
  p.setValue("decoys", "ignore");
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

START_SECTION(([EXTRA] Closed search with c/z ions toggled - ETD-style fragmentation))
{
  // ProSE can score c/z fragment ions (e.g. ETD/ECD data) via the
  // ions:add_c_ions / ions:add_z_ions toggles. Build spectra that contain ONLY
  // c/z ions and confirm a c/z-enabled search identifies the peptides, while a
  // default (b/y) search on the same spectra does not -- the c/z peaks are
  // shifted ~16-17 Da from b/y and cannot be matched as b/y.
  vector<FASTAFile::FASTAEntry> fasta_db = {
    {"P01", "Test", "MSDEREKVLGFHQRMPNASTICYWDLKEGFVRTHQPSANLDIKCMYKWTE"
                    "RHASGDFLKPIVEQNCTMYRGWSADELKHPFNQGTICMSYREWDAVLKPH"},
  };

  // TheoreticalSpectrumGenerator configured to emit c/z ions only
  TheoreticalSpectrumGenerator tsg;
  Param tsg_param = tsg.getParameters();
  tsg_param.setValue("add_b_ions", "false");
  tsg_param.setValue("add_y_ions", "false");
  tsg_param.setValue("add_c_ions", "true");
  tsg_param.setValue("add_z_ions", "true");
  tsg_param.setValue("add_metainfo", "true");
  tsg.setParameters(tsg_param);

  // fully-tryptic peptides of the protein with no C/M (no fixed/variable mods)
  vector<string> test_seqs = { "VLGFHQR", "THQPSANLDIK" };

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
    spec.setNativeID("spectrum=" + StringUtils::toStr(spectra.size()));
    spectra.addSpectrum(std::move(spec));
  }

  auto run_search = [&](bool enable_cz) {
    ProSEAlgorithm algo;
    Param p = algo.getParameters();
    p.setValue("precursor:mass_tolerance_lower", 10.0);
    p.setValue("precursor:mass_tolerance_upper", 10.0);
    p.setValue("precursor:mass_tolerance_unit", "ppm");
    p.setValue("fragment:mass_tolerance", 20.0);
    p.setValue("fragment:mass_tolerance_unit", "ppm");
    p.setValue("modifications:fixed", vector<string>{});
    p.setValue("modifications:variable", vector<string>{});
    p.setValue("decoys", "ignore");
    p.setValue("peptide:min_size", 7);
    p.setValue("peptide:max_size", 40);
    p.setValue("peptide:missed_cleavages", 1);
    if (enable_cz)
    {
      p.setValue("ions:add_b_ions", "false");
      p.setValue("ions:add_y_ions", "false");
      p.setValue("ions:add_c_ions", "true");
      p.setValue("ions:add_z_ions", "true");
    }
    algo.setParameters(p);
    vector<ProteinIdentification> prot_ids;
    PeptideIdentificationList pep_ids;
    auto ec = algo.search(spectra, fasta_db, prot_ids, pep_ids);
    TEST_EQUAL(ec == ProSEAlgorithm::ExitCodes::EXECUTION_OK, true)
    return pep_ids;
  };

  // (1) c/z-enabled search identifies the peptides from the c/z spectra
  PeptideIdentificationList cz_ids = run_search(true);
  TEST_TRUE(cz_ids.size() > 0)
  std::set<std::string> found;
  for (const auto& pid : cz_ids)
    for (const auto& hit : pid.getHits())
      found.insert(hit.getSequence().toUnmodifiedString());
  TEST_EQUAL(found.count("VLGFHQR") + found.count("THQPSANLDIK") > 0, true)

  // (2) a default (b/y) search on the same c/z spectra matches nothing
  PeptideIdentificationList by_ids = run_search(false);
  Size by_hits = 0;
  for (const auto& pid : by_ids) by_hits += pid.getHits().size();
  TEST_EQUAL(by_hits, 0)
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
    spec.setNativeID("spectrum=" + StringUtils::toStr(spectra.size()));

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
  p.setValue("decoys", "ignore");
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
  TEST_STRING_EQUAL(StringUtils::toStr(prot_ids[0].getMetaValue(Constants::UserParam::IM)), "1/K0")
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
    p.setValue("decoys", "ignore");
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

START_SECTION((ExitCodes search(const std::string &, const std::string &, std::vector<ProteinIdentification> &, PeptideIdentificationList &) const))
{
  // The single-file (file-path) search applies protein-level picked FDR, because a single
  // input file IS a complete experiment. This locks the valid single-file protein-FDR path
  // that the ProSE TOPP tool relies on for 1-input runs (see the single-file block in ProSE.cpp).
  std::vector<FASTAFile::FASTAEntry> fasta_db;
  PeakMap spectra;
  buildSyntheticProteinFDRData(fasta_db, spectra);
  TEST_TRUE(spectra.size() > 500)

  std::string tmp_mzml;
  NEW_TMP_FILE(tmp_mzml)
  tmp_mzml += ".mzML";
  FileHandler().storeExperiment(tmp_mzml, spectra, {FileTypes::MZML});
  std::string tmp_fasta;
  NEW_TMP_FILE(tmp_fasta)
  tmp_fasta += ".fasta";
  FASTAFile().store(tmp_fasta, fasta_db);

  ProSEAlgorithm algo;
  Param p = algo.getParameters();
  p.setValue("precursor:mass_tolerance_lower", 500.0);
  p.setValue("precursor:mass_tolerance_upper", 500.0);
  p.setValue("precursor:mass_tolerance_unit", "Da");
  p.setValue("fragment:mass_tolerance", 20.0);
  p.setValue("fragment:mass_tolerance_unit", "ppm");
  p.setValue("modifications:fixed", std::vector<std::string>{"Carbamidomethyl (C)"});
  p.setValue("decoys", "generate");
  p.setValue("FDR:PSM", 0.05);
  p.setValue("FDR:protein", 0.5);   // lenient: keep proteins but exercise the picked-FDR path
  algo.setParameters(p);

  std::vector<ProteinIdentification> prot_ids;
  PeptideIdentificationList pep_ids;
  auto ec = algo.search(tmp_mzml, tmp_fasta, prot_ids, pep_ids);
  TEST_EQUAL(ec == ProSEAlgorithm::ExitCodes::EXECUTION_OK, true)
  TEST_EQUAL(prot_ids.size(), 1)
  TEST_TRUE(prot_ids[0].getHits().size() > 0)

  // Protein FDR ran: picked-protein FDR + cleanup removes the decoy proteins from the report.
  Size decoy_proteins = 0;
  for (const auto& ph : prot_ids[0].getHits())
  {
    if (ph.getAccession().rfind("DECOY_", 0) == 0) { ++decoy_proteins; }
  }
  TEST_EQUAL(decoy_proteins, 0)

  // The FDR-filtered result must be a valid idXML: storing throws on dangling protein
  // references (groups or peptide evidence pointing at removed decoy proteins).
  std::string tmp_out;
  NEW_TMP_FILE(tmp_out)
  tmp_out += ".idXML";
  FileHandler().storeIdentifications(tmp_out, prot_ids, pep_ids, {FileTypes::IDXML});
  std::vector<ProteinIdentification> rprot;
  PeptideIdentificationList rpep;
  FileHandler().loadIdentifications(tmp_out, rprot, rpep, {FileTypes::IDXML});
  TEST_EQUAL(rprot.size(), 1)
  Size reloaded_decoys = 0;
  for (const auto& ph : rprot[0].getHits()) { if (ph.getAccession().rfind("DECOY_", 0) == 0) { ++reloaded_decoys; } }
  TEST_EQUAL(reloaded_decoys, 0)
}
END_SECTION

START_SECTION(([EXTRA] file-based single-file search retains decoys when protein FDR is OFF))
{
  // Decoy reporting is tied to protein-level FDR, NOT to PSM-level FDR. With FDR:protein==0
  // the single-file (file-path) search must RETAIN decoys after PSM filtering: they are the
  // intermediate evidence a later/global protein FDR or cross-file merge needs. (FDR:protein>0
  // finalizes and removes them — see the section above.) This pins the decoupling of PSM-level
  // FDR from decoy removal.
  std::vector<FASTAFile::FASTAEntry> fasta_db;
  PeakMap spectra;
  buildSyntheticProteinFDRData(fasta_db, spectra);

  std::string tmp_mzml;
  NEW_TMP_FILE(tmp_mzml)
  tmp_mzml += ".mzML";
  FileHandler().storeExperiment(tmp_mzml, spectra, {FileTypes::MZML});
  std::string tmp_fasta;
  NEW_TMP_FILE(tmp_fasta)
  tmp_fasta += ".fasta";
  FASTAFile().store(tmp_fasta, fasta_db);

  ProSEAlgorithm algo;
  Param p = algo.getParameters();
  p.setValue("precursor:mass_tolerance_lower", 500.0);
  p.setValue("precursor:mass_tolerance_upper", 500.0);
  p.setValue("precursor:mass_tolerance_unit", "Da");
  p.setValue("fragment:mass_tolerance", 20.0);
  p.setValue("fragment:mass_tolerance_unit", "ppm");
  p.setValue("modifications:fixed", std::vector<std::string>{"Carbamidomethyl (C)"});
  p.setValue("decoys", "generate");
  p.setValue("FDR:PSM", 0.5);       // PSM filtering ON (lenient, so decoys survive the q-value cut) ...
  p.setValue("FDR:protein", 0.0);   // ... but protein FDR OFF -> decoys must be retained
  algo.setParameters(p);

  std::vector<ProteinIdentification> prot_ids;
  PeptideIdentificationList pep_ids;
  auto ec = algo.search(tmp_mzml, tmp_fasta, prot_ids, pep_ids);
  TEST_EQUAL(ec == ProSEAlgorithm::ExitCodes::EXECUTION_OK, true)
  TEST_EQUAL(prot_ids.size(), 1)

  // Decoy proteins survive (no protein-FDR finalization happened).
  Size decoy_proteins = 0;
  for (const auto& ph : prot_ids[0].getHits())
  {
    if (ph.getAccession().rfind("DECOY_", 0) == 0) { ++decoy_proteins; }
  }
  TEST_TRUE(decoy_proteins > 0)

  // Decoy PSMs survive PSM-level FDR filtering (PSM FDR annotates + filters, never strips decoys).
  Size decoy_psms = 0;
  for (const auto& pid : pep_ids)
  {
    for (const auto& hit : pid.getHits())
    {
      if (hit.metaValueExists("target_decoy")
          && hit.getMetaValue("target_decoy").toString().find("decoy") != std::string::npos) { ++decoy_psms; }
    }
  }
  TEST_TRUE(decoy_psms > 0)

  // The decoy-retaining result is still valid idXML (stores + reloads).
  std::string tmp_out;
  NEW_TMP_FILE(tmp_out)
  tmp_out += ".idXML";
  FileHandler().storeIdentifications(tmp_out, prot_ids, pep_ids, {FileTypes::IDXML});
  std::vector<ProteinIdentification> rprot;
  PeptideIdentificationList rpep;
  FileHandler().loadIdentifications(tmp_out, rprot, rpep, {FileTypes::IDXML});
  TEST_EQUAL(rprot.size(), 1)
}
END_SECTION

START_SECTION(([EXTRA] in-memory search applies PSM-level FDR only, never protein FDR))
{
  // Per-file / multi-file searches must NOT apply protein FDR: FDR does not compose across
  // runs, so picked-protein FDR is valid only on a COMPLETE set (a single file, or the merged
  // aggregate). This pins the "PSM-only" contract of the in-memory search() overload used by
  // the multi-file wrapper — applying protein FDR per file would inflate the combined FDR.
  std::vector<FASTAFile::FASTAEntry> fasta_db;
  PeakMap spectra;
  buildSyntheticProteinFDRData(fasta_db, spectra);

  ProSEAlgorithm algo;
  Param p = algo.getParameters();
  p.setValue("precursor:mass_tolerance_lower", 500.0);
  p.setValue("precursor:mass_tolerance_upper", 500.0);
  p.setValue("precursor:mass_tolerance_unit", "Da");
  p.setValue("fragment:mass_tolerance", 20.0);
  p.setValue("fragment:mass_tolerance_unit", "ppm");
  p.setValue("modifications:fixed", std::vector<std::string>{"Carbamidomethyl (C)"});
  p.setValue("decoys", "generate");
  p.setValue("FDR:PSM", 0.0);       // no PSM filtering, so decoys are retained...
  p.setValue("FDR:protein", 0.5);   // ...and this overload must NOT remove them via protein FDR
  algo.setParameters(p);

  std::vector<ProteinIdentification> prot_ids;
  PeptideIdentificationList pep_ids;
  auto ec = algo.search(spectra, fasta_db, prot_ids, pep_ids);
  TEST_EQUAL(ec == ProSEAlgorithm::ExitCodes::EXECUTION_OK, true)
  TEST_EQUAL(prot_ids.size(), 1)

  // Protein FDR was NOT applied by this overload: decoy proteins survive (picked-protein FDR
  // would have removed them). That is the multi-file/per-file path's intended contract.
  Size decoy_proteins = 0;
  for (const auto& ph : prot_ids[0].getHits())
  {
    if (ph.getAccession().rfind("DECOY_", 0) == 0) { ++decoy_proteins; }
  }
  TEST_TRUE(decoy_proteins > 0)
}
END_SECTION

START_SECTION(([EXTRA] in-memory search retains decoys after PSM-level FDR filtering))
{
  // PSM-level FDR must NOT remove decoys (decoupled from decoy removal): the in-memory search()
  // overload produces per-file/multi-file results that a later protein FDR or cross-file merge
  // relies on having decoys for. With FDR:PSM>0 and FDR:protein==0, decoy PSMs that pass the
  // q-value threshold are retained (previously they were stripped here).
  std::vector<FASTAFile::FASTAEntry> fasta_db;
  PeakMap spectra;
  buildSyntheticProteinFDRData(fasta_db, spectra);

  ProSEAlgorithm algo;
  Param p = algo.getParameters();
  p.setValue("precursor:mass_tolerance_lower", 500.0);
  p.setValue("precursor:mass_tolerance_upper", 500.0);
  p.setValue("precursor:mass_tolerance_unit", "Da");
  p.setValue("fragment:mass_tolerance", 20.0);
  p.setValue("fragment:mass_tolerance_unit", "ppm");
  p.setValue("modifications:fixed", std::vector<std::string>{"Carbamidomethyl (C)"});
  p.setValue("decoys", "generate");
  p.setValue("FDR:PSM", 0.5);       // PSM filtering ON (lenient) ...
  p.setValue("FDR:protein", 0.0);   // ... protein FDR OFF -> decoys retained
  algo.setParameters(p);

  std::vector<ProteinIdentification> prot_ids;
  PeptideIdentificationList pep_ids;
  auto ec = algo.search(spectra, fasta_db, prot_ids, pep_ids);
  TEST_EQUAL(ec == ProSEAlgorithm::ExitCodes::EXECUTION_OK, true)
  TEST_EQUAL(prot_ids.size(), 1)

  // Decoy PSMs survive PSM-level FDR (the contract this overload now pins).
  Size decoy_psms = 0;
  for (const auto& pid : pep_ids)
  {
    for (const auto& hit : pid.getHits())
    {
      if (hit.metaValueExists("target_decoy")
          && hit.getMetaValue("target_decoy").toString().find("decoy") != std::string::npos) { ++decoy_psms; }
    }
  }
  TEST_TRUE(decoy_psms > 0)
}
END_SECTION

START_SECTION((SearchResult searchWithModificationAnalysis(const std::string &, const std::string &, const std::string &) const))
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
      spec.setNativeID("spectrum=" + StringUtils::toStr(spectra.size()));

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
  p.setValue("decoys", "ignore");
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

START_SECTION((MultiFileSearchResult searchWithModificationAnalysis(const std::vector<std::string>&, const std::vector<FASTAFile::FASTAEntry>&, const std::vector<std::string>&, const std::string&) const))
{
  // Verify the multi-file in-memory FASTA overload validates input list lengths.
  ProSEAlgorithm algo;
  Param p = algo.getParameters();
  p.setValue("decoys", "ignore");
  algo.setParameters(p);

  vector<FASTAFile::FASTAEntry> fasta_db = {{"P01", "Test", "MSDEREKVLGFHQRMPNASTICYWDLK"}};
  vector<std::string> in_files = {"a.mzML", "b.mzML"};
  vector<std::string> mismatched_base_names = {"a"}; // wrong size

  TEST_EXCEPTION(Exception::InvalidParameter,
                 algo.searchWithModificationAnalysis(in_files, fasta_db, mismatched_base_names, ""))

  // Empty input file list returns INPUT_FILE_EMPTY (no exception).
  auto empty_res = algo.searchWithModificationAnalysis(std::vector<std::string>{}, fasta_db, std::vector<std::string>{}, "");
  TEST_EQUAL(empty_res.per_file.empty(), true)
  TEST_EQUAL(empty_res.aggregate.exit_code == ProSEAlgorithm::ExitCodes::INPUT_FILE_EMPTY, true)
}
END_SECTION

START_SECTION((MultiFileSearchResult searchWithModificationAnalysis(const std::vector<std::string>&, const std::string&, const std::vector<std::string>&, const std::string&) const))
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
  TEST_EQUAL(hit.metaValueExists(Constants::UserParam::MATCHED_PREFIX_IONS), true)
  TEST_EQUAL(hit.metaValueExists(Constants::UserParam::MATCHED_SUFFIX_IONS), true)
  TEST_EQUAL(hit.metaValueExists(Constants::UserParam::LONGEST_PEPTIDE_ION_SEQUENCE), true)

  int num_matched = hit.getMetaValue(Constants::UserParam::NUM_MATCHED_PEAKS);
  int prefix_ions = hit.getMetaValue(Constants::UserParam::MATCHED_PREFIX_IONS);
  int suffix_ions = hit.getMetaValue(Constants::UserParam::MATCHED_SUFFIX_IONS);
  int longest_run = hit.getMetaValue(Constants::UserParam::LONGEST_PEPTIDE_ION_SEQUENCE);

  TEST_EQUAL(num_matched > 0, true)
  TEST_EQUAL(num_matched, prefix_ions + suffix_ions)
  TEST_EQUAL(prefix_ions > 0, true)
  TEST_EQUAL(suffix_ions > 0, true)

  // Perfect match: longest run should be substantial (peptide length - 1 for one series)
  TEST_EQUAL(longest_run >= 3, true)

  // Delta score is emitted on every retained hit. With a single database peptide,
  // there is no competing candidate, so delta = full score (same "no competition
  // = maximum delta" convention as Sage/MSFragger).
  TEST_EQUAL(hit.metaValueExists(Constants::UserParam::DELTA_SCORE), true)
  double delta = hit.getMetaValue(Constants::UserParam::DELTA_SCORE);
  TEST_REAL_SIMILAR(delta, hit.getScore())

  // MIC = sum of experimental intensities over matched peaks. The synthetic
  // spectrum copies theoretical peaks into the experimental one, so MIC should
  // equal the sum of theoretical peak intensities. Verifies the MIC code
  // accumulates exactly once per matched peak (no double-counting).
  TEST_EQUAL(hit.metaValueExists(Constants::UserParam::MATCHED_ION_CURRENT), true)
  double expected_mic = 0.0;
  for (const auto& p : theo) { expected_mic += p.getIntensity(); }
  double mic = hit.getMetaValue(Constants::UserParam::MATCHED_ION_CURRENT);
  TEST_REAL_SIMILAR(mic, expected_mic)

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

  // Positive bias => cal_lower > cal_upper. Under the (lower, upper) convention
  // signed error e = observed - theoretical lies in [-cal_upper, +cal_lower]; a
  // strictly-positive bias means the +99.5% quantile exceeds the |-0.5% quantile|,
  // so cal_lower (= max positive error) must exceed cal_upper (= |max negative error|).
  // A regression that swapped the endpoints would flip the ordering.
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

START_SECTION(([EXTRA] preprocessSpectra_ never aborts; gates deisotoping on the Deisotoper limit (OpenMS#9619)))
{
  // Regression for OpenMS#9619: preprocessSpectra_ must never let Deisotoper throw
  // inside its OpenMP region (an escaping exception calls std::terminate). It gates
  // the Deisotoper call on Deisotoper::isToleranceSupported(), so even
  // deisotope_requested=true with a low-resolution (out-of-range) tolerance is a
  // safe no-op rather than an abort. Mode resolution (auto/true/false) is covered
  // via the param in the next section.
  auto make_exp = []()
  {
    PeakMap exp;
    MSSpectrum s;
    s.setMSLevel(2);
    s.setRT(1.0);
    Precursor prec;
    prec.setMZ(500.0);
    prec.setCharge(2);
    s.getPrecursors().push_back(prec);
    for (double mz : {110.07, 120.08, 130.10, 200.10, 201.10, 300.20, 350.25, 500.30})
    {
      Peak1D p;
      p.setMZ(mz);
      p.setIntensity(1000.0f);
      s.push_back(p);
    }
    exp.addSpectrum(s);
    return exp;
  };

  // Low-resolution tolerance: requested true OR false -> never throws (deisotoping
  // is skipped because the tolerance is out of the Deisotoper's supported range).
  {
    PeakMap exp = make_exp();
    ProSEAlgorithm_test::preprocessSpectra_(exp, 0.5, false, true, 0, 20);
    TEST_EQUAL(exp.size(), 1)
    TEST_EQUAL(exp[0].empty(), false)
  }
  {
    PeakMap exp = make_exp();
    ProSEAlgorithm_test::preprocessSpectra_(exp, 150.0, true, true, 0, 20);
    TEST_EQUAL(exp.size(), 1)
  }
  {
    PeakMap exp = make_exp();
    ProSEAlgorithm_test::preprocessSpectra_(exp, 0.5, false, false, 0, 20);
    TEST_EQUAL(exp.size(), 1)
  }

  // High-resolution tolerance: requested true -> deisotoping path runs (no throw);
  // requested false -> skipped.
  {
    PeakMap exp = make_exp();
    ProSEAlgorithm_test::preprocessSpectra_(exp, 0.05, false, true, 0, 20);
    TEST_EQUAL(exp.size(), 1)
  }
  {
    PeakMap exp = make_exp();
    ProSEAlgorithm_test::preprocessSpectra_(exp, 20.0, true, false, 0, 20);
    TEST_EQUAL(exp.size(), 1)
  }
}
END_SECTION

START_SECTION(([EXTRA] auto peak retention (peaks:keep_n=0) is resolution-aware))
{
  // Low-resolution fragment tolerances admit many spurious matches; auto retention keeps far
  // fewer peaks at low-res than at high-res (where behavior is unchanged). A dense spectrum so
  // the cap actually bites.
  auto dense = []() {
    PeakMap exp; MSSpectrum s; s.setMSLevel(2);
    Precursor p; p.setMZ(800.0); p.setCharge(2); s.setPrecursors({p}); s.setRT(1.0);
    for (int i = 0; i < 500; ++i) { Peak1D pk; pk.setMZ(150.0 + i * 3.0); pk.setIntensity(1.0 + (i % 50)); s.push_back(pk); }
    s.sortByPosition(); exp.addSpectrum(s); return exp;
  };
  PeakMap hi = dense();  // high-res (0.02 Da, within deisotoper range) -> legacy cap (400)
  ProSEAlgorithm_test::preprocessSpectra_(hi, 0.02, false, false, 0, 20);
  PeakMap lo = dense();  // low-res (0.5 Da) -> auto formula -> ~80
  ProSEAlgorithm_test::preprocessSpectra_(lo, 0.5, false, false, 0, 20);
  TEST_TRUE(lo[0].size() < hi[0].size())   // low-res retains strictly fewer peaks
  TEST_TRUE(lo[0].size() <= 90)            // auto cap at 0.5 Da is 80 (+ headroom)
  TEST_TRUE(lo[0].size() >= 60)            // clamp floor
  PeakMap ov = dense();                    // explicit value overrides auto, any resolution
  ProSEAlgorithm_test::preprocessSpectra_(ov, 0.5, false, false, 50, 20);
  TEST_TRUE(ov[0].size() <= 50)
}
END_SECTION

START_SECTION(([EXTRA] fragment:deisotope parameter + validation (OpenMS#9619)))
{
  ProSEAlgorithm_test algo;
  // Default is the instrument-aware "auto".
  TEST_EQUAL(algo.getParameters().getValue("fragment:deisotope").toString(), "auto")

  // deisotope=true with a low-resolution (Da > 0.1) tolerance is rejected up front,
  // rather than aborting later inside the search.
  Param p = algo.getParameters();
  p.setValue("fragment:deisotope", "true");
  p.setValue("fragment:mass_tolerance", 0.5);
  p.setValue("fragment:mass_tolerance_unit", "Da");
  TEST_EXCEPTION(Exception::InvalidParameter, algo.setParameters(p))

  // deisotope=true with a high-resolution tolerance is accepted.
  p.setValue("fragment:mass_tolerance", 0.02);
  p.setValue("fragment:mass_tolerance_unit", "Da");
  algo.setParameters(p);
  TEST_EQUAL(algo.getParameters().getValue("fragment:deisotope").toString(), "true")

  // "auto" and "false" accept any tolerance (incl. low-res).
  p.setValue("fragment:deisotope", "auto");
  p.setValue("fragment:mass_tolerance", 0.5);
  p.setValue("fragment:mass_tolerance_unit", "Da");
  algo.setParameters(p);
  p.setValue("fragment:deisotope", "false");
  algo.setParameters(p);
  TEST_EQUAL(algo.getParameters().getValue("fragment:deisotope").toString(), "false")
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
