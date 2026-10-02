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
#include <OpenMS/CHEMISTRY/DecoyGenerator.h>
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
#include <map>
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
  using ProSEAlgorithm::scoringMaxCharge_;
  using ProSEAlgorithm::last_mod_match_tolerance_used_;
  using ProSEAlgorithm::CalibrationResult_;
  using ProSEAlgorithm::preprocessSpectra_;
  using ProSEAlgorithm::filterLocalPeaks_;
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

// etd_ions: c/z+1 instead of b/y fragments
static PeakMap build_calibration_spectra_(const vector<double>& ppm_shifts, bool etd_ions = false)
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
  if (etd_ions)
  {
    tsg_param.setValue("add_b_ions", "false");
    tsg_param.setValue("add_y_ions", "false");
    tsg_param.setValue("add_c_ions", "true");
    tsg_param.setValue("add_zp1_ions", "true");
  }
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
// Appends @p per_protein spectra per protein of decoy peptides: the pseudo-reversed protein
// that ProSE generates as decoy (DecoyGenerator::reversePeptides() with trypsin), digested
// without missed cleavages and with Carbamidomethyl (C), so that a search with decoys
// yields decoy PSMs by construction. The other synthetic spectra are noise-free and
// explained by their targets: without these spectra, no decoy is ever a top hit.
// ---------------------------------------------------------------------------
void addDecoySpectra(PeakMap& spectra, const std::vector<FASTAFile::FASTAEntry>& fasta_db, Size per_protein, double& rt)
{
  ProteaseDigestion digester;
  digester.setEnzyme("Trypsin");
  digester.setMissedCleavages(0);
  const ModifiedPeptideGenerator::MapToResidueType fixed_mods =
    ModifiedPeptideGenerator::getModifications({"Carbamidomethyl (C)"});
  TheoreticalSpectrumGenerator tsg;
  Param tsg_param = tsg.getParameters();
  tsg_param.setValue("add_first_prefix_ion", "true");
  tsg_param.setValue("add_metainfo", "true");
  tsg.setParameters(tsg_param);
  DecoyGenerator decoy_generator;
  for (const auto& entry : fasta_db)
  {
    // Match ProSE's default: preserve the initial Met when generating decoys.
    const bool preserve_met = entry.sequence.size() > 1 && entry.sequence[0] == 'M';
    std::string decoy_sequence
      = decoy_generator.reversePeptides(AASequence::fromString(preserve_met ? entry.sequence.substr(1) : entry.sequence), "Trypsin").toString();
    if (preserve_met) { decoy_sequence.insert(decoy_sequence.begin(), 'M'); }
    const AASequence decoy_protein = AASequence::fromString(decoy_sequence);
    std::vector<AASequence> peptides;
    digester.digest(decoy_protein, peptides, 8, 40);
    peptides.resize(std::min(peptides.size(), per_protein));
    for (AASequence& pep : peptides)
    {
      ModifiedPeptideGenerator::applyFixedModifications(fixed_mods, pep);
      MSSpectrum spec;
      tsg.getSpectrum(spec, pep, 1, 1);
      spec.sortByPosition();
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
}

// ---------------------------------------------------------------------------
// Shared synthetic search problem for the protein-FDR contract tests below.
// 10 proteins + many modified-precursor spectra under a wide precursor window, plus
// spectra of decoy peptides, so ProSEAlgorithm with decoys=true reliably produces BOTH
// target and decoy protein hits — the prerequisite for exercising picked-protein FDR.
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
  addDecoySpectra(spectra, fasta_db, 2, rt);
}

// A run with one ETD spectrum (c/z+1 ions) and one HCD spectrum, each tagged with its activation
// method as in mzML, for the ions:by_activation tests below. The HCD spectrum holds c/z+1 peaks
// besides its b/y peaks: scored with b/y ions only, as HCD spectra are, they stay unannotated.
static PeakMap build_etd_hcd_spectra_()
{
  PeakMap spectra;
  auto add_spectrum = [&spectra](const std::string& seq_str, bool etd)
  {
    TheoreticalSpectrumGenerator tsg;
    Param tsg_param = tsg.getParameters();
    tsg_param.setValue("add_b_ions", etd ? "false" : "true");
    tsg_param.setValue("add_y_ions", etd ? "false" : "true");
    tsg_param.setValue("add_c_ions", "true");
    tsg_param.setValue("add_zp1_ions", "true");
    tsg.setParameters(tsg_param);
    const AASequence seq = AASequence::fromString(seq_str);
    MSSpectrum spec;
    tsg.getSpectrum(spec, seq, 1, 1);
    spec.sortByPosition();
    spec.setMSLevel(2);
    spec.setRT(100.0 + spectra.size());
    Precursor prec;
    prec.setMZ(seq.getMZ(2));
    prec.setCharge(2);
    prec.setActivationMethods({etd ? Precursor::ActivationMethod::ETD : Precursor::ActivationMethod::HCD});
    spec.setPrecursors({prec});
    spec.setNativeID("spectrum=" + StringUtils::toStr(spectra.size()));
    spectra.addSpectrum(std::move(spec));
  };
  add_spectrum("VLGFHQR", true);
  add_spectrum("THQPSANLDIK", false);
  return spectra;
}

// HCD spectrum with the y ions of THQPSANLDIK and, by chance, the c and z+1 ions of NDSIQLHTAPK, a
// peptide of the same composition and hence the same precursor mass. Matched against c and z+1 ions
// as well, NDSIQLHTAPK explains more peaks than THQPSANLDIK; matched against b/y ions, it explains
// hardly any.
static MSSpectrum build_displacement_hcd_spectrum_()
{
  auto ions = [](const std::string& seq_str, bool y, bool c_zp1)
  {
    TheoreticalSpectrumGenerator tsg;
    Param tsg_param = tsg.getParameters();
    tsg_param.setValue("add_b_ions", "false");
    tsg_param.setValue("add_y_ions", y ? "true" : "false");
    tsg_param.setValue("add_c_ions", c_zp1 ? "true" : "false");
    tsg_param.setValue("add_zp1_ions", c_zp1 ? "true" : "false");
    tsg.setParameters(tsg_param);
    MSSpectrum ion_spectrum;
    tsg.getSpectrum(ion_spectrum, AASequence::fromString(seq_str), 1, 1);
    return ion_spectrum;
  };
  MSSpectrum spec = ions("THQPSANLDIK", true, false);
  for (const Peak1D& peak : ions("NDSIQLHTAPK", false, true)) spec.push_back(peak);
  spec.sortByPosition();
  spec.setMSLevel(2);
  spec.setRT(200.0);
  Precursor prec;
  prec.setMZ(AASequence::fromString("THQPSANLDIK").getMZ(2));
  prec.setCharge(2);
  prec.setActivationMethods({Precursor::ActivationMethod::HCD});
  spec.setPrecursors({prec});
  spec.setNativeID("spectrum=hcd");
  return spec;
}

static void configure_by_activation_params_(ProSEAlgorithm& algo, bool by_activation)
{
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
  p.setValue("ions:by_activation", by_activation ? "true" : "false");
  algo.setParameters(p);
}

// top hit per spectrum, keyed by native ID
static std::map<std::string, PeptideHit> top_hits_by_spectrum_(PeptideIdentificationList& pep_ids)
{
  std::map<std::string, PeptideHit> top_hits;
  for (PeptideIdentification& pid : pep_ids)
  {
    if (pid.getHits().empty()) continue;
    pid.sort();
    top_hits[pid.getSpectrumReference()] = pid.getHits()[0];
  }
  return top_hits;
}

static Size count_annotations_(const PeptideHit& hit, const std::string& prefix)
{
  Size n = 0;
  for (const auto& pa : hit.getPeakAnnotations())
  {
    if (StringUtils::hasPrefix(pa.annotation, prefix)) ++n;
  }
  return n;
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

  // --- stop codons: trailing ones are removed, an entry of stop codons only is dropped ---
  // The emptied entry used to reach DecoyGenerator::reversePeptides(), which crashes on a
  // protein without residues.
  {
    ProSEAlgorithm_test algo;
    Param p = algo.getParameters();
    p.setValue("decoys", "generate");
    algo.setParameters(p);
    const std::vector<FASTAFile::FASTAEntry> stop_db = {
      FASTAFile::FASTAEntry("sp|P1|A", "", "PEPTIDEKAAR*"),
      FASTAFile::FASTAEntry("sp|P2|B", "", "**") };
    ProSEAlgorithm_test::DecoyStrategy_ s = algo.resolveDecoyStrategy_(stop_db);
    std::vector<FASTAFile::FASTAEntry> built = algo.buildDecoyAugmentedDB_(stop_db, s);
    TEST_EQUAL(built.size(), 2)             // P1 and its decoy
    for (const auto& e : built)
    {
      TEST_EQUAL(e.sequence.size(), 11)
      TEST_EQUAL(e.identifier.find("P2"), std::string::npos)
    }
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
  addDecoySpectra(spectra, fasta_db, 2, rt); // decoy PSMs for the FDR filter to remove
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

START_SECTION(([EXTRA] Closed search with c/z+1 ions - ETD fragmentation))
{
  // ETD, EThcD and ETciD spectra are dominated by c and z+1 (z-dot) ions. z+1 ions are one
  // hydrogen atom heavier than the z ions of ions:add_z_ions (y - NH3), so only
  // ions:add_zp1_ions matches the C-terminal fragments of such spectra.
  vector<FASTAFile::FASTAEntry> fasta_db = {
    {"P01", "Test", "MSDEREKVLGFHQRMPNASTICYWDLKEGFVRTHQPSANLDIKCMYKWTE"
                    "RHASGDFLKPIVEQNCTMYRGWSADELKHPFNQGTICMSYREWDAVLKPH"},
  };

  TheoreticalSpectrumGenerator tsg;
  Param tsg_param = tsg.getParameters();
  tsg_param.setValue("add_b_ions", "false");
  tsg_param.setValue("add_y_ions", "false");
  tsg_param.setValue("add_c_ions", "true");
  tsg_param.setValue("add_zp1_ions", "true");
  tsg.setParameters(tsg_param);

  const vector<string> test_seqs = { "VLGFHQR", "THQPSANLDIK" };
  PeakMap spectra;
  for (const auto& seq_str : test_seqs)
  {
    const AASequence seq = AASequence::fromString(seq_str);
    MSSpectrum spec;
    tsg.getSpectrum(spec, seq, 1, 1);
    spec.sortByPosition();
    spec.setMSLevel(2);
    spec.setRT(100.0 + spectra.size());
    Precursor prec;
    prec.setMZ(seq.getMZ(2));
    prec.setCharge(2);
    spec.setPrecursors({prec});
    spec.setNativeID("spectrum=" + StringUtils::toStr(spectra.size()));
    spectra.addSpectrum(std::move(spec));
  }

  // top hit per spectrum, keyed by sequence
  auto run_search = [&](const std::string& z_ion_series) {
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
    p.setValue("ions:add_b_ions", "false");
    p.setValue("ions:add_y_ions", "false");
    p.setValue("ions:add_c_ions", "true");
    p.setValue(z_ion_series, "true");
    algo.setParameters(p);
    vector<ProteinIdentification> prot_ids;
    PeptideIdentificationList pep_ids;
    auto ec = algo.search(spectra, fasta_db, prot_ids, pep_ids);
    TEST_EQUAL(ec == ProSEAlgorithm::ExitCodes::EXECUTION_OK, true)
    std::map<std::string, PeptideHit> top_hits;
    for (PeptideIdentification& pid : pep_ids)
    {
      if (pid.getHits().empty()) continue;
      pid.sort();
      top_hits[pid.getHits()[0].getSequence().toUnmodifiedString()] = pid.getHits()[0];
    }
    return top_hits;
  };

  // (1) with z+1 ions, both peptides are identified and all their z+1 ions are matched.
  // z+1 ordinals ("z.3+") count for the longest ion series: z1..z(n-1) is one run longer
  // than c2..c(n-1) (the spectra have no c1).
  std::map<std::string, PeptideHit> zp1_hits = run_search("ions:add_zp1_ions");
  TEST_EQUAL(zp1_hits.size(), test_seqs.size())
  for (const auto& seq_str : test_seqs)
  {
    ABORT_IF(zp1_hits.count(seq_str) != 1)
    const PeptideHit& hit = zp1_hits[seq_str];
    const int n_suffix = static_cast<int>(seq_str.size()) - 1;
    TEST_EQUAL(static_cast<int>(hit.getMetaValue(Constants::UserParam::MATCHED_SUFFIX_IONS)), n_suffix)
    TEST_EQUAL(static_cast<int>(hit.getMetaValue(Constants::UserParam::LONGEST_PEPTIDE_ION_SEQUENCE)), n_suffix)
    Size zp1_annotations = 0;
    for (const auto& pa : hit.getPeakAnnotations())
    {
      if (StringUtils::hasPrefix(pa.annotation, "z.")) ++zp1_annotations;
    }
    TEST_EQUAL(zp1_annotations, static_cast<Size>(n_suffix))
  }

  // (2) ProSE's z ions (y - NH3) miss every z+1 peak; the c ions alone still identify a peptide
  std::map<std::string, PeptideHit> z_hits = run_search("ions:add_z_ions");
  TEST_FALSE(z_hits.empty())
  for (const auto& [seq_str, hit] : z_hits)
  {
    TEST_EQUAL(static_cast<int>(hit.getMetaValue(Constants::UserParam::MATCHED_SUFFIX_IONS)), 0)
  }
}
END_SECTION

START_SECTION(([EXTRA] ions:by_activation scores electron-activated spectra with c/z+1 ions))
{
  // With the default ion series (b/y), ions:by_activation (default on) adds c and z+1 ions for
  // the ETD spectrum only; the HCD spectrum is scored with b/y ions alone, so its c/z+1 peaks
  // stay unannotated.
  vector<FASTAFile::FASTAEntry> fasta_db = {
    {"P01", "Test", "MSDEREKVLGFHQRMPNASTICYWDLKEGFVRTHQPSANLDIKCMYKWTE"
                    "RHASGDFLKPIVEQNCTMYRGWSADELKHPFNQGTICMSYREWDAVLKPH"},
  };
  PeakMap spectra = build_etd_hcd_spectra_();

  ProSEAlgorithm algo;
  configure_by_activation_params_(algo, true);
  vector<ProteinIdentification> prot_ids;
  PeptideIdentificationList pep_ids;
  TEST_EQUAL(algo.search(spectra, fasta_db, prot_ids, pep_ids) == ProSEAlgorithm::ExitCodes::EXECUTION_OK, true)
  std::map<std::string, PeptideHit> hits = top_hits_by_spectrum_(pep_ids);
  ABORT_IF(hits.count("spectrum=0") != 1 || hits.count("spectrum=1") != 1)

  const PeptideHit& etd_hit = hits["spectrum=0"];
  TEST_STRING_EQUAL(etd_hit.getSequence().toUnmodifiedString(), "VLGFHQR")
  TEST_EQUAL(static_cast<int>(etd_hit.getMetaValue(Constants::UserParam::MATCHED_SUFFIX_IONS)), 6)
  TEST_EQUAL(count_annotations_(etd_hit, "z."), 6)

  const PeptideHit& hcd_hit = hits["spectrum=1"];
  TEST_STRING_EQUAL(hcd_hit.getSequence().toUnmodifiedString(), "THQPSANLDIK")
  TEST_EQUAL(count_annotations_(hcd_hit, "y") > 0, true)
  TEST_EQUAL(count_annotations_(hcd_hit, "c"), 0)
  TEST_EQUAL(count_annotations_(hcd_hit, "z."), 0)

  // switched off, the ETD spectrum is scored with b/y ions only and none of its peaks match
  PeakMap spectra_off = build_etd_hcd_spectra_();
  ProSEAlgorithm algo_off;
  configure_by_activation_params_(algo_off, false);
  vector<ProteinIdentification> prot_ids_off;
  PeptideIdentificationList pep_ids_off;
  TEST_EQUAL(algo_off.search(spectra_off, fasta_db, prot_ids_off, pep_ids_off) == ProSEAlgorithm::ExitCodes::EXECUTION_OK, true)
  std::map<std::string, PeptideHit> hits_off = top_hits_by_spectrum_(pep_ids_off);
  TEST_EQUAL(hits_off.count("spectrum=0"), 0)
  TEST_EQUAL(hits_off.count("spectrum=1"), 1)
}
END_SECTION

START_SECTION(([EXTRA] ions:by_activation selects the candidates of other spectra with their ion series alone))
{
  // Only one candidate per spectrum is scored. For the HCD spectrum it must be THQPSANLDIK, as when
  // the spectrum is searched alone, also when an ETD spectrum in the same run makes the index hold
  // c and z+1 ions: against those, NDSIQLHTAPK matches more peaks.
  vector<FASTAFile::FASTAEntry> fasta_db = {
    {"P01", "Test", "MSDEREKVLGFHQRMPNASTICYWDLKEGFVRTHQPSANLDIKCMYKWTE"
                    "RHASGDFLKPIVEQNCTMYRGWSADELKHPFNQGTICMSYREWDAVLKPH"},
    {"P02", "Competitor", "MSGRNDSIQLHTAPKWEAGR"},
  };
  ProSEAlgorithm algo;
  configure_by_activation_params_(algo, true);
  Param p = algo.getParameters();
  p.setValue("scoring:max_candidates_per_spectrum", 1);
  algo.setParameters(p);

  PeakMap hcd_alone;
  hcd_alone.addSpectrum(build_displacement_hcd_spectrum_());
  vector<ProteinIdentification> prot_ids;
  PeptideIdentificationList pep_ids;
  TEST_EQUAL(algo.search(hcd_alone, fasta_db, prot_ids, pep_ids) == ProSEAlgorithm::ExitCodes::EXECUTION_OK, true)
  std::map<std::string, PeptideHit> alone = top_hits_by_spectrum_(pep_ids);
  ABORT_IF(alone.count("spectrum=hcd") != 1)
  TEST_STRING_EQUAL(alone["spectrum=hcd"].getSequence().toUnmodifiedString(), "THQPSANLDIK")

  PeakMap mixed = build_etd_hcd_spectra_();
  PeakMap run;
  run.addSpectrum(mixed[0]); // ETD spectrum of VLGFHQR
  run.addSpectrum(build_displacement_hcd_spectrum_());
  vector<ProteinIdentification> prot_ids_mixed;
  PeptideIdentificationList pep_ids_mixed;
  TEST_EQUAL(algo.search(run, fasta_db, prot_ids_mixed, pep_ids_mixed) == ProSEAlgorithm::ExitCodes::EXECUTION_OK, true)
  std::map<std::string, PeptideHit> hits = top_hits_by_spectrum_(pep_ids_mixed);
  ABORT_IF(hits.count("spectrum=0") != 1)
  TEST_STRING_EQUAL(hits["spectrum=0"].getSequence().toUnmodifiedString(), "VLGFHQR")
  TEST_EQUAL(hits.count("spectrum=hcd"), 1)
  ABORT_IF(hits.count("spectrum=hcd") != 1)
  TEST_STRING_EQUAL(hits["spectrum=hcd"].getSequence().toUnmodifiedString(), "THQPSANLDIK")
  TEST_REAL_SIMILAR(hits["spectrum=hcd"].getScore(), alone["spectrum=hcd"].getScore())
}
END_SECTION

START_SECTION(([EXTRA] ions:by_activation leaves a prepared context unchanged))
{
  // prepareContext(fasta_db) does not know the spectra, so its index holds the configured b/y ions
  // only. search() must not change the context (concurrent searches may share it): for the ETD
  // spectrum it builds a temporary index with c and z+1 ions. prepareContext(fasta_db, true)
  // prepares a context that has them.
  vector<FASTAFile::FASTAEntry> fasta_db = {
    {"P01", "Test", "MSDEREKVLGFHQRMPNASTICYWDLKEGFVRTHQPSANLDIKCMYKWTE"
                    "RHASGDFLKPIVEQNCTMYRGWSADELKHPFNQGTICMSYREWDAVLKPH"},
  };
  ProSEAlgorithm algo;
  configure_by_activation_params_(algo, true);
  ProSEAlgorithm::SearchContext ctx = algo.prepareContext(fasta_db);
  TEST_EQUAL(ctx.electron_ions, false)
  const Size fragments_before = ctx.fragment_index.getNumFragments();

  PeakMap spectra = build_etd_hcd_spectra_();
  vector<ProteinIdentification> prot_ids;
  PeptideIdentificationList pep_ids;
  TEST_EQUAL(algo.search(spectra, ctx, prot_ids, pep_ids) == ProSEAlgorithm::ExitCodes::EXECUTION_OK, true)
  TEST_EQUAL(ctx.electron_ions, false)
  TEST_EQUAL(ctx.fragment_index.getNumFragments(), fragments_before)
  std::map<std::string, PeptideHit> hits = top_hits_by_spectrum_(pep_ids);
  ABORT_IF(hits.count("spectrum=0") != 1)
  TEST_STRING_EQUAL(hits["spectrum=0"].getSequence().toUnmodifiedString(), "VLGFHQR")
  TEST_EQUAL(count_annotations_(hits["spectrum=0"], "z."), 6)

  ProSEAlgorithm::SearchContext electron_ctx = algo.prepareContext(fasta_db, true);
  TEST_EQUAL(electron_ctx.electron_ions, true)
  TEST_EQUAL(electron_ctx.fragment_index.getNumFragments() > fragments_before, true)
}
END_SECTION

START_SECTION(([EXTRA] ions:by_activation gives each file of a multi-file search the results it gets alone))
{
  // An HCD file and an ETD file searched together, in both orders, without and with chunks
  // (database:chunk_size). The ETD file makes the index hold c and z+1 ions. The HCD spectrum is the
  // one above: matched against these ions, NDSIQLHTAPK (in the same chunk) would displace
  // THQPSANLDIK. Each file must get the top hit it gets alone.
  vector<FASTAFile::FASTAEntry> fasta_db = {
    {"P01", "Test", "MSDEREKTHQPSANLDIKCMYKWTERNDSIQLHTAPKWEAGR"},
    {"P02", "Test", "MSDEREKVLGFHQRMPNASTICYWDLKEGFVR"},
  };
  PeakMap etd_spectra, hcd_spectra;
  etd_spectra.addSpectrum(build_etd_hcd_spectra_()[0]);
  hcd_spectra.addSpectrum(build_displacement_hcd_spectrum_());
  std::string etd_file, hcd_file;
  NEW_TMP_FILE(etd_file)
  etd_file += ".mzML";
  NEW_TMP_FILE(hcd_file)
  hcd_file += ".mzML";
  FileHandler().storeExperiment(etd_file, etd_spectra, {FileTypes::MZML});
  FileHandler().storeExperiment(hcd_file, hcd_spectra, {FileTypes::MZML});

  ProSEAlgorithm algo;
  configure_by_activation_params_(algo, true);
  Param p = algo.getParameters();
  p.setValue("scoring:max_candidates_per_spectrum", 1);
  algo.setParameters(p);
  const ProSEAlgorithm::SearchContext electron_ctx = algo.prepareContext(fasta_db, true);
  const Size electron_fragments = electron_ctx.fragment_index.getNumFragments();
  Size chunked_peptides = 0, chunked_fragments = 0;
  for (const auto& protein : fasta_db)
  {
    const auto chunk = algo.prepareContext(vector<FASTAFile::FASTAEntry> {protein}, true);
    chunked_peptides += chunk.fragment_index.getPeptides().size();
    chunked_fragments += chunk.fragment_index.getNumFragments();
  }

  // the top hit of each file searched alone
  auto search_alone = [&algo, &fasta_db](PeakMap alone)
  {
    vector<ProteinIdentification> prot_ids;
    PeptideIdentificationList pep_ids;
    algo.search(alone, fasta_db, prot_ids, pep_ids);
    return top_hits_by_spectrum_(pep_ids);
  };
  const std::map<std::string, PeptideHit> hcd_alone = search_alone(hcd_spectra);
  const std::map<std::string, PeptideHit> etd_alone = search_alone(etd_spectra);
  ABORT_IF(hcd_alone.size() != 1 || etd_alone.size() != 1)
  TEST_STRING_EQUAL(hcd_alone.begin()->second.getSequence().toUnmodifiedString(), "THQPSANLDIK")
  TEST_STRING_EQUAL(etd_alone.begin()->second.getSequence().toUnmodifiedString(), "VLGFHQR")
  TEST_EQUAL(count_annotations_(etd_alone.begin()->second, "z."), 6)

  auto test_same_top_hit = [](PeptideIdentificationList& pep_ids, const std::map<std::string, PeptideHit>& alone)
  {
    std::map<std::string, PeptideHit> hits = top_hits_by_spectrum_(pep_ids);
    TEST_EQUAL(hits.size(), 1)
    if (hits.size() != 1) return;
    TEST_EQUAL(hits.begin()->second.getSequence(), alone.begin()->second.getSequence())
    TEST_REAL_SIMILAR(hits.begin()->second.getScore(), alone.begin()->second.getScore())
  };

  for (int chunk_size : {0, 1})
  {
    p.setValue("database:chunk_size", chunk_size);
    algo.setParameters(p);
    auto hcd_first = algo.searchWithModificationAnalysis(vector<std::string>{hcd_file, etd_file}, fasta_db, vector<std::string>{}, "", false);
    auto etd_first = algo.searchWithModificationAnalysis(vector<std::string>{etd_file, hcd_file}, fasta_db, vector<std::string>{}, "", false);
    ABORT_IF(hcd_first.per_file.size() != 2 || etd_first.per_file.size() != 2)
    for (const auto* res : {&hcd_first, &etd_first})
    {
      TEST_EQUAL(res->shared.chunked, chunk_size > 0)
      // Peptides shared between proteins are indexed once per index, hence once per chunk.
      TEST_EQUAL(res->shared.indexed_peptides, chunk_size > 0 ? chunked_peptides : electron_ctx.fragment_index.getPeptides().size())
      TEST_EQUAL(res->shared.indexed_fragments, chunk_size > 0 ? chunked_fragments : electron_fragments)
    }
    test_same_top_hit(hcd_first.per_file[0].peptide_ids, hcd_alone);
    test_same_top_hit(etd_first.per_file[1].peptide_ids, hcd_alone);
    test_same_top_hit(hcd_first.per_file[1].peptide_ids, etd_alone);
    test_same_top_hit(etd_first.per_file[0].peptide_ids, etd_alone);
  }
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

START_SECTION(([EXTRA] Stop codons in the database: a trailing one is removed, an inner one does not abort the search))
{
  // Sequences translated from genomes (e.g. SGD's yeast database) end with a stop codon ('*'), and
  // a few contain one. P02's VLGFHQ*R has the precursor mass and fragments of VLGFHQR: it used to
  // be indexed, and scoring it aborted the search, because AASequence parses '*' as a weightless
  // X. DIVSAGSLYL, the C-terminal peptide of P03, is only searchable without the stop codon.
  vector<FASTAFile::FASTAEntry> fasta_db = {
    {"P01", "Test", "MSDEREKVLGFHQRMPNASTICYWDLKEGFVRTHQPSANLDIK*"},
    {"P02", "Test", "MSTEKVLGFHQ*RGWSADEK*"},
    {"P03", "Test", "MDSTEKLIHRDIVSAGSLYL*"},
  };

  TheoreticalSpectrumGenerator tsg;
  Param tsg_param = tsg.getParameters();
  tsg_param.setValue("add_first_prefix_ion", "true");
  tsg.setParameters(tsg_param);

  PeakMap spectra;
  for (const std::string seq_str : {"VLGFHQR", "DIVSAGSLYL"})
  {
    const AASequence seq = AASequence::fromString(seq_str);
    MSSpectrum spec;
    tsg.getSpectrum(spec, seq, 1, 1);
    spec.sortByPosition();
    spec.setMSLevel(2);
    spec.setRT(100.0 + spectra.size());
    Precursor prec;
    prec.setMZ(seq.getMZ(2));
    prec.setCharge(2);
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

  std::map<std::string, const PeptideHit*> top_hits;
  for (PeptideIdentification& pid : pep_ids)
  {
    pid.sort();
    for (const PeptideHit& hit : pid.getHits())
    {
      TEST_EQUAL(hit.getSequence().toString().find('X'), std::string::npos)
    }
    if (!pid.getHits().empty()) top_hits[pid.getHits()[0].getSequence().toString()] = &pid.getHits()[0];
  }
  TEST_EQUAL(top_hits.size(), 2)
  TEST_EQUAL(top_hits.count("VLGFHQR"), 1)
  ABORT_IF(top_hits.count("DIVSAGSLYL") != 1)
  const std::vector<PeptideEvidence>& evidences = top_hits["DIVSAGSLYL"]->getPeptideEvidences();
  ABORT_IF(evidences.size() != 1)
  TEST_STRING_EQUAL(evidences[0].getProteinAccession(), "P03")
  TEST_EQUAL(evidences[0].getAAAfter(), PeptideEvidence::C_TERMINAL_AA)
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

START_SECTION((MultiFileSearchResult searchWithModificationAnalysis(const std::vector<std::string>&, const std::vector<FASTAFile::FASTAEntry>&, const std::vector<std::string>&, const std::string&, bool) const))
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

START_SECTION((MultiFileSearchResult searchWithModificationAnalysis(const std::vector<std::string>&, const std::string&, const std::vector<std::string>&, const std::string&, bool) const))
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

START_SECTION(([EXTRA] calibration scores candidates with the configured ion series))
{
  // Same fixture as above, but with ETD-type spectra (c/z+1 ions) and a search for c/z+1
  // ions only. The calibration pass must score its candidates with these ions; a
  // generator left at the b/y defaults matches none of their peaks and calibration fails.
  const vector<double> ppm_shifts = {
    0.0, 2.0, 4.0, 5.0, 6.0, 7.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0
  };
  PeakMap spectra = build_calibration_spectra_(ppm_shifts, /*etd_ions*/ true);
  auto fasta_db = calibration_fasta_db_();

  ProSEAlgorithm_test algo;
  configure_calibration_params_(algo, /*lower_ppm*/ 20.0, /*upper_ppm*/ 30.0,
                                /*min_psms*/ 3);
  Param p = algo.getParameters();
  p.setValue("ions:add_b_ions", "false");
  p.setValue("ions:add_y_ions", "false");
  p.setValue("ions:add_c_ions", "true");
  p.setValue("ions:add_zp1_ions", "true");
  algo.setParameters(p);

  vector<ProteinIdentification> prot_ids;
  PeptideIdentificationList pep_ids;
  auto ec = algo.search(spectra, fasta_db, prot_ids, pep_ids);
  TEST_EQUAL(ec == ProSEAlgorithm::ExitCodes::EXECUTION_OK, true)

  const auto& cal = algo.last_calibration_result_;
  TEST_EQUAL(cal.success, true)
  TEST_EQUAL(cal.extreme_bias, false)
  TEST_EQUAL(cal.precursor_shift > 0.0, true)
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

START_SECTION(([EXTRA] deisotoping keeps a fragment ion that has a small peak one isotope spacing below it))
{
  // b2 of TMTpro-YMATQLLAK-TMTpro (eclipse_tmtpro_10855, scan 30241): the TMTpro label's isotope
  // impurity puts a 5% peak 1.00335 Da below the ion. That peak must not become the monoisotopic
  // peak of the envelope; previously the ion and its +1 isotope were removed as its isotopes.
  // Regular envelopes are still deisotoped: a 1+ envelope keeps its monoisotopic peak only, and
  // a 2+ envelope is converted to its singly charged m/z.
  PeakMap exp;
  MSSpectrum s;
  s.setMSLevel(2);
  s.setRT(1.0);
  Precursor prec;
  prec.setMZ(600.0);
  prec.setCharge(3);
  s.getPrecursors().push_back(prec);
  const std::vector<std::pair<double, float>> peaks = {
    {450.2500, 0.50f}, {450.7517, 0.20f}, {451.2534, 0.05f},     // 2+ envelope
    {598.3157, 0.04f}, {599.3193, 0.74f}, {600.3250, 0.10f},     // shadow peak, ion, +1 isotope
    {700.4000, 1.00f}, {701.4034, 0.35f}, {702.4067, 0.08f}};    // 1+ envelope
  for (const auto& [mz, intensity] : peaks) s.emplace_back(mz, intensity);
  exp.addSpectrum(s);
  ProSEAlgorithm_test::preprocessSpectra_(exp, 20.0, true, true, 0, 20);

  auto has_peak = [&exp](double mz)
  {
    return std::any_of(exp[0].begin(), exp[0].end(), [mz](const Peak1D& p) { return std::fabs(p.getMZ() - mz) <= 20e-6 * mz; });
  };
  TEST_EQUAL(has_peak(599.3193), true)   // the fragment ion survives
  TEST_EQUAL(has_peak(700.4000), true)   // 1+ envelope: monoisotopic peak kept ...
  TEST_EQUAL(has_peak(701.4034), false)  // ... isotopes removed
  TEST_EQUAL(has_peak(702.4067), false)
  TEST_EQUAL(has_peak(450.2500), false)  // 2+ envelope converted ...
  TEST_EQUAL(has_peak(450.2500 * 2.0 - Constants::PROTON_MASS_U), true)  // ... to its 1+ m/z
}
END_SECTION

START_SECTION(([EXTRA] peptidoform deduplication preserves protein evidence and candidate statistics across chunks))
{
  // Eight target proteins and one decoy hold the same peptide. Without deduplication its copies fill
  // the candidate cap; with it, the cap holds three distinct sequences in every search path.
  const AASequence peptide = AASequence::fromString("THQPSANLDIK");
  vector<FASTAFile::FASTAEntry> db;
  for (int i = 0; i < 8; ++i)
  {
    db.push_back({"P" + std::to_string(i), "", peptide.toString()});
  }
  db.push_back({"DECOY_shared", "", peptide.toString()});
  db.push_back({"P_variant", "", "THQPSALNDIK"});
  db.push_back({"DECOY_variant", "", "THQPSADNLIK"});

  MSSpectrum spectrum;
  TheoreticalSpectrumGenerator generator;
  Param gp = generator.getParameters();
  gp.setValue("add_first_prefix_ion", "true");
  generator.setParameters(gp);
  generator.getSpectrum(spectrum, peptide, 1, 1);
  spectrum.setMSLevel(2);
  spectrum.setNativeID("scan=1");
  Precursor precursor;
  precursor.setMZ(peptide.getMZ(2));
  precursor.setCharge(2);
  spectrum.setPrecursors({precursor});
  PeakMap input;
  input.addSpectrum(spectrum);
  std::string input_file;
  NEW_TMP_FILE(input_file)
  FileHandler().storeExperiment(input_file, input, {FileTypes::MZML});

  ProSEAlgorithm algo;
  Param p = algo.getParameters();
  p.setValue("fragment:mass_tolerance", 20.0);
  p.setValue("fragment:mass_tolerance_unit", "ppm");
  p.setValue("fragment:deisotope", "false");
  p.setValue("fragment:min_mz", 0);
  p.setValue("fragment:min_ion_index", 0);
  p.setValue("fragment:min_matched_ions", 3);
  p.setValue("peptide:missed_cleavages", 0);
  p.setValue("scoring:max_candidates_per_spectrum", 3);
  p.setValue("report:top_hits", 10);
  p.setValue("calibration:enabled", "false");
  p.setValue("FDR:PSM", 0.0);
  p.setValue("FDR:protein", 0.0);
  p.setValue("modifications:fixed", vector<string> {});
  p.setValue("modifications:variable", vector<string> {});
  p.setValue("annotate:PSM", vector<string> {"ALL"});
  p.setValue("peptide:deduplicate", "false");
  algo.setParameters(p);
  PeakMap spectra = input;
  vector<ProteinIdentification> proteins;
  PeptideIdentificationList legacy;
  algo.search(spectra, db, proteins, legacy);
  ABORT_IF(legacy.size() != 1)
  TEST_EQUAL(legacy[0].getHits().size(), 3)
  for (const auto& hit : legacy[0].getHits())
  {
    TEST_EQUAL(hit.getSequence(), peptide) // Protein copies exhaust the candidate cap.
  }

  p.setValue("peptide:deduplicate", "true");
  vector<PeptideHit> reference;
  for (Int chunk_size : {0, 1, 4})
  {
    p.setValue("database:chunk_size", chunk_size);
    algo.setParameters(p);
    spectra = input;
    PeptideIdentificationList ids;
    algo.search(spectra, db, proteins, ids);
    ABORT_IF(ids.size() != 1)
    const auto& hits = ids[0].getHits();
    TEST_EQUAL(hits.size(), 3)
    TEST_EQUAL(hits[0].getSequence(), peptide)
    TEST_EQUAL(hits[0].extractProteinAccessionsSet().size(), 9)
    TEST_EQUAL(hits[0].getMetaValue("target_decoy").toString(), "target+decoy")
    TEST_REAL_SIMILAR(static_cast<double>(hits[0].getMetaValue(Constants::UserParam::LN_NUM_CANDIDATES)), std::log1p(3.0))
    TEST_TRUE(static_cast<double>(hits[0].getMetaValue(Constants::UserParam::DELTA_SCORE)) > 0.0)
    set<string> sequences;
    for (const auto& hit : hits)
    {
      sequences.insert(hit.getSequence().toString());
      if (hit.getSequence().toString() == "THQPSADNLIK") { TEST_EQUAL(hit.getMetaValue("target_decoy").toString(), "decoy") }
    }
    TEST_EQUAL(sequences.size(), 3)
    if (chunk_size == 0) reference = hits;
    else
      TEST_TRUE(hits == reference)

    const auto multi = algo.searchWithModificationAnalysis(vector<string> {input_file, input_file}, db, vector<string> {}, "", false);
    TEST_EQUAL(multi.per_file.size(), 2)
    for (const auto& file : multi.per_file)
    {
      ABORT_IF(file.peptide_ids.size() != 1)
      TEST_TRUE(file.peptide_ids[0].getHits() == reference)
    }
  }
}
END_SECTION

START_SECTION(([EXTRA] peptidoform deduplication keeps separate charge and isotope hypotheses))
{
  const AASequence peptide = AASequence::fromString("THQPSANLDIK");
  const vector<FASTAFile::FASTAEntry> db = {{"P1", "", peptide.toString()}, {"P2", "", peptide.toString()}};
  MSSpectrum spectrum;
  TheoreticalSpectrumGenerator().getSpectrum(spectrum, peptide, 1, 1);
  spectrum.setMSLevel(2);
  spectrum.setNativeID("scan=1");
  Precursor precursor;
  spectrum.setPrecursors({precursor});
  ProSEAlgorithm algo;
  Param p = algo.getParameters();
  p.setValue("peptide:deduplicate", "true");
  p.setValue("fragment:mass_tolerance", 20.0);
  p.setValue("fragment:mass_tolerance_unit", "ppm");
  p.setValue("fragment:deisotope", "false");
  p.setValue("fragment:min_matched_ions", 3);
  p.setValue("precursor:mass_tolerance_unit", "Da");
  p.setValue("precursor:min_charge", 2);
  p.setValue("precursor:max_charge", 3);
  p.setValue("precursor:isotope_error_min", 0);
  p.setValue("precursor:isotope_error_max", 1);
  p.setValue("report:top_hits", 10);
  p.setValue("calibration:enabled", "false");
  p.setValue("FDR:PSM", 0.0);
  p.setValue("decoys", "ignore");
  p.setValue("modifications:fixed", vector<string> {});
  p.setValue("modifications:variable", vector<string> {});
  p.setValue("annotate:PSM", vector<string> {"ALL"});
  for (bool unknown_charge : {false, true})
  {
    // Closed windows overlap two isotope assignments; the open window admits
    // two precursor charges but (by design) no separate isotope hypotheses.
    const double tolerance = unknown_charge ? 1000.0 : 0.75;
    p.setValue("precursor:mass_tolerance_lower", tolerance);
    p.setValue("precursor:mass_tolerance_upper", tolerance);
    spectrum.getPrecursors()[0].setMZ(peptide.getMZ(2) - (unknown_charge ? 0.0 : 0.25));
    spectrum.getPrecursors()[0].setCharge(unknown_charge ? 0 : 2);
    for (Int chunk_size : {0, 1})
    {
      p.setValue("database:chunk_size", chunk_size);
      algo.setParameters(p);
      PeakMap spectra;
      spectra.addSpectrum(spectrum);
      vector<ProteinIdentification> proteins;
      PeptideIdentificationList ids;
      algo.search(spectra, db, proteins, ids);
      ABORT_IF(ids.size() != 1)
      TEST_EQUAL(ids[0].getHits().size(), 2)
      set<int> hypotheses;
      for (const auto& hit : ids[0].getHits())
      {
        TEST_EQUAL(hit.getSequence(), peptide)
        TEST_REAL_SIMILAR(static_cast<double>(hit.getMetaValue(Constants::UserParam::LN_NUM_CANDIDATES)), std::log1p(2.0))
        hypotheses.insert(unknown_charge ? hit.getCharge() : static_cast<int>(hit.getMetaValue(Constants::UserParam::ISOTOPE_ERROR)));
      }
      TEST_TRUE(hypotheses == (unknown_charge ? set<int> {2, 3} : set<int> {0, 1}))
    }
  }
}
END_SECTION

START_SECTION(([EXTRA] scoring:fragment_charges scores multiply charged fragments of spectra that are not deisotoped))
{
  // Only doubly charged fragments of a 3+ precursor, as ion-trap CID spectra often hold them.
  const AASequence peptide = AASequence::fromString("THQPSANLDIK");
  const vector<FASTAFile::FASTAEntry> fasta_db = {{"P01", "Test", peptide.toString()}};
  TheoreticalSpectrumGenerator tsg;
  MSSpectrum spec;
  tsg.getSpectrum(spec, peptide, 2, 2);
  spec.setMSLevel(2);
  spec.setRT(100.0);
  spec.setNativeID("scan=1");
  Precursor prec;
  prec.setMZ(peptide.getMZ(3));
  prec.setCharge(3);
  prec.setActivationMethods({Precursor::ActivationMethod::CID});
  spec.setPrecursors({prec});

  ProSEAlgorithm algo;
  Param p = algo.getParameters();
  TEST_EQUAL(p.getValue("scoring:fragment_charges").toString(), "auto")
  p.setValue("fragment:mass_tolerance", 0.01);
  p.setValue("fragment:mass_tolerance_unit", "Da");
  p.setValue("fragment:deisotope", "false");
  p.setValue("fragment:min_ion_index", 0);
  p.setValue("fragment:min_matched_ions", 3);
  p.setValue("fragment:min_mz", 0);
  p.setValue("decoys", "ignore");
  p.setValue("calibration:enabled", "false");
  p.setValue("modifications:fixed", vector<string>{});
  p.setValue("modifications:variable", vector<string>{});
  p.setValue("annotate:PSM", vector<string>{"ALL"});
  p.setValue("FDR:PSM", 0.0);

  auto search = [&](const Param& params, int charge, vector<ProteinIdentification>& proteins)
  {
    algo.setParameters(params);
    PeakMap spectra;
    MSSpectrum input = spec;
    input.getPrecursors()[0].setCharge(charge);
    input.getPrecursors()[0].setMZ(peptide.getMZ(charge));
    spectra.addSpectrum(input);
    PeptideIdentificationList peptides;
    const auto result = algo.search(spectra, fasta_db, proteins, peptides);
    TEST_TRUE(result == ProSEAlgorithm::ExitCodes::EXECUTION_OK)
    return peptides;
  };
  vector<ProteinIdentification> proteins;
  p.setValue("scoring:fragment_charges", "single");
  TEST_TRUE(search(p, 3, proteins).empty()) // 1+ theory cannot explain the doubly charged peaks.
  ABORT_IF(proteins.size() != 1)
  TEST_EQUAL(proteins[0].getSearchParameters().getMetaValue("scoring:fragment_charges_resolved").toString(), "single")

  p.setValue("scoring:fragment_charges", "multiple");
  const auto multiple = search(p, 3, proteins);
  ABORT_IF(multiple.size() != 1 || multiple[0].getHits().empty())
  const PeptideHit& hit = multiple[0].getHits()[0];
  TEST_EQUAL(hit.getSequence(), peptide)
  TEST_TRUE(hit.getScore() > 0.0)
  TEST_TRUE(static_cast<int>(hit.getMetaValue(Constants::UserParam::NUM_MATCHED_PEAKS)) >= 10)
  const double matched_prefix = hit.getMetaValue(Constants::UserParam::MATCHED_PREFIX_IONS);
  TEST_REAL_SIMILAR(static_cast<double>(hit.getMetaValue(Constants::UserParam::MATCHED_PREFIX_IONS_FRACTION)), matched_prefix / (2.0 * peptide.size()))
  TEST_FALSE(hit.getPeakAnnotations().empty())
  for (const auto& annotation : hit.getPeakAnnotations())
  {
    TEST_EQUAL(annotation.charge, 2)
  }
  TEST_EQUAL(proteins[0].getSearchParameters().getMetaValue("scoring:fragment_charges").toString(), "multiple")
  TEST_EQUAL(proteins[0].getSearchParameters().getMetaValue("scoring:fragment_charges_resolved").toString(), "multiple")
  TEST_TRUE(search(p, 2, proteins).empty()) // A 2+ precursor gets no 2+ fragments.
  p.setValue("fragment:max_charge", 1);
  TEST_TRUE(search(p, 3, proteins).empty()) // The fragment charge cap applies.
  p.setValue("fragment:max_charge", 2);

  // 'auto' scores multiple charges exactly when the spectra are not deisotoped: with
  // fragment:deisotope=false, and at an ion-trap tolerance the deisotoper does not support.
  const auto explicit_multiple = search(p, 3, proteins);
  p.setValue("scoring:fragment_charges", "auto");
  for (int chunk_size : {0, 1})
  {
    p.setValue("database:chunk_size", chunk_size);
    const auto automatic = search(p, 3, proteins);
    ABORT_IF(automatic.size() != 1)
    TEST_TRUE(automatic[0].getHits() == explicit_multiple[0].getHits())
    TEST_EQUAL(proteins[0].getSearchParameters().getMetaValue("scoring:fragment_charges").toString(), "auto")
    TEST_EQUAL(proteins[0].getSearchParameters().getMetaValue("scoring:fragment_charges_resolved").toString(), "multiple")
  }
  p.setValue("database:chunk_size", 0);
  p.setValue("fragment:deisotope", "auto");
  p.setValue("fragment:mass_tolerance", 0.5);
  TEST_EQUAL(search(p, 3, proteins).size(), 1)
  TEST_EQUAL(proteins[0].getSearchParameters().getMetaValue("scoring:fragment_charges_resolved").toString(), "multiple")
  // High-resolution spectra are deisotoped to charge 1, so 'auto' keeps single charges.
  p.setValue("fragment:mass_tolerance", 20.0);
  p.setValue("fragment:mass_tolerance_unit", "ppm");
  search(p, 3, proteins);
  ABORT_IF(proteins.size() != 1)
  TEST_EQUAL(proteins[0].getSearchParameters().getMetaValue("scoring:fragment_charges_resolved").toString(), "single")

  // The precursor-calibration pass scores with the same fragment charges.
  p.setValue("fragment:deisotope", "false");
  p.setValue("fragment:mass_tolerance", 0.01);
  p.setValue("fragment:mass_tolerance_unit", "Da");
  p.setValue("calibration:enabled", "true");
  p.setValue("calibration:subset_ratio", 1.0);
  p.setValue("calibration:min_psms", 1);
  ProSEAlgorithm_test calibrated;
  calibrated.setParameters(p);
  TEST_EQUAL(calibrated.scoringMaxCharge_(3), 2)
  TEST_EQUAL(calibrated.scoringMaxCharge_(2), 1)
  PeakMap spectra;
  spectra.addSpectrum(spec);
  PeptideIdentificationList peptide_ids;
  TEST_TRUE(calibrated.search(spectra, fasta_db, proteins, peptide_ids) == ProSEAlgorithm::ExitCodes::EXECUTION_OK)
  TEST_TRUE(calibrated.last_calibration_result_.success)
  TEST_EQUAL(peptide_ids.size(), 1)
}
END_SECTION

START_SECTION(([EXTRA] high resolution local filtering preserves short final windows and aligned peak data))
{
  auto filter = [](double tolerance, bool ppm, const std::string& mode) {
    PeakMap exp;
    MSSpectrum spectrum;
    spectrum.setMSLevel(2);
    spectrum.setNativeID("scan=17");
    const std::vector<double> mz {100.0, 110.0, 120.0, 200.0, 201.0};
    const std::vector<float> intensity {5.0f, 5.0f, 5.0f, 7.0f, 2.0f};
    spectrum.getFloatDataArrays().emplace_back();
    spectrum.getFloatDataArrays().back().setName("ion_mobility");
    spectrum.getStringDataArrays().emplace_back();
    spectrum.getStringDataArrays().back().setName("annotation");
    spectrum.getIntegerDataArrays().emplace_back();
    spectrum.getIntegerDataArrays().back().setName("original_index");
    for (Size i = 0; i < mz.size(); ++i)
    {
      Peak1D peak;
      peak.setMZ(mz[i]);
      peak.setIntensity(intensity[i]);
      spectrum.push_back(peak);
      spectrum.getIntegerDataArrays().back().push_back(static_cast<Int>(i));
      spectrum.getFloatDataArrays().back().push_back(static_cast<float>(i) / 10.0f);
      spectrum.getStringDataArrays().back().push_back(std::to_string(i));
    }
    exp.addSpectrum(spectrum);
    ProSEAlgorithm_test::preprocessSpectra_(exp, tolerance, ppm, false, 400, 2, mode);
    return exp[0];
  };

  const MSSpectrum full = filter(20.0, true, "auto");
  TEST_EQUAL(full.size(), 4)
  TEST_EQUAL(full.getNativeID(), "scan=17")
  TEST_REAL_SIMILAR(full[0].getMZ(), 100.0)
  TEST_REAL_SIMILAR(full[1].getMZ(), 110.0)
  TEST_REAL_SIMILAR(full[2].getMZ(), 200.0)
  TEST_REAL_SIMILAR(full[3].getMZ(), 201.0)
  TEST_EQUAL(full.getIntegerDataArrays()[0][0], 0)
  TEST_EQUAL(full.getIntegerDataArrays()[0][1], 1)
  TEST_EQUAL(full.getIntegerDataArrays()[0][2], 3)
  TEST_EQUAL(full.getIntegerDataArrays()[0][3], 4)
  TEST_EQUAL(full.getFloatDataArrays()[0].getName(), "ion_mobility")
  TEST_EQUAL(full.getStringDataArrays()[0].getName(), "annotation")
  for (Size i = 0; i < full.size(); ++i)
  {
    const Int original_index = full.getIntegerDataArrays()[0][i];
    TEST_REAL_SIMILAR(full.getFloatDataArrays()[0][i], original_index / 10.0)
    TEST_EQUAL(full.getStringDataArrays()[0][i], std::to_string(original_index))
  }
  TEST_TRUE(full == filter(20.0, true, "jump_full"))
  TEST_TRUE(full == filter(100.0, true, "auto"))
  TEST_TRUE(full == filter(0.1, false, "auto"))
  TEST_TRUE(full == filter(0.5, false, "jump_full"))

  // Legacy jump filtering rounds the short final window's quota down to zero.
  const MSSpectrum legacy = filter(20.0, true, "jump");
  TEST_EQUAL(legacy.size(), 2)
  TEST_TRUE(legacy == filter(0.5, false, "auto"))
  TEST_TRUE(legacy == filter(101.0, true, "auto"))
  TEST_TRUE(legacy == filter(0.1001, false, "auto"))

  PeakMap empty;
  empty.addSpectrum(MSSpectrum());
  ProSEAlgorithm_test::preprocessSpectra_(empty, 20.0, true, false, 400, 20);
  TEST_TRUE(empty[0].empty())

  PeakMap singleton;
  MSSpectrum spectrum;
  Peak1D peak;
  peak.setMZ(1000.0);
  peak.setIntensity(1000.0);
  spectrum.push_back(peak);
  singleton.addSpectrum(spectrum);
  ProSEAlgorithm_test::preprocessSpectra_(singleton, 20.0, true, false, 400, 20);
  TEST_EQUAL(singleton[0].size(), 1)
  TEST_REAL_SIMILAR(singleton[0][0].getMZ(), 1000.0)
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

START_SECTION(([EXTRA] dense spectra keep peaks:dense_window_top peaks per window))
{
  // One 100 Da window (0.5 Da spacing): 20 strong peaks and n weak ones. With 80 weak peaks the top-20 filter removes 80 of
  // 280 intensity units (29%), so the spectrum is dense; with 10 weak peaks it removes 10 of 210 (5%).
  auto window = [](Size weak)
  {
    MSSpectrum s;
    s.setMSLevel(2);
    for (Size i = 0; i < 20 + weak; ++i) { s.emplace_back(200.0 + i * 0.5, i < 20 ? 10.0f : 1.0f); }
    s.sortByPosition();
    return s;
  };

  MSSpectrum dense = window(80);
  TEST_EQUAL(ProSEAlgorithm_test::filterLocalPeaks_(dense, 20, 100, 0.2), true)
  TEST_EQUAL(dense.size(), 100)

  MSSpectrum sparse = window(10);
  TEST_EQUAL(ProSEAlgorithm_test::filterLocalPeaks_(sparse, 20, 100, 0.2), false)
  TEST_EQUAL(sparse.size(), 20)
  for (const Peak1D& p : sparse) { TEST_REAL_SIMILAR(p.getIntensity(), 10.0) }

  // A larger allowed loss, a dense quota of 0 or one not above the regular quota keep the regular quota.
  MSSpectrum tolerant = window(80);
  TEST_EQUAL(ProSEAlgorithm_test::filterLocalPeaks_(tolerant, 20, 100, 0.3), false)
  TEST_EQUAL(tolerant.size(), 20)
  MSSpectrum disabled = window(80);
  TEST_EQUAL(ProSEAlgorithm_test::filterLocalPeaks_(disabled, 20, 0, 0.2), false)
  TEST_EQUAL(disabled.size(), 20)
  MSSpectrum not_larger = window(80);
  TEST_EQUAL(ProSEAlgorithm_test::filterLocalPeaks_(not_larger, 20, 20, 0.2), false)
  TEST_EQUAL(not_larger.size(), 20)

  // The dense quota still caps each window: 150 weak peaks, 100 of 170 peaks kept, the strongest first.
  MSSpectrum capped = window(150);
  TEST_EQUAL(ProSEAlgorithm_test::filterLocalPeaks_(capped, 20, 100, 0.2), true)
  TEST_EQUAL(capped.size(), 100)
  TEST_EQUAL(std::count_if(capped.begin(), capped.end(), [](const Peak1D& p) { return p.getIntensity() == 10.0f; }), 20)

  // preprocessSpectra_ applies it only where the local filter keeps full quotas (high resolution under 'auto',
  // or 'jump_full') and reports the number of dense spectra. The legacy 'jump' filter is unchanged.
  auto run = [&window](double tolerance, bool ppm, const std::string& type, Size dense_top)
  {
    PeakMap exp;
    exp.addSpectrum(window(80));
    exp.addSpectrum(window(10));
    const Size n_dense = ProSEAlgorithm_test::preprocessSpectra_(exp, tolerance, ppm, false, 400, 20, type, dense_top, 0.2);
    return std::make_tuple(n_dense, exp[0].size(), exp[1].size());
  };
  TEST_TRUE(run(20.0, true, "auto", 100) == std::make_tuple(Size(1), Size(100), Size(20)))
  TEST_TRUE(run(0.5, false, "jump_full", 100) == std::make_tuple(Size(1), Size(100), Size(20)))
  TEST_TRUE(run(20.0, true, "auto", 0) == std::make_tuple(Size(0), Size(20), Size(20)))
  TEST_EQUAL(std::get<0>(run(0.5, false, "auto", 100)), 0)
  TEST_TRUE(run(0.5, false, "auto", 100) == run(0.5, false, "auto", 0))
  TEST_EQUAL(std::get<0>(run(20.0, true, "jump", 100)), 0)
  TEST_TRUE(run(20.0, true, "jump", 100) == run(20.0, true, "jump", 0))

  ProSEAlgorithm_test algo;
  TEST_EQUAL(int(algo.getParameters().getValue("peaks:dense_window_top")), 100)
  TEST_REAL_SIMILAR(double(algo.getParameters().getValue("peaks:dense_intensity_loss")), 0.2)
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

START_SECTION(([EXTRA] preprocessSpectra_ removes only peaks without intensity, so results do not depend on the intensity scale))
{
  // Zero-intensity peaks are dropped; every positive intensity survives, however small on an
  // absolute scale. The ThresholdMower default of 0.05 absolute, applied before normalization,
  // deleted real peaks from intensity-scaled spectra.
  {
    PeakMap exp;
    MSSpectrum s;
    s.setMSLevel(2);
    s.setRT(1.0);
    Precursor prec;
    prec.setMZ(500.0);
    prec.setCharge(2);
    s.getPrecursors().push_back(prec);
    const std::vector<std::pair<double, float>> peaks = {
      {110.07, 0.0f}, {120.08, 1e-4f}, {130.10, 0.01f}, {200.10, 0.04f}, {300.20, 0.5f}, {350.25, 1.0f}};
    for (const auto& [mz, intensity] : peaks) s.emplace_back(mz, intensity);
    exp.addSpectrum(s);
    ProSEAlgorithm_test::preprocessSpectra_(exp, 20.0, true, false, 0, 20);
    ABORT_IF(exp.size() != 1)
    TEST_EQUAL(exp[0].size(), 5) // the empty peak is gone, the weak ones stay
    TEST_REAL_SIMILAR(exp[0][0].getMZ(), 120.08)
    TEST_TRUE(exp[0][0].getIntensity() > 0.0f)
  }

  // Spectra scaled to a base peak of 1e-3 (every peak below the old cutoff) give the same
  // identifications and scores as the unscaled spectra.
  const vector<double> no_shift(12, 0.0);
  PeakMap spectra = build_calibration_spectra_(no_shift);
  PeakMap scaled = spectra;
  for (MSSpectrum& spectrum : scaled)
  {
    for (Peak1D& peak : spectrum) peak.setIntensity(peak.getIntensity() * 1e-3f);
  }
  auto fasta_db = calibration_fasta_db_();

  ProSEAlgorithm algo;
  configure_calibration_params_(algo, 20.0, 30.0, 3);
  Param p = algo.getParameters();
  p.setValue("calibration:enabled", "false");
  algo.setParameters(p);

  vector<ProteinIdentification> prot_ids, scaled_prot_ids;
  PeptideIdentificationList pep_ids, scaled_pep_ids;
  TEST_EQUAL(algo.search(spectra, fasta_db, prot_ids, pep_ids) == ProSEAlgorithm::ExitCodes::EXECUTION_OK, true)
  TEST_EQUAL(algo.search(scaled, fasta_db, scaled_prot_ids, scaled_pep_ids) == ProSEAlgorithm::ExitCodes::EXECUTION_OK, true)
  const std::map<std::string, PeptideHit> hits = top_hits_by_spectrum_(pep_ids);
  const std::map<std::string, PeptideHit> scaled_hits = top_hits_by_spectrum_(scaled_pep_ids);
  TEST_TRUE(hits.size() > 0)
  TEST_EQUAL(scaled_hits.size(), hits.size())
  for (const auto& [reference, hit] : hits)
  {
    const auto it = scaled_hits.find(reference);
    TEST_EQUAL(it != scaled_hits.end(), true)
    if (it == scaled_hits.end()) continue;
    TEST_EQUAL(it->second.getSequence(), hit.getSequence())
    TEST_EQUAL(it->second.getCharge(), hit.getCharge())
    TEST_REAL_SIMILAR(it->second.getScore(), hit.getScore())
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
START_SECTION(([EXTRA] ProSE searches initial methionine loss by default in single and chunked multi - file searches))
{
  const vector<FASTAFile::FASTAEntry> database {{"target", "", "MPEPTIDER"}, {"DECOY_control", "", "MPEPTIDEK"}};
  vector<AASequence> sequences {AASequence::fromString("PEPTIDER"), AASequence::fromString("PEPTIDER"), AASequence::fromString("PEPTIDEK"),
                                AASequence::fromString("MPEPTIDER")};
  sequences[1].setNTerminalModification("Acetyl (Protein N-term)");
  PeakMap original;
  for (Size i = 0; i < sequences.size(); ++i)
  {
    MSSpectrum spectrum;
    TheoreticalSpectrumGenerator().getSpectrum(spectrum, sequences[i], 1, 1);
    spectrum.sortByPosition();
    spectrum.setMSLevel(2);
    spectrum.setRT(i + 1);
    spectrum.setNativeID("scan=" + std::to_string(i + 1));
    Precursor precursor;
    precursor.setCharge(2);
    precursor.setMZ(sequences[i].getMZ(2));
    spectrum.setPrecursors({precursor});
    original.addSpectrum(spectrum);
  }
  std::string spectrum_file;
  NEW_TMP_FILE(spectrum_file)
  spectrum_file += ".mzML";
  FileHandler().storeExperiment(spectrum_file, original, {FileTypes::MZML});
  ProSEAlgorithm algo;
  auto p = algo.getParameters();
  TEST_TRUE(p.getValue("peptide:clip_nterm_methionine").toBool())
  p.setValue("calibration:enabled", "false");
  p.setValue("fragment:deisotope", "false");
  p.setValue("fragment:min_ion_index", 0);
  p.setValue("precursor:isotope_error_min", 0);
  p.setValue("precursor:isotope_error_max", 0);
  p.setValue("modifications:fixed", StringList {});
  p.setValue("modifications:variable", StringList {"Acetyl (Protein N-term)"});
  p.setValue("modifications:variable_max_per_peptide", 1);
  p.setValue("peptide:missed_cleavages", 0);
  algo.setParameters(p);

  auto check = [&](const vector<ProteinIdentification>& proteins, const PeptideIdentificationList& peptides, bool enabled) {
    TEST_EQUAL(peptides.size(), enabled ? 4 : 1)
    TEST_EQUAL(proteins.size(), 1)
    if (! proteins.empty()) { TEST_EQUAL(proteins[0].getSearchParameters().getMetaValue("peptide:clip_nterm_methionine").toBool(), enabled) }
    for (const auto& id : peptides)
    {
      const Size i = static_cast<Size>(std::stoi(id.getSpectrumReference().substr(5)) - 1);
      TEST_EQUAL(id.getHits().size(), 1)
      if (id.getHits().empty() || i >= sequences.size()) { continue; }
      const auto& hit = id.getHits()[0];
      TEST_EQUAL(hit.getSequence(), sequences[i])
      TEST_EQUAL(hit.getMetaValue("target_decoy").toString(), i == 2 ? "decoy" : "target")
      TEST_TRUE(hit.getScore() > 0)
      TEST_EQUAL(hit.getPeptideEvidences().size(), 1)
      for (const auto& evidence : hit.getPeptideEvidences())
      {
        TEST_EQUAL(evidence.getStart(), i == 3 ? 0 : 1)
        TEST_EQUAL(evidence.getEnd(), 8)
        if (i != 3) { TEST_EQUAL(evidence.getAABefore(), 'M') }
      }
    }
  };
  vector<ProteinIdentification> proteins;
  PeptideIdentificationList peptides;
  PeakMap spectra = original;
  TEST_TRUE(algo.search(spectra, database, proteins, peptides) == ProSEAlgorithm::ExitCodes::EXECUTION_OK)
  check(proteins, peptides, true);

  for (bool enabled : {true, false})
  {
    p.setValue("peptide:clip_nterm_methionine", enabled ? "true" : "false");
    for (int chunk_size : {0, 1})
    {
      p.setValue("database:chunk_size", chunk_size);
      algo.setParameters(p);
      auto result = algo.searchWithModificationAnalysis(vector<string> {spectrum_file, spectrum_file}, database, {}, "", false);
      TEST_EQUAL(result.per_file.size(), 2)
      TEST_EQUAL(result.shared.chunked, chunk_size > 0)
      for (const auto& file : result.per_file)
      {
        check(file.protein_ids, file.peptide_ids, enabled);
      }
    }
  }
}
END_SECTION

START_SECTION(([EXTRA] generated decoys preserve initial methionine when clipping is enabled))
{
  const vector<FASTAFile::FASTAEntry> database {{"protein", "", "MACDEKAGHILR"}, {"non_m", "", "ACDEKAGHILR"}, {"single_m", "", "M"}};
  for (bool enabled : {false, true})
  {
    for (const string specificity : {"full", "none"})
    {
      ProSEAlgorithm_test algo;
      auto p = algo.getParameters();
      p.setValue("decoys", "generate");
      p.setValue("peptide:clip_nterm_methionine", enabled ? "true" : "false");
      p.setValue("peptide:enzyme_specificity", specificity);
      algo.setParameters(p);
      auto strategy = algo.resolveDecoyStrategy_(database);
      auto result = algo.buildDecoyAugmentedDB_(database, strategy);
      TEST_EQUAL(result.size(), 6)
      DecoyGenerator generator;
      for (const auto& original : database)
      {
        auto expected = original.sequence;
        const bool preserve_met = enabled && expected.size() > 1 && expected[0] == 'M';
        if (preserve_met) { expected.erase(0, 1); }
        const auto seq = AASequence::fromString(expected);
        expected = (specificity == "none" ? generator.reverseProtein(seq) : generator.reversePeptides(seq, "Trypsin")).toString();
        if (preserve_met) { expected.insert(expected.begin(), 'M'); }
        for (const auto& entry : result)
        {
          if (entry.identifier == original.identifier) { TEST_EQUAL(entry.sequence, original.sequence) }
          if (entry.identifier == "DECOY_" + original.identifier) { TEST_EQUAL(entry.sequence, expected) }
        }
      }
    }
  }
}
END_SECTION

END_TEST
