// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Raphael Förster $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////////
#include <OpenMS/ANALYSIS/ID/FragmentIndex.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <chrono>
#include <OpenMS/KERNEL/Peak1D.h>
#include <limits>

/*
  FragmentIndex tests

  This suite verifies:
  - build(): digestion and peptide generation across enzyme, length/mass limits, missed cleavages, and modifications; asserts ordering invariants for
  peptides/fragments.
  - clear(): resets index state.
  - querySpectrum(): candidate generation across precursor charges with and without known precursor charge.
  - isotope_error: precursor m/z isotope offsets map to expected peptide subsequences.
  - tolerance: fragment and precursor tolerance handling using small deterministic m/z jitter.

  Invariants validated by helper methods:
  - Peptides are sorted by precursor_mz_ (non-decreasing).
  - Fragments are bucketed and within each bucket sorted by peptide_idx_.
*/
using namespace OpenMS;
using namespace std;

// Helper test subclass exposing internal invariants (fi_peptides_, fi_fragments_, bucketsize_).
// Only used in tests to assert ordering and to craft white-box expectations.
class FragmentIndex_test : public FragmentIndex
{
public:
  // Verifies that the generated peptide set matches the expected set exactly
  // (by subsequence window and modification index).
  bool testDigestion(const std::vector<FragmentIndex::Peptide>& expected)
  {
    if (expected.size() != fi_peptides_.size()) return false;
    for (const auto& exp : expected)
    {
      bool found = false;
      for (const auto& act : fi_peptides_)
      {
        if ((exp.sequence_ == act.sequence_) && (exp.modification_idx_ == act.modification_idx_))
        {
          found = true;
          break;
        }
      }
      if (! found) return false;
    }
    return true;
  }
  // Checks non-decreasing order of precursor_mz_ across all peptides (invariant of build()).
  bool peptidesSorted()
  {
    float last_mz = std::numeric_limits<float>::lowest();
    for (const auto& pep : fi_peptides_)
    {
      if (pep.precursor_mz_ >= last_mz) { last_mz = pep.precursor_mz_; }
      else { return false; }
    }
    return true;
  }

  // Validates that within each fragment bucket, peptide_idx_ is non-decreasing.
  // This captures the two-dimensional ordering constraint of the index.
  bool fragmentsSorted()
  {
    for (size_t fi_idx = 0; fi_idx < fi_fragment_mzs_.size(); fi_idx += bucketsize_)
    {
      UInt32 last_idx = 0;
      const size_t end = std::min(fi_idx + bucketsize_, fi_fragment_mzs_.size());
      for (size_t bucket_idx = fi_idx; bucket_idx < end; ++bucket_idx)
      {
        if (fi_fragment_peptide_idxs_[bucket_idx] < last_idx) return false;
        last_idx = fi_fragment_peptide_idxs_[bucket_idx];
      }
    }
    return true;
  }

  bool testQuery(const UInt32 charge, const bool precursor_mz_known, const std::vector<FASTAFile::FASTAEntry>& entries)
  {
    // fetch parameters for modification generation
    auto params = getParameters();
    const StringList modifications_fixed_ = ListUtils::toStringList<std::string>(params.getValue("modifications:fixed"));
    const StringList modifications_variable_ = ListUtils::toStringList<std::string>(params.getValue("modifications:variable"));
    const ModifiedPeptideGenerator::MapToResidueType fixed_modifications = ModifiedPeptideGenerator::getModifications(modifications_fixed_);
    const ModifiedPeptideGenerator::MapToResidueType variable_modifications = ModifiedPeptideGenerator::getModifications(modifications_variable_);

    // Create theoretical spectra for different charges
    TheoreticalSpectrumGenerator tsg;
    std::vector<AASequence> mod_peptides;
    PeakSpectrum b_y_ions;
    MSSpectrum spec_theo;
    Precursor prec_theo;

    const std::vector<FragmentIndex::Peptide>& peptides = getPeptides();
    bool test = true;

    // Create different ms/ms spectra with different charges

    size_t peptide_idx = 0; // use size_t to match SpectrumMatch::peptide_idx_ type
    // For each peptide that was created, we now generate a theoretical spectra for the given charge
    // Each peptide should hit its own entry in the db. In this case the test returns true
    for (const auto& pep : peptides)
    {
      FragmentIndex::SpectrumMatchesTopN sms;
      b_y_ions.clear(true);
      mod_peptides.clear();
      spec_theo.clear(true);

      prec_theo.clearMetaInfo();
      const AASequence unmod_peptide = AASequence::fromString(entries[0].sequence.substr(pep.sequence_.first, pep.sequence_.second));
      AASequence mod_peptide = AASequence(unmod_peptide); // copy the peptide
      ModifiedPeptideGenerator::applyFixedModifications(fixed_modifications, mod_peptide);
      ModifiedPeptideGenerator::applyVariableModifications(variable_modifications, mod_peptide,
                                                           params.getValue("modifications:variable_max_per_peptide"), mod_peptides);
      mod_peptide = mod_peptides[pep.modification_idx_];
      tsg.getSpectrum(b_y_ions, mod_peptide, charge, charge);
      prec_theo.setMZ(mod_peptide.getMZ(charge));
      if (precursor_mz_known) { prec_theo.setCharge(charge); }
      spec_theo.setMSLevel(2);
      spec_theo.setPrecursors({prec_theo});
      for (const auto& ion : b_y_ions)
      {
        spec_theo.push_back(ion);
      }

      querySpectrum(spec_theo, sms);
      bool found = false;

      // iterate candidates and check matching count for the exact peptide/charge
      for (const auto& s : sms.hits_)
      {
        if ((s.peptide_idx_ == peptide_idx) && (s.precursor_charge_ == charge))
        {
          // All generated peaks must be matched and the correct precursor charge identified
          found = (s.num_matched_ >= spec_theo.size());
        }
      }
      test = test && found;
      peptide_idx++;
    }
    return test;
  }
};

// Helper: build a FragmentIndex with given params and protein sequence
FragmentIndex_test buildTestIndex(const std::string& protein_seq,
                                  double frag_tol = 0.05,
                                  const std::string& frag_unit = "Da",
                                  double prec_tol = 10.0,
                                  const std::string& prec_unit = "ppm",
                                  int missed_cleavages = 1,
                                  int min_size = 5,
                                  int max_size = 40)
{
  FragmentIndex_test fi;
  Param p = fi.getParameters();
  p.setValue("ions:add_b_ions", "true");
  p.setValue("ions:add_y_ions", "true");
  p.setValue("peptide:min_size", min_size);
  p.setValue("peptide:max_size", max_size);
  p.setValue("peptide:missed_cleavages", missed_cleavages);
  p.setValue("peptide:min_mass", 0);
  p.setValue("peptide:max_mass", 9000);
  p.setValue("fragment:min_mz", 0);
  p.setValue("fragment:max_mz", 90000);
  p.setValue("fragment:mass_tolerance", frag_tol);
  p.setValue("fragment:mass_tolerance_unit", frag_unit);
  p.setValue("precursor:mass_tolerance", prec_tol);
  p.setValue("precursor:mass_tolerance_unit", prec_unit);
  p.setValue("precursor:min_charge", 1);
  p.setValue("precursor:max_charge", 3);
  fi.setParameters(p);

  std::vector<FASTAFile::FASTAEntry> entries;
  entries.push_back({"sp|TEST", "test protein", protein_seq});
  fi.build(entries);
  return fi;
}

// Helper: compare SIMD vs scalar query results for every peak in a spectrum
bool simdScalarMatch(FragmentIndex_test& fi, const MSSpectrum& spec,
                     const std::pair<size_t,size_t>& range, uint16_t charge,
                     size_t& simd_count, size_t& scalar_count)
{
  bool match = true;
  for (const auto& peak : spec)
  {
    auto simd_hits = fi.query(peak, range, charge);
    auto scalar_hits = fi.queryScalar(peak, range, charge);
    simd_count += simd_hits.size();
    scalar_count += scalar_hits.size();
    if (simd_hits.size() != scalar_hits.size()) match = false;
  }
  return match;
}

//////////////////////////////
START_TEST(FragmentIndex, "$Id")

//////////////////////////////

/// Test the build for peptides
START_SECTION(build())
{
  // Test proteins used to generate expected peptides for multiple parameterizations
  /*
    Format of expected peptide descriptors below and their mapping to FragmentIndex::Peptide fields:
      { protein_idx, modification_idx_, { start, length }, precursor_mz_ }

    Where:
    - protein_idx: 0-based index into the FASTA entries vector passed to build(); selects the source protein.
    - modification_idx_: index into mod_peptides returned by ModifiedPeptideGenerator for the given unmodified subsequence
                        (0 = unmodified; higher values enumerate concrete variable-mod combinations).
    - start: 0-based start offset within the selected protein sequence.
    - length: number of residues for the peptide (used as std::string::substr(start, length)).
    - precursor_mz_: mono-isotopic m/z at charge 1 (M+H)+. In these tests we often use a dummy value, as only ordering
                    invariants on peptides/fragment buckets are asserted.

    Note: testDigestion() compares expected vs. built peptides only by {sequence_, modification_idx_}.
  */
  const std::vector<FASTAFile::FASTAEntry> entries0 {{"t", "t", "ARGEPADSSRKDFDMDMDM"}, {"t2", "t2", "HALLORTSCHSM"}};
  // Expected peptides when enabling fixed Carbamidomethyl (C) and variable Oxidation (M)
  std::vector<FragmentIndex::Peptide> peptides_we_should_hit_mod {{0, 0, {2, 8}, 5},  {0, 0, {11, 8}, 5}, {0, 1, {11, 8}, 5}, {0, 2, {11, 8}, 5},
                                                                  {0, 3, {11, 8}, 5}, {0, 4, {11, 8}, 5}, {0, 5, {11, 8}, 5}, {0, 6, {11, 8}, 5},
                                                                  {1, 0, {0, 6}, 5},  {1, 0, {6, 6}, 5},  {1, 1, {6, 6}, 5}

  };
  // Expected peptides without min/max size constraints (no missed cleavages, no modifications)
  std::vector<FragmentIndex::Peptide> peptides_unmod_no_minmax {{0, 0, {0, 2}, 5},  {0, 0, {2, 8}, 5}, {0, 0, {10, 1}, 5},
                                                                {0, 0, {11, 8}, 5}, {1, 0, {0, 6}, 5}, {1, 0, {6, 6}, 5}};

  // Expected peptides with size in [min_size, max_size] only
  std::vector<FragmentIndex::Peptide> peptides_unmod_minmax {{0, 0, {0, 2}, 5}, {1, 0, {0, 6}, 5}, {1, 0, {6, 6}, 5}};
  // Expected peptides with one missed cleavage allowed
  std::vector<FragmentIndex::Peptide> peptides_unmod_minmax_missed_cleavage {{0, 0, {0, 2}, 5},  {0, 0, {2, 8}, 5}, {0, 0, {11, 8}, 5},
                                                                            {0, 0, {0, 10}, 5}, {0, 0, {2, 9}, 5}, {0, 0, {10, 9}, 5},
                                                                            {1, 0, {0, 6}, 5},  {1, 0, {6, 6}, 5}, {1, 0, {0, 12}, 5}};


  FragmentIndex_test buildTest;
  auto params = buildTest.getParameters();
  params.setValue("enzyme", "Trypsin");
  params.setValue("peptide:missed_cleavages", 0);
  params.setValue("peptide:min_mass", 0);
  params.setValue("peptide:min_size", 0);
  params.setValue("peptide:max_mass", 5000);
  params.setValue("modifications:variable", std::vector<std::string> {});
  params.setValue("modifications:fixed", std::vector<std::string> {});
  buildTest.setParameters(params);

  buildTest.build(entries0);
  TEST_TRUE(buildTest.testDigestion(peptides_unmod_no_minmax))
  TEST_TRUE(buildTest.peptidesSorted())
  TEST_TRUE(buildTest.fragmentsSorted())

  buildTest.clear();
  params.setValue("peptide:min_size", 2);
  params.setValue("peptide:max_size", 6);
  buildTest.setParameters(params);
  buildTest.build(entries0);
  TEST_TRUE(buildTest.testDigestion(peptides_unmod_minmax))
  TEST_TRUE(buildTest.peptidesSorted())
  TEST_TRUE(buildTest.fragmentsSorted())

  buildTest.clear();
  params.setValue("peptide:max_size", 100);
  params.setValue("peptide:missed_cleavages", 1);
  buildTest.setParameters(params);
  buildTest.build(entries0);
  TEST_TRUE(buildTest.testDigestion(peptides_unmod_minmax_missed_cleavage))
  TEST_TRUE(buildTest.peptidesSorted())
  TEST_TRUE(buildTest.fragmentsSorted())

  buildTest.clear();
  params.setValue("enzyme", "Trypsin");
  params.setValue("peptide:missed_cleavages", 0);
  params.setValue("peptide:min_mass", 0);
  params.setValue("peptide:min_size", 6);
  params.setValue("modifications:variable", std::vector<std::string> {"Oxidation (M)"});
  params.setValue("modifications:fixed", std::vector<std::string> {"Carbamidomethyl (C)"});
  buildTest.setParameters(params);
  buildTest.build(entries0);
  TEST_TRUE(buildTest.testDigestion(peptides_we_should_hit_mod))
  TEST_TRUE(buildTest.peptidesSorted())
  TEST_TRUE(buildTest.fragmentsSorted())
}
END_SECTION

// Verify that clear() resets the internal peptide container.
START_SECTION(clear())
{
  const std::vector<FASTAFile::FASTAEntry> entries0 {{"t", "t", "ARGEPADSSRKDFDMDMDM"}, {"t2", "t2", "HALLORTSCHS"}};
  FragmentIndex clearTest;
  clearTest.build(entries0);
  clearTest.clear();

  TEST_TRUE(clearTest.getPeptides().empty())
}
END_SECTION


////TEST Different Charges of the query Spectrum ////
// For each charge (1..4), a peptide's own theoretical spectrum should self-hit,
// with and without explicitly setting the precursor charge.
START_SECTION(void querySpectrum(const MSSpectrum& spectrum, SpectrumMatchesTopN& sms))
{
  const std::vector<FASTAFile::FASTAEntry> entries {
    {"test1", "test1",
    "MSDEREVAEAATGEDASSPPPKTEAASDPQHPAASEGAAAAAASPPLLRCLVLTGFGGYDKVKLQSRPAAPPAPGPGQLTLRLRACGLNFADLMARQGLYDRLPPLPVTPGMEGAGVVIAVGEGVSDRKAGDRVMVLNRSGMWQE"
    "EVTVPSVQTFLIPEAMTFEEAAALLVNYITAYMVLFDFGNLQPGHSVLVHMAAGGVGMAAVQLCRTVENVTVFGTASASKHEALKENGVTHPIDYHTTDYVDEIKKISPKGVDIVMDPLGGSDTAKGYNLLKPMGKVVTYGMANL"
    "LTGPKRNLMALARTWWNQFSVTALQLLQANRAVCGFHLGYLDGEVELVSGVVARLLALYNQGHIKPHIDSVWPFEKVADAMKQMQEKKNVGKVLLVPGPEKEN"}};

  FragmentIndex_test queryTest;

  auto params = queryTest.getParameters();
  params.setValue("fragment:max_charge", 4);
  params.setValue("precursor:min_charge", 1);
  params.setValue("precursor:max_charge", 4);
  params.setValue("fragment:min_mz", 0);
  // ensure all peptides/fragments are generated for exhaustive self-hit checks
  params.setValue("fragment:max_mz", 5000000);
  queryTest.setParameters(params);

  queryTest.build(entries);

  // Create different ms/ms spectra with different charges

  for (uint16_t charge = 1; charge <= 4; ++charge)
  {
    TEST_TRUE(queryTest.testQuery(charge, false, entries))
    TEST_TRUE(queryTest.testQuery(charge, true, entries))
  }
}
END_SECTION

// Shift the precursor by integer isotope errors [-3..3] and expect stable peptide window mapping.
START_SECTION(isotope_error)
{
  const std::vector<FASTAFile::FASTAEntry> entries {
    {"test1", "test1",
    "MSDEREVAEAATGEDASSPPPKTEAASDPQHPAASEGAAAAAASPPLLRCLVLTGFGGYDKVKLQSRPAAPPAPGPGQLTLRLRACGLNFADLMARQGLYDRLPPLPVTPGMEGAGVVIAVGEGVSDRKAGDRVMVLNRSGMWQE"
    "EVTVPSVQTFLIPEAMTFEEAAALLVNYITAYMVLFDFGNLQPGHSVLVHMAAGGVGMAAVQLCRTVENVTVFGTASASKHEALKENGVTHPIDYHTTDYVDEIKKISPKGVDIVMDPLGGSDTAKGYNLLKPMGKVVTYGMANL"
    "LTGPKRNLMALARTWWNQFSVTALQLLQANRAVCGFHLGYLDGEVELVSGVVARLLALYNQGHIKPHIDSVWPFEKVADAMKQMQEKKNVGKVLLVPGPEKEN"}};

  FragmentIndex_test isoTest;

  // Configure parameters before building the index (isotope error and fragment m/z bounds)
  auto params = isoTest.getParameters();
  params.setValue("precursor:isotope_error_min", -3);
  params.setValue("precursor:isotope_error_max", 3);
  params.setValue("fragment:min_mz", 0);
  params.setValue("fragment:max_mz", 90000);
  params.setValue("modifications:variable", std::vector<std::string> {});
  params.setValue("modifications:fixed", std::vector<std::string> {});
  isoTest.setParameters(params);

  // build after parameterization
  isoTest.build(entries);

  TheoreticalSpectrumGenerator tsg;
  PeakSpectrum b_y_ions;
  AASequence peptide = AASequence::fromString("EVAEAATGEDASSPPPK");
  tsg.getSpectrum(b_y_ions, peptide, 1, 1);
  MSSpectrum theo_spec;
  Precursor theo_prec;
  theo_prec.setCharge(1);
  theo_spec.setMSLevel(2);

  for (const auto& peak : b_y_ions)
  {
    theo_spec.push_back(peak);
  }

  for (int iso = -3; iso <= 3; ++iso)
  {
    theo_prec.setMZ(peptide.getMZ(1) + iso * Constants::C13C12_MASSDIFF_U);
    theo_spec.setPrecursors({theo_prec});
    FragmentIndex::SpectrumMatchesTopN sms;
    isoTest.querySpectrum(theo_spec, sms);
    bool found = false;

    for (const auto& hit : sms.hits_)
    {
      auto result = isoTest.getPeptides()[hit.peptide_idx_];
      auto psize = peptide.size();
      TEST_EQUAL(result.sequence_.first, 5)
      TEST_EQUAL(result.sequence_.second, psize)
      found = true;
    }
    TEST_TRUE(found);
  }
}
END_SECTION

// Apply small deterministic fragment m/z jitter and a precursor offset within tolerances;
// expect the correct peptide hit and zero isotope error.
START_SECTION(tolerance)
{
  const std::vector<FASTAFile::FASTAEntry> entries {
    {"test1", "test1",
     "MSDEREVAEAATGEDASSPPPKTEAASDPQHPAASEGAAAAAASPPLLRCLVLTGFGGYDKVKLQSRPAAPPAPGPGQLTLRLRACGLNFADLMARQGLYDRLPPLPVTPGMEGAGVVIAVGEGVSDRKAGDRVMVLNRSGMWQE"
     "EVTVPSVQTFLIPEAMTFEEAAALLVNYITAYMVLFDFGNLQPGHSVLVHMAAGGVGMAAVQLCRTVENVTVFGTASASKHEALKENGVTHPIDYHTTDYVDEIKKISPKGVDIVMDPLGGSDTAKGYNLLKPMGKVVTYGMANL"
     "LTGPKRNLMALARTWWNQFSVTALQLLQANRAVCGFHLGYLDGEVELVSGVVARLLALYNQGHIKPHIDSVWPFEKVADAMKQMQEKKNVGKVLLVPGPEKEN"}};

  FragmentIndex_test tolTest;

  auto params = tolTest.getParameters();
  params.setValue("fragment:min_mz", 0);
  params.setValue("fragment:max_mz", 90000);
  params.setValue("fragment:mass_tolerance", 0.05);
  params.setValue("fragment:mass_tolerance_unit", "Da");
  params.setValue("precursor:mass_tolerance", 2.0);
  params.setValue("precursor:mass_tolerance_unit", "Da");
  params.setValue("modifications:variable", std::vector<std::string> {});
  params.setValue("modifications:fixed", std::vector<std::string> {});
  tolTest.setParameters(params);

  tolTest.build(entries);

  TheoreticalSpectrumGenerator tsg;
  PeakSpectrum b_y_ions;

  AASequence peptide = AASequence::fromString("EVAEAATGEDASSPPPK");

  tsg.getSpectrum(b_y_ions, peptide, 1, 1);

  MSSpectrum theo_spec;
  Precursor theo_prec;
  theo_prec.setCharge(1);
  theo_prec.setMZ(peptide.getMZ(1) + 1.9);
  theo_spec.setMSLevel(2);
  theo_spec.setPrecursors({theo_prec});
  // Deterministic, small m/z jitter within ±0.045 Da to exercise tolerance handling
  constexpr float kJitterStep = 0.001f;
  constexpr int kJitterHalfWidth = 45;
  size_t i = 0;

  for (auto& peak : b_y_ions)
  {
    const float factor = (static_cast<int>(i % (2 * kJitterHalfWidth + 1)) - kJitterHalfWidth) * kJitterStep;
    peak.setMZ(peak.getMZ() + factor);
    theo_spec.push_back(peak);
    ++i;
  }

  FragmentIndex::SpectrumMatchesTopN sms;
  tolTest.querySpectrum(theo_spec, sms);
  bool found = false;
  for (const auto& hit : sms.hits_)
  {
    auto sequence = tolTest.getPeptides()[hit.peptide_idx_].sequence_;
    if ((sequence.first == 5) && (sequence.second == peptide.size()) && (hit.isotope_error_ == 0))
    {
      found = true;
      TEST_TRUE(hit.num_matched_ >= theo_spec.size());
    }
  }
  TEST_TRUE(found);
}
END_SECTION

START_SECTION(simd_scalar_equivalence_basic)
{
  // Test SIMD vs scalar equivalence with a standard protein
  auto fi = buildTestIndex(
    "EVAEAATGEDASSPPPKMDAILRGLNFEQLEEVIQTLHNAFSKENEPQLYPAIERIHFHTYRQVKAFPAYQLDKK"
    "MGKDYTRIVFVEDPTFHQIITSLENYRSMKARASNFKSARLAKYGVDLEKELVARFAQHTAIKNPNFTKLQGSR");

  TheoreticalSpectrumGenerator tsg;
  auto tp = tsg.getParameters();
  tp.setValue("add_b_ions", "true");
  tp.setValue("add_y_ions", "true");
  tsg.setParameters(tp);

  // Test with multiple peptides at different charges
  std::vector<std::string> test_peptides = {
    "EVAEAATGEDASSPPPK", "MDAILR", "GLNFEQLEEVIQTLHNAFSKENEPQLYPAIER",
    "AFPAYQLDKK", "DYTRIVFVEDPTFHQIITSLENYR"
  };

  for (const auto& pep_str : test_peptides)
  {
    AASequence pep = AASequence::fromString(pep_str);
    for (int charge = 1; charge <= 3; ++charge)
    {
      MSSpectrum spec;
      tsg.getSpectrum(spec, pep, charge, charge);
      float prec_mz = pep.getMZ(charge);
      auto range = fi.getPeptidesInPrecursorRange(prec_mz, {0.0f, 0.0f});
      size_t simd_n = 0, scalar_n = 0;
      TEST_TRUE(simdScalarMatch(fi, spec, range, charge, simd_n, scalar_n))
      TEST_EQUAL(simd_n, scalar_n)
    }
  }
}
END_SECTION

START_SECTION(simd_scalar_equivalence_real_fasta)
{
  // Use the SSE/Comet test FASTA (real protein sequences)
  auto fi = buildTestIndex(
    // proteins.fasta test sequence
    "GSMTVDMQEIGSTEMPYEVPTQPNATSASAGRGWFDGPSFKVPSVPTRPSGIFRRPSRIKPEFSFKEKVSELVSPAVYTFGLFVQNASESLTSDDPSDVPTQRTFKSDFQSV"
    "GSMTVDMQEIGSTEMPYEVPTQPNATSASAGRGWFDGPSFKVPSVPTRPSGIFRRPSRIKPEFSFKEKVSELVSPAVYTFGLFVQNASESLTSDDPSDVPTQRTFKSDFQSV"
    "AXXSTFDFYQRRLVTLAESPRAPSP"
    "GSMTVDMQEIGSTEMPYEVPTQPNATSASAGRGWFDGPSFKVPSVPTRPSGIFRRPSRIKPEFSFKEKVSELVSPAVYTFGLFVQNASESLTSDDPSDVPTQRTFKSDFQSV",
    0.05, "Da", 10.0, "ppm", 2); // 2 missed cleavages for more peptides

  TheoreticalSpectrumGenerator tsg;
  auto tp = tsg.getParameters();
  tp.setValue("add_b_ions", "true");
  tp.setValue("add_y_ions", "true");
  tsg.setParameters(tp);

  // Query with multiple known peptides from this sequence
  std::vector<std::string> peps = {"GSMTVDMQEIGSTEMPYEVPTQPNATSA", "GWFDGPSFK", "VPSVPTR",
                                    "PSGIFR", "EKVSELVSPAVYTFGLFVQNASESLTSDDPSDVPTQR"};
  for (const auto& p : peps)
  {
    AASequence pep = AASequence::fromString(p);
    MSSpectrum spec;
    tsg.getSpectrum(spec, pep, 1, 1);
    float prec_mz = pep.getMZ(1);
    auto range = fi.getPeptidesInPrecursorRange(prec_mz, {0.0f, 0.0f});
    size_t simd_n = 0, scalar_n = 0;
    TEST_TRUE(simdScalarMatch(fi, spec, range, 1, simd_n, scalar_n))
    TEST_EQUAL(simd_n, scalar_n)
  }
}
END_SECTION

START_SECTION(simd_scalar_edge_empty_range)
{
  // Edge case: empty peptide range [k, k) — should produce 0 hits in both
  auto fi = buildTestIndex("EVAEAATGEDASSPPPKMDAILR");

  Peak1D peak;
  peak.setMZ(500.0);
  peak.setIntensity(1.0);

  // Use a precursor mass that doesn't match anything
  auto range = fi.getPeptidesInPrecursorRange(99999.0f, {0.0f, 0.0f});
  TEST_EQUAL(range.first, range.second) // empty range

  auto simd_hits = fi.query(peak, range, 1);
  auto scalar_hits = fi.queryScalar(peak, range, 1);
  TEST_EQUAL(simd_hits.size(), 0)
  TEST_EQUAL(scalar_hits.size(), 0)
}
END_SECTION

START_SECTION(simd_scalar_edge_single_fragment_bucket)
{
  // Edge case: very short protein producing < 4 fragments (tests scalar remainder path only)
  auto fi = buildTestIndex("AAAAAKR", 0.05, "Da", 10.0, "ppm", 0, 5, 7);

  Peak1D peak;
  peak.setMZ(300.0);
  peak.setIntensity(1.0);
  auto range = fi.getPeptidesInPrecursorRange(300.0f, {-300.0f, 300.0f});

  auto simd_hits = fi.query(peak, range, 1);
  auto scalar_hits = fi.queryScalar(peak, range, 1);
  TEST_EQUAL(simd_hits.size(), scalar_hits.size())
}
END_SECTION

START_SECTION(simd_scalar_edge_ppm_tolerance)
{
  // Edge case: ppm tolerance instead of Da
  auto fi = buildTestIndex(
    "EVAEAATGEDASSPPPKMDAILRGLNFEQLEEVIQTLHNAFSKENEPQLYPAIER",
    20.0, "ppm"); // 20 ppm fragment tolerance

  TheoreticalSpectrumGenerator tsg;
  auto tp = tsg.getParameters();
  tp.setValue("add_b_ions", "true");
  tp.setValue("add_y_ions", "true");
  tsg.setParameters(tp);

  AASequence pep = AASequence::fromString("EVAEAATGEDASSPPPK");
  MSSpectrum spec;
  tsg.getSpectrum(spec, pep, 2, 2);
  float prec_mz = pep.getMZ(2);
  auto range = fi.getPeptidesInPrecursorRange(prec_mz, {0.0f, 0.0f});

  size_t simd_n = 0, scalar_n = 0;
  TEST_TRUE(simdScalarMatch(fi, spec, range, 2, simd_n, scalar_n))
  TEST_EQUAL(simd_n, scalar_n)
}
END_SECTION

START_SECTION(simd_scalar_edge_tolerance_boundary)
{
  // Edge case: peaks right at tolerance boundary
  auto fi = buildTestIndex("EVAEAATGEDASSPPPKMDAILR", 0.01, "Da");

  TheoreticalSpectrumGenerator tsg;
  auto tp = tsg.getParameters();
  tp.setValue("add_b_ions", "true");
  tp.setValue("add_y_ions", "true");
  tsg.setParameters(tp);

  AASequence pep = AASequence::fromString("EVAEAATGEDASSPPPK");
  MSSpectrum spec;
  tsg.getSpectrum(spec, pep, 1, 1);

  // Shift all peaks by exactly the tolerance (0.01 Da) — should still match
  MSSpectrum shifted_spec;
  for (const auto& peak : spec)
  {
    Peak1D p;
    p.setMZ(peak.getMZ() + 0.009); // just inside tolerance
    p.setIntensity(peak.getIntensity());
    shifted_spec.push_back(p);
  }

  float prec_mz = pep.getMZ(1);
  auto range = fi.getPeptidesInPrecursorRange(prec_mz, {0.0f, 0.0f});

  size_t simd_n = 0, scalar_n = 0;
  TEST_TRUE(simdScalarMatch(fi, shifted_spec, range, 1, simd_n, scalar_n))
  TEST_EQUAL(simd_n, scalar_n)

  // Shift peaks just outside tolerance — should produce 0 hits
  MSSpectrum outside_spec;
  for (const auto& peak : spec)
  {
    Peak1D p;
    p.setMZ(peak.getMZ() + 0.02); // outside 0.01 Da tolerance
    p.setIntensity(peak.getIntensity());
    outside_spec.push_back(p);
  }

  size_t simd_n2 = 0, scalar_n2 = 0;
  TEST_TRUE(simdScalarMatch(fi, outside_spec, range, 1, simd_n2, scalar_n2))
  TEST_EQUAL(simd_n2, scalar_n2)
  TEST_EQUAL(simd_n2, 0)
}
END_SECTION

START_SECTION(simd_scalar_edge_multi_charge)
{
  // Edge case: multiple fragment charges
  auto fi = buildTestIndex(
    "EVAEAATGEDASSPPPKMDAILRGLNFEQLEEVIQTLHNAFSKENEPQLYPAIER"
    "IHFHTYRQVKAFPAYQLDKKMGKDYTRIVFVEDPTFHQIITSLENYR");

  TheoreticalSpectrumGenerator tsg;
  auto tp = tsg.getParameters();
  tp.setValue("add_b_ions", "true");
  tp.setValue("add_y_ions", "true");
  tsg.setParameters(tp);

  AASequence pep = AASequence::fromString("GLNFEQLEEVIQTLHNAFSKENEPQLYPAIER");
  MSSpectrum spec;
  tsg.getSpectrum(spec, pep, 1, 3); // generate ions at charges 1-3

  float prec_mz = pep.getMZ(3);
  auto range = fi.getPeptidesInPrecursorRange(prec_mz, {0.0f, 0.0f});

  // Test with fragment charges 1, 2, and 3
  for (uint16_t fc = 1; fc <= 3; ++fc)
  {
    size_t simd_n = 0, scalar_n = 0;
    TEST_TRUE(simdScalarMatch(fi, spec, range, fc, simd_n, scalar_n))
    TEST_EQUAL(simd_n, scalar_n)
  }
}
END_SECTION

START_SECTION(simd_scalar_edge_wide_precursor_window)
{
  // Edge case: wide precursor window (many candidates spanning multiple buckets)
  auto fi = buildTestIndex(
    "EVAEAATGEDASSPPPKMDAILRGLNFEQLEEVIQTLHNAFSKENEPQLYPAIER"
    "IHFHTYRQVKAFPAYQLDKKMGKDYTRIVFVEDPTFHQIITSLENYR"
    "SMKARASNFKSARLAKYGVDLEKELVARFAQHTAIKNPNFTKLQGSR",
    0.05, "Da", 500.0, "Da"); // very wide precursor tolerance

  TheoreticalSpectrumGenerator tsg;
  auto tp = tsg.getParameters();
  tp.setValue("add_b_ions", "true");
  tp.setValue("add_y_ions", "true");
  tsg.setParameters(tp);

  AASequence pep = AASequence::fromString("EVAEAATGEDASSPPPK");
  MSSpectrum spec;
  tsg.getSpectrum(spec, pep, 1, 1);

  float prec_mz = pep.getMZ(1);
  auto range = fi.getPeptidesInPrecursorRange(prec_mz, {0.0f, 0.0f});

  size_t simd_n = 0, scalar_n = 0;
  TEST_TRUE(simdScalarMatch(fi, spec, range, 1, simd_n, scalar_n))
  TEST_EQUAL(simd_n, scalar_n)
  TEST_TRUE(simd_n > 0) // wide window should produce hits
}
END_SECTION

START_SECTION(simd_scalar_edge_soa_invariants)
{
  // Verify SoA buffer alignment invariant (via fragmentsSorted which accesses both arrays)
  auto fi = buildTestIndex("EVAEAATGEDASSPPPKMDAILRGLNFEQLEEVIQTLHNAFSKENEPQLYPAIER");
  TEST_TRUE(fi.fragmentsSorted())
  TEST_TRUE(fi.peptidesSorted())
}
END_SECTION

START_SECTION(benchmark_simd_vs_scalar)
{
  auto fi = buildTestIndex(
    "EVAEAATGEDASSPPPKMDAILRGLNFEQLEEVIQTLHNAFSKENEPQLYPAIERIHFHTYRQVKAFPAYQLDKK"
    "MGKDYTRIVFVEDPTFHQIITSLENYRSMKARASNFKSARLAKYGVDLEKELVARFAQHTAIKNPNFTKLQGSR");

  TheoreticalSpectrumGenerator tsg;
  auto tp = tsg.getParameters();
  tp.setValue("add_b_ions", "true");
  tp.setValue("add_y_ions", "true");
  tsg.setParameters(tp);

  AASequence pep = AASequence::fromString("GLNFEQLEEVIQTLHNAFSKENEPQLYPAIER");
  MSSpectrum spec;
  tsg.getSpectrum(spec, pep, 1, 1);
  float prec_mz = pep.getMZ(1);
  auto range = fi.getPeptidesInPrecursorRange(prec_mz, {0.0f, 0.0f});

  // Warm up
  for (int w = 0; w < 10; ++w)
    for (const auto& peak : spec)
    {
      fi.query(peak, range, 1);
      fi.queryScalar(peak, range, 1);
    }

  const int ITERATIONS = 1000;

  auto t0 = std::chrono::high_resolution_clock::now();
  size_t simd_total = 0;
  for (int iter = 0; iter < ITERATIONS; ++iter)
    for (const auto& peak : spec)
      simd_total += fi.query(peak, range, 1).size();
  auto t1 = std::chrono::high_resolution_clock::now();

  size_t scalar_total = 0;
  auto t2 = std::chrono::high_resolution_clock::now();
  for (int iter = 0; iter < ITERATIONS; ++iter)
    for (const auto& peak : spec)
      scalar_total += fi.queryScalar(peak, range, 1).size();
  auto t3 = std::chrono::high_resolution_clock::now();

  double simd_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
  double scalar_ms = std::chrono::duration<double, std::milli>(t3 - t2).count();

  std::cout << "\n=== FragmentIndex Query Benchmark ===" << std::endl;
  std::cout << "Iterations: " << ITERATIONS << " x " << spec.size() << " peaks" << std::endl;
  std::cout << "SIMD:   " << simd_ms << " ms (hits: " << simd_total << ")" << std::endl;
  std::cout << "Scalar: " << scalar_ms << " ms (hits: " << scalar_total << ")" << std::endl;
  std::cout << "Speedup: " << (scalar_ms / simd_ms) << "x" << std::endl;
  std::cout << "======================================\n" << std::endl;

  TEST_EQUAL(simd_total, scalar_total)
}
END_SECTION

END_TEST
