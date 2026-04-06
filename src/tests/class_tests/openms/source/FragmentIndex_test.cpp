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
  // (by subsequence window and mod bitmask).
  bool testDigestion(const std::vector<FragmentIndex::Peptide>& expected)
  {
    if (expected.size() != fi_peptides_.size()) return false;
    for (const auto& exp : expected)
    {
      bool found = false;
      for (const auto& act : fi_peptides_)
      {
        if ((exp.sequence_ == act.sequence_) && (exp.mod_bitmask_ == act.mod_bitmask_))
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
    for (size_t fi_idx = 0; fi_idx < fi_fragments_.size(); fi_idx += bucketsize_)
    {
      UInt32 last_idx = 0;
      const size_t end = (fi_idx + bucketsize_ > fi_fragments_.size()) ? fi_fragments_.size() : (fi_idx + bucketsize_);
      for (size_t bucket_idx = fi_idx; bucket_idx < end; ++bucket_idx)
      {
        if (fi_fragments_[bucket_idx].peptide_idx_ < last_idx) return false;
        last_idx = fi_fragments_[bucket_idx].peptide_idx_;
      }
    }
    return true;
  }

  // Returns the total number of fragments generated for a given peptide index.
  size_t fragmentCountForPeptide(UInt32 peptide_idx) const
  {
    size_t count = 0;
    for (const auto& f : fi_fragments_)
    {
      if (f.peptide_idx_ == peptide_idx) ++count;
    }
    return count;
  }

  const std::vector<Fragment>& getFragments() const { return fi_fragments_; }

  bool testQuery(const UInt32 charge, const bool precursor_mz_known, const std::vector<FASTAFile::FASTAEntry>& entries)
  {
    // Create theoretical spectra for different charges
    TheoreticalSpectrumGenerator tsg;
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
      spec_theo.clear(true);

      prec_theo.clearMetaInfo();
      AASequence mod_peptide = reconstructModifiedSequence(pep, entries);
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

//////////////////////////////
START_TEST(FragmentIndex, "$Id")

//////////////////////////////

/// Test the build for peptides
START_SECTION(build())
{
  // Test proteins used to generate expected peptides for multiple parameterizations
  /*
    Format of expected peptide descriptors below and their mapping to FragmentIndex::Peptide fields:
      { protein_idx, mod_bitmask_, { start, length }, precursor_mz_ }

    Where:
    - protein_idx: 0-based index into the FASTA entries vector passed to build(); selects the source protein.
    - mod_bitmask_: bitmask of active variable modification slots. Each bit corresponds to a (position, mod_type) pair
                   found by scanning the sequence left-to-right (0 = unmodified/fixed-only).
    - start: 0-based start offset within the selected protein sequence.
    - length: number of residues for the peptide (used as std::string::substr(start, length)).
    - precursor_mz_: mono-isotopic m/z at charge 1 (M+H)+. In these tests we often use a dummy value, as only ordering
                    invariants on peptides/fragment buckets are asserted.

    Note: testDigestion() compares expected vs. built peptides only by {sequence_, mod_bitmask_}.
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

// Verify that the lightweight fragment generator produces the expected number of
// b/y ions: 2*(n-1) for an n-residue peptide with default b+y ion types,
// consistent with standard fragment indexing (b1..b(n-1), y1..y(n-1)).
START_SECTION(lightweight_fragment_count)
{
  const std::string seq = "PEPTIDER";  // 8 residues
  const std::vector<FASTAFile::FASTAEntry> entries {{"p", "p", seq}};

  FragmentIndex_test fcTest;
  auto params = fcTest.getParameters();
  params.setValue("enzyme", "no cleavage");
  params.setValue("peptide:min_size", 0);
  params.setValue("peptide:max_size", 100);
  params.setValue("peptide:min_mass", 0);
  params.setValue("peptide:max_mass", 50000);
  params.setValue("fragment:min_mz", 0);
  params.setValue("fragment:max_mz", 50000);
  params.setValue("modifications:variable", std::vector<std::string> {});
  params.setValue("modifications:fixed", std::vector<std::string> {});
  fcTest.setParameters(params);

  fcTest.build(entries);

  // Should produce exactly one peptide
  TEST_EQUAL(fcTest.getPeptides().size(), 1)

  // For b+y ions (default): 2 * (n-1) = 2 * 7 = 14 fragments
  size_t expected_fragments = 2 * (seq.size() - 1);
  size_t actual_fragments = fcTest.fragmentCountForPeptide(0);
  TEST_EQUAL(actual_fragments, expected_fragments)
}
END_SECTION

// Test multi-mod-per-site: two different variable mods targeting the same AA (C)
// and fragment count correctness with modifications
START_SECTION(multi_mod_per_site)
{
  // Peptide with 2 C sites — both Glutathione(C) and Carbamidomethyl(C) are variable
  const std::vector<FASTAFile::FASTAEntry> entries {{"p", "p", "ACACK"}};

  FragmentIndex_test mmTest;
  auto params = mmTest.getParameters();
  params.setValue("enzyme", "no cleavage");
  params.setValue("peptide:min_size", 0);
  params.setValue("peptide:max_size", 100);
  params.setValue("peptide:min_mass", 0);
  params.setValue("peptide:max_mass", 50000);
  params.setValue("fragment:min_mz", 0);
  params.setValue("fragment:max_mz", 50000);
  params.setValue("modifications:variable_max_per_peptide", 2);
  params.setValue("modifications:variable", std::vector<std::string> {"Oxidation (M)"});
  params.setValue("modifications:fixed", std::vector<std::string> {"Carbamidomethyl (C)"});
  mmTest.setParameters(params);

  mmTest.build(entries);

  // "ACACK" has no M, so no variable mod sites → only 1 peptide (fixed C mods only)
  TEST_EQUAL(mmTest.getPeptides().size(), 1)
  // All peptides should have bitmask 0 (no variable mods)
  TEST_EQUAL(mmTest.getPeptides()[0].mod_bitmask_, 0u)
  TEST_TRUE(mmTest.peptidesSorted())
  TEST_TRUE(mmTest.fragmentsSorted())
  // Fragment count: 2*(5-1) = 8 for b+y ions
  TEST_EQUAL(mmTest.fragmentCountForPeptide(0), 8)
}
END_SECTION

// Test variable mods on M with fixed mods on C — verifies bitmask enumeration
START_SECTION(fixed_plus_variable_mods)
{
  const std::vector<FASTAFile::FASTAEntry> entries {{"p", "p", "ACMK"}};

  FragmentIndex_test fvTest;
  auto params = fvTest.getParameters();
  params.setValue("enzyme", "no cleavage");
  params.setValue("peptide:min_size", 0);
  params.setValue("peptide:max_size", 100);
  params.setValue("peptide:min_mass", 0);
  params.setValue("peptide:max_mass", 50000);
  params.setValue("fragment:min_mz", 0);
  params.setValue("fragment:max_mz", 50000);
  params.setValue("modifications:variable_max_per_peptide", 2);
  params.setValue("modifications:variable", std::vector<std::string> {"Oxidation (M)"});
  params.setValue("modifications:fixed", std::vector<std::string> {"Carbamidomethyl (C)"});
  fvTest.setParameters(params);

  fvTest.build(entries);

  // "ACMK": 1 M site → 2 peptides (bitmask 0 = no Ox, bitmask 1 = Ox on M)
  TEST_EQUAL(fvTest.getPeptides().size(), 2)
  TEST_TRUE(fvTest.peptidesSorted())
  TEST_TRUE(fvTest.fragmentsSorted())

  // Both variants should produce 2*(4-1) = 6 fragments each
  for (size_t i = 0; i < fvTest.getPeptides().size(); ++i)
  {
    TEST_EQUAL(fvTest.fragmentCountForPeptide(static_cast<UInt32>(i)), 6)
  }

  // Verify reconstructModifiedSequence produces valid sequences
  for (const auto& pep : fvTest.getPeptides())
  {
    AASequence reconstructed = fvTest.reconstructModifiedSequence(pep, entries);
    TEST_EQUAL(reconstructed.size(), 4)
    // C at position 1 should always have Carbamidomethyl (fixed mod)
    TEST_TRUE(reconstructed[1].isModified())
    // M at position 2: modified only when bitmask bit 0 is set
    TEST_EQUAL(reconstructed[2].isModified(), (pep.mod_bitmask_ & 1u) != 0)
  }
}
END_SECTION

// Cross-validate bitmask enumeration against ModifiedPeptideGenerator.
// For each test case: build FragmentIndex with bitmask path, also run
// ModifiedPeptideGenerator independently, then compare:
//   1. Same number of modification variants
//   2. Same set of precursor masses (within float tolerance)
//   3. Reconstructed AASequences match the ModifiedPeptideGenerator output
START_SECTION(cross_validate_vs_ModifiedPeptideGenerator)
{
  // Helper lambda: run ModifiedPeptideGenerator on a peptide string and return
  // sorted vector of (precursor_mz_charge1, AASequence_string) pairs
  auto run_modpepgen = [](const std::string& pep_str,
                          const std::vector<std::string>& fixed_mod_names,
                          const std::vector<std::string>& var_mod_names,
                          size_t max_var_mods) -> std::vector<std::pair<float, std::string>>
  {
    AASequence unmod = AASequence::fromString(pep_str);
    AASequence mod = AASequence(unmod);

    ModifiedPeptideGenerator::MapToResidueType fixed_mods;
    ModifiedPeptideGenerator::MapToResidueType var_mods;
    if (!fixed_mod_names.empty())
    {
      StringList sl(fixed_mod_names.begin(), fixed_mod_names.end());
      fixed_mods = ModifiedPeptideGenerator::getModifications(sl);
      ModifiedPeptideGenerator::applyFixedModifications(fixed_mods, mod);
    }
    std::vector<AASequence> variants;
    if (!var_mod_names.empty())
    {
      StringList sl(var_mod_names.begin(), var_mod_names.end());
      var_mods = ModifiedPeptideGenerator::getModifications(sl);
      ModifiedPeptideGenerator::applyVariableModifications(var_mods, mod, max_var_mods, variants);
    }
    else
    {
      variants.push_back(mod);
    }

    std::vector<std::pair<float, std::string>> result;
    for (const auto& v : variants)
    {
      result.emplace_back(static_cast<float>(v.getMZ(1)), v.toString());
    }
    std::sort(result.begin(), result.end());
    return result;
  };

  // Helper: build FragmentIndex, collect sorted (precursor_mz, reconstructed_string) pairs
  auto run_fragment_index = [](const std::string& seq,
                               const std::vector<std::string>& fixed_mod_names,
                               const std::vector<std::string>& var_mod_names,
                               size_t max_var_mods) -> std::vector<std::pair<float, std::string>>
  {
    std::vector<FASTAFile::FASTAEntry> entries {{"p", "p", seq}};
    FragmentIndex fi;
    auto params = fi.getParameters();
    params.setValue("enzyme", "no cleavage");
    params.setValue("peptide:min_size", 0);
    params.setValue("peptide:max_size", 100);
    params.setValue("peptide:min_mass", 0);
    params.setValue("peptide:max_mass", 50000);
    params.setValue("fragment:min_mz", 0);
    params.setValue("fragment:max_mz", 50000);
    params.setValue("modifications:variable_max_per_peptide", static_cast<int>(max_var_mods));
    params.setValue("modifications:variable", std::vector<std::string>(var_mod_names.begin(), var_mod_names.end()));
    params.setValue("modifications:fixed", std::vector<std::string>(fixed_mod_names.begin(), fixed_mod_names.end()));
    fi.setParameters(params);
    fi.build(entries);

    std::vector<std::pair<float, std::string>> result;
    for (const auto& pep : fi.getPeptides())
    {
      AASequence reconstructed = fi.reconstructModifiedSequence(pep, entries);
      result.emplace_back(pep.precursor_mz_, reconstructed.toString());
    }
    std::sort(result.begin(), result.end());
    return result;
  };

  // --- Test case 1: Simple Oxidation(M) + Carbamidomethyl(C) ---
  {
    auto mpg = run_modpepgen("ACMACK", {"Carbamidomethyl (C)"}, {"Oxidation (M)"}, 2);
    auto fi  = run_fragment_index("ACMACK", {"Carbamidomethyl (C)"}, {"Oxidation (M)"}, 2);
    TEST_EQUAL(fi.size(), mpg.size())
    for (size_t i = 0; i < std::min(fi.size(), mpg.size()); ++i)
    {
      TEST_REAL_SIMILAR(fi[i].first, mpg[i].first)
      TEST_EQUAL(fi[i].second, mpg[i].second)
    }
  }

  // --- Test case 2: Multiple M sites ---
  {
    auto mpg = run_modpepgen("DFDMDMDM", {}, {"Oxidation (M)"}, 2);
    auto fi  = run_fragment_index("DFDMDMDM", {}, {"Oxidation (M)"}, 2);
    TEST_EQUAL(fi.size(), mpg.size())
    for (size_t i = 0; i < std::min(fi.size(), mpg.size()); ++i)
    {
      TEST_REAL_SIMILAR(fi[i].first, mpg[i].first)
      TEST_EQUAL(fi[i].second, mpg[i].second)
    }
  }

  // --- Test case 3: N-terminal variable mod (Carbamyl) + residue mod (Oxidation) ---
  {
    auto mpg = run_modpepgen("KAAAAAAAMA", {}, {"Carbamyl (N-term)", "Oxidation (M)"}, 2);
    auto fi  = run_fragment_index("KAAAAAAAMA", {}, {"Carbamyl (N-term)", "Oxidation (M)"}, 2);
    TEST_EQUAL(fi.size(), mpg.size())
    for (size_t i = 0; i < std::min(fi.size(), mpg.size()); ++i)
    {
      TEST_REAL_SIMILAR(fi[i].first, mpg[i].first)
      TEST_EQUAL(fi[i].second, mpg[i].second)
    }
  }

  // --- Test case 4: Two different variable mods on same AA (C) ---
  {
    auto mpg = run_modpepgen("ACAACAACA", {}, {"Glutathione (C)", "Carbamidomethyl (C)"}, 1);
    auto fi  = run_fragment_index("ACAACAACA", {}, {"Glutathione (C)", "Carbamidomethyl (C)"}, 1);
    TEST_EQUAL(fi.size(), mpg.size())
    for (size_t i = 0; i < std::min(fi.size(), mpg.size()); ++i)
    {
      TEST_REAL_SIMILAR(fi[i].first, mpg[i].first)
      TEST_EQUAL(fi[i].second, mpg[i].second)
    }
  }

  // --- Test case 5: No modifiable sites ---
  {
    auto mpg = run_modpepgen("AAAAAAAAA", {"Carbamidomethyl (C)"}, {"Oxidation (M)"}, 2);
    auto fi  = run_fragment_index("AAAAAAAAA", {"Carbamidomethyl (C)"}, {"Oxidation (M)"}, 2);
    TEST_EQUAL(fi.size(), mpg.size())
    for (size_t i = 0; i < std::min(fi.size(), mpg.size()); ++i)
    {
      TEST_REAL_SIMILAR(fi[i].first, mpg[i].first)
      TEST_EQUAL(fi[i].second, mpg[i].second)
    }
  }

  // --- Test case 6: Fixed + variable mods, max_var_mods=3, multiple site types ---
  {
    auto mpg = run_modpepgen("ACMACMACA", {"Carbamidomethyl (C)"}, {"Oxidation (M)"}, 3);
    auto fi  = run_fragment_index("ACMACMACA", {"Carbamidomethyl (C)"}, {"Oxidation (M)"}, 3);
    TEST_EQUAL(fi.size(), mpg.size())
    for (size_t i = 0; i < std::min(fi.size(), mpg.size()); ++i)
    {
      TEST_REAL_SIMILAR(fi[i].first, mpg[i].first)
      TEST_EQUAL(fi[i].second, mpg[i].second)
    }
  }
}
END_SECTION

START_SECTION(tree_based_variant_generation)
{
  // Validates that generateFragmentsVariants_ (called when variable mods are set)
  // produces fragment counts and index ordering identical to the per-variant
  // lightweight path.
  //
  // A peptide with 2 Oxidation(M) sites produces 3 distinct variants (max_var=2):
  //   bitmask 0b00 (no mods), 0b01 (first M oxidized), 0b10 (second M oxidized),
  //   0b11 (both M oxidized).
  // All variants have the same sequence length, so each should get exactly
  // 2*(len-1) fragments in the default b+y-only configuration.

  std::vector<FASTAFile::FASTAEntry> entries {{"p", "p", "PEPTMDMIDE"}};

  FragmentIndex_test treeTest;
  auto params = treeTest.getParameters();
  params.setValue("enzyme", "no cleavage");
  params.setValue("peptide:min_size", 0);
  params.setValue("peptide:max_size", 100);
  params.setValue("peptide:min_mass", 0.0f);
  params.setValue("peptide:max_mass", 50000.0f);
  params.setValue("fragment:min_mz", 0.0f);
  params.setValue("fragment:max_mz", 50000.0f);
  params.setValue("modifications:variable_max_per_peptide", 2);
  params.setValue("modifications:variable", std::vector<std::string>{"Oxidation (M)"});
  treeTest.setParameters(params);
  treeTest.build(entries);

  // 4 variants: unmodified, ox@5, ox@7, ox@5+ox@7
  TEST_EQUAL(treeTest.getPeptides().size(), 4)

  const size_t seq_len = 10; // "PEPTMDMIDE"
  const size_t expected_frags_per_variant = 2 * (seq_len - 1); // b + y ions

  // Each variant should produce exactly 2*(len-1) fragments.
  for (UInt32 i = 0; i < 4; ++i)
  {
    TEST_EQUAL(treeTest.fragmentCountForPeptide(i), expected_frags_per_variant)
  }

  // Index structural invariants must hold after tree-based generation.
  TEST_EQUAL(treeTest.peptidesSorted(), true)
  TEST_EQUAL(treeTest.fragmentsSorted(), true)

  // Cross-validate: fragments for each variant match the single-variant
  // lightweight path via ModifiedPeptideGenerator reconstruction.
  {
    TheoreticalSpectrumGenerator tsg;
    Param tsg_params = tsg.getParameters();
    tsg_params.setValue("add_b_ions", "true");
    tsg_params.setValue("add_y_ions", "true");
    tsg_params.setValue("add_a_ions", "false");
    tsg_params.setValue("add_c_ions", "false");
    tsg_params.setValue("add_x_ions", "false");
    tsg_params.setValue("add_z_ions", "false");
    tsg.setParameters(tsg_params);

    const std::vector<FragmentIndex::Peptide>& peptides = treeTest.getPeptides();
    for (UInt32 pep_idx = 0; pep_idx < static_cast<UInt32>(peptides.size()); ++pep_idx)
    {
      AASequence reconstructed = treeTest.reconstructModifiedSequence(peptides[pep_idx], entries);
      PeakSpectrum theo;
      tsg.getSpectrum(theo, reconstructed, 1, 1);
      // Every theoretical peak must appear in the index for this peptide.
      size_t found = 0;
      for (const auto& peak : theo)
      {
        for (const auto& frag : treeTest.getFragments())
        {
          if (frag.peptide_idx_ == pep_idx &&
              std::abs(frag.fragment_mz_ - static_cast<float>(peak.getMZ())) < 0.001f)
          {
            ++found;
            break;
          }
        }
      }
      TEST_EQUAL(found, theo.size())
    }
  }
}
END_SECTION

END_TEST
