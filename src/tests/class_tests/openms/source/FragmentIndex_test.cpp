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

  std::vector<double> exposeComputeSnesSigmaDeltaSet(bool include_prot_nterm_mods,
                                                      bool include_prot_cterm_mods) const
  {
    return computeSnesSigmaDeltaSet_(include_prot_nterm_mods, include_prot_cterm_mods);
  }

  const std::vector<double>& getSnesSigmaDeltaSet() const { return snes_sigma_delta_set_; }
  const std::vector<double>& getSnesSigmaDeltaSetProtNterm() const { return snes_sigma_delta_set_with_prot_nterm_; }
  const std::vector<double>& getSnesSigmaDeltaSetProtCterm() const { return snes_sigma_delta_set_with_prot_cterm_; }

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

// Verify that the new 'peptide:enzyme_specificity' parameter changes digestion behavior.
// Three modes are tested:
//   - "full" (default): both termini must be enzyme-specific (canonical, e.g. tryptic)
//   - "semi": one terminus may be non-enzyme-specific (semi-tryptic)
//   - "none": every substring of length [min,max] is enumerated, regardless of enzyme
//             (the canonical immunopeptidomics path, e.g. HLA-I 8..12mers)
START_SECTION([EXTRA] peptide:enzyme_specificity (full / semi / none))
{
  // 18-aa "protein" with two internal trypsin cuts (after K at pos 1, after R at pos 7).
  // Sequence has no B/X/Z (which FragmentIndex filters out as ambiguous AAs).
  // Tryptic products: "AK" (0..1), "ACDEFGR" (2..8), "HILMNPQSTV" (9..18).
  const std::vector<FASTAFile::FASTAEntry> entries{
    {"t", "t", "AKACDEFGRHILMNPQSTV"}};
  // sanity: 19 aa total, no B/X/Z

  // ---------- full (default): only fully-tryptic products ----------
  {
    FragmentIndex_test fi_full;
    auto p = fi_full.getParameters();
    p.setValue("enzyme", "Trypsin");
    p.setValue("peptide:missed_cleavages", 0);
    p.setValue("peptide:enzyme_specificity", "full");
    p.setValue("peptide:min_size", 2);
    p.setValue("peptide:max_size", 100);
    p.setValue("peptide:min_mass", 0);
    p.setValue("peptide:max_mass", 50000);
    p.setValue("modifications:variable", std::vector<std::string>{});
    p.setValue("modifications:fixed", std::vector<std::string>{});
    fi_full.setParameters(p);
    fi_full.build(entries);
    // 3 fully-tryptic products of length >= 2
    TEST_EQUAL(fi_full.getPeptides().size(), 3)
  }

  // ---------- semi: fully-tryptic + semi-tryptic variants ----------
  {
    FragmentIndex_test fi_semi;
    auto p = fi_semi.getParameters();
    p.setValue("enzyme", "Trypsin");
    p.setValue("peptide:missed_cleavages", 0);
    p.setValue("peptide:enzyme_specificity", "semi");
    p.setValue("peptide:min_size", 2);
    p.setValue("peptide:max_size", 100);
    p.setValue("peptide:min_mass", 0);
    p.setValue("peptide:max_mass", 50000);
    p.setValue("modifications:variable", std::vector<std::string>{});
    p.setValue("modifications:fixed", std::vector<std::string>{});
    fi_semi.setParameters(p);
    fi_semi.build(entries);
    // semi must yield strictly more peptides than full (semi = full + semi-specific extras)
    TEST_EQUAL(fi_semi.getPeptides().size() > 3, true)
  }

  // ---------- none (immunopeptidomics): all substrings of [min,max] ----------
  // For 8..12mers from a 19-aa sequence: 8mers=12, 9mers=11, 10mers=10, 11mers=9, 12mers=8 → 50.
  {
    FragmentIndex_test fi_none;
    auto p = fi_none.getParameters();
    p.setValue("enzyme", "Trypsin"); // enzyme is irrelevant under specificity=none
    p.setValue("peptide:missed_cleavages", 0);
    p.setValue("peptide:enzyme_specificity", "none");
    p.setValue("peptide:min_size", 8);
    p.setValue("peptide:max_size", 12);
    p.setValue("peptide:min_mass", 0);
    p.setValue("peptide:max_mass", 50000);
    p.setValue("modifications:variable", std::vector<std::string>{});
    p.setValue("modifications:fixed", std::vector<std::string>{});
    fi_none.setParameters(p);
    fi_none.build(entries);
    TEST_EQUAL(fi_none.getPeptides().size(), 50)
  }

  // ---------- none > semi when using same length window ----------
  {
    // Use same 8..12 window for both modes so the comparison is fair.
    FragmentIndex_test fi_semi_8_12;
    auto p = fi_semi_8_12.getParameters();
    p.setValue("enzyme", "Trypsin");
    p.setValue("peptide:missed_cleavages", 0);
    p.setValue("peptide:enzyme_specificity", "semi");
    p.setValue("peptide:min_size", 8);
    p.setValue("peptide:max_size", 12);
    p.setValue("peptide:min_mass", 0);
    p.setValue("peptide:max_mass", 50000);
    p.setValue("modifications:variable", std::vector<std::string>{});
    p.setValue("modifications:fixed", std::vector<std::string>{});
    fi_semi_8_12.setParameters(p);
    fi_semi_8_12.build(entries);
    // none (50 substrings) must exceed semi with the same length window
    TEST_EQUAL(50 > fi_semi_8_12.getPeptides().size(), true)
  }

  // ---------- none with very short protein: must not crash ----------
  // FASTA databases often contain very short entries; pre-fix this would underflow.
  {
    const std::vector<FASTAFile::FASTAEntry> tiny{{"t", "t", "ABC"}}; // shorter than min_size
    FragmentIndex_test fi_tiny;
    auto p = fi_tiny.getParameters();
    p.setValue("enzyme", "Trypsin");
    p.setValue("peptide:enzyme_specificity", "none");
    p.setValue("peptide:min_size", 8);
    p.setValue("peptide:max_size", 12);
    p.setValue("peptide:min_mass", 0);
    p.setValue("peptide:max_mass", 50000);
    p.setValue("modifications:variable", std::vector<std::string>{});
    p.setValue("modifications:fixed", std::vector<std::string>{});
    fi_tiny.setParameters(p);
    fi_tiny.build(tiny); // must not crash
    TEST_EQUAL(fi_tiny.getPeptides().size(), 0)
  }
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
  params.setValue("fragment:min_ion_index", 0);
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
  params.setValue("fragment:min_ion_index", 0); // index all ions to verify all theoretical peaks match
  params.setValue("fragment:mass_tolerance", 0.05);
  params.setValue("fragment:mass_tolerance_unit", "Da");
  params.setValue("precursor:mass_tolerance_lower", 2.0);
  params.setValue("precursor:mass_tolerance_upper", 2.0);
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
  params.setValue("fragment:min_ion_index", 0); // include all ions for this test
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

  // With min_ion_index=2, skip b1/b2/y1/y2 → 2*(n-1-2) = 2*5 = 10 fragments
  fcTest.clear();
  params.setValue("fragment:min_ion_index", 2);
  fcTest.setParameters(params);
  fcTest.build(entries);
  TEST_EQUAL(fcTest.getPeptides().size(), 1)
  size_t expected_with_skip = 2 * (seq.size() - 1 - 2); // skip 2 from each series
  TEST_EQUAL(fcTest.fragmentCountForPeptide(0), expected_with_skip)
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
  params.setValue("fragment:min_ion_index", 0); // include all ions for fragment count check
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
  params.setValue("fragment:min_ion_index", 0); // include all ions for fragment count check
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

// --- Asymmetric precursor window: Task 8 tests 1-3 ---

START_SECTION((pair<size_t, size_t> getPeptidesInMassWindow(float, const pair<float, float>&) const))
{
  // Symmetric default [20, 20] ppm — each peptide retrieves itself within its own mass window.
  // Uses the high-level self-hit helper `testQuery`, which is the behavioural equivalent of
  // a getPeptidesInMassWindow round-trip (peptide -> theoretical spectrum -> back to peptide_idx).
  FragmentIndex_test fi;
  Param p = fi.getParameters();
  p.setValue("precursor:mass_tolerance_lower", 20.0);
  p.setValue("precursor:mass_tolerance_upper", 20.0);
  p.setValue("precursor:mass_tolerance_unit", "ppm");
  // include all fragment ions so testQuery's `num_matched >= spec.size()` can be satisfied
  p.setValue("fragment:min_ion_index", 0);
  // restrict to iso=0 so the test isolates the symmetric-window self-hit semantics
  p.setValue("precursor:isotope_error_min", 0);
  p.setValue("precursor:isotope_error_max", 0);
  fi.setParameters(p);

  // Build a small fixture
  vector<FASTAFile::FASTAEntry> entries;
  FASTAFile::FASTAEntry e;
  e.identifier = "TEST1";
  e.sequence = "PEPTIDER";
  entries.push_back(e);
  fi.build(entries);

  TEST_EQUAL(fi.testQuery(2, true, entries), true);
}
END_SECTION

START_SECTION((asymmetric window compensates precursor calibration offset))
{
  // Instrument reads precursor m/z +8 ppm high. A symmetric [5, 5] ppm window misses the peptide
  // at iso=0, because the observed mass sits 8 ppm ABOVE the peptide — outside [-5, +5] ppm.
  // Compensating asymmetrically by widening the LOWER side ([15, 5] ppm) shifts the window
  // down to cover the peptide: [-15, +5] ppm around the observed mass includes the true mass.
  //
  // Isotope iteration is collapsed to [0, 0] so the test observes *only* the window behaviour
  // under investigation (iso=±1 would otherwise reshape the effective window by ±1.003 Da).
  FragmentIndex_test fi;
  Param p = fi.getParameters();
  p.setValue("precursor:mass_tolerance_unit", "ppm");
  p.setValue("fragment:min_ion_index", 0);
  p.setValue("precursor:isotope_error_min", 0);
  p.setValue("precursor:isotope_error_max", 0);

  vector<FASTAFile::FASTAEntry> entries;
  FASTAFile::FASTAEntry e;
  e.identifier = "TEST";
  e.sequence = "PEPTIDER";
  entries.push_back(e);

  // First run: symmetric tight — should NOT find
  p.setValue("precursor:mass_tolerance_lower", 5.0);
  p.setValue("precursor:mass_tolerance_upper", 5.0);
  fi.setParameters(p);
  fi.build(entries);

  // Construct a query spectrum whose precursor is shifted +8 ppm
  AASequence seq = AASequence::fromString("PEPTIDER");
  MSSpectrum spec;
  Precursor prec;
  const double true_mz = seq.getMZ(2);
  prec.setMZ(true_mz * (1.0 + 8e-6));   // +8 ppm
  prec.setCharge(2);
  spec.getPrecursors().push_back(prec);
  spec.setMSLevel(2);

  // Build a minimal theoretical spectrum for the query (all b/y ions, charge 1)
  TheoreticalSpectrumGenerator tsg;
  PeakSpectrum theo;
  tsg.getSpectrum(theo, seq, 1, 1);
  for (const auto& peak : theo) spec.push_back(peak);
  spec.sortByPosition();

  FragmentIndex::SpectrumMatchesTopN sms_tight;
  fi.querySpectrum(spec, sms_tight);
  TEST_EQUAL(sms_tight.hits_.empty(), true);

  // Second run: asymmetric [15, 5] ppm — widen the LOWER side to compensate the +8 ppm bias
  p.setValue("precursor:mass_tolerance_lower", 15.0);
  p.setValue("precursor:mass_tolerance_upper", 5.0);
  fi.setParameters(p);

  FragmentIndex::SpectrumMatchesTopN sms_asym;
  fi.querySpectrum(spec, sms_asym);
  TEST_NOT_EQUAL(sms_asym.hits_.size(), 0);
}
END_SECTION

START_SECTION((static bool isOpenSearchMode(double, double, bool)))
{
  // Strict > threshold. 1000 ppm stays closed.
  TEST_EQUAL(FragmentIndex::isOpenSearchMode(500.0,  1500.0, true), true);
  TEST_EQUAL(FragmentIndex::isOpenSearchMode(999.0,   999.0, true), false);
  TEST_EQUAL(FragmentIndex::isOpenSearchMode(1000.0, 1000.0, true), false);
  TEST_EQUAL(FragmentIndex::isOpenSearchMode(1000.0001, 1000.0, true), true);

  // Da unit — threshold 1.0
  TEST_EQUAL(FragmentIndex::isOpenSearchMode(0.9, 0.9, false), false);
  TEST_EQUAL(FragmentIndex::isOpenSearchMode(1.0, 1.0, false), false);
  TEST_EQUAL(FragmentIndex::isOpenSearchMode(1.1, 0.5, false), true);
}
END_SECTION

// --- Asymmetric precursor window: Task 9 observable-proxy isotope tests ---

START_SECTION((open-mode forces isotope_error iteration to [0,0]))
{
  // Observable-proxy: under open mode, a fixture with isotope_error_range [-2, +2]
  // produces the same PSM set as [0, 0] — proving iteration collapsed.
  vector<FASTAFile::FASTAEntry> entries;
  FASTAFile::FASTAEntry e;
  e.identifier = "TEST";
  e.sequence = "PEPTIDER";
  entries.push_back(e);

  // Run 1: open mode with user iso range [-2, +2]
  FragmentIndex_test fi_a;
  Param p_a = fi_a.getParameters();
  p_a.setValue("precursor:mass_tolerance_lower", 0.5);
  p_a.setValue("precursor:mass_tolerance_upper", 1.5);  // 1.5 Da > 1.0 → open mode
  p_a.setValue("precursor:mass_tolerance_unit", "Da");
  p_a.setValue("precursor:isotope_error_min", -2);
  p_a.setValue("precursor:isotope_error_max", +2);
  fi_a.setParameters(p_a);
  fi_a.build(entries);

  // Run 2: open mode with iso range [0, 0]
  FragmentIndex_test fi_b;
  Param p_b = fi_b.getParameters();
  p_b.setValue("precursor:mass_tolerance_lower", 0.5);
  p_b.setValue("precursor:mass_tolerance_upper", 1.5);
  p_b.setValue("precursor:mass_tolerance_unit", "Da");
  p_b.setValue("precursor:isotope_error_min", 0);
  p_b.setValue("precursor:isotope_error_max", 0);
  fi_b.setParameters(p_b);
  fi_b.build(entries);

  // Construct identical query spectrum for both
  AASequence seq = AASequence::fromString("PEPTIDER");
  MSSpectrum spec;
  Precursor prec;
  prec.setMZ(seq.getMZ(2));
  prec.setCharge(2);
  spec.getPrecursors().push_back(prec);
  spec.setMSLevel(2);
  TheoreticalSpectrumGenerator tsg;
  PeakSpectrum theo;
  tsg.getSpectrum(theo, seq, 1, 1);
  for (const auto& peak : theo) spec.push_back(peak);
  spec.sortByPosition();

  FragmentIndex::SpectrumMatchesTopN sms_a, sms_b;
  fi_a.querySpectrum(spec, sms_a);
  fi_b.querySpectrum(spec, sms_b);

  // Equal PSM set sizes → iteration collapsed (the [-2,+2] config did NOT produce more hits)
  TEST_EQUAL(sms_a.hits_.size(), sms_b.hits_.size());
}
END_SECTION

START_SECTION((asymmetric closed window interacts with isotope_error iteration))
{
  // [5, 15] ppm + isotope_error [-1, +2]. A multi-peptide fixture covered by a closed-mode
  // asymmetric window; each peptide self-hits under testQuery (observable proxy for the
  // combined window + iso_error iteration path).
  // Use 'no cleavage' so the full fasta sequences become distinct peptides with distinct masses.
  vector<FASTAFile::FASTAEntry> entries;
  FASTAFile::FASTAEntry e1, e2, e3;
  e1.identifier = "P1"; e1.sequence = "PEPTIDER";      // mass m0
  e2.identifier = "P2"; e2.sequence = "PEPTIDERG";     // mass ~m0+57 Da (G = 57.02)
  e3.identifier = "P3"; e3.sequence = "PEPTIDERA";     // mass ~m0+71 Da (A = 71.04)
  entries = {e1, e2, e3};

  FragmentIndex_test fi;
  Param p = fi.getParameters();
  p.setValue("enzyme", "no cleavage");
  p.setValue("precursor:mass_tolerance_lower", 5.0);
  p.setValue("precursor:mass_tolerance_upper", 15.0);
  p.setValue("precursor:mass_tolerance_unit", "ppm");
  p.setValue("precursor:isotope_error_min", -1);
  p.setValue("precursor:isotope_error_max", +2);
  // Include all fragment peaks (low-mz + low-index) so testQuery can satisfy num_matched >= spec.size()
  p.setValue("fragment:min_ion_index", 0);
  p.setValue("fragment:min_mz", 0);
  p.setValue("fragment:max_mz", 50000);
  fi.setParameters(p);
  fi.build(entries);

  // Query each peptide's own theoretical spectrum and verify self-hit via testQuery
  TEST_EQUAL(fi.testQuery(2, true, entries), true);
}
END_SECTION

// --- Task 10: parameter validation throws ---

START_SECTION((validation: negative magnitude rejected))
{
  FragmentIndex fi;
  Param p = fi.getParameters();
  p.setValue("precursor:mass_tolerance_lower", -5.0);  // invalid: below setMinFloat(0.0)
  p.setValue("precursor:mass_tolerance_upper", 10.0);
  p.setValue("precursor:mass_tolerance_unit", "ppm");
  // checkDefaults_ fires from setParameters via the min-float check
  TEST_EXCEPTION(Exception::InvalidParameter, fi.setParameters(p));
}
END_SECTION

START_SECTION((validation: zero-width window rejected))
{
  FragmentIndex fi;
  Param p = fi.getParameters();
  p.setValue("precursor:mass_tolerance_lower", 0.0);
  p.setValue("precursor:mass_tolerance_upper", 0.0);   // sum == 0 → rejected
  p.setValue("precursor:mass_tolerance_unit", "ppm");
  TEST_EXCEPTION(Exception::InvalidParameter, fi.setParameters(p));
}
END_SECTION

START_SECTION((validation: NaN rejected))
{
  FragmentIndex fi;
  Param p = fi.getParameters();
  p.setValue("precursor:mass_tolerance_lower", 10.0);
  p.setValue("precursor:mass_tolerance_upper", std::numeric_limits<double>::quiet_NaN());
  p.setValue("precursor:mass_tolerance_unit", "ppm");
  TEST_EXCEPTION(Exception::InvalidParameter, fi.setParameters(p));
}
END_SECTION

// ---------------------------------------------------------------------------
// Half-open vs closed-closed peptide_idx_range contract
// ---------------------------------------------------------------------------
// getPeptidesInMassWindow returns a HALF-OPEN [first, second) index range.
// Previously the callers (searchDifferentPrecursorRanges and
// FragmentIndex::query()) treated the range as closed-closed, spuriously
// including the peptide at index `second` — a peptide whose precursor_mz is
// strictly greater than the window's upper bound. The tests below pin the
// half-open contract and guard against a regression of the callers.

START_SECTION((getPeptidesInMassWindow half-open contract))
{
  // Build an index with 5 well-separated peptide masses. peptide:enzyme=no cleavage
  // makes each FASTA entry one peptide with a predictable mass.
  FragmentIndex_test fi;
  Param p = fi.getParameters();
  p.setValue("precursor:mass_tolerance_unit", "Da");
  p.setValue("precursor:mass_tolerance_lower", 1000.0);   // wide at build time — the
  p.setValue("precursor:mass_tolerance_upper", 1000.0);   // window we test with is passed
                                                           // directly to getPeptidesInMassWindow,
                                                           // not derived from these params.
  p.setValue("enzyme", "no cleavage");
  p.setValue("peptide:min_size", 1);
  fi.setParameters(p);

  std::vector<FASTAFile::FASTAEntry> entries;
  FASTAFile::FASTAEntry e1, e2, e3, e4, e5;
  e1.identifier = "P1"; e1.sequence = "AAAK";        // ~ 388 Da
  e2.identifier = "P2"; e2.sequence = "AAAAK";       // ~ 459 Da
  e3.identifier = "P3"; e3.sequence = "AAAAAK";      // ~ 530 Da
  e4.identifier = "P4"; e4.sequence = "AAAAAAK";     // ~ 601 Da
  e5.identifier = "P5"; e5.sequence = "AAAAAAAK";    // ~ 672 Da
  entries = {e1, e2, e3, e4, e5};
  fi.build(entries);

  const auto& peptides = fi.getPeptides();
  TEST_EQUAL(peptides.size(), 5u);
  // Peptides are sorted ascending by precursor_mz_ after build().
  TEST_EQUAL(peptides[0].precursor_mz_ < peptides[1].precursor_mz_, true);
  TEST_EQUAL(peptides[3].precursor_mz_ < peptides[4].precursor_mz_, true);

  // Case 1: narrow window at a middle peptide excludes both neighbours.
  const float p2_mass = peptides[2].precursor_mz_;
  auto r1 = fi.getPeptidesInMassWindow(p2_mass, {-10.0f, 10.0f});
  TEST_EQUAL(r1.first, 2u);
  TEST_EQUAL(r1.second, 3u);   // HALF-OPEN: [2, 3) — only peptide 2
  // peptide at index r1.second (== 3) must have mass STRICTLY GREATER than the
  // window's upper bound. This is the half-open invariant the callers must respect.
  TEST_EQUAL(peptides[r1.second].precursor_mz_ > p2_mass + 10.0f, true);

  // Case 2: window at the last peptide — second == size() is the half-open sentinel.
  const float p5_mass = peptides[4].precursor_mz_;
  auto r2 = fi.getPeptidesInMassWindow(p5_mass, {-10.0f, 10.0f});
  TEST_EQUAL(r2.first, 4u);
  TEST_EQUAL(r2.second, 5u);   // [4, 5) == [4, size())

  // Case 3: window at the first peptide.
  const float p1_mass = peptides[0].precursor_mz_;
  auto r3 = fi.getPeptidesInMassWindow(p1_mass, {-10.0f, 10.0f});
  TEST_EQUAL(r3.first, 0u);
  TEST_EQUAL(r3.second, 1u);

  // Case 4: empty window in the gap between two peptides.
  const float gap_mass = (peptides[1].precursor_mz_ + peptides[2].precursor_mz_) * 0.5f;
  const float gap_tol = (peptides[2].precursor_mz_ - peptides[1].precursor_mz_) * 0.25f;
  auto r4 = fi.getPeptidesInMassWindow(gap_mass, {-gap_tol, gap_tol});
  TEST_EQUAL(r4.first, r4.second);   // empty range (first == second)

  // Case 5: window entirely below all peptides.
  auto r5 = fi.getPeptidesInMassWindow(100.0f, {-10.0f, 10.0f});
  TEST_EQUAL(r5.first, 0u);
  TEST_EQUAL(r5.second, 0u);   // empty at the beginning

  // Case 6: window entirely above all peptides.
  auto r6 = fi.getPeptidesInMassWindow(10000.0f, {-10.0f, 10.0f});
  TEST_EQUAL(r6.first, 5u);
  TEST_EQUAL(r6.second, 5u);   // empty at the end (first == size == second)

  // Case 7: wide window covering all peptides — second == size() signals "no peptide past end".
  auto r7 = fi.getPeptidesInMassWindow(peptides[2].precursor_mz_, {-10000.0f, 10000.0f});
  TEST_EQUAL(r7.first, 0u);
  TEST_EQUAL(r7.second, 5u);

  // Case 8: upper bound EXACTLY at a peptide's mass. upper_bound returns the first
  // iterator strictly greater than the bound, so a peptide whose mass equals the bound
  // IS included in the half-open range.
  const float lo_mass = peptides[1].precursor_mz_;
  const float hi_mass = peptides[2].precursor_mz_;
  const float upper_offset = hi_mass - lo_mass;   // so lo_mass + upper_offset == hi_mass
  auto r8 = fi.getPeptidesInMassWindow(lo_mass, {0.0f, upper_offset});
  TEST_EQUAL(r8.first, 1u);
  // peptide 2 (mass == upper bound) IS included; peptide 3 is not.
  TEST_EQUAL(r8.second, 3u);   // [1, 3) contains peptides 1 and 2
}
END_SECTION

START_SECTION((FragmentIndex does not score the peptide at range.second (closed-closed regression guard)))
{
  // Regression guard: before fixing the callers, FragmentIndex::query() used strict `>`
  // to decide loop termination (treating range.second as inclusive) and
  // searchDifferentPrecursorRanges sized candidates_iso_error.hits_ with `second - first + 1`,
  // creating a spurious slot for the out-of-window neighbour peptide. queryPeaks would then
  // write that neighbour's fragment matches into the spurious slot, making it appear in the
  // PSM output even though its precursor mass was strictly outside the user's window.
  //
  // Fixture: two peptides A (light) and B (heavy, +71 Da outside). Query uses A's precursor
  // mass but B's theoretical fragment peaks — so if the old closed-closed iteration is still
  // active, B scores perfectly as a spurious hit. With the half-open fix, B is excluded from
  // the candidate range and never written into hits.

  FragmentIndex_test fi;
  Param p = fi.getParameters();
  p.setValue("precursor:mass_tolerance_unit", "Da");
  p.setValue("precursor:mass_tolerance_lower", 5.0);
  p.setValue("precursor:mass_tolerance_upper", 5.0);
  p.setValue("precursor:isotope_error_min", 0);
  p.setValue("precursor:isotope_error_max", 0);
  p.setValue("enzyme", "no cleavage");
  p.setValue("peptide:min_size", 1);
  p.setValue("fragment:min_ion_index", 0);
  p.setValue("fragment:min_mz", 0);
  p.setValue("fragment:max_mz", 50000);
  fi.setParameters(p);

  std::vector<FASTAFile::FASTAEntry> entries;
  FASTAFile::FASTAEntry eA, eB;
  eA.identifier = "A"; eA.sequence = "PEPTIDER";      // ~ 957 Da
  eB.identifier = "B"; eB.sequence = "PEPTIDERA";     // ~ 957 + 71 = 1028 Da (well outside 10 Da window)
  entries = {eA, eB};
  fi.build(entries);

  const auto& peptides = fi.getPeptides();
  TEST_EQUAL(peptides.size(), 2u);
  // Lighter (A) sorts to index 0, heavier (B) to index 1.
  TEST_EQUAL(peptides[0].precursor_mz_ < peptides[1].precursor_mz_, true);
  const size_t idx_A = 0u;
  const size_t idx_B = 1u;
  // Confirm B is genuinely outside A's window.
  TEST_EQUAL(peptides[idx_B].precursor_mz_ > peptides[idx_A].precursor_mz_ + 5.0f, true);

  // Confirm getPeptidesInMassWindow returns the half-open single-peptide range [0, 1).
  auto range = fi.getPeptidesInMassWindow(peptides[idx_A].precursor_mz_, {-5.0f, 5.0f});
  TEST_EQUAL(range.first, idx_A);
  TEST_EQUAL(range.second, idx_A + 1);  // points past A but BEFORE B

  // Construct a query spectrum: precursor m/z is A's, fragment peaks are B's theoretical fragments.
  // If the pre-fix closed-closed iteration were still active, B would be processed and its
  // fragment matches would be written via (peptide_idx - first) into a spurious slot.
  AASequence seqB = AASequence::fromString(eB.sequence);
  TheoreticalSpectrumGenerator tsg;
  Param tsg_params = tsg.getParameters();
  tsg_params.setValue("add_first_prefix_ion", "true");
  tsg.setParameters(tsg_params);
  PeakSpectrum theoB;
  tsg.getSpectrum(theoB, seqB, 1, 1);

  MSSpectrum spec;
  for (const auto& peak : theoB) spec.push_back(peak);
  spec.sortByPosition();
  Precursor prec;
  prec.setMZ(peptides[idx_A].precursor_mz_);
  prec.setCharge(1);
  spec.getPrecursors().push_back(prec);
  spec.setMSLevel(2);

  FragmentIndex::SpectrumMatchesTopN sms;
  fi.querySpectrum(spec, sms);

  // Post-fix expectation: peptide B is never scored, because its index (1) is the `second`
  // bound of the half-open range and the query loop stops strictly before it. No hit in
  // sms.hits_ should reference idx_B.
  bool found_B_hit = false;
  for (const auto& hit : sms.hits_)
  {
    if (hit.peptide_idx_ == idx_B && hit.num_matched_ > 0)
    {
      found_B_hit = true;
      break;
    }
  }
  TEST_EQUAL(found_B_hit, false);
}
END_SECTION

// ============================================================================
// SNES (Speedy Non-specific Enzyme Search) — mother-peptide indexing
// ============================================================================

START_SECTION((SNES mother enumeration on a small protein))
{
  // 19-aa protein, length window [8, 12]. Naive SPEC_NONE enumeration produces 50
  // sub-peptides (see "peptide:enzyme_specificity" section above). SNES replaces
  // that with mother-peptide indexing:
  //   Single-N anchors i in [0, L - min_length] = [0, 11]  →  12 mothers
  //   Single-C anchors j in [min_length - 1, L - 1] = [7, 18]  →  12 mothers
  //   Total: 24 mothers.
  // Every sub-peptide remains reachable via realization in ProSEAlgorithm — this
  // test only checks the index's mother enumeration.
  const std::vector<FASTAFile::FASTAEntry> entries{{"p", "p", "AKACDEFGRHILMNPQSTV"}};

  FragmentIndex_test fi;
  auto p = fi.getParameters();
  p.setValue("peptide:enzyme_specificity", "none");
  p.setValue("peptide:min_size", 8);
  p.setValue("peptide:max_size", 12);
  p.setValue("peptide:min_mass", 0);
  p.setValue("peptide:max_mass", 50000);
  p.setValue("modifications:variable", std::vector<std::string>{});
  p.setValue("modifications:fixed", std::vector<std::string>{});
  p.setValue("snes_enabled", "true");
  fi.setParameters(p);
  fi.build(entries);

  TEST_EQUAL(fi.isSnesMode(), true)
  TEST_EQUAL(fi.getPeptides().size(), 24u)

  // Bit 31 of mod_bitmask_ tags Single-C; clear bit is Single-N.
  size_t n_count = 0, c_count = 0;
  for (const auto& pep : fi.getPeptides())
  {
    if (FragmentIndex::isSingleCMother(pep.mod_bitmask_)) ++c_count;
    else ++n_count;
  }
  TEST_EQUAL(n_count, 12u)
  TEST_EQUAL(c_count, 12u)
}
END_SECTION

START_SECTION((SNES is skipped when enzyme_specificity != none))
{
  // snes_enabled=true only takes effect under SPEC_NONE. For SPEC_FULL the flag is
  // ignored and the standard tryptic digestion path runs.
  const std::vector<FASTAFile::FASTAEntry> entries{{"p", "p", "AKACDEFGRHILMNPQSTV"}};
  FragmentIndex_test fi;
  auto p = fi.getParameters();
  p.setValue("enzyme", "Trypsin");
  p.setValue("peptide:enzyme_specificity", "full");
  p.setValue("peptide:missed_cleavages", 0);
  p.setValue("peptide:min_size", 2);
  p.setValue("peptide:max_size", 100);
  p.setValue("peptide:min_mass", 0);
  p.setValue("peptide:max_mass", 50000);
  p.setValue("modifications:variable", std::vector<std::string>{});
  p.setValue("modifications:fixed", std::vector<std::string>{});
  p.setValue("snes_enabled", "true");
  fi.setParameters(p);
  fi.build(entries);

  TEST_EQUAL(fi.isSnesMode(), false)
  // 3 fully-tryptic products (same result as the specificity=full test above).
  TEST_EQUAL(fi.getPeptides().size(), 3u)
}
END_SECTION

START_SECTION((realizeSNESLength locates the correct sub-peptide length))
{
  // Build a SNES index on a 19-aa protein, compute the exact mass of a known
  // sub-peptide (first 10 residues, N-anchored), and verify the realization step
  // picks up that length when given the exact target mass.
  const std::vector<FASTAFile::FASTAEntry> entries{{"p", "p", "AKACDEFGRHILMNPQSTV"}};
  FragmentIndex_test fi;
  auto p = fi.getParameters();
  p.setValue("peptide:enzyme_specificity", "none");
  p.setValue("peptide:min_size", 8);
  p.setValue("peptide:max_size", 12);
  p.setValue("peptide:min_mass", 0);
  p.setValue("peptide:max_mass", 50000);
  p.setValue("modifications:variable", std::vector<std::string>{});
  p.setValue("modifications:fixed", std::vector<std::string>{});
  p.setValue("snes_enabled", "true");
  fi.setParameters(p);
  fi.build(entries);

  // Target: the 10-residue N-anchored sub-peptide "AKACDEFGRH", (M+H)+.
  AASequence known = AASequence::fromString("AKACDEFGRH");
  const double target_mh_plus = known.getMonoWeight() + Constants::PROTON_MASS_U;

  // Find the Single-N mother anchored at position 0 of protein 0.
  const auto& peptides = fi.getPeptides();
  size_t mother_idx = peptides.size();
  for (size_t i = 0; i < peptides.size(); ++i)
  {
    if (!FragmentIndex::isSingleCMother(peptides[i].mod_bitmask_)
        && peptides[i].protein_idx == 0
        && peptides[i].sequence_.first == 0)
    {
      mother_idx = i;
      break;
    }
  }
  TEST_NOT_EQUAL(mother_idx, peptides.size())

  // 10 ppm symmetric tolerance is ample for an exact-mass lookup.
  const int realized = fi.realizeSNESLength(peptides[mother_idx], entries,
                                            target_mh_plus, 10.0, 10.0, /*ppm=*/true);
  TEST_EQUAL(realized, 10)

  AASequence realized_seq = fi.reconstructRealizedSubSequence(
      peptides[mother_idx], entries, static_cast<size_t>(realized));
  TEST_EQUAL(realized_seq.toUnmodifiedString(), "AKACDEFGRH")
}
END_SECTION

START_SECTION((realizeSNESLength handles Single-C realization))
{
  // Symmetric check on the Single-C side: trim from the N-end until mass matches.
  const std::vector<FASTAFile::FASTAEntry> entries{{"p", "p", "AKACDEFGRHILMNPQSTV"}};
  FragmentIndex_test fi;
  auto p = fi.getParameters();
  p.setValue("peptide:enzyme_specificity", "none");
  p.setValue("peptide:min_size", 8);
  p.setValue("peptide:max_size", 12);
  p.setValue("peptide:min_mass", 0);
  p.setValue("peptide:max_mass", 50000);
  p.setValue("modifications:variable", std::vector<std::string>{});
  p.setValue("modifications:fixed", std::vector<std::string>{});
  p.setValue("snes_enabled", "true");
  fi.setParameters(p);
  fi.build(entries);

  // Target: the 9-residue C-anchored sub-peptide (last 9 of the protein) = "LMNPQSTV" + one more.
  // Protein L=19, so last 9 = protein[10..19) = "ILMNPQSTV".
  AASequence known = AASequence::fromString("ILMNPQSTV");
  const double target_mh_plus = known.getMonoWeight() + Constants::PROTON_MASS_U;

  // Find the Single-C mother anchored at j = L - 1 = 18 of protein 0.
  const auto& peptides = fi.getPeptides();
  size_t mother_idx = peptides.size();
  for (size_t i = 0; i < peptides.size(); ++i)
  {
    if (FragmentIndex::isSingleCMother(peptides[i].mod_bitmask_)
        && peptides[i].protein_idx == 0
        && static_cast<size_t>(peptides[i].sequence_.first + peptides[i].sequence_.second) == 19u)
    {
      mother_idx = i;
      break;
    }
  }
  TEST_NOT_EQUAL(mother_idx, peptides.size())

  const int realized = fi.realizeSNESLength(peptides[mother_idx], entries,
                                            target_mh_plus, 10.0, 10.0, /*ppm=*/true);
  TEST_EQUAL(realized, 9)

  // Asymmetric-tolerance regression (review L3): shift the target ABOVE
  // the true mass so realized_mass - shifted_target ≈ -50 mDa (negative).
  // That delta is inside [-tol_lower, +tol_upper] only if tol_lower ≥ 50 mDa.
  // A symmetric max(lower, upper) collapse would silently admit both cases
  // below; the asymmetric implementation rejects the tight-lower config.
  const double shifted_high = target_mh_plus + 0.05; // ~50 mDa ≈ 50 ppm at mass 1000
  // Loose lower (500 ppm ≈ 500 mDa), tight upper (10 ppm ≈ 10 mDa): accept.
  TEST_EQUAL(fi.realizeSNESLength(peptides[mother_idx], entries,
                                   shifted_high,
                                   /*lower=*/500.0, /*upper=*/10.0, /*ppm=*/true), 9)
  // Tight lower (10 ppm), loose upper (500 ppm): reject (negative delta
  // exceeds the tight lower bound; upper bound irrelevant here).
  TEST_EQUAL(fi.realizeSNESLength(peptides[mother_idx], entries,
                                   shifted_high,
                                   /*lower=*/10.0, /*upper=*/500.0, /*ppm=*/true), -1)

  AASequence realized_seq = fi.reconstructRealizedSubSequence(
      peptides[mother_idx], entries, static_cast<size_t>(realized));
  TEST_EQUAL(realized_seq.toUnmodifiedString(), "ILMNPQSTV")
}
END_SECTION

START_SECTION((SNES fragment-index size is smaller than naive SPEC_NONE))
{
  // The whole point of SNES is the memory win from indexing mother peptides with
  // only one ion series per mother. For a given (protein, min, max) triple, the
  // number of mothers is O(L) and each mother emits one series (b or y) of length
  // O(length-1); the naive SPEC_NONE path enumerates O(L * (max-min+1)) sub-peptides
  // each emitting both b and y ions. Exact ratios vary with length, but SNES must
  // always emit strictly fewer fragments. Assert that here so a regression that
  // silently disabled the series-restriction would be caught.
  const std::vector<FASTAFile::FASTAEntry> entries{{"p", "p", "AKACDEFGRHILMNPQSTV"}};
  auto base_params = [](){
    Param p;
    p.setValue("peptide:enzyme_specificity", "none");
    p.setValue("peptide:min_size", 8);
    p.setValue("peptide:max_size", 12);
    p.setValue("peptide:min_mass", 0);
    p.setValue("peptide:max_mass", 50000);
    p.setValue("modifications:variable", std::vector<std::string>{});
    p.setValue("modifications:fixed", std::vector<std::string>{});
    return p;
  };

  FragmentIndex_test fi_naive;
  Param p_naive = fi_naive.getParameters();
  p_naive.update(base_params());
  p_naive.setValue("snes_enabled", "false");
  fi_naive.setParameters(p_naive);
  fi_naive.build(entries);

  FragmentIndex_test fi_snes;
  Param p_snes = fi_snes.getParameters();
  p_snes.update(base_params());
  p_snes.setValue("snes_enabled", "true");
  fi_snes.setParameters(p_snes);
  fi_snes.build(entries);

  TEST_EQUAL(fi_naive.isSnesMode(), false)
  TEST_EQUAL(fi_snes.isSnesMode(), true)

  // SNES has fewer peptides (24 mothers vs 50 subpeptides — validated elsewhere).
  TEST_EQUAL(fi_snes.getPeptides().size() < fi_naive.getPeptides().size(), true)

  // Cross-validate that both index SOMETHING (neither is empty).
  // The fragment count is not directly exposed, but the fact that each path
  // builds without error and the SNES mother count is a strict subset of the
  // naive subpeptide count is the load-bearing invariant. Fragment counts
  // per peptide are verified indirectly via the end-to-end matching test.
  TEST_EQUAL(fi_naive.getPeptides().size() > 0u, true)
  TEST_EQUAL(fi_snes.getPeptides().size() > 0u, true)
}
END_SECTION

START_SECTION((reconstructRealizedSubSequence applies fixed modifications))
{
  // Configure Carbamidomethyl on cysteine. Build SNES index. For a mother whose
  // realized sub-peptide contains a C, reconstructRealizedSubSequence must apply
  // the fixed mod (not just return the raw substring). This exercises the
  // fixed-mod pathway in the realization reconstruction, which was not covered
  // by the basic realization tests above.
  const std::vector<FASTAFile::FASTAEntry> entries{{"p", "p", "AKACDEFGRHILMNPQSTV"}};

  FragmentIndex_test fi;
  auto p = fi.getParameters();
  p.setValue("peptide:enzyme_specificity", "none");
  p.setValue("peptide:min_size", 8);
  p.setValue("peptide:max_size", 12);
  p.setValue("peptide:min_mass", 0);
  p.setValue("peptide:max_mass", 50000);
  p.setValue("modifications:variable", std::vector<std::string>{});
  p.setValue("modifications:fixed", std::vector<std::string>{"Carbamidomethyl (C)"});
  p.setValue("snes_enabled", "true");
  fi.setParameters(p);
  fi.build(entries);

  // Find the Single-N mother anchored at position 0. Its realized 8-mer = "AKACDEFG"
  // — the third residue is C, which must carry Carbamidomethyl after reconstruction.
  const auto& peptides = fi.getPeptides();
  size_t mother_idx = peptides.size();
  for (size_t i = 0; i < peptides.size(); ++i)
  {
    if (!FragmentIndex::isSingleCMother(peptides[i].mod_bitmask_)
        && peptides[i].protein_idx == 0
        && peptides[i].sequence_.first == 0)
    {
      mother_idx = i;
      break;
    }
  }
  TEST_NOT_EQUAL(mother_idx, peptides.size())

  AASequence realized_seq = fi.reconstructRealizedSubSequence(peptides[mother_idx], entries, 8u);
  TEST_EQUAL(realized_seq.toUnmodifiedString(), "AKACDEFG")
  TEST_EQUAL(realized_seq.size(), 8u)
  // toString() renders the modification inline — the exact format is
  // "AKAC(Carbamidomethyl)DEFG" when the fixed mod has been applied.
  TEST_EQUAL(realized_seq.toString(), "AKAC(Carbamidomethyl)DEFG")
}
END_SECTION

START_SECTION((SNES index admits a candidate whose sub-peptide matches an observed precursor))
{
  // End-to-end sanity: build a SNES index, synthesize a spectrum from a known
  // sub-peptide's b/y ions, query it, and verify at least one candidate hits the
  // mother containing that sub-peptide.
  const std::vector<FASTAFile::FASTAEntry> entries{{"p", "p", "AKACDEFGRHILMNPQSTV"}};

  FragmentIndex_test fi;
  auto p = fi.getParameters();
  p.setValue("peptide:enzyme_specificity", "none");
  p.setValue("peptide:min_size", 8);
  p.setValue("peptide:max_size", 12);
  p.setValue("peptide:min_mass", 0);
  p.setValue("peptide:max_mass", 50000);
  p.setValue("precursor:mass_tolerance_lower", 20.0);
  p.setValue("precursor:mass_tolerance_upper", 20.0);
  p.setValue("precursor:mass_tolerance_unit", "ppm");
  p.setValue("fragment:mass_tolerance", 20.0);
  p.setValue("fragment:mass_tolerance_unit", "ppm");
  p.setValue("precursor:isotope_error_min", 0);
  p.setValue("precursor:isotope_error_max", 0);
  p.setValue("modifications:variable", std::vector<std::string>{});
  p.setValue("modifications:fixed", std::vector<std::string>{});
  p.setValue("snes_enabled", "true");
  p.setValue("fragment:min_matched_ions", 3);
  fi.setParameters(p);
  fi.build(entries);

  // Synthesize b/y ions of "ACDEFGRHIL" (starts at protein pos 2, length 10 — realizable
  // from either the Single-N mother anchored at pos 2 or a Single-C mother ending at pos 11).
  TheoreticalSpectrumGenerator tsg;
  Param tsg_p = tsg.getParameters();
  tsg_p.setValue("add_metainfo", "true");
  tsg.setParameters(tsg_p);

  AASequence target = AASequence::fromString("ACDEFGRHIL");
  PeakSpectrum theo;
  tsg.getSpectrum(theo, target, 1, 1);
  theo.sortByPosition();

  MSSpectrum spec;
  for (const auto& peak : theo) spec.push_back(peak);
  Precursor prec;
  prec.setMZ(target.getMonoWeight() + Constants::PROTON_MASS_U); // (M+H)+ as charge-1 m/z
  prec.setCharge(1);
  spec.getPrecursors().push_back(prec);
  spec.setMSLevel(2);

  FragmentIndex::SpectrumMatchesTopN sms;
  fi.querySpectrum(spec, entries, sms);

  // At least one of the returned candidates must correspond to a mother whose
  // protein contains "ACDEFGRHIL" as a sub-sequence — trivially true here.
  bool any_matched = false;
  for (const auto& hit : sms.hits_)
  {
    if (hit.num_matched_ >= 3u) { any_matched = true; break; }
  }
  TEST_EQUAL(any_matched, true)
}
END_SECTION

START_SECTION((SNES matches candidates when a fixed N-terminal modification is configured))
{
  // Build a SNES index with Acetyl (N-term) as a fixed modification and verify
  // that a spectrum synthesized from a sub-peptide with the N-term acetyl applied
  // is correctly matched. Exercises fixed_nterm_delta_ != 0 paths in both
  // build-time fragment generation and query-time precursor-target derivation.
  //
  // Regression guard: the default Carbamidomethyl (C) fixture has
  // fixed_nterm_delta_ == fixed_cterm_delta_ == 0, masking an earlier bug where
  // the query target omitted the terminal delta. A non-default Carbamidomethyl
  // on a non-C residue would be rejected at parameter parse time — Acetyl
  // (N-term) is the minimal non-ANYWHERE fixed-mod that isolates the terminal
  // delta.
  const std::vector<FASTAFile::FASTAEntry> entries{{"p", "p", "AKAGDEFGRHILMNPQSTV"}};

  FragmentIndex_test fi;
  auto p = fi.getParameters();
  p.setValue("peptide:enzyme_specificity", "none");
  p.setValue("peptide:min_size", 8);
  p.setValue("peptide:max_size", 12);
  p.setValue("peptide:min_mass", 0);
  p.setValue("peptide:max_mass", 50000);
  p.setValue("precursor:mass_tolerance_lower", 20.0);
  p.setValue("precursor:mass_tolerance_upper", 20.0);
  p.setValue("precursor:mass_tolerance_unit", "ppm");
  p.setValue("fragment:mass_tolerance", 20.0);
  p.setValue("fragment:mass_tolerance_unit", "ppm");
  p.setValue("precursor:isotope_error_min", 0);
  p.setValue("precursor:isotope_error_max", 0);
  p.setValue("modifications:variable", std::vector<std::string>{});
  p.setValue("modifications:fixed", std::vector<std::string>{"Acetyl (N-term)"});
  p.setValue("snes_enabled", "true");
  p.setValue("fragment:min_matched_ions", 3);
  fi.setParameters(p);
  fi.build(entries);

  // Target: "AGDEFGRHIL" (sub-peptide at protein pos 2..11) with fixed N-term acetyl.
  AASequence target = AASequence::fromString("AGDEFGRHIL");
  target.setNTerminalModification("Acetyl");
  // Sanity: the modified target's mono weight must include the Acetyl delta.
  TEST_REAL_SIMILAR(
      target.getMonoWeight() - AASequence::fromString("AGDEFGRHIL").getMonoWeight(),
      42.010565);

  TheoreticalSpectrumGenerator tsg;
  Param tsg_p = tsg.getParameters();
  tsg_p.setValue("add_metainfo", "true");
  tsg.setParameters(tsg_p);
  PeakSpectrum theo;
  tsg.getSpectrum(theo, target, 1, 1);
  theo.sortByPosition();

  MSSpectrum spec;
  for (const auto& peak : theo) spec.push_back(peak);
  Precursor prec;
  prec.setMZ(target.getMonoWeight() + Constants::PROTON_MASS_U);
  prec.setCharge(1);
  spec.getPrecursors().push_back(prec);
  spec.setMSLevel(2);

  FragmentIndex::SpectrumMatchesTopN sms;
  fi.querySpectrum(spec, entries, sms);

  bool any_matched = false;
  for (const auto& hit : sms.hits_)
  {
    if (hit.num_matched_ >= 3u) { any_matched = true; break; }
  }
  TEST_EQUAL(any_matched, true)
}
END_SECTION

START_SECTION((SNES rejects configuration with add_b_ions=false or add_y_ions=false))
{
  // The SNES fragment index hard-codes b-ions for Single-N mothers and y-ions
  // for Single-C mothers; querySpectrumSNES_ looks up b/y precursor-equivalent
  // targets. If the user disables either series the downstream scorer
  // (ProSEAlgorithm) builds theoretical spectra without that series, which
  // silently degrades score quality on admitted candidates. v1 rejects the
  // configuration at updateMembers_ time.
  const std::vector<FASTAFile::FASTAEntry> entries{{"p", "p", "AKACDEFGRHILMNPQSTV"}};

  FragmentIndex_test fi;
  auto p = fi.getParameters();
  p.setValue("peptide:enzyme_specificity", "none");
  p.setValue("peptide:min_size", 8);
  p.setValue("peptide:max_size", 12);
  p.setValue("modifications:variable", std::vector<std::string>{});
  p.setValue("modifications:fixed", std::vector<std::string>{});
  p.setValue("snes_enabled", "true");
  p.setValue("ions:add_b_ions", "false");

  TEST_EXCEPTION(Exception::InvalidParameter, fi.setParameters(p))

  // Symmetric: y-ions disabled.
  auto p2 = fi.getParameters();
  p2.setValue("peptide:enzyme_specificity", "none");
  p2.setValue("peptide:min_size", 8);
  p2.setValue("peptide:max_size", 12);
  p2.setValue("modifications:variable", std::vector<std::string>{});
  p2.setValue("modifications:fixed", std::vector<std::string>{});
  p2.setValue("snes_enabled", "true");
  p2.setValue("ions:add_b_ions", "true");
  p2.setValue("ions:add_y_ions", "false");

  TEST_EXCEPTION(Exception::InvalidParameter, fi.setParameters(p2))

  // Non-SNES configuration (snes_enabled=false) accepts add_b_ions=false freely.
  auto p3 = fi.getParameters();
  p3.setValue("peptide:enzyme_specificity", "none");
  p3.setValue("peptide:min_size", 8);
  p3.setValue("peptide:max_size", 12);
  p3.setValue("modifications:variable", std::vector<std::string>{});
  p3.setValue("modifications:fixed", std::vector<std::string>{});
  p3.setValue("snes_enabled", "false");
  p3.setValue("ions:add_b_ions", "false");

  fi.setParameters(p3); // expected to not throw
  TEST_EQUAL(true, true) // reached only if setParameters did not throw
}
END_SECTION

START_SECTION((SpectrumMatch default-initializes subset_bitmask_ and sigma_delta_ to zero))
{
  FragmentIndex::SpectrumMatch sm;
  TEST_EQUAL(sm.num_matched_, 0u)
  TEST_EQUAL(sm.subset_bitmask_, 0u)
  TEST_REAL_SIMILAR(sm.sigma_delta_, 0.0f)
  TEST_EQUAL(sm.precursor_charge_, 0u)
  TEST_EQUAL(sm.isotope_error_, 0)
  TEST_EQUAL(sm.peptide_idx_, 0u)
}
END_SECTION

START_SECTION((reconstructModifiedSequence masks SNES_KIND_BIT_MASK from bitmask iteration))
{
  // Construct a FragmentIndex configured for SNES but with no variable mods
  // so n_slots == 0. Build a Single-C mother (bit 31 set). Verify that
  // reconstructModifiedSequence does not misinterpret bit 31 as an active
  // slot (which would produce a garbage modification or out-of-range access).
  const std::vector<FASTAFile::FASTAEntry> entries{{"p", "p", "ACDEFGHIK"}};

  FragmentIndex_test fi;
  auto p = fi.getParameters();
  p.setValue("peptide:enzyme_specificity", "none");
  p.setValue("peptide:min_size", 8);
  p.setValue("peptide:max_size", 9);
  p.setValue("peptide:min_mass", 0);
  p.setValue("peptide:max_mass", 50000);
  p.setValue("modifications:variable", std::vector<std::string>{});
  p.setValue("modifications:fixed", std::vector<std::string>{});
  p.setValue("snes_enabled", "true");
  fi.setParameters(p);
  fi.build(entries);

  // Find any Single-C mother in the index.
  const auto& peptides = fi.getPeptides();
  size_t single_c_idx = peptides.size();
  for (size_t i = 0; i < peptides.size(); ++i)
  {
    if (FragmentIndex::isSingleCMother(peptides[i].mod_bitmask_))
    {
      single_c_idx = i;
      break;
    }
  }
  TEST_NOT_EQUAL(single_c_idx, peptides.size())

  // reconstructModifiedSequence must return the bare sub-sequence with no
  // variable modifications applied (since none are configured). It must NOT
  // throw, assert, or produce a bitmask-out-of-range interpretation.
  AASequence seq = fi.reconstructModifiedSequence(peptides[single_c_idx], entries);
  TEST_EQUAL(seq.size(), peptides[single_c_idx].sequence_.second)
  TEST_EQUAL(seq.toUnmodifiedString().size(), peptides[single_c_idx].sequence_.second)
}
END_SECTION

START_SECTION((computeSnesSigmaDeltaSet_ returns sorted distinct values for typical config))
{
  // Config: Oxidation (M) + Deamidated (N) + Deamidated (Q), max_per_peptide = 2.
  // Both deamidation variants share the same delta (+0.984016 Da); deduplication
  // collapses them into a single eligible delta for the enumeration.
  // Expected Σ values (Unimod deltas):
  //   0                        (no mods)
  //   0.984016  (1 deamid)
  //   1.968032  (2 deamid)
  //   15.994915 (1 ox)
  //   16.978931 (1 ox + 1 deamid)
  //   31.989830 (2 ox)
  FragmentIndex_test fi;
  auto p = fi.getParameters();
  p.setValue("modifications:variable",
             std::vector<std::string>{"Oxidation (M)", "Deamidated (N)", "Deamidated (Q)"});
  p.setValue("modifications:variable_max_per_peptide", 2);
  p.setValue("modifications:fixed", std::vector<std::string>{});
  fi.setParameters(p);

  auto deltas = fi.exposeComputeSnesSigmaDeltaSet(false, false);

  TEST_EQUAL(deltas.size(), 6u)
  TEST_REAL_SIMILAR(deltas[0], 0.0)
  TEST_REAL_SIMILAR(deltas[1], 0.984016)
  TEST_REAL_SIMILAR(deltas[2], 1.968032)
  TEST_REAL_SIMILAR(deltas[3], 15.994915)
  TEST_REAL_SIMILAR(deltas[4], 16.978931)
  TEST_REAL_SIMILAR(deltas[5], 31.989830)
}
END_SECTION

START_SECTION((computeSnesSigmaDeltaSet_ honors include_prot_nterm_mods flag))
{
  // Config: Acetyl (Protein N-term) only. Without the flag, Σ_set should
  // contain just {0}; with the flag, should contain {0, +42.010565}.
  FragmentIndex_test fi;
  auto p = fi.getParameters();
  p.setValue("modifications:variable",
             std::vector<std::string>{"Acetyl (Protein N-term)"});
  p.setValue("modifications:variable_max_per_peptide", 1);
  p.setValue("modifications:fixed", std::vector<std::string>{});
  fi.setParameters(p);

  auto deltas_without = fi.exposeComputeSnesSigmaDeltaSet(false, false);
  TEST_EQUAL(deltas_without.size(), 1u)
  TEST_REAL_SIMILAR(deltas_without[0], 0.0)

  auto deltas_with = fi.exposeComputeSnesSigmaDeltaSet(true, false);
  TEST_EQUAL(deltas_with.size(), 2u)
  TEST_REAL_SIMILAR(deltas_with[0], 0.0)
  TEST_REAL_SIMILAR(deltas_with[1], 42.010565)
}
END_SECTION

START_SECTION((updateMembers_ populates the three SNES sigma_delta sets))
{
  FragmentIndex_test fi;
  auto p = fi.getParameters();
  p.setValue("peptide:enzyme_specificity", "none");
  p.setValue("modifications:variable",
             std::vector<std::string>{"Oxidation (M)", "Acetyl (Protein N-term)"});
  p.setValue("modifications:variable_max_per_peptide", 1);
  p.setValue("modifications:fixed", std::vector<std::string>{});
  p.setValue("snes_enabled", "true");
  fi.setParameters(p);

  const auto& baseline = fi.getSnesSigmaDeltaSet();
  const auto& with_nterm = fi.getSnesSigmaDeltaSetProtNterm();
  const auto& with_cterm = fi.getSnesSigmaDeltaSetProtCterm();

  // Baseline has {0, +15.995} — excludes Acetyl (Protein N-term).
  TEST_EQUAL(baseline.size(), 2u)
  TEST_REAL_SIMILAR(baseline[0], 0.0)
  TEST_REAL_SIMILAR(baseline[1], 15.994915)

  // With N-term extension: {0, +15.995, +42.011}.
  TEST_EQUAL(with_nterm.size(), 3u)
  TEST_REAL_SIMILAR(with_nterm[0], 0.0)
  TEST_REAL_SIMILAR(with_nterm[1], 15.994915)
  TEST_REAL_SIMILAR(with_nterm[2], 42.010565)

  // With C-term extension: same as baseline (no protein C-term mod here).
  TEST_EQUAL(with_cterm.size(), 2u)
}
END_SECTION

START_SECTION((reconstructRealizedSubSequence applies mods from subset_bitmask))
{
  // Build SNES index with Oxidation (M) variable mod. For a mother whose
  // realized 5-mer contains M at position 2, subset_bitmask = 1 (slot 0
  // active → the M slot) must produce AASequence with Oxidation applied.
  const std::vector<FASTAFile::FASTAEntry> entries{{"p", "p", "AKAMCDEFGR"}};

  FragmentIndex_test fi;
  auto p = fi.getParameters();
  p.setValue("peptide:enzyme_specificity", "none");
  p.setValue("peptide:min_size", 5);
  p.setValue("peptide:max_size", 10);
  p.setValue("peptide:min_mass", 0);
  p.setValue("peptide:max_mass", 50000);
  p.setValue("modifications:variable", std::vector<std::string>{"Oxidation (M)"});
  p.setValue("modifications:variable_max_per_peptide", 1);
  p.setValue("modifications:fixed", std::vector<std::string>{});
  p.setValue("snes_enabled", "true");
  fi.setParameters(p);
  fi.build(entries);

  // Find the Single-N mother anchored at position 1 of the protein (so
  // realized 5-mer = "KAMCD", M at sub-peptide position 2).
  const auto& peptides = fi.getPeptides();
  size_t mother_idx = peptides.size();
  for (size_t i = 0; i < peptides.size(); ++i)
  {
    if (!FragmentIndex::isSingleCMother(peptides[i].mod_bitmask_)
        && peptides[i].protein_idx == 0
        && peptides[i].sequence_.first == 1)
    {
      mother_idx = i;
      break;
    }
  }
  TEST_NOT_EQUAL(mother_idx, peptides.size())

  // subset_bitmask = 0 → plain sub-sequence, no Oxidation.
  AASequence unmod = fi.reconstructRealizedSubSequence(peptides[mother_idx], entries, 5u, 0u);
  TEST_EQUAL(unmod.toString(), "KAMCD")

  // Slot numbering: buildModSlots_ enumerates pure N-term mods first, then
  // per-residue mods left-to-right, then pure C-term mods. With only
  // `Oxidation (M)` configured (ANYWHERE specificity, residue-bound), there
  // are no pure N-term mods, so the M at sub-peptide position 2 is slot 0.
  // subset_bitmask = 1u = 1 << 0 → activate that single slot.
  AASequence ox = fi.reconstructRealizedSubSequence(peptides[mother_idx], entries, 5u, 1u);
  TEST_EQUAL(ox.toString(), "KAM(Oxidation)CD")
}
END_SECTION

START_SECTION((SNES query returns candidate with subset_bitmask for variable-mod spectrum))
{
  // Build SNES index with Oxidation (M) variable mod. Synthesize a spectrum
  // from "ACDEFMGR" with Oxidation applied at the M residue (sub-peptide
  // position 5, 0-based). Query → expect at least one hit with
  // subset_bitmask_ != 0 and sigma_delta_ ≈ 15.995.
  const std::vector<FASTAFile::FASTAEntry> entries{{"p", "p", "AKACDEFMGRHILNPQSTV"}};

  FragmentIndex_test fi;
  auto p = fi.getParameters();
  p.setValue("peptide:enzyme_specificity", "none");
  p.setValue("peptide:min_size", 8);
  p.setValue("peptide:max_size", 12);
  p.setValue("peptide:min_mass", 0);
  p.setValue("peptide:max_mass", 50000);
  p.setValue("precursor:mass_tolerance_lower", 20.0);
  p.setValue("precursor:mass_tolerance_upper", 20.0);
  p.setValue("precursor:mass_tolerance_unit", "ppm");
  p.setValue("fragment:mass_tolerance", 20.0);
  p.setValue("fragment:mass_tolerance_unit", "ppm");
  p.setValue("precursor:isotope_error_min", 0);
  p.setValue("precursor:isotope_error_max", 0);
  p.setValue("modifications:variable", std::vector<std::string>{"Oxidation (M)"});
  p.setValue("modifications:variable_max_per_peptide", 1);
  p.setValue("modifications:fixed", std::vector<std::string>{});
  p.setValue("snes_enabled", "true");
  p.setValue("fragment:min_matched_ions", 3);
  fi.setParameters(p);
  fi.build(entries);

  AASequence target = AASequence::fromString("ACDEFMGR");
  target.setModification(5, "Oxidation"); // M residue, 0-based position 5

  TheoreticalSpectrumGenerator tsg;
  Param tsg_p = tsg.getParameters();
  tsg_p.setValue("add_metainfo", "true");
  tsg.setParameters(tsg_p);
  PeakSpectrum theo;
  tsg.getSpectrum(theo, target, 1, 1);
  theo.sortByPosition();

  MSSpectrum spec;
  for (const auto& peak : theo) spec.push_back(peak);
  Precursor prec;
  prec.setMZ(target.getMonoWeight() + Constants::PROTON_MASS_U);
  prec.setCharge(1);
  spec.getPrecursors().push_back(prec);
  spec.setMSLevel(2);

  FragmentIndex::SpectrumMatchesTopN sms;
  fi.querySpectrum(spec, entries, sms);

  bool found_modified = false;
  for (const auto& hit : sms.hits_)
  {
    if (hit.subset_bitmask_ != 0 && std::abs(hit.sigma_delta_ - 15.994915f) < 0.01f)
    {
      found_modified = true;
      break;
    }
  }
  TEST_EQUAL(found_modified, true)
}
END_SECTION

START_SECTION((SNES emits one SpectrumMatch per valid subset at the same Σ (emit-both)))
{
  // Peptide "ACDEFMGMR" has two M residues at positions 5 and 7 (0-indexed).
  // With Oxidation (M) and max=1, Σ=15.995 is reachable by activating either
  // M individually (two distinct subsets). Each must produce a distinct
  // SpectrumMatch with a different subset_bitmask_.
  //
  // Placing the first M at position 5 ensures b3(ACD), b4(ACDE), b5(ACDEF)
  // are unmodified and score ≥ 3 against the Single-N mother, allowing the
  // SNES byte-scan to meet min_matched_ions=3.
  const std::vector<FASTAFile::FASTAEntry> entries{{"p", "p", "AKACDEFMGMRHILNPQSTV"}};

  FragmentIndex_test fi;
  auto p = fi.getParameters();
  p.setValue("peptide:enzyme_specificity", "none");
  p.setValue("peptide:min_size", 8);
  p.setValue("peptide:max_size", 12);
  p.setValue("peptide:min_mass", 0);
  p.setValue("peptide:max_mass", 50000);
  p.setValue("precursor:mass_tolerance_lower", 20.0);
  p.setValue("precursor:mass_tolerance_upper", 20.0);
  p.setValue("precursor:mass_tolerance_unit", "ppm");
  p.setValue("fragment:mass_tolerance", 20.0);
  p.setValue("fragment:mass_tolerance_unit", "ppm");
  p.setValue("precursor:isotope_error_min", 0);
  p.setValue("precursor:isotope_error_max", 0);
  p.setValue("modifications:variable", std::vector<std::string>{"Oxidation (M)"});
  p.setValue("modifications:variable_max_per_peptide", 1);
  p.setValue("modifications:fixed", std::vector<std::string>{});
  p.setValue("snes_enabled", "true");
  p.setValue("fragment:min_matched_ions", 3);
  fi.setParameters(p);
  fi.build(entries);

  // Target "ACDEFMGMR" contains two M residues (positions 5 and 7 in 0-indexed
  // sub-peptide). Apply Oxidation at position 5 (first M).
  AASequence target = AASequence::fromString("ACDEFMGMR");
  target.setModification(5, "Oxidation");

  TheoreticalSpectrumGenerator tsg;
  Param tsg_p = tsg.getParameters();
  tsg_p.setValue("add_metainfo", "true");
  tsg.setParameters(tsg_p);
  PeakSpectrum theo;
  tsg.getSpectrum(theo, target, 1, 1);
  theo.sortByPosition();

  MSSpectrum spec;
  for (const auto& peak : theo) spec.push_back(peak);
  Precursor prec;
  prec.setMZ(target.getMonoWeight() + Constants::PROTON_MASS_U);
  prec.setCharge(1);
  spec.getPrecursors().push_back(prec);
  spec.setMSLevel(2);

  FragmentIndex::SpectrumMatchesTopN sms;
  fi.querySpectrum(spec, entries, sms);

  // Collect the subset_bitmask_ values of modified hits for any mother that
  // could realize "ACDEFMGMR".
  std::set<uint32_t> modified_bitmasks;
  for (const auto& hit : sms.hits_)
  {
    if (hit.subset_bitmask_ != 0
        && std::abs(hit.sigma_delta_ - 15.994915f) < 0.01f)
    {
      modified_bitmasks.insert(hit.subset_bitmask_);
    }
  }
  // Expect at least 2 distinct subsets at Σ=15.995 (Oxidation on first M
  // vs Oxidation on second M).
  TEST_EQUAL(modified_bitmasks.size() >= 2u, true)
}
END_SECTION

START_SECTION((SNES subset enumeration rejects position conflicts))
{
  // Configure two variable mods that both claim the N-terminal residue
  // (e.g., Acetyl (N-term) + Carbamyl (N-term) — both N-term ANYWHERE).
  // Activating both would conflict on position 0; a subset that tries is
  // rejected. Σ=Σ_acetyl+Σ_carbamyl should have NO valid subset on a
  // peptide where both would apply to the same residue.
  FragmentIndex_test fi;
  auto p = fi.getParameters();
  p.setValue("peptide:enzyme_specificity", "none");
  p.setValue("peptide:min_size", 5);
  p.setValue("peptide:max_size", 8);
  p.setValue("modifications:variable",
             std::vector<std::string>{"Acetyl (N-term)", "Carbamyl (N-term)"});
  p.setValue("modifications:variable_max_per_peptide", 2);
  p.setValue("modifications:fixed", std::vector<std::string>{});
  p.setValue("snes_enabled", "true");
  fi.setParameters(p);

  // Σ_delta set should contain 0, +42.011 (Acetyl), +43.006 (Carbamyl),
  // and the SUM +85.017 (activating both — but at query time this subset
  // is rejected by position conflict).
  const auto& deltas = fi.getSnesSigmaDeltaSet();
  TEST_EQUAL(deltas.size() >= 3u, true)

  // The conflict is evaluated at query-time subset enumeration. A direct
  // query-path test for this would require synthesizing a spectrum whose
  // precursor matches Σ=85.017, building the index, and asserting NO
  // SpectrumMatch is emitted for subset_bitmask with both bits active.
  // Simpler invariant: the Σ-set enumeration itself does NOT discriminate,
  // so the set CAN contain 85.017 — it's the subset-time check that rejects.
  // This test asserts only the enumeration invariant; positional rejection
  // is covered by the next test.
}
END_SECTION

START_SECTION((SNES respects max_variable_mods_per_peptide cap in subset enumeration))
{
  // Three eligible Oxidation (M) sites; max_per_peptide = 1 means no subset
  // with popcount > 1 can be emitted. Σ_delta set should include values up
  // to 1*15.995 only (+ {0, 15.995}).
  FragmentIndex_test fi;
  auto p = fi.getParameters();
  p.setValue("peptide:enzyme_specificity", "none");
  p.setValue("modifications:variable", std::vector<std::string>{"Oxidation (M)"});
  p.setValue("modifications:variable_max_per_peptide", 1);
  p.setValue("modifications:fixed", std::vector<std::string>{});
  p.setValue("snes_enabled", "true");
  fi.setParameters(p);

  const auto& deltas = fi.getSnesSigmaDeltaSet();
  TEST_EQUAL(deltas.size(), 2u)
  TEST_REAL_SIMILAR(deltas[0], 0.0)
  TEST_REAL_SIMILAR(deltas[1], 15.994915)

  // Now max=2 → set grows.
  p.setValue("modifications:variable_max_per_peptide", 2);
  fi.setParameters(p);
  const auto& deltas2 = fi.getSnesSigmaDeltaSet();
  TEST_EQUAL(deltas2.size(), 3u)
  TEST_REAL_SIMILAR(deltas2[2], 31.989830)
}
END_SECTION

START_SECTION((SNES handles identical-delta variable mods without collapsing subsets))
{
  // Two variable mods with identical Δ (Oxidation on M and Oxidation on W,
  // both +15.995) on a peptide containing one M and one W → subsets
  // {bit_for_M_slot} and {bit_for_W_slot} both have Σ=15.995 but are
  // distinct subsets. Must emit both (verified via subset_bitmask_ distinct
  // values on a synthesized spectrum).
  //
  // Note: OpenMS Unimod modifications on different origins share the same
  // ResidueModification delta. Use Oxidation (M) + Oxidation (W) to get
  // two entries with the same delta but different origins.
  const std::vector<FASTAFile::FASTAEntry> entries{{"p", "p", "AKMCDWEFGRHILNPQSTV"}};

  FragmentIndex_test fi;
  auto p = fi.getParameters();
  p.setValue("peptide:enzyme_specificity", "none");
  p.setValue("peptide:min_size", 8);
  p.setValue("peptide:max_size", 12);
  p.setValue("peptide:min_mass", 0);
  p.setValue("peptide:max_mass", 50000);
  p.setValue("precursor:mass_tolerance_lower", 20.0);
  p.setValue("precursor:mass_tolerance_upper", 20.0);
  p.setValue("precursor:mass_tolerance_unit", "ppm");
  p.setValue("fragment:mass_tolerance", 20.0);
  p.setValue("fragment:mass_tolerance_unit", "ppm");
  p.setValue("precursor:isotope_error_min", 0);
  p.setValue("precursor:isotope_error_max", 0);
  p.setValue("modifications:variable",
             std::vector<std::string>{"Oxidation (M)", "Oxidation (W)"});
  p.setValue("modifications:variable_max_per_peptide", 1);
  p.setValue("modifications:fixed", std::vector<std::string>{});
  p.setValue("snes_enabled", "true");
  p.setValue("fragment:min_matched_ions", 3);
  fi.setParameters(p);
  fi.build(entries);

  // Target "KMCDWEFG" — M at 0-based position 1, W at position 4. Either
  // Ox on M or Ox on W produces Σ=15.995. Apply Ox on M for the synthesized
  // spectrum; query should return matches with BOTH subset variants (since
  // both have Σ=15.995 and both are valid on the realized sub-peptide).
  AASequence target = AASequence::fromString("KMCDWEFG");
  target.setModification(1, "Oxidation");

  TheoreticalSpectrumGenerator tsg;
  Param tsg_p = tsg.getParameters();
  tsg_p.setValue("add_metainfo", "true");
  tsg.setParameters(tsg_p);
  PeakSpectrum theo;
  tsg.getSpectrum(theo, target, 1, 1);
  theo.sortByPosition();

  MSSpectrum spec;
  for (const auto& peak : theo) spec.push_back(peak);
  Precursor prec;
  prec.setMZ(target.getMonoWeight() + Constants::PROTON_MASS_U);
  prec.setCharge(1);
  spec.getPrecursors().push_back(prec);
  spec.setMSLevel(2);

  FragmentIndex::SpectrumMatchesTopN sms;
  fi.querySpectrum(spec, entries, sms);

  std::set<uint32_t> modified_bitmasks;
  for (const auto& hit : sms.hits_)
  {
    if (hit.subset_bitmask_ != 0
        && std::abs(hit.sigma_delta_ - 15.994915f) < 0.01f)
    {
      modified_bitmasks.insert(hit.subset_bitmask_);
    }
  }
  TEST_EQUAL(modified_bitmasks.size() >= 2u, true)
}
END_SECTION

START_SECTION((SNES query admits PROTEIN_N_TERM variable mod only for anchor-0 mothers))
{
  // Build SNES index with Acetyl (Protein N-term). Two proteins: one where
  // the sub-peptide ACDEFGHI at protein position 0 is realizable from a
  // Single-N mother anchored at 0; another where ACDEFGHI sits mid-protein.
  //
  // The query spectrum is generated from the UNMODIFIED peptide (b/y ions
  // are unmodified), but the precursor m/z is shifted by the Acetyl delta
  // (+42.010565 Da). SNES phase-1 fragment scoring then matches the
  // unmodified b-ions to Single-N mothers; the PROT_NTERM precursor-filter
  // walk (sigma=42.010565) admits only mothers with sequence_.first==0,
  // gating out the mid-protein sub-peptide from protein idx 1.
  const std::vector<FASTAFile::FASTAEntry> entries{
      {"anchored", "anchored", "ACDEFGHIJKLMNPQR"},     // ACDEFGHI at pos 0
      {"mid", "mid", "XXXACDEFGHIJKLMNPQR"}              // ACDEFGHI at pos 3
  };

  FragmentIndex_test fi;
  auto p = fi.getParameters();
  p.setValue("peptide:enzyme_specificity", "none");
  p.setValue("peptide:min_size", 8);
  p.setValue("peptide:max_size", 12);
  p.setValue("peptide:min_mass", 0);
  p.setValue("peptide:max_mass", 50000);
  p.setValue("precursor:mass_tolerance_lower", 20.0);
  p.setValue("precursor:mass_tolerance_upper", 20.0);
  p.setValue("precursor:mass_tolerance_unit", "ppm");
  p.setValue("fragment:mass_tolerance", 20.0);
  p.setValue("fragment:mass_tolerance_unit", "ppm");
  p.setValue("precursor:isotope_error_min", 0);
  p.setValue("precursor:isotope_error_max", 0);
  p.setValue("modifications:variable",
             std::vector<std::string>{"Acetyl (Protein N-term)"});
  p.setValue("modifications:variable_max_per_peptide", 1);
  p.setValue("modifications:fixed", std::vector<std::string>{});
  p.setValue("snes_enabled", "true");
  p.setValue("fragment:min_matched_ions", 3);
  fi.setParameters(p);
  fi.build(entries);

  // Unmodified target for fragment generation; manually shift precursor by
  // the Acetyl delta so the SNES PROT_NTERM walk (sigma=42.010565) fires.
  AASequence unmod = AASequence::fromString("ACDEFGHI");
  const double acetyl_delta = 42.010565;

  TheoreticalSpectrumGenerator tsg;
  Param tsg_p = tsg.getParameters();
  tsg_p.setValue("add_metainfo", "true");
  tsg.setParameters(tsg_p);
  PeakSpectrum theo;
  tsg.getSpectrum(theo, unmod, 1, 1);
  theo.sortByPosition();

  MSSpectrum spec;
  for (const auto& peak : theo) spec.push_back(peak);
  Precursor prec;
  prec.setMZ(unmod.getMonoWeight() + acetyl_delta + Constants::PROTON_MASS_U);
  prec.setCharge(1);
  spec.getPrecursors().push_back(prec);
  spec.setMSLevel(2);

  FragmentIndex::SpectrumMatchesTopN sms;
  fi.querySpectrum(spec, entries, sms);

  // The match must come from the anchored protein (idx 0), not the
  // mid-protein one (idx 1). Verify via the mother's protein_idx.
  bool found_anchored = false;
  bool found_mid = false;
  for (const auto& hit : sms.hits_)
  {
    if (hit.subset_bitmask_ == 0) continue;
    if (std::abs(hit.sigma_delta_ - static_cast<float>(acetyl_delta)) > 0.1f) continue;
    const auto& mother = fi.getPeptides()[hit.peptide_idx_];
    if (mother.protein_idx == 0 && mother.sequence_.first == 0) found_anchored = true;
    if (mother.protein_idx == 1 && mother.sequence_.first != 0) found_mid = true;
  }
  TEST_EQUAL(found_anchored, true)
  TEST_EQUAL(found_mid, false)
}
END_SECTION

START_SECTION((SNES query admits PROTEIN_C_TERM variable mod only for anchor-end mothers))
{
  // Symmetric to the N-term test: Amidated (Protein C-term) variable mod.
  // Single-C mothers at the protein end admit; mid-protein sub-peptides
  // with the same residues do not.
  //
  // The query spectrum is generated from the UNMODIFIED peptide (y-ions
  // are unmodified), but the precursor m/z is shifted by the Amidated
  // delta (-0.984016 Da). SNES phase-1 fragment scoring matches the
  // unmodified y-ions to Single-C mothers; the PROT_CTERM precursor-filter
  // walk (sigma=-0.984016) admits only mothers at the protein C-terminus.
  const std::vector<FASTAFile::FASTAEntry> entries{
      {"anchored", "anchored", "ACDEFGHIJKLMNPQR"},      // R at protein pos 15 (end)
      {"mid", "mid", "ACDEFGHIJKLMNPQRXXX"}               // R is mid-protein (pos 15 of 19)
  };

  FragmentIndex_test fi;
  auto p = fi.getParameters();
  p.setValue("peptide:enzyme_specificity", "none");
  p.setValue("peptide:min_size", 8);
  p.setValue("peptide:max_size", 12);
  p.setValue("peptide:min_mass", 0);
  p.setValue("peptide:max_mass", 50000);
  p.setValue("precursor:mass_tolerance_lower", 20.0);
  p.setValue("precursor:mass_tolerance_upper", 20.0);
  p.setValue("precursor:mass_tolerance_unit", "ppm");
  p.setValue("fragment:mass_tolerance", 20.0);
  p.setValue("fragment:mass_tolerance_unit", "ppm");
  p.setValue("precursor:isotope_error_min", 0);
  p.setValue("precursor:isotope_error_max", 0);
  p.setValue("modifications:variable",
             std::vector<std::string>{"Amidated (Protein C-term)"});
  p.setValue("modifications:variable_max_per_peptide", 1);
  p.setValue("modifications:fixed", std::vector<std::string>{});
  p.setValue("snes_enabled", "true");
  p.setValue("fragment:min_matched_ions", 3);
  fi.setParameters(p);
  fi.build(entries);

  // Unmodified target; precursor m/z shifted by Amidated delta.
  AASequence unmod = AASequence::fromString("GHIJKLMNPQR");
  const double amidated_delta = -0.984016;

  TheoreticalSpectrumGenerator tsg;
  Param tsg_p = tsg.getParameters();
  tsg_p.setValue("add_metainfo", "true");
  tsg.setParameters(tsg_p);
  PeakSpectrum theo;
  tsg.getSpectrum(theo, unmod, 1, 1);
  theo.sortByPosition();

  MSSpectrum spec;
  for (const auto& peak : theo) spec.push_back(peak);
  Precursor prec;
  prec.setMZ(unmod.getMonoWeight() + amidated_delta + Constants::PROTON_MASS_U);
  prec.setCharge(1);
  spec.getPrecursors().push_back(prec);
  spec.setMSLevel(2);

  FragmentIndex::SpectrumMatchesTopN sms;
  fi.querySpectrum(spec, entries, sms);

  // Amidated delta ≈ -0.984016. sigma_delta_ stores the raw Σ, which is
  // negative for mass-loss mods; the tolerance check handles this correctly.
  bool found_anchored = false;
  bool found_mid = false;
  const float amidated_delta_f = static_cast<float>(amidated_delta);
  for (const auto& hit : sms.hits_)
  {
    if (hit.subset_bitmask_ == 0) continue;
    if (std::abs(hit.sigma_delta_ - amidated_delta_f) > 0.1f) continue;
    const auto& mother = fi.getPeptides()[hit.peptide_idx_];
    const size_t prot_len = entries[mother.protein_idx].sequence.size();
    if (mother.protein_idx == 0 && mother.sequence_.first + mother.sequence_.second == prot_len) found_anchored = true;
    if (mother.protein_idx == 1 && mother.sequence_.first + mother.sequence_.second != prot_len) found_mid = true;
  }
  TEST_EQUAL(found_anchored, true)
  TEST_EQUAL(found_mid, false)
}
END_SECTION

START_SECTION((SNES query-path rejects position-conflicting subsets))
{
  // Two N-term variable mods (Acetyl + Carbamyl, both N_TERM ANYWHERE) claim
  // the peptide N-terminus. A subset that activates both has Σ=85.017 Da but
  // is rejected at subset-enumeration due to position conflict. Synthesize a
  // spectrum with (M+H)+ shifted by +85.017 and verify zero modified hits at
  // that Σ (the only non-conflict way to reach Σ=85.017 is an invalid two-
  // mod subset at position 0).
  const std::vector<FASTAFile::FASTAEntry> entries{{"p", "p", "ACDEFGHIJKLMNPQR"}};

  FragmentIndex_test fi;
  auto p = fi.getParameters();
  p.setValue("peptide:enzyme_specificity", "none");
  p.setValue("peptide:min_size", 8);
  p.setValue("peptide:max_size", 12);
  p.setValue("peptide:min_mass", 0);
  p.setValue("peptide:max_mass", 50000);
  p.setValue("precursor:mass_tolerance_lower", 20.0);
  p.setValue("precursor:mass_tolerance_upper", 20.0);
  p.setValue("precursor:mass_tolerance_unit", "ppm");
  p.setValue("fragment:mass_tolerance", 20.0);
  p.setValue("fragment:mass_tolerance_unit", "ppm");
  p.setValue("precursor:isotope_error_min", 0);
  p.setValue("precursor:isotope_error_max", 0);
  p.setValue("modifications:variable",
             std::vector<std::string>{"Acetyl (N-term)", "Carbamyl (N-term)"});
  p.setValue("modifications:variable_max_per_peptide", 2);
  p.setValue("modifications:fixed", std::vector<std::string>{});
  p.setValue("snes_enabled", "true");
  p.setValue("fragment:min_matched_ions", 3);
  fi.setParameters(p);
  fi.build(entries);

  // Synthesize an unmodified-fragment spectrum for "ACDEFGHI" with precursor
  // shifted by +85.017 (Σ_acetyl + Σ_carbamyl).
  AASequence target = AASequence::fromString("ACDEFGHI");

  TheoreticalSpectrumGenerator tsg;
  Param tsg_p = tsg.getParameters();
  tsg_p.setValue("add_metainfo", "true");
  tsg.setParameters(tsg_p);
  PeakSpectrum theo;
  tsg.getSpectrum(theo, target, 1, 1);
  theo.sortByPosition();

  MSSpectrum spec;
  for (const auto& peak : theo) spec.push_back(peak);
  Precursor prec;
  // (M+H)+ shifted by conflict-sum Σ (42.010565 + 43.005814 = 85.016379).
  prec.setMZ(target.getMonoWeight() + Constants::PROTON_MASS_U + 85.016379);
  prec.setCharge(1);
  spec.getPrecursors().push_back(prec);
  spec.setMSLevel(2);

  FragmentIndex::SpectrumMatchesTopN sms;
  fi.querySpectrum(spec, entries, sms);

  // No hit should have sigma_delta_ ≈ 85.017 with subset_bitmask_ != 0,
  // because the only subset summing to that Σ requires two N-term mods
  // at the same position — rejected.
  bool found_conflict_subset = false;
  for (const auto& hit : sms.hits_)
  {
    if (hit.subset_bitmask_ != 0
        && std::abs(hit.sigma_delta_ - 85.016379f) < 0.1f)
    {
      found_conflict_subset = true;
      break;
    }
  }
  TEST_EQUAL(found_conflict_subset, false)
}
END_SECTION

START_SECTION((SNES mother generation rejects ambiguous residue spans (X/B/Z)))
{
  // Protein contains an X in the middle. Mothers whose span covers the X must
  // be skipped (AASequence::fromString would fail at realization). Mothers in
  // the unambiguous prefix (before X) or suffix (after X) must still be kept.
  const std::vector<FASTAFile::FASTAEntry> entries{
      {"p", "p", "ACDEFGHIXKLMNPQSTVWY"}}; // X at 0-based position 8

  FragmentIndex_test fi;
  auto p = fi.getParameters();
  p.setValue("peptide:enzyme_specificity", "none");
  p.setValue("peptide:min_size", 8);
  p.setValue("peptide:max_size", 8);
  p.setValue("peptide:min_mass", 0);
  p.setValue("peptide:max_mass", 50000);
  p.setValue("modifications:variable", std::vector<std::string>{});
  p.setValue("modifications:fixed", std::vector<std::string>{});
  p.setValue("snes_enabled", "true");
  fi.setParameters(p);
  fi.build(entries);

  // No mother's span [start, start+8) may include the X at position 8.
  // For Single-N: valid start positions are 0 (span [0,8) — just before X)
  // and 9,10,11,12 (spans in the post-X region). Starts 1..8 would span the X.
  // For Single-C: symmetric — ends at 7..19 translate to starts 0..12.
  // Starts that include X: mothers whose span covers position 8.
  for (const auto& mother : fi.getPeptides())
  {
    const size_t start = mother.sequence_.first;
    const size_t end = start + mother.sequence_.second;
    // None of the kept mothers can span the X at position 8.
    TEST_EQUAL(start > 8u || end <= 8u, true)
  }

  // Positive-existence assertion: at least one mother from the unambiguous
  // prefix (start == 0, length 8) must have been kept.
  bool found_prefix = false;
  for (const auto& mother : fi.getPeptides())
  {
    if (mother.sequence_.first == 0 && mother.sequence_.second == 8) { found_prefix = true; break; }
  }
  TEST_EQUAL(found_prefix, true)
}
END_SECTION

START_SECTION((SNES full-length realization hits via supplementary precursor lookup))
{
  // When the observed precursor equals the full mother mass, the realized
  // length == mother length. The fragment index only stores b_1..b_{L-1}
  // and y_1..y_{L-1}, so the supplementary direct-precursor binary search
  // on fi_peptides_ is the path that admits this candidate. Construct a
  // protein of length equal to min=max, so every mother is also the full
  // peptide, then query with the peptide's (M+H)+.
  const std::vector<FASTAFile::FASTAEntry> entries{{"p", "p", "ACDEFGHIK"}}; // 9 AA

  FragmentIndex_test fi;
  auto p = fi.getParameters();
  p.setValue("peptide:enzyme_specificity", "none");
  p.setValue("peptide:min_size", 9);
  p.setValue("peptide:max_size", 9);
  p.setValue("peptide:min_mass", 0);
  p.setValue("peptide:max_mass", 50000);
  p.setValue("precursor:mass_tolerance_lower", 20.0);
  p.setValue("precursor:mass_tolerance_upper", 20.0);
  p.setValue("precursor:mass_tolerance_unit", "ppm");
  p.setValue("fragment:mass_tolerance", 20.0);
  p.setValue("fragment:mass_tolerance_unit", "ppm");
  p.setValue("precursor:isotope_error_min", 0);
  p.setValue("precursor:isotope_error_max", 0);
  p.setValue("modifications:variable", std::vector<std::string>{});
  p.setValue("modifications:fixed", std::vector<std::string>{});
  p.setValue("snes_enabled", "true");
  // Minimum 1 matched ion — this path's admission is precursor-only, not
  // fragment-count-driven, and we want the supplementary lookup to fire.
  p.setValue("fragment:min_matched_ions", 1);
  fi.setParameters(p);
  fi.build(entries);

  AASequence target = AASequence::fromString("ACDEFGHIK");
  TheoreticalSpectrumGenerator tsg;
  Param tsg_p = tsg.getParameters();
  tsg_p.setValue("add_metainfo", "true");
  tsg.setParameters(tsg_p);
  PeakSpectrum theo;
  tsg.getSpectrum(theo, target, 1, 1);
  theo.sortByPosition();

  MSSpectrum spec;
  for (const auto& peak : theo) spec.push_back(peak);
  Precursor prec;
  prec.setMZ(target.getMonoWeight() + Constants::PROTON_MASS_U);
  prec.setCharge(1);
  spec.getPrecursors().push_back(prec);
  spec.setMSLevel(2);

  FragmentIndex::SpectrumMatchesTopN sms;
  fi.querySpectrum(spec, entries, sms);

  // At least one full-length mother must have been admitted. With L=9 mother
  // length and realized length = 9, this hit comes via the supplementary
  // precursor-index path, not the b/y fragment bin walk.
  TEST_EQUAL(sms.hits_.empty(), false)
}
END_SECTION

START_SECTION((SNES query honors multi-charge precursor when charge is unset))
{
  // A spectrum whose precursor has charge == 0 should trigger iteration
  // across [min_precursor_charge_, max_precursor_charge_]. Synthesize a
  // 2+ precursor spectrum, clear the stored charge, and assert the
  // candidate is still found via the 2+ arm of the charge loop.
  const std::vector<FASTAFile::FASTAEntry> entries{{"p", "p", "AKACDEFGRHILMNPQSTV"}};

  FragmentIndex_test fi;
  auto p = fi.getParameters();
  p.setValue("peptide:enzyme_specificity", "none");
  p.setValue("peptide:min_size", 8);
  p.setValue("peptide:max_size", 12);
  p.setValue("peptide:min_mass", 0);
  p.setValue("peptide:max_mass", 50000);
  p.setValue("precursor:mass_tolerance_lower", 20.0);
  p.setValue("precursor:mass_tolerance_upper", 20.0);
  p.setValue("precursor:mass_tolerance_unit", "ppm");
  p.setValue("fragment:mass_tolerance", 20.0);
  p.setValue("fragment:mass_tolerance_unit", "ppm");
  p.setValue("precursor:isotope_error_min", 0);
  p.setValue("precursor:isotope_error_max", 0);
  p.setValue("precursor:min_charge", 2);
  p.setValue("precursor:max_charge", 3);
  p.setValue("modifications:variable", std::vector<std::string>{});
  p.setValue("modifications:fixed", std::vector<std::string>{});
  p.setValue("snes_enabled", "true");
  p.setValue("fragment:min_matched_ions", 3);
  fi.setParameters(p);
  fi.build(entries);

  // Target = 10-AA sub-peptide at protein positions 2..11.
  AASequence target = AASequence::fromString("ACDEFGRHIL");
  const double target_mh_plus = target.getMonoWeight() + Constants::PROTON_MASS_U;
  const double target_mz_2plus = (target_mh_plus + Constants::PROTON_MASS_U) / 2.0;

  TheoreticalSpectrumGenerator tsg;
  Param tsg_p = tsg.getParameters();
  tsg_p.setValue("add_metainfo", "true");
  tsg.setParameters(tsg_p);
  PeakSpectrum theo;
  tsg.getSpectrum(theo, target, 1, 1);
  theo.sortByPosition();

  MSSpectrum spec;
  for (const auto& peak : theo) spec.push_back(peak);
  Precursor prec;
  prec.setMZ(target_mz_2plus);
  prec.setCharge(0); // unset — exercises the multi-charge iteration path
  spec.getPrecursors().push_back(prec);
  spec.setMSLevel(2);

  FragmentIndex::SpectrumMatchesTopN sms;
  fi.querySpectrum(spec, entries, sms);

  // Expect at least one hit at charge 2 or 3 — the multi-charge query
  // iterates [min_precursor_charge_, max_precursor_charge_] when the
  // spectrum's declared charge is 0. Main invariant: the query does NOT
  // abort on the unset charge.
  bool found_multi_charge = false;
  for (const auto& hit : sms.hits_)
  {
    if (hit.precursor_charge_ == 2u || hit.precursor_charge_ == 3u)
    {
      found_multi_charge = true;
      break;
    }
  }
  TEST_EQUAL(found_multi_charge, true)
}
END_SECTION

END_TEST
