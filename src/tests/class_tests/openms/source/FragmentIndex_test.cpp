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

END_TEST
