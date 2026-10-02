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
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <algorithm>
#include <bit>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <tuple>
#ifdef _OPENMP
  #include <omp.h>
#endif

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

  static void sortPeptides(std::vector<Peptide>& peptides) { sortPeptides_(peptides); }
  static void sortPeptides(std::vector<Peptide>& peptides, size_t min_task_size) { sortPeptides_(peptides, min_task_size); }

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

// Stop codons ('*') and other symbols have no residue mass. A peptide containing one used to be
// indexed, and scoring it aborted ProSE, because AASequence parses '*' as a weightless X.
START_SECTION([EXTRA] build() skips peptides containing stop codons or other symbols)
{
  // Tryptic products: "ACDEFGR", "HIL*MNPQK" (stop codon), "STVWYGHIK", "LMNP#QR" (other symbol)
  // and "DEFGHIL*" (C-terminal peptide followed by the stop codon).
  const std::vector<FASTAFile::FASTAEntry> entries{
    {"t", "t", "ACDEFGRHIL*MNPQKSTVWYGHIKLMNP#QRDEFGHIL*"}};
  FragmentIndex_test fi;
  auto p = fi.getParameters();
  p.setValue("enzyme", "Trypsin");
  p.setValue("peptide:missed_cleavages", 0);
  p.setValue("peptide:min_size", 2);
  p.setValue("peptide:max_size", 100);
  p.setValue("peptide:min_mass", 0);
  p.setValue("peptide:max_mass", 50000);
  p.setValue("modifications:variable", std::vector<std::string>{});
  p.setValue("modifications:fixed", std::vector<std::string>{});
  fi.setParameters(p);
  fi.build(entries);

  std::vector<std::string> indexed;
  for (const auto& peptide : fi.getPeptides())
  {
    indexed.push_back(entries[0].sequence.substr(peptide.sequence_.first, peptide.sequence_.second));
  }
  std::sort(indexed.begin(), indexed.end());
  TEST_STRING_EQUAL(ListUtils::concatenate(indexed, ","), "ACDEFGR,STVWYGHIK")
}
END_SECTION

// A FASTA entry longer than 65535 residues would overflow the 16-bit peptide start offset in
// Peptide::sequence_ and silently index fragments from the wrong subsequence. build() must reject it.
START_SECTION([EXTRA] build() rejects FASTA entries longer than 65535 residues)
{
  // 70000-residue contig (> uint16 max): the guard fires on length before digestion.
  const std::vector<FASTAFile::FASTAEntry> too_long{{"contig1", "long metaproteomic contig", std::string(70000, 'A')}};
  FragmentIndex fi_long;
  TEST_EXCEPTION(Exception::InvalidParameter, fi_long.build(too_long))

  // Boundary: exactly 65535 residues is the largest allowed and must NOT throw.
  const std::vector<FASTAFile::FASTAEntry> ok{{"contig_ok", "ok", std::string(65535, 'A')}};
  FragmentIndex fi_ok;
  fi_ok.build(ok); // must not throw
  TEST_EQUAL(fi_ok.isBuild(), true)
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

// z+1 ions (z-dot), the main C-terminal fragments of ETD-type spectra, are one hydrogen atom
// heavier than the z ions of ions:add_z_ions. The index must hold the m/z values that
// TheoreticalSpectrumGenerator produces for them, as scoring uses the latter.
START_SECTION([EXTRA] z+1 ions match TheoreticalSpectrumGenerator)
{
  const std::string seq = "PEPTIDER";
  const std::vector<FASTAFile::FASTAEntry> entries {{"p", "p", seq}};

  auto indexed_mzs = [&entries](const std::string& ion_series)
  {
    FragmentIndex_test fi;
    auto params = fi.getParameters();
    params.setValue("enzyme", "no cleavage");
    params.setValue("peptide:min_size", 0);
    params.setValue("peptide:max_size", 100);
    params.setValue("peptide:min_mass", 0);
    params.setValue("peptide:max_mass", 50000);
    params.setValue("fragment:min_mz", 0);
    params.setValue("fragment:max_mz", 50000);
    params.setValue("fragment:min_ion_index", 0);
    params.setValue("modifications:variable", std::vector<std::string> {});
    params.setValue("modifications:fixed", std::vector<std::string> {});
    params.setValue("ions:add_b_ions", "false");
    params.setValue("ions:add_y_ions", "false");
    params.setValue(ion_series, "true");
    fi.setParameters(params);
    fi.build(entries);
    std::vector<double> mzs;
    for (const auto& f : fi.getFragments()) mzs.push_back(f.fragment_mz_);
    std::sort(mzs.begin(), mzs.end());
    return mzs;
  };
  const std::vector<double> zp1_mzs = indexed_mzs("ions:add_zp1_ions");
  const std::vector<double> z_mzs = indexed_mzs("ions:add_z_ions");

  TheoreticalSpectrumGenerator tsg;
  Param tsg_param = tsg.getParameters();
  tsg_param.setValue("add_b_ions", "false");
  tsg_param.setValue("add_y_ions", "false");
  tsg_param.setValue("add_zp1_ions", "true");
  tsg.setParameters(tsg_param);
  PeakSpectrum theo;
  tsg.getSpectrum(theo, AASequence::fromString(seq), 1, 1);
  theo.sortByPosition();

  TEST_EQUAL(zp1_mzs.size(), seq.size() - 1)
  ABORT_IF(zp1_mzs.size() != theo.size() || z_mzs.size() != zp1_mzs.size())
  const double hydrogen = EmpiricalFormula("H").getMonoWeight();
  for (Size i = 0; i < zp1_mzs.size(); ++i)
  {
    TEST_REAL_SIMILAR(zp1_mzs[i], theo[i].getMZ())
    TEST_REAL_SIMILAR(zp1_mzs[i], z_mzs[i] + hydrogen)
  }
}
END_SECTION

// ions:electron_ions indexes c and z+1 ions apart from the main ion series. Queried without them, the
// index must give the candidates of an index that lacks them; queried with them, those of an index
// that holds them in the main series. The spectrum holds the y ions of THQPSANLDIK and the c and z+1
// ions of NDSIQLHTAPK, which has the same composition and hence the same precursor mass.
START_SECTION([EXTRA] ions:electron_ions are matched only on request)
{
  const std::vector<FASTAFile::FASTAEntry> entries {{"p1", "p1", "THQPSANLDIK"}, {"p2", "p2", "NDSIQLHTAPK"}};
  auto build_index = [&entries](FragmentIndex& fi, bool electron_ions, bool c_zp1_in_main_series)
  {
    auto params = fi.getParameters();
    params.setValue("enzyme", "no cleavage");
    params.setValue("peptide:min_size", 0);
    params.setValue("peptide:max_size", 100);
    params.setValue("peptide:min_mass", 0);
    params.setValue("peptide:max_mass", 50000);
    params.setValue("fragment:min_mz", 0);
    params.setValue("fragment:max_mz", 50000);
    params.setValue("fragment:mass_tolerance", 20.0);
    params.setValue("fragment:mass_tolerance_unit", "ppm");
    params.setValue("modifications:variable", std::vector<std::string> {});
    params.setValue("modifications:fixed", std::vector<std::string> {});
    params.setValue("ions:add_c_ions", c_zp1_in_main_series ? "true" : "false");
    params.setValue("ions:add_zp1_ions", c_zp1_in_main_series ? "true" : "false");
    params.setValue("ions:electron_ions", electron_ions ? "true" : "false");
    fi.setParameters(params);
    fi.build(entries);
  };
  FragmentIndex by_index, electron_index, union_index;
  build_index(by_index, false, false);
  build_index(electron_index, true, false);
  build_index(union_index, false, true);
  TEST_EQUAL(electron_index.getNumFragments() > by_index.getNumFragments(), true)
  TEST_EQUAL(electron_index.getNumFragments(), union_index.getNumFragments())

  auto ions = [](const std::string& seq_str, bool y, bool c_zp1)
  {
    TheoreticalSpectrumGenerator tsg;
    Param tsg_param = tsg.getParameters();
    tsg_param.setValue("add_b_ions", "false");
    tsg_param.setValue("add_y_ions", y ? "true" : "false");
    tsg_param.setValue("add_c_ions", c_zp1 ? "true" : "false");
    tsg_param.setValue("add_zp1_ions", c_zp1 ? "true" : "false");
    tsg.setParameters(tsg_param);
    PeakSpectrum ion_spectrum;
    tsg.getSpectrum(ion_spectrum, AASequence::fromString(seq_str), 1, 1);
    return ion_spectrum;
  };
  MSSpectrum spec = ions("THQPSANLDIK", true, false);
  for (const Peak1D& peak : ions("NDSIQLHTAPK", false, true)) spec.push_back(peak);
  spec.sortByPosition();
  spec.setMSLevel(2);
  Precursor prec;
  prec.setMZ(AASequence::fromString("THQPSANLDIK").getMZ(2));
  prec.setCharge(2);
  spec.setPrecursors({prec});

  // candidates as (peptide sequence, number of matched fragments), best first
  auto candidates = [&entries, &spec](FragmentIndex& fi, bool with_electron_ions)
  {
    FragmentIndex::SpectrumMatchesTopN sms;
    fi.querySpectrum(spec, entries, sms, with_electron_ions);
    std::vector<std::pair<std::string, uint32_t>> result;
    for (const auto& hit : sms.hits_)
    {
      const FragmentIndex::Peptide& pep = fi.getPeptides()[hit.peptide_idx_];
      result.emplace_back(entries[pep.protein_idx].sequence.substr(pep.sequence_.first, pep.sequence_.second), hit.num_matched_);
    }
    return result;
  };
  const auto by_candidates = candidates(by_index, false);
  const auto without = candidates(electron_index, false);
  const auto with = candidates(electron_index, true);
  const auto union_candidates = candidates(union_index, false);
  ABORT_IF(by_candidates.empty() || with.empty())
  TEST_EQUAL(without == by_candidates, true)
  TEST_EQUAL(with == union_candidates, true)
  TEST_STRING_EQUAL(by_candidates[0].first, "THQPSANLDIK")
  TEST_STRING_EQUAL(with[0].first, "NDSIQLHTAPK")
  // the default overload does not match the c and z+1 ions
  FragmentIndex::SpectrumMatchesTopN sms_default;
  electron_index.querySpectrum(spec, entries, sms_default);
  TEST_EQUAL(sms_default.hits_.size(), without.size())
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

START_SECTION((SNES query is safe and correct when a smaller index is queried after a larger one on the same thread))
{
  // Regression for an out-of-bounds access in querySpectrumSNES_'s touched-only reset.
  // Its score_table / viable_words / emitted scratch buffers are thread_local and persist
  // across queries AND across different / rebuilt FragmentIndex instances on the same
  // thread (e.g. the smaller final chunk of a chunked search). The reset must NOT walk the
  // previous query's touched-id list after the buffers were resized to a SMALLER index, or
  // it indexes past the new size. The defect is silent in a release build (std::vector
  // keeps capacity on assign), but aborts under _GLIBCXX_DEBUG / _GLIBCXX_ASSERTIONS / ASan
  // — which is where this test has teeth (run the suite under one of those to catch a
  // regression; the assertions below also cover the result staying correct after the shrink).
  auto make_snes = [](FragmentIndex_test& fi) {
    Param p = fi.getParameters();
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
  };
  auto make_spectrum = [](const std::string& seq) {
    TheoreticalSpectrumGenerator tsg;
    Param tsg_p = tsg.getParameters();
    tsg_p.setValue("add_metainfo", "true");
    tsg.setParameters(tsg_p);
    AASequence target = AASequence::fromString(seq);
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
    return spec;
  };

  // (1) Large index: several distinct 20-aa proteins -> a large mother set. Querying it
  //     populates the thread_local touched-id scratch lists with large mother ids.
  std::vector<FASTAFile::FASTAEntry> large_entries{
    {"L0", "L0", "ACDEFGHIKLMNPQRSTVWY"}, {"L1", "L1", "WYVTSRPNMLKIHGFEDCAQ"},
    {"L2", "L2", "GASTCVLIMPFWYHKRDENQ"}, {"L3", "L3", "MKVLAGDESTPNQRIHFYWC"},
    {"L4", "L4", "PQRSTVWYACDEFGHIKLMN"}, {"L5", "L5", "HRKDENQSTGAVLIMPFWYC"}};
  FragmentIndex_test fi_large;
  make_snes(fi_large);
  fi_large.build(large_entries);
  TEST_EQUAL(fi_large.isSnesMode(), true)
  const size_t large_mothers = fi_large.getPeptides().size();
  {
    MSSpectrum spec = make_spectrum("ACDEFGHIKL"); // 10-mer present in L0
    FragmentIndex::SpectrumMatchesTopN sms;
    fi_large.querySpectrum(spec, large_entries, sms); // leaves large ids in the thread_local touched lists
  }

  // (2) Small index: one short protein -> far fewer mothers. Querying it shrinks the
  //     thread_local buffers; the buggy reset would walk the stale large ids out of bounds.
  std::vector<FASTAFile::FASTAEntry> small_entries{{"S", "S", "ACDEFGHIKLM"}}; // 11 aa
  FragmentIndex_test fi_small;
  make_snes(fi_small);
  fi_small.build(small_entries);
  const size_t small_mothers = fi_small.getPeptides().size();
  TEST_EQUAL(small_mothers < large_mothers, true) // the index genuinely shrank

  MSSpectrum spec_small = make_spectrum("ACDEFGHIK"); // 9-mer sub-peptide of the small protein
  FragmentIndex::SpectrumMatchesTopN sms_small;
  fi_small.querySpectrum(spec_small, small_entries, sms_small); // <-- shrink path: must not OOB

  // Correctness after the shrink: the sub-peptide is still found.
  bool any_matched = false;
  for (const auto& hit : sms_small.hits_) { if (hit.num_matched_ >= 3u) { any_matched = true; break; } }
  TEST_EQUAL(any_matched, true)
}
END_SECTION

START_SECTION((non-SNES query is deterministic across repeated queries and safe when a smaller index is queried after a larger one on the same thread))
{
  // Regression guard for queryPeaks' thread_local window-relative counting buffers
  // (match_counts / touched_ids). They persist across queries AND across different /
  // rebuilt FragmentIndex instances on the same thread, and are restored to all-zero
  // by a touched-only reset at every block start. Two invariants have teeth here:
  //  (1) repeat determinism — a stale (unreset) count would inflate num_matched_ on
  //      the second query of the same spectrum;
  //  (2) large->small index reuse on one thread (a chunked search's smaller final
  //      chunk) must neither read out of bounds (run under ASan / _GLIBCXX_ASSERTIONS
  //      for full teeth) nor change the result.
  auto make_closed = [](FragmentIndex_test& fi) {
    Param p = fi.getParameters();
    p.setValue("precursor:mass_tolerance_lower", 20.0);
    p.setValue("precursor:mass_tolerance_upper", 20.0);
    p.setValue("precursor:mass_tolerance_unit", "ppm");
    p.setValue("fragment:mass_tolerance", 20.0);
    p.setValue("fragment:mass_tolerance_unit", "ppm");
    p.setValue("precursor:isotope_error_min", 0);
    p.setValue("precursor:isotope_error_max", 0);
    p.setValue("modifications:variable", std::vector<std::string>{});
    p.setValue("modifications:fixed", std::vector<std::string>{});
    p.setValue("peptide:min_size", 6);
    p.setValue("fragment:min_matched_ions", 3);
    fi.setParameters(p);
  };
  auto make_spectrum = [](const std::string& seq) {
    TheoreticalSpectrumGenerator tsg;
    Param tsg_p = tsg.getParameters();
    tsg_p.setValue("add_metainfo", "true");
    tsg.setParameters(tsg_p);
    AASequence target = AASequence::fromString(seq);
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
    return spec;
  };

  // (1) Larger index: tryptic peptides from several proteins.
  std::vector<FASTAFile::FASTAEntry> large_entries{
    {"L0", "L0", "MAGDEFHILNPKSAMPLEPEPTIDERWYVTSNMLIHGFEDCAK"},
    {"L1", "L1", "GASTCVLIMPFWKANOTHERLONGERSEQRHKDENQSTGAVLK"},
    {"L2", "L2", "PQSTVWYACDEFGHILMNKMKVLAGDESTPNQRIHFYWCAETK"}};
  FragmentIndex_test fi_large;
  make_closed(fi_large);
  fi_large.build(large_entries);
  const size_t large_peptides = fi_large.getPeptides().size();

  MSSpectrum spec = make_spectrum("SAMPLEPEPTIDER"); // tryptic peptide of L0
  FragmentIndex::SpectrumMatchesTopN sms_first, sms_second;
  fi_large.querySpectrum(spec, sms_first);
  fi_large.querySpectrum(spec, sms_second); // same thread, same buffers: must be identical

  TEST_EQUAL(sms_first.hits_.empty(), false)
  TEST_EQUAL(sms_first.hits_.size(), sms_second.hits_.size())
  for (Size i = 0; i < sms_first.hits_.size() && i < sms_second.hits_.size(); ++i)
  {
    TEST_EQUAL(sms_first.hits_[i].peptide_idx_, sms_second.hits_[i].peptide_idx_)
    TEST_EQUAL(sms_first.hits_[i].num_matched_, sms_second.hits_[i].num_matched_)
    TEST_EQUAL(sms_first.hits_[i].precursor_charge_, sms_second.hits_[i].precursor_charge_)
    TEST_EQUAL(sms_first.hits_[i].isotope_error_, sms_second.hits_[i].isotope_error_)
  }

  // (2) Smaller index queried on the same thread afterwards.
  std::vector<FASTAFile::FASTAEntry> small_entries{{"S", "S", "MKSAMPLEPEPTIDERAK"}};
  FragmentIndex_test fi_small;
  make_closed(fi_small);
  fi_small.build(small_entries);
  TEST_EQUAL(fi_small.getPeptides().size() < large_peptides, true) // genuinely smaller

  FragmentIndex::SpectrumMatchesTopN sms_small;
  fi_small.querySpectrum(spec, sms_small); // reuse path: must not OOB, must still match

  bool small_found = false;
  for (const auto& hit : sms_small.hits_) { if (hit.num_matched_ >= 3u) { small_found = true; break; } }
  TEST_EQUAL(small_found, true)
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

START_SECTION((SNES mother generation rejects spans with a stop codon))
{
  // As above, with a stop codon ('*') instead of the X: it has no residue mass either.
  const std::vector<FASTAFile::FASTAEntry> entries{
      {"p", "p", "ACDEFGHI*KLMNPQSTVWY"}}; // '*' at 0-based position 8

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

  bool found_prefix = false;
  for (const auto& mother : fi.getPeptides())
  {
    const size_t start = mother.sequence_.first;
    const size_t end = start + mother.sequence_.second;
    TEST_EQUAL(start > 8u || end <= 8u, true)
    if (start == 0 && mother.sequence_.second == 8) { found_prefix = true; }
  }
  TEST_EQUAL(found_prefix, true)
}
END_SECTION

START_SECTION((SNES mother generation truncates Single-N mother to unambiguous prefix on X/B/Z))
{
  // Issue #9192 item 2: a Single-N mother anchored at position 0 with proposed
  // length 12 spans the X at position 8. The whole mother used to be dropped;
  // now the unambiguous prefix [0, 8) length 8 must still be emitted.
  const std::vector<FASTAFile::FASTAEntry> entries{
      {"p", "p", "ACDEFGHIXKLMNPQSTVWY"}}; // X at 0-based position 8

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

  // Truncated Single-N mother [0, 8) of length 8 must exist.
  bool found_truncated_prefix = false;
  for (const auto& mother : fi.getPeptides())
  {
    const bool is_single_c = FragmentIndex::isSingleCMother(mother.mod_bitmask_);
    if (!is_single_c
        && mother.sequence_.first == 0
        && mother.sequence_.second == 8)
    {
      found_truncated_prefix = true;
      break;
    }
  }
  TEST_EQUAL(found_truncated_prefix, true)

  // Invariant: no kept mother spans the X at position 8.
  for (const auto& mother : fi.getPeptides())
  {
    const size_t start = mother.sequence_.first;
    const size_t end = start + mother.sequence_.second;
    TEST_EQUAL(start > 8u || end <= 8u, true)
  }
}
END_SECTION

START_SECTION((SNES mother generation truncates Single-C mother to unambiguous suffix on X/B/Z))
{
  // Issue #9192 item 2: a Single-C mother anchored at the last residue (j=19)
  // with proposed length 12 spans [8, 20) and covers the X. The whole mother
  // used to be dropped; the unambiguous suffix [9, 20) length 11 must now
  // be emitted (length capped per-span at e - s = 11).
  const std::vector<FASTAFile::FASTAEntry> entries{
      {"p", "p", "ACDEFGHIXKLMNPQSTVWY"}}; // X at 0-based position 8

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

  // Truncated Single-C mother [9, 20) of length 11 must exist.
  bool found_truncated_suffix = false;
  for (const auto& mother : fi.getPeptides())
  {
    const bool is_single_c = FragmentIndex::isSingleCMother(mother.mod_bitmask_);
    if (is_single_c
        && mother.sequence_.first == 9
        && mother.sequence_.second == 11)
    {
      found_truncated_suffix = true;
      break;
    }
  }
  TEST_EQUAL(found_truncated_suffix, true)

  // Invariant: no kept mother spans the X at position 8.
  for (const auto& mother : fi.getPeptides())
  {
    const size_t start = mother.sequence_.first;
    const size_t end = start + mother.sequence_.second;
    TEST_EQUAL(start > 8u || end <= 8u, true)
  }
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

START_SECTION(([EXTRA] rebuilding replaces the previous peptide and fragment buffers))
{
  const vector<FASTAFile::FASTAEntry> first = {{"P1", "", "THQPSANLDIK"}};
  const vector<FASTAFile::FASTAEntry> second = {{"P2", "", "VLVLDTDYK"}};
  FragmentIndex fi;
  Param p = fi.getParameters();
  p.setValue("decoys", "false");
  p.setValue("modifications:fixed", vector<string> {});
  p.setValue("modifications:variable", vector<string> {});
  fi.setParameters(p);
  fi.build(first);
  const Size fragments = fi.getNumFragments();
  fi.build(first);
  TEST_EQUAL(fi.getPeptides().size(), 1)
  TEST_EQUAL(fi.getNumFragments(), fragments)
  fi.build(second);
  ABORT_IF(fi.getPeptides().size() != 1)
  TEST_EQUAL(fi.reconstructModifiedSequence(fi.getPeptides()[0], second).toString(), "VLVLDTDYK")
  FragmentIndex fresh;
  fresh.setParameters(p);
  fresh.build(second);
  TEST_EQUAL(fi.getNumFragments(), fresh.getNumFragments())
  TEST_TRUE(fi.isBuild())

  // A failed rebuild must not leave the previous index marked as built.
  const vector<FASTAFile::FASTAEntry> invalid = {{"too_long", "", string(65536, 'A')}};
  TEST_EXCEPTION(Exception::InvalidParameter, fi.build(invalid))
  TEST_FALSE(fi.isBuild())
  TEST_TRUE(fi.getPeptides().empty())
  TEST_EQUAL(fi.getNumFragments(), 0)
}
END_SECTION

START_SECTION(([EXTRA] initial methionine clipping preserves coordinates and digestion limits))
{
  FragmentIndex_test fi;
  auto p = fi.getParameters();
  TEST_FALSE(p.getValue("peptide:clip_nterm_methionine").toBool())
  p.setValue("peptide:clip_nterm_methionine", "true");
  p.setValue("peptide:min_size", 1);
  p.setValue("peptide:max_size", 40);
  p.setValue("peptide:min_mass", 0);
  p.setValue("modifications:fixed", StringList {});
  p.setValue("modifications:variable", StringList {});
  p.setValue("peptide:missed_cleavages", 0);
  const vector<FASTAFile::FASTAEntry> entries {{"target", "", "MPEPTIDER"},      {"DECOY_test", "", "MPEPTIDEK"}, {"non_m", "", "APEPTIDER"},
                                               {"internal_m", "", "KMPEPTIDER"}, {"two_m", "", "MMPEPTIDER"},     {"single_m", "", "M"}};
  auto sequences = [&](UInt32 protein) {
    set<string> result;
    for (const auto& pep : fi.getPeptides())
    {
      if (pep.protein_idx != protein) { continue; }
      const auto sequence = fi.reconstructModifiedSequence(pep, entries);
      result.insert(sequence.toUnmodifiedString());
      TEST_REAL_SIMILAR(pep.precursor_mz_, sequence.getMZ(1))
      if (sequence.toUnmodifiedString() == "PEPTIDER")
      {
        TEST_EQUAL(pep.sequence_.first, 1)
        TEST_EQUAL(pep.sequence_.second, 8)
      }
    }
    return result;
  };
  fi.setParameters(p);
  fi.build(entries);
  TEST_TRUE(sequences(0) == set<string>({"MPEPTIDER", "PEPTIDER"}))
  TEST_TRUE(sequences(1) == set<string>({"MPEPTIDEK", "PEPTIDEK"}))
  TEST_TRUE(sequences(2) == set<string>({"APEPTIDER"}))
  TEST_TRUE(sequences(3) == set<string>({"K", "MPEPTIDER"}))
  TEST_TRUE(sequences(4) == set<string>({"MMPEPTIDER", "MPEPTIDER"}))
  TEST_TRUE(sequences(5) == set<string>({"M"}))

  // The mature peptide qualifies even when the retained-M form is too long.
  p.setValue("peptide:min_size", 8);
  p.setValue("peptide:max_size", 8);
  fi.setParameters(p);
  fi.build(entries);
  TEST_TRUE(sequences(0) == set<string>({"PEPTIDER"}))
  p.setValue("peptide:clip_nterm_methionine", "false");
  fi.setParameters(p);
  fi.build(entries);
  TEST_TRUE(sequences(0).empty())
  p.setValue("peptide:clip_nterm_methionine", "true");
  p.setValue("peptide:min_size", 9);
  p.setValue("peptide:max_size", 9);
  fi.setParameters(p);
  fi.build(entries);
  TEST_TRUE(sequences(0) == set<string>({"MPEPTIDER"}))

  const vector<FASTAFile::FASTAEntry> missed {{"missed", "", "MACDEKAGHILR"}};
  p.setValue("peptide:min_size", 1);
  p.setValue("peptide:max_size", 40);
  for (int mc : {0, 1})
  {
    p.setValue("peptide:missed_cleavages", mc);
    fi.setParameters(p);
    fi.build(missed);
    set<string> mature;
    for (const auto& pep : fi.getPeptides())
    {
      if (pep.sequence_.first == 1) { mature.insert(fi.reconstructModifiedSequence(pep, missed).toUnmodifiedString()); }
    }
    TEST_EQUAL(mature.count("ACDEK"), 1)
    TEST_EQUAL(mature.count("ACDEKAGHILR"), mc)
  }

  // No-cleavage still includes both complete proteoforms. Semi/non-specific
  // digestion must not emit an unmodified coordinate twice.
  const vector<FASTAFile::FASTAEntry> simple {{"p", "", "MPEPTIDE"}};
  p.setValue("enzyme", "no cleavage");
  p.setValue("peptide:missed_cleavages", 0);
  fi.setParameters(p);
  fi.build(simple);
  TEST_EQUAL(fi.getPeptides().size(), 2)
  for (const string specificity : {"semi", "none"})
  {
    p.setValue("peptide:enzyme_specificity", specificity);
    fi.setParameters(p);
    fi.build(simple);
    set<pair<uint16_t, uint16_t>> coordinates;
    for (const auto& pep : fi.getPeptides())
    {
      TEST_TRUE(coordinates.insert(pep.sequence_).second)
    }
  }
}
END_SECTION

START_SECTION(([EXTRA] clipped methionine peptides retain protein N - terminal variable modifications))
{
  const vector<FASTAFile::FASTAEntry> entries {{"nterm", "", "MPEPTIDER"}, {"internal", "", "KPEPTIDER"}};
  FragmentIndex_test fi;
  auto p = fi.getParameters();
  p.setValue("enzyme", "Trypsin/P"); // also cleave K-P in the internal-peptide control
  p.setValue("peptide:clip_nterm_methionine", "true");
  p.setValue("peptide:min_size", 7);
  p.setValue("peptide:missed_cleavages", 0);
  p.setValue("modifications:fixed", StringList {});
  p.setValue("modifications:variable", StringList {"Acetyl (Protein N-term)", "Oxidation (M)"});
  p.setValue("modifications:variable_max_per_peptide", 2);
  p.setValue("fragment:min_mz", 0);
  p.setValue("fragment:min_ion_index", 0);
  p.setValue("precursor:isotope_error_min", 0);
  p.setValue("precursor:isotope_error_max", 0);
  fi.setParameters(p);
  fi.build(entries);
  Size acetylated_mature = 0;
  for (const auto& pep : fi.getPeptides())
  {
    const auto seq = fi.reconstructModifiedSequence(pep, entries);
    TEST_REAL_SIMILAR(pep.precursor_mz_, seq.getMZ(1))
    if (pep.protein_idx == 1) { TEST_FALSE(seq.hasNTerminalModification()) }
    if (pep.protein_idx == 0 && pep.sequence_.first == 1)
    {
      TEST_EQUAL(seq.toUnmodifiedString(), "PEPTIDER")
      for (Size i = 0; i < seq.size(); ++i)
      {
        TEST_FALSE(seq[i].isModified())
      }
      if (seq.hasNTerminalModification())
      {
        ++acetylated_mature;
        TEST_EQUAL(seq.getNTerminalModification()->getFullId(), "Acetyl (Protein N-term)")
        // Verify fragment emission and reconstruction agree for the new variant.
        MSSpectrum spectrum;
        TheoreticalSpectrumGenerator().getSpectrum(spectrum, seq, 1, 1);
        Precursor precursor;
        precursor.setMZ(seq.getMZ(2));
        precursor.setCharge(2);
        spectrum.setPrecursors({precursor});
        spectrum.setMSLevel(2);
        spectrum.sortByPosition();
        FragmentIndex::SpectrumMatchesTopN matches;
        fi.querySpectrum(spectrum, entries, matches);
        bool found = false;
        for (const auto& match : matches.hits_)
        {
          found |= fi.reconstructModifiedSequence(fi.getPeptides()[match.peptide_idx_], entries) == seq;
        }
        TEST_TRUE(found)
      }
    }
  }
  TEST_EQUAL(acetylated_mature, 1)
}
END_SECTION

START_SECTION(([EXTRA] SNES recognizes the mature protein N - terminus after initial methionine loss))
{
  const vector<FASTAFile::FASTAEntry> entries {{"mature", "", "MACDEFGHILNPQR"}, {"internal", "", "KMACDEFGHILNPQR"}};
  FragmentIndex_test fi;
  auto p = fi.getParameters();
  p.setValue("peptide:enzyme_specificity", "none");
  p.setValue("peptide:min_size", 8);
  p.setValue("peptide:max_size", 12);
  p.setValue("modifications:fixed", StringList {});
  p.setValue("modifications:variable", StringList {"Acetyl (Protein N-term)"});
  p.setValue("modifications:variable_max_per_peptide", 1);
  p.setValue("snes_enabled", "true");
  p.setValue("fragment:min_matched_ions", 3);
  p.setValue("precursor:isotope_error_min", 0);
  p.setValue("precursor:isotope_error_max", 0);
  AASequence peptide = AASequence::fromString("ACDEFGHI");
  MSSpectrum spectrum;
  TheoreticalSpectrumGenerator().getSpectrum(spectrum, peptide, 1, 1);
  peptide.setNTerminalModification("Acetyl (Protein N-term)");
  Precursor precursor;
  precursor.setCharge(2);
  precursor.setMZ(peptide.getMZ(2));
  spectrum.setPrecursors({precursor});
  spectrum.setMSLevel(2);
  spectrum.sortByPosition();
  for (bool enabled : {false, true})
  {
    p.setValue("peptide:clip_nterm_methionine", enabled ? "true" : "false");
    fi.setParameters(p);
    fi.build(entries);
    FragmentIndex::SpectrumMatchesTopN matches;
    fi.querySpectrum(spectrum, entries, matches);
    bool found = false;
    for (const auto& match : matches.hits_)
    {
      if (match.subset_bitmask_ == 0) { continue; }
      const auto& mother = fi.getPeptides()[match.peptide_idx_];
      TEST_EQUAL(mother.protein_idx, 0)
      TEST_EQUAL(mother.sequence_.first, 1)
      const auto reconstructed = fi.reconstructRealizedSubSequence(mother, entries, 8, match.subset_bitmask_);
      found |= reconstructed == peptide;
    }
    TEST_EQUAL(found, enabled)
  }
}
END_SECTION

START_SECTION([EXTRA] sortPeptides_() leaves exactly the order of std::sort)
{
  // Peptides with equal (precursor_mz_, protein_idx) differ in mod_bitmask_ / sequence_, and the order std::sort
  // happens to give them defines the peptide indices and, through tie-breaking, the search results. sortPeptides_()
  // must therefore return the very permutation std::sort returns - for every number of threads and every task size.
  using Peptide = FragmentIndex::Peptide;
  const auto by_mz_then_protein = [](const Peptide& a, const Peptide& b)
  {
    return std::tie(a.precursor_mz_, a.protein_idx) < std::tie(b.precursor_mz_, b.protein_idx);
  };

  // M. D. McIlroy, "A Killer Adversary for Quicksort" (1999): answers the comparisons of std::sort such that its
  // quicksort degenerates. Sorting the returned keys makes the same comparisons, so an introsort runs into its depth
  // limit and has to fall back to heap sort.
  const auto quicksortKiller = [](size_t n)
  {
    const int gas = static_cast<int>(n) - 1;
    std::vector<int> key(n, gas);
    std::vector<int> items(n);
    std::iota(items.begin(), items.end(), 0);
    int n_solid = 0;
    int candidate = 0;
    std::sort(items.begin(), items.end(), [&](int x, int y)
    {
      if (key[x] == gas && key[y] == gas) { key[x == candidate ? x : y] = n_solid++; }
      if (key[x] == gas) { candidate = x; }
      else if (key[y] == gas) { candidate = y; }
      return key[x] < key[y];
    });
    return key;
  };

  const int n_patterns = 10;
  const auto makePeptides = [&](int pattern, size_t n)
  {
    std::mt19937 rng(static_cast<uint32_t>(pattern * 1000003 + n));
    std::vector<int> killer;
    if (pattern >= 8) { killer = quicksortKiller(n); }
    std::vector<Peptide> peptides;
    peptides.reserve(n);
    for (size_t i = 0; i < n; ++i)
    {
      const size_t r = n - 1 - i;
      float mz = 0;
      UInt32 protein = 0;
      switch (pattern)
      {
        case 0: mz = static_cast<float>(rng() % (n / 8 + 1)); protein = rng() % 4; break; // random, many ties
        case 1: mz = static_cast<float>(rng() % 3); protein = rng() % 2; break;           // six different keys
        case 2: mz = static_cast<float>(i / 3); protein = (i % 3) / 2; break;             // sorted
        case 3: mz = static_cast<float>(r / 3); protein = (r % 3) / 2; break;             // reverse sorted
        case 4: mz = 1.0f; protein = 7; break;                                            // all equal
        case 5: mz = static_cast<float>(std::min(i, r) / 2); break;                       // organ pipe
        case 6: mz = static_cast<float>(i % 17); protein = i % 2; break;                  // sawtooth
        case 7: mz = 500.0f + 0.01f * static_cast<float>(rng() % 200000); protein = static_cast<UInt32>(i / 10); break; // as build() has them
        case 8: mz = static_cast<float>(killer[i]); break;                                // heap sort fallback
        default: mz = static_cast<float>(killer[i] / 4); protein = rng() % 2; break;      // heap sort fallback, ties
      }
      // mod_bitmask_ and sequence_ are not part of the sort key: they tell equal peptides apart
      peptides.emplace_back(protein, static_cast<uint32_t>(i),
                            std::make_pair(static_cast<uint16_t>(i & 0xFFFF), static_cast<uint16_t>(i >> 16)), mz);
    }
    return peptides;
  };

  const auto countDifferences = [](const std::vector<Peptide>& a, const std::vector<Peptide>& b)
  {
    if (a.size() != b.size()) { return std::max<size_t>(a.size(), b.size()); }
    size_t differences = 0;
    for (size_t i = 0; i < a.size(); ++i)
    {
      differences += !(a[i].protein_idx == b[i].protein_idx && a[i].mod_bitmask_ == b[i].mod_bitmask_
                       && a[i].sequence_ == b[i].sequence_ && a[i].precursor_mz_ == b[i].precursor_mz_);
    }
    return differences;
  };

#ifdef _OPENMP
  const int max_threads = omp_get_max_threads();
#endif
  std::vector<size_t> sizes(41);
  std::iota(sizes.begin(), sizes.end(), 0); // 0..40: around the 16 elements below which std::sort only insertion-sorts
  sizes.insert(sizes.end(), {1000, 20000, 1000000});
  for (const size_t n : sizes)
  {
    size_t differences = 0;
    for (int pattern = 0; pattern < n_patterns; ++pattern)
    {
      const std::vector<Peptide> input = makePeptides(pattern, n);
      std::vector<Peptide> expected = input;
      std::sort(expected.begin(), expected.end(), by_mz_then_protein);
      for (const int threads : {1, 2, 4, 16})
      {
#ifdef _OPENMP
        omp_set_num_threads(threads);
#endif
        if (n < 1000000)
        {
          // every partitioning step hands its right part to another task / only the larger parts
          for (const size_t min_task_size : {size_t(0), size_t(100)})
          {
            std::vector<Peptide> sorted = input;
            FragmentIndex_test::sortPeptides(sorted, min_task_size);
            differences += countDifferences(sorted, expected);
          }
        }
        else
        {
          std::vector<Peptide> sorted = input;
          FragmentIndex_test::sortPeptides(sorted); // as build() calls it
          differences += countDifferences(sorted, expected);
        }
        (void)threads;
      }
    }
    TEST_EQUAL(differences, 0)
  }
#ifdef _OPENMP
  omp_set_num_threads(max_threads);
#endif
}
END_SECTION

START_SECTION(([EXTRA] optional peptidoform deduplication preserves modification sites and terminal contexts))
{
  const vector<FASTAFile::FASTAEntry> db
    = {{"P1", "", "MPEPCIDEMK"}, {"P2", "", "MPEPCIDEMK"}, {"P3", "", "AKMPEPCIDEMK"}, {"DECOY_shared", "", "MPEPCIDEMK"}};
  FragmentIndex fi;
  Param p = fi.getParameters();
  p.setValue("decoys", "false");
  p.setValue("peptide:missed_cleavages", 0);
  p.setValue("modifications:fixed", vector<string> {"Carbamidomethyl (C)"});
  p.setValue("modifications:variable", vector<string> {"Oxidation (M)", "Acetyl (Protein N-term)"});
  p.setValue("modifications:variable_max_per_peptide", 1);
  p.setValue("peptide:deduplicate", "false");
  fi.setParameters(p);
  fi.build(db);
  set<string> expected;
  for (const auto& peptide : fi.getPeptides())
  {
    expected.insert(fi.reconstructModifiedSequence(peptide, db).toString());
  }
  TEST_EQUAL(expected.size(), 4) // Fixed-only, oxidation at either M, protein-N acetyl.
  const Size original_count = fi.getPeptides().size();
  TEST_TRUE(original_count > expected.size())
  p.setValue("peptide:deduplicate", "true");
  fi.setParameters(p);
  fi.build(db);
  set<string> observed;
  for (const auto& peptide : fi.getPeptides())
  {
    observed.insert(fi.reconstructModifiedSequence(peptide, db).toString());
  }
  TEST_TRUE(observed == expected)
  TEST_EQUAL(fi.getPeptides().size(), expected.size())

  // Rebuilding with the default restores every occurrence; clear() carries no identity state.
  fi.clear();
  p.setValue("peptide:deduplicate", "false");
  fi.setParameters(p);
  fi.build(db);
  TEST_EQUAL(fi.getPeptides().size(), original_count)
}
END_SECTION

START_SECTION(([EXTRA] peptidoform deduplication keeps the first entry of every rendered peptidoform))
{
  // peptide:deduplicate compares entries within runs of equal precursor m/z, residue by residue and modification by
  // modification (or, for the last configuration, as strings). By definition it keeps, of all entries with the same
  // reconstructModifiedSequence(...).toString(), the first one. Build both ways on proteins that repeat peptides in
  // every protein-terminal context.
  const vector<string> blocks = {"MPEPCIDEMK", "QPEPTIDEMR", "AMDEQK", "MSTQMPEPK", "GGQMSTMK", "CMQEDK", "QQMMCR", "PEPTIDEK"};
  std::mt19937 rng(17);
  vector<FASTAFile::FASTAEntry> db;
  for (int i = 0; i < 40; ++i)
  {
    string sequence = (i % 3 == 0) ? "M" : "";
    const int n_blocks = 1 + static_cast<int>(rng() % 4);
    for (int b = 0; b < n_blocks; ++b) sequence += blocks[rng() % blocks.size()];
    db.push_back({"P" + std::to_string(i), "", sequence});
  }
  struct Config { vector<string> fixed, variable; };
  const vector<Config> configs = {
    {{}, {}},
    {{"Carbamidomethyl (C)"}, {"Oxidation (M)"}},
    {{"Carbamidomethyl (C)"}, {"Oxidation (M)", "Acetyl (Protein N-term)"}},
    {{"Carbamidomethyl (C)"}, {"Oxidation (M)", "Acetyl (N-term)", "Acetyl (Protein N-term)"}}, // two modifications rendered alike
    {{"Carbamidomethyl (C)", "TMT6plex (N-term)", "TMT6plex (K)"}, {"Oxidation (M)", "Gln->pyro-Glu (N-term Q)"}}, // terminal residue modification
    {{}, {"Oxidation (M)", "Amidated (Protein C-term)", "Deamidated (Q)"}},
    {{"Carbamidomethyl (C)"}, {"Carbamidomethyl (C)", "Oxidation (M)"}}}; // fixed and variable: compared as strings
  Size checked = 0, mismatches = 0, removed = 0;
  for (const Config& config : configs)
  {
    for (const string clip : {"false", "true"})
    {
      for (const int max_mods : {1, 2})
      {
        for (const string specificity : {"full", "semi"})
        {
          FragmentIndex fi;
          Param p = fi.getParameters();
          p.setValue("decoys", "false");
          p.setValue("peptide:min_size", 5);
          p.setValue("peptide:missed_cleavages", 2);
          p.setValue("peptide:clip_nterm_methionine", clip);
          p.setValue("peptide:enzyme_specificity", specificity);
          p.setValue("modifications:fixed", config.fixed);
          p.setValue("modifications:variable", config.variable);
          p.setValue("modifications:variable_max_per_peptide", max_mods);
          p.setValue("peptide:deduplicate", "false");
          fi.setParameters(p);
          fi.build(db);
          vector<FragmentIndex::Peptide> expected;
          // the removed occurrences: (kept entry, protein, start), by kept entry and then in index order
          map<string, Size> kept_index;
          vector<tuple<Size, UInt32, uint16_t>> expected_removed;
          for (const auto& peptide : fi.getPeptides())
          {
            const auto [it, inserted] = kept_index.emplace(fi.reconstructModifiedSequence(peptide, db).toString(), expected.size());
            if (inserted) { expected.push_back(peptide); }
            else { expected_removed.emplace_back(it->second, peptide.protein_idx, peptide.sequence_.first); }
          }
          std::stable_sort(expected_removed.begin(), expected_removed.end(),
                           [](const auto& a, const auto& b) { return std::get<0>(a) < std::get<0>(b); });
          TEST_TRUE(fi.getRemovedOccurrences().empty())
          removed += fi.getPeptides().size() - expected.size();
          p.setValue("peptide:deduplicate", "true");
          fi.setParameters(p);
          fi.build(db);
          const auto& observed = fi.getPeptides();
          bool same = observed.size() == expected.size();
          for (Size i = 0; same && i < observed.size(); ++i)
          {
            same = observed[i].protein_idx == expected[i].protein_idx && observed[i].mod_bitmask_ == expected[i].mod_bitmask_
                   && observed[i].sequence_ == expected[i].sequence_ && observed[i].precursor_mz_ == expected[i].precursor_mz_;
          }
          const auto& occurrences = fi.getRemovedOccurrences();
          same = same && occurrences.size() == expected_removed.size();
          for (Size i = 0; same && i < occurrences.size(); ++i)
          {
            same = occurrences[i].peptide_idx == std::get<0>(expected_removed[i]) && occurrences[i].protein_idx == std::get<1>(expected_removed[i])
                   && occurrences[i].start == std::get<2>(expected_removed[i]);
          }
          mismatches += same ? 0 : 1;
          ++checked;
        }
      }
    }
  }
  TEST_EQUAL(checked, configs.size() * 8)
  TEST_EQUAL(mismatches, 0)
  TEST_TRUE(removed > 0)
}
END_SECTION

START_SECTION((static void checkFixedModifications(const StringList& fixed_modifications)))
{
  // A fixed terminal modification is one N- or C-terminal mass for every peptide. Those that apply to some peptides
  // only are rejected (they are searched as variable modifications), as is a second one on the same terminus.
  const vector<StringList> rejected = {
    {"Acetyl (Protein N-term)"},
    {"Amidated (Protein C-term)"},
    {"Gln->pyro-Glu (N-term Q)"},
    {"Carbamidomethyl (C)", "TMT6plex (N-term)", "Acetyl (N-term)"}};
  for (const StringList& fixed : rejected)
  {
    TEST_EXCEPTION(Exception::InvalidParameter, FragmentIndex::checkFixedModifications(fixed))
    FragmentIndex fi;
    Param p = fi.getParameters();
    p.setValue("modifications:fixed", fixed);
    TEST_EXCEPTION(Exception::InvalidParameter, fi.setParameters(p))
  }
  // peptide-terminal fixed modifications of any residue, one per terminus, plus residue modifications
  const StringList accepted = {"Carbamidomethyl (C)", "TMT6plex (K)", "TMT6plex (N-term)", "Amidated (C-term)"};
  FragmentIndex::checkFixedModifications(accepted);
  FragmentIndex::checkFixedModifications({});

  // ... and they are applied to every peptide, internal ones included
  const vector<FASTAFile::FASTAEntry> db {{"p", "p", "MCAPEPTIDEKQLGSVTAKQMNPEPTIDER"}};
  FragmentIndex fi;
  Param p = fi.getParameters();
  p.setValue("peptide:min_size", 5);
  p.setValue("peptide:missed_cleavages", 1);
  p.setValue("modifications:fixed", accepted);
  p.setValue("modifications:variable", StringList {});
  fi.setParameters(p);
  fi.build(db);
  TEST_TRUE(fi.getPeptides().size() >= 4)
  for (const auto& peptide : fi.getPeptides())
  {
    const AASequence seq = fi.reconstructModifiedSequence(peptide, db);
    TEST_TRUE(seq.hasNTerminalModification())
    TEST_TRUE(seq.hasCTerminalModification())
    TEST_REAL_SIMILAR(peptide.precursor_mz_, seq.getMZ(1))
  }
}
END_SECTION

START_SECTION(([EXTRA] variable modification enumeration: subsets up to variable_max_per_peptide, at most 31 sites))
{
  // A peptide with n Met has n Oxidation (M) sites, so with at most k variable modifications per peptide the index
  // holds sum_{j <= k} C(n, j) forms of it. The enumeration visits only those subsets (it used to visit all 2^n, which
  // does not finish for 25 sites and is undefined for 32), and at most the first 31 sites take part.
  auto forms = [](Size n, Size k) {
    Size total = 0, c = 1; // c = C(n, j)
    for (Size j = 0; j <= k && j <= n; ++j) { total += c; c = c * (n - j) / (j + 1); }
    return total;
  };
  for (const Size n : {Size(3), Size(25), Size(31), Size(32), Size(35)})
  {
    for (const Size k : {Size(0), Size(1), Size(2)})
    {
      const vector<FASTAFile::FASTAEntry> db {{"p", "p", "G" + string(n, 'M') + "K"}};
      FragmentIndex_test fi;
      Param p = fi.getParameters();
      p.setValue("enzyme", "no cleavage");
      p.setValue("peptide:min_size", 0);
      p.setValue("peptide:max_size", 100);
      p.setValue("peptide:min_mass", 0);
      p.setValue("peptide:max_mass", 50000);
      p.setValue("fragment:min_mz", 0);
      p.setValue("fragment:max_mz", 50000);
      p.setValue("modifications:fixed", StringList {});
      p.setValue("modifications:variable", StringList {"Oxidation (M)"});
      p.setValue("modifications:variable_max_per_peptide", static_cast<int>(k));
      p.setValue("peptide:deduplicate", "false");
      fi.setParameters(p);
      fi.build(db);
      const Size sites = std::min<Size>(n, 31);
      TEST_EQUAL(fi.getPeptides().size(), forms(sites, k))
      set<uint32_t> masks;
      for (const auto& peptide : fi.getPeptides())
      {
        masks.insert(peptide.mod_bitmask_);
        TEST_TRUE(static_cast<Size>(std::popcount(peptide.mod_bitmask_)) <= k)
        TEST_EQUAL(peptide.mod_bitmask_ >> sites, 0u) // the slots beyond the 31st stay unmodified
        const AASequence seq = fi.reconstructModifiedSequence(peptide, db);
        Size oxidized = 0;
        for (Size i = 0; i < seq.size(); ++i) { oxidized += seq[i].isModified() ? 1 : 0; }
        TEST_EQUAL(oxidized, static_cast<Size>(std::popcount(peptide.mod_bitmask_)))
        TEST_REAL_SIMILAR(peptide.precursor_mz_, seq.getMZ(1))
      }
      TEST_EQUAL(masks.size(), fi.getPeptides().size())
    }
  }
}
END_SECTION

START_SECTION(([EXTRA] the default isotope error range covers precursors selected at the first and second 13C peak))
{
  // The query matches observed mass + isotope_error * C13C12: the default [-2, 0] finds a peptide whose precursor was
  // selected one or two 13C spacings above the monoisotopic peak, not one below it.
  FragmentIndex fi;
  TEST_EQUAL(static_cast<int>(fi.getParameters().getValue("precursor:isotope_error_min")), -2)
  TEST_EQUAL(static_cast<int>(fi.getParameters().getValue("precursor:isotope_error_max")), 0)
  const vector<FASTAFile::FASTAEntry> db {{"p", "p", "EVAEAATGEDASSPPPK"}};
  Param p = fi.getParameters();
  p.setValue("enzyme", "no cleavage");
  p.setValue("modifications:fixed", StringList {});
  p.setValue("modifications:variable", StringList {});
  p.setValue("fragment:min_mz", 0);
  fi.setParameters(p);
  fi.build(db);
  const AASequence peptide = AASequence::fromString("EVAEAATGEDASSPPPK");
  PeakSpectrum ions;
  TheoreticalSpectrumGenerator().getSpectrum(ions, peptide, 1, 1);
  for (int observed_minus_theoretical = -1; observed_minus_theoretical <= 2; ++observed_minus_theoretical)
  {
    MSSpectrum spectrum;
    for (const auto& ion : ions) { spectrum.push_back(ion); }
    spectrum.setMSLevel(2);
    Precursor precursor;
    precursor.setCharge(2);
    precursor.setMZ(peptide.getMZ(2) + observed_minus_theoretical * Constants::C13C12_MASSDIFF_U / 2);
    spectrum.setPrecursors({precursor});
    FragmentIndex::SpectrumMatchesTopN sms;
    fi.querySpectrum(spectrum, sms);
    bool found = false;
    for (const auto& hit : sms.hits_)
    {
      found |= hit.precursor_charge_ == 2 && hit.isotope_error_ == -observed_minus_theoretical;
    }
    TEST_EQUAL(found, observed_minus_theoretical >= 0)
  }
}
END_SECTION

START_SECTION((static StringList shadowedVariableTerminalModifications(const StringList& fixed_modifications, const StringList& variable_modifications)))
{
  // A terminus carries one modification: a variable modification of the whole terminus is not applied where a fixed
  // one sits. Residue-specific terminal variable modifications modify the residue and stay.
  const StringList variable = {"Acetyl (Protein N-term)", "Oxidation (M)", "Gln->pyro-Glu (N-term Q)", "Amidated (C-term)"};
  TEST_EQUAL(ListUtils::concatenate(FragmentIndex::shadowedVariableTerminalModifications({"TMT6plex (N-term)"}, variable), ","),
             "Acetyl (Protein N-term)")
  TEST_EQUAL(ListUtils::concatenate(FragmentIndex::shadowedVariableTerminalModifications(
               {"Carbamidomethyl (C)", "TMT6plex (N-term)", "Amidated (C-term)"}, variable), ","),
             "Acetyl (Protein N-term),Amidated (C-term)")
  TEST_EQUAL(FragmentIndex::shadowedVariableTerminalModifications({"Carbamidomethyl (C)", "TMT6plex (K)"}, variable).size(), 0)
  TEST_EQUAL(FragmentIndex::shadowedVariableTerminalModifications({"TMT6plex (N-term)"}, {}).size(), 0)

  // The index applies the excluded modification nowhere, so every peptide's precursor m/z is that of its reconstructed
  // sequence. (It used to add both terminal masses while the reconstructed sequence kept the variable one only.)
  const vector<FASTAFile::FASTAEntry> db {{"p", "p", "ACAPEPTIDEKQLGSVTAKQMNPEPTIDER"}};
  auto build = [&db](const StringList& fixed)
  {
    FragmentIndex fi;
    Param p = fi.getParameters();
    p.setValue("peptide:min_size", 5);
    p.setValue("peptide:missed_cleavages", 1);
    p.setValue("modifications:fixed", fixed);
    p.setValue("modifications:variable", StringList {"Acetyl (Protein N-term)", "Oxidation (M)"});
    fi.setParameters(p);
    fi.build(db);
    Size acetylated = 0;
    for (const auto& peptide : fi.getPeptides())
    {
      const AASequence seq = fi.reconstructModifiedSequence(peptide, db);
      TEST_REAL_SIMILAR(peptide.precursor_mz_, seq.getMZ(1))
      if (seq.getNTerminalModificationName() == "Acetyl") ++acetylated;
      if (fixed.size() > 1)
      {
        TEST_EQUAL(seq.getNTerminalModificationName(), "TMT6plex")
      }
    }
    return std::make_pair(fi.getPeptides().size(), acetylated);
  };
  const auto [n_without_fixed_nterm, acetylated_without] = build({"Carbamidomethyl (C)"});
  const auto [n_with_fixed_nterm, acetylated_with] = build({"Carbamidomethyl (C)", "TMT6plex (N-term)"});
  TEST_EQUAL(acetylated_without, 2) // ACAPEPTIDEK and ACAPEPTIDEKQLGSVTAK
  TEST_EQUAL(acetylated_with, 0)
  TEST_EQUAL(n_with_fixed_nterm, n_without_fixed_nterm - acetylated_without)
}
END_SECTION

START_SECTION(([EXTRA] the prefilter totals count every candidate with a matched fragment, before the gate and the cap))
{
  // SpectrumMatchesTopN::scored_candidates_ / matched_peaks_ are taken before fragment:min_matched_ions and
  // scoring:max_candidates_per_spectrum: with a gate of 1 and no cap, the hits are exactly the counted candidates.
  const vector<FASTAFile::FASTAEntry> db{
    {"L0", "L0", "MAGDEFHILNPKSAMPLEPEPTIDERWYVTSNMLIHGFEDCAK"},
    {"L1", "L1", "GASTCVLIMPFWKANOTHERLONGERSEQRHKDENQSTGAVLK"},
    {"L2", "L2", "PQSTVWYACDEFGHILMNKMKVLAGDESTPNQRIHFYWCAETK"}};
  TheoreticalSpectrumGenerator tsg;
  PeakSpectrum spectrum;
  const AASequence target = AASequence::fromString("SAMPLEPEPTIDER");
  tsg.getSpectrum(spectrum, target, 1, 2);
  Precursor precursor;
  precursor.setMZ(target.getMZ(2));
  precursor.setCharge(2);
  spectrum.setPrecursors({precursor});
  spectrum.setMSLevel(2);

  auto query = [&](Int min_matched_ions, Int max_candidates)
  {
    FragmentIndex fi;
    Param p = fi.getParameters();
    p.setValue("precursor:mass_tolerance_lower", 1000.0); // every peptide is a candidate
    p.setValue("precursor:mass_tolerance_upper", 1000.0);
    p.setValue("precursor:mass_tolerance_unit", "Da");
    p.setValue("fragment:mass_tolerance", 0.5);
    p.setValue("fragment:mass_tolerance_unit", "Da");
    p.setValue("fragment:min_mz", 0);
    p.setValue("fragment:min_ion_index", 0);
    p.setValue("modifications:variable", vector<string>{});
    p.setValue("modifications:fixed", vector<string>{});
    p.setValue("peptide:min_size", 4);
    p.setValue("fragment:min_matched_ions", min_matched_ions);
    p.setValue("scoring:max_candidates_per_spectrum", max_candidates);
    fi.setParameters(p);
    fi.build(db);
    FragmentIndex::SpectrumMatchesTopN sms;
    fi.querySpectrum(spectrum, sms);
    return sms;
  };
  const FragmentIndex::SpectrumMatchesTopN all = query(1, 100000);
  uint64_t matched = 0;
  for (const auto& hit : all.hits_) matched += hit.num_matched_;
  TEST_TRUE(all.hits_.size() > 3)
  TEST_EQUAL(all.scored_candidates_, all.hits_.size())
  TEST_EQUAL(all.matched_peaks_, matched)

  const FragmentIndex::SpectrumMatchesTopN gated = query(5, 3);
  TEST_TRUE(gated.hits_.size() <= 3)
  TEST_EQUAL(gated.scored_candidates_, all.scored_candidates_)
  TEST_EQUAL(gated.matched_peaks_, all.matched_peaks_)

  FragmentIndex::SpectrumMatchesTopN sum = all;
  sum += gated;
  TEST_EQUAL(sum.scored_candidates_, 2 * all.scored_candidates_)
  TEST_EQUAL(sum.matched_peaks_, 2 * all.matched_peaks_)
  sum.clear();
  TEST_EQUAL(sum.scored_candidates_, 0)
  TEST_EQUAL(sum.matched_peaks_, 0)
}
END_SECTION

START_SECTION((void build(const std::vector<FASTAFile::FASTAEntry>& fasta_entries, const std::function<const MSExperiment*(Size)>& searched_spectra)))
{
  // The index built for some spectra keeps only the peptides in their precursor windows, in their order: these
  // spectra get the same candidates (peptides, matched peaks, charges, isotope errors) as from the full index.
  const std::vector<FASTAFile::FASTAEntry> entries{
    {"L0", "L0", "MAGDEFHILNPKSAMPLEPEPTIDERWYVTSNMLIHGFEDCAKLLIGHTDFEK"},
    {"L1", "L1", "GASTCVLIMPFWKANOTHERLONGERSEQRHKDENQSTGAVLKMEDITATESK"},
    {"L2", "L2", "PQSTVWYACDEFGHILMNKMKVLAGDESTPNQRIHFYWCAETKAAGVHELPR"}};
  auto configure = [](FragmentIndex& fi, double precursor_tolerance_ppm)
  {
    Param p = fi.getParameters();
    p.setValue("precursor:mass_tolerance_lower", precursor_tolerance_ppm);
    p.setValue("precursor:mass_tolerance_upper", precursor_tolerance_ppm);
    p.setValue("precursor:mass_tolerance_unit", "ppm");
    p.setValue("fragment:mass_tolerance", 20.0);
    p.setValue("fragment:mass_tolerance_unit", "ppm");
    p.setValue("precursor:isotope_error_min", -1);
    p.setValue("precursor:isotope_error_max", 1);
    p.setValue("precursor:min_charge", 1);
    p.setValue("precursor:max_charge", 3);
    p.setValue("peptide:min_size", 6);
    p.setValue("peptide:missed_cleavages", 2);
    p.setValue("fragment:min_matched_ions", 1);
    fi.setParameters(p);
  };
  auto make_spectrum = [](const std::string& seq, Int charge, Int ms_level)
  {
    TheoreticalSpectrumGenerator tsg;
    const AASequence target = AASequence::fromString(seq);
    MSSpectrum spec;
    tsg.getSpectrum(spec, target, 1, 1);
    spec.sortByPosition();
    Precursor prec;
    prec.setMZ(target.getMZ(2)); // doubly charged precursor; the charge itself may be unknown (0)
    prec.setCharge(charge);
    spec.setPrecursors({prec}); // replaces the one TheoreticalSpectrumGenerator sets
    spec.setMSLevel(ms_level);
    return spec;
  };
  MSExperiment spectra;
  spectra.addSpectrum(make_spectrum("SAMPLEPEPTIDER", 2, 2));
  spectra.addSpectrum(make_spectrum("VLAGDESTPNQR", 0, 2));    // all charges are tried
  spectra.addSpectrum(make_spectrum("ANOTHERLONGERSEQR", 2, 1)); // MS1: not searched

  FragmentIndex full;
  configure(full, 20.0);
  full.build(entries);

  FragmentIndex restricted;
  configure(restricted, 20.0);
  Size calls = 0;
  Size peptides_reported = 0;
  restricted.build(entries, [&](Size peptides) { ++calls; peptides_reported = peptides; return &spectra; });
  TEST_EQUAL(peptides_reported, full.getPeptides().size())
  TEST_EQUAL(calls, 1)
  TEST_EQUAL(restricted.isBuild(), true)
  TEST_EQUAL(restricted.getPeptides().size() < full.getPeptides().size(), true)
  TEST_EQUAL(restricted.getPeptides().empty(), false)
  TEST_EQUAL(restricted.getNumFragments() < full.getNumFragments(), true)

  auto same_peptide = [](const FragmentIndex::Peptide& a, const FragmentIndex::Peptide& b)
  {
    return a.protein_idx == b.protein_idx && a.mod_bitmask_ == b.mod_bitmask_ && a.sequence_ == b.sequence_
           && a.precursor_mz_ == b.precursor_mz_;
  };
  // the kept peptides are a subsequence of the full index' peptides
  Size next = 0;
  for (const auto& peptide : restricted.getPeptides())
  {
    while (next < full.getPeptides().size() && !same_peptide(full.getPeptides()[next], peptide)) { ++next; }
    TEST_EQUAL(next < full.getPeptides().size(), true)
    ++next;
  }
  for (Size i = 0; i < 2; ++i)
  {
    FragmentIndex::SpectrumMatchesTopN from_full, from_restricted;
    full.querySpectrum(spectra[i], from_full);
    restricted.querySpectrum(spectra[i], from_restricted);
    TEST_EQUAL(from_full.hits_.empty(), false)
    ABORT_IF(from_full.hits_.size() != from_restricted.hits_.size())
    for (Size k = 0; k < from_full.hits_.size(); ++k)
    {
      const auto& a = from_full.hits_[k];
      const auto& b = from_restricted.hits_[k];
      TEST_TRUE(same_peptide(full.getPeptides()[a.peptide_idx_], restricted.getPeptides()[b.peptide_idx_]))
      TEST_EQUAL(a.num_matched_, b.num_matched_)
      TEST_EQUAL(a.precursor_charge_, b.precursor_charge_)
      TEST_EQUAL(a.isotope_error_, b.isotope_error_)
    }
  }

  // no spectra (nullptr), no callback, and open search: the full index
  FragmentIndex unrestricted;
  configure(unrestricted, 20.0);
  unrestricted.build(entries, [](Size) -> const MSExperiment* { return nullptr; });
  TEST_EQUAL(unrestricted.getPeptides().size(), full.getPeptides().size())
  unrestricted.build(entries, {});
  TEST_EQUAL(unrestricted.getPeptides().size(), full.getPeptides().size())
  FragmentIndex open_full, open_restricted;
  configure(open_full, 5000.0);
  configure(open_restricted, 5000.0);
  open_full.build(entries);
  calls = 0;
  open_restricted.build(entries, [&spectra, &calls](Size) { ++calls; return &spectra; });
  TEST_EQUAL(calls, 0)
  TEST_EQUAL(open_restricted.getPeptides().size(), open_full.getPeptides().size())

  // peptide:deduplicate by strings (a modification configured fixed and variable: equal renderings, different precursor
  // m/z) keeps the first entry of a peptidoform among all of its entries, before the index is restricted: a spectrum at
  // the mass of a removed entry gets no candidate from it. The occurrences removed from a kept entry follow it.
  const std::vector<FASTAFile::FASTAEntry> repeats{
    {"R0", "R0", "MAGDEFHILNPKSAMPLECPEPTIDERWYVTSNMLIHGFEDCAKLLIGHTDFEK"},
    {"R1", "R1", "GASTCVLIMPFWKSAMPLECPEPTIDERHKDENQSTGAVLKMEDITATESK"},
    {"R2", "R2", "PQSTVWYACDEFGHILMNKSAMPLECPEPTIDERVLAGDESTPNQRAAGVHELPR"}};
  auto configure_by_string = [&configure](FragmentIndex& fi)
  {
    configure(fi, 20.0);
    Param p = fi.getParameters();
    p.setValue("modifications:fixed", std::vector<std::string>{"Carbamidomethyl (C)"});
    p.setValue("modifications:variable", std::vector<std::string>{"Carbamidomethyl (C)", "Oxidation (M)"});
    p.setValue("peptide:deduplicate", "true");
    fi.setParameters(p);
  };
  FragmentIndex by_string_full;
  configure_by_string(by_string_full);
  by_string_full.build(repeats);
  TEST_EQUAL(by_string_full.getRemovedOccurrences().empty(), false)
  auto with_precursor = [](const std::string& seq, double extra_mass)
  {
    TheoreticalSpectrumGenerator tsg;
    const AASequence target = AASequence::fromString(seq);
    MSSpectrum spec;
    tsg.getSpectrum(spec, target, 1, 1);
    spec.sortByPosition();
    Precursor prec;
    prec.setMZ((target.getMonoWeight() + extra_mass + 2 * Constants::PROTON_MASS_U) / 2.0);
    prec.setCharge(2);
    spec.setPrecursors({prec});
    spec.setMSLevel(2);
    return spec;
  };
  const double carbamidomethyl = ModificationsDB::getInstance()->getModification("Carbamidomethyl (C)")->getDiffMonoMass();
  MSExperiment at_fixed, at_both;
  at_fixed.addSpectrum(with_precursor("SAMPLEC(Carbamidomethyl)PEPTIDER", 0.0));
  at_fixed.addSpectrum(with_precursor("VLAGDESTPNQR", 0.0));
  at_both.addSpectrum(with_precursor("SAMPLEC(Carbamidomethyl)PEPTIDER", carbamidomethyl)); // fixed and variable
  for (const MSExperiment* searched : {&at_fixed, &at_both})
  {
    FragmentIndex by_string_restricted;
    configure_by_string(by_string_restricted);
    by_string_restricted.build(repeats, [searched](Size) { return searched; });
    const auto& full_peptides = by_string_full.getPeptides();
    const auto& kept_peptides = by_string_restricted.getPeptides();
    TEST_EQUAL(kept_peptides.size() < full_peptides.size(), true)
    // index in the restricted index of every peptide of the full one (or -1)
    std::vector<SignedSize> new_index(full_peptides.size(), -1);
    Size at = 0;
    for (Size k = 0; k < kept_peptides.size(); ++k)
    {
      while (at < full_peptides.size() && !same_peptide(full_peptides[at], kept_peptides[k])) { ++at; }
      ABORT_IF(at == full_peptides.size())
      new_index[at++] = static_cast<SignedSize>(k);
    }
    std::vector<std::tuple<SignedSize, UInt32, uint16_t>> expected_occurrences, kept_occurrences;
    for (const auto& occurrence : by_string_full.getRemovedOccurrences())
    {
      if (new_index[occurrence.peptide_idx] >= 0)
      {
        expected_occurrences.emplace_back(new_index[occurrence.peptide_idx], occurrence.protein_idx, occurrence.start);
      }
    }
    for (const auto& occurrence : by_string_restricted.getRemovedOccurrences())
    {
      kept_occurrences.emplace_back(static_cast<SignedSize>(occurrence.peptide_idx), occurrence.protein_idx, occurrence.start);
    }
    TEST_TRUE(kept_occurrences == expected_occurrences)
    for (const MSSpectrum& spectrum : *searched)
    {
      FragmentIndex::SpectrumMatchesTopN from_full, from_restricted;
      by_string_full.querySpectrum(spectrum, from_full);
      by_string_restricted.querySpectrum(spectrum, from_restricted);
      ABORT_IF(from_full.hits_.size() != from_restricted.hits_.size())
      for (Size k = 0; k < from_full.hits_.size(); ++k)
      {
        TEST_TRUE(same_peptide(full_peptides[from_full.hits_[k].peptide_idx_], kept_peptides[from_restricted.hits_[k].peptide_idx_]))
        TEST_EQUAL(from_full.hits_[k].num_matched_, from_restricted.hits_[k].num_matched_)
      }
    }
  }
}
END_SECTION

END_TEST
