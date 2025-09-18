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
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <limits>


/*
  FragmentIndex tests

  This suite verifies:
  - build(): digestion and peptide generation across enzyme, length/mass limits, missed cleavages, and modifications; asserts ordering invariants for peptides/fragments.
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
class FragmentIndex_test: public FragmentIndex
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
      if (!found) return false;
    }
    return true;
  }
  // Checks non-decreasing order of precursor_mz_ across all peptides (invariant of build()).
  bool peptidesSorted()
  {
    float last_mz = std::numeric_limits<float>::lowest();
    for (const auto& pep : fi_peptides_)
    {
      if (pep.precursor_mz_ >= last_mz)
      {
        last_mz = pep.precursor_mz_;
      }
      else
      {
        return false;
      }
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

  bool testQuery(const UInt32 charge, const bool precursor_mz_known, const std::vector<FASTAFile::FASTAEntry>& entries)
  {
    // fetch parameters for modification generation
    auto params = getParameters();
    const StringList modifications_fixed_ = ListUtils::toStringList<std::string>(params.getValue("modifications:fixed"));
    const StringList modifications_variable_ = ListUtils::toStringList<std::string>(params.getValue("modifications:variable"));
    const ModifiedPeptideGenerator::MapToResidueType fixed_modifications = ModifiedPeptideGenerator::getModifications(modifications_fixed_);
    const ModifiedPeptideGenerator::MapToResidueType variable_modifications = ModifiedPeptideGenerator::getModifications(modifications_variable_);

    //Create theoretical spectra for different charges
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
        ModifiedPeptideGenerator::applyVariableModifications(variable_modifications, mod_peptide, params.getValue("modifications:variable_max_per_peptide"), mod_peptides);
        mod_peptide = mod_peptides[pep.modification_idx_];
        tsg.getSpectrum(b_y_ions, mod_peptide, charge, charge);
        prec_theo.setMZ(mod_peptide.getMZ(charge));
        if (precursor_mz_known)
        {
          prec_theo.setCharge(charge);
        }
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
          if ((s.peptide_idx_ == peptide_idx) && (s.precursor_charge_ ==  charge))
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
// Test proteins used to generate expected peptides for multiple parameterizations
// Format of expected peptide descriptors below:
//   {protein_idx, modification_idx, {start, length}, dummy precursor_mz}
const std::vector<FASTAFile::FASTAEntry> entries0{{"t", "t", "ARGEPADSSRKDFDMDMDM"},
                                             {"t2", "t2", "HALLORTSCHSM"}};
// Expected peptides when enabling fixed Carbamidomethyl (C) and variable Oxidation (M)
std::vector<FragmentIndex::Peptide> peptides_we_should_hit_mod{{0,0,{2,8},5},
                                                            {0,0,{11,8},5},
                                                            {0,1,{11,8},5},
                                                            {0,2,{11,8},5},
                                                            {0,3,{11,8},5},
                                                            {0,4,{11,8},5},
                                                            {0,5,{11,8},5},
                                                            {0,6,{11,8},5},
                                                            {1,0 ,{0,6}, 5},
                                                            {1, 0,{6, 6}, 5},
                                                            {1, 1,{6, 6}, 5}

};
// Expected peptides without min/max size constraints (no missed cleavages, no modifications)
std::vector<FragmentIndex::Peptide> peptides_unmod_no_minmax{{0,0,{0,2},5},
                                                              {0,0,{2,8},5},
                                                              {0,0,{10,1},5},
                                                              {0,0,{11,8},5},
                                                              {1,0,{0,6},5},
                                                              {1,0,{6,6},5}};

// Expected peptides with size in [min_size, max_size] only
std::vector<FragmentIndex::Peptide> peptides_unmod_minmax{{0,0,{0,2},5},
                                                           {1,0,{0,6},5},
                                                           {1,0,{6,6},5}};
// Expected peptides with one missed cleavage allowed
std::vector<FragmentIndex::Peptide> peptides_unmod_minmax_missed_cleavage{{0,0,{0,2},5},
                                                                           {0,0,{2,8},5},
                                                                           {0,0,{11,8},5},
                                                                           {0,0,{0,10},5},
                                                                           {0,0,{2,9},5},
                                                                           {0,0,{10,9},5},
                                                                           {1,0,{0,6},5},
                                                                           {1,0,{6,6},5},
                                                                           {1,0,{0,12},5}};


FragmentIndex_test buildTest;
auto params = buildTest.getParameters();
params.setValue("enzyme", "Trypsin");
params.setValue("peptide:missed_cleavages", 0);
params.setValue("peptide:min_mass", 0);
params.setValue("peptide:min_size", 0);
params.setValue("peptide:max_mass", 5000);
params.setValue("modifications:variable", std::vector<std::string>{});
params.setValue("modifications:fixed", std::vector<std::string>{});
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
params.setValue("modifications:variable", std::vector<std::string>{"Oxidation (M)"});
params.setValue("modifications:fixed", std::vector<std::string>{"Carbamidomethyl (C)"});
buildTest.setParameters(params);
buildTest.build(entries0);
TEST_TRUE(buildTest.testDigestion(peptides_we_should_hit_mod))
TEST_TRUE(buildTest.peptidesSorted())
TEST_TRUE(buildTest.fragmentsSorted())

END_SECTION

// Verify that clear() resets the internal peptide container.
START_SECTION(clear())
const std::vector<FASTAFile::FASTAEntry> entries0{{"t", "t", "ARGEPADSSRKDFDMDMDM"},
                                             {"t2", "t2", "HALLORTSCHS"}};
FragmentIndex clearTest;
clearTest.build(entries0);
clearTest.clear();

TEST_TRUE(clearTest.getPeptides().empty())

END_SECTION



////TEST Different Charges of the query Spectrum ////
 // For each charge (1..4), a peptide's own theoretical spectrum should self-hit,
 // with and without explicitly setting the precursor charge.
START_SECTION(void querySpectrum(const MSSpectrum& spectrum, SpectrumMatchesTopN& sms))



const std::vector<FASTAFile::FASTAEntry> entries{{"test1", "test1","MSDEREVAEAATGEDASSPPPKTEAASDPQHPAASEGAAAAAASPPLLRCLVLTGFGGYDKVKLQSRPAAPPAPGPGQLTLRLRACGLNFADLMARQGLYDRLPPLPVTPGMEGAGVVIAVGEGVSDRKAGDRVMVLNRSGMWQEEVTVPSVQTFLIPEAMTFEEAAALLVNYITAYMVLFDFGNLQPGHSVLVHMAAGGVGMAAVQLCRTVENVTVFGTASASKHEALKENGVTHPIDYHTTDYVDEIKKISPKGVDIVMDPLGGSDTAKGYNLLKPMGKVVTYGMANLLTGPKRNLMALARTWWNQFSVTALQLLQANRAVCGFHLGYLDGEVELVSGVVARLLALYNQGHIKPHIDSVWPFEKVADAMKQMQEKKNVGKVLLVPGPEKEN"}};

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

END_SECTION

// Shift the precursor by integer isotope errors [-3..3] and expect stable peptide window mapping.
START_SECTION(isotope_error)
const std::vector<FASTAFile::FASTAEntry> entries{{"test1", "test1","MSDEREVAEAATGEDASSPPPKTEAASDPQHPAASEGAAAAAASPPLLRCLVLTGFGGYDKVKLQSRPAAPPAPGPGQLTLRLRACGLNFADLMARQGLYDRLPPLPVTPGMEGAGVVIAVGEGVSDRKAGDRVMVLNRSGMWQEEVTVPSVQTFLIPEAMTFEEAAALLVNYITAYMVLFDFGNLQPGHSVLVHMAAGGVGMAAVQLCRTVENVTVFGTASASKHEALKENGVTHPIDYHTTDYVDEIKKISPKGVDIVMDPLGGSDTAKGYNLLKPMGKVVTYGMANLLTGPKRNLMALARTWWNQFSVTALQLLQANRAVCGFHLGYLDGEVELVSGVVARLLALYNQGHIKPHIDSVWPFEKVADAMKQMQEKKNVGKVLLVPGPEKEN"}};

FragmentIndex_test isoTest;

// Configure parameters before building the index (isotope error and fragment m/z bounds)
auto params = isoTest.getParameters();
params.setValue("min_isotope_error", -3);
params.setValue("max_isotope_error", 3);
params.setValue("fragment:min_mz", 0);
params.setValue("fragment:max_mz", 90000);
params.setValue("modifications:variable", std::vector<std::string>{});
params.setValue("modifications:fixed", std::vector<std::string>{});
isoTest.setParameters(params);

// build after parameterization
isoTest.build(entries);

TheoreticalSpectrumGenerator tsg;
PeakSpectrum b_y_ions;
AASequence peptide = AASequence::fromString("EVAEAATGEDASSPPPK");
tsg.getSpectrum(b_y_ions, peptide,1,1);
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
  for (const auto& hit : sms.hits_){
      auto result = isoTest.getPeptides()[hit.peptide_idx_];
      auto psize = peptide.size();
      TEST_EQUAL(result.sequence_.first, 5)
      TEST_EQUAL(result.sequence_.second, psize)
      found = true;
  }
  TEST_TRUE(found);
}

END_SECTION

// Apply small deterministic fragment m/z jitter and a precursor offset within tolerances;
// expect the correct peptide hit and zero isotope error.
START_SECTION(tolerance)
const std::vector<FASTAFile::FASTAEntry> entries{{"test1", "test1","MSDEREVAEAATGEDASSPPPKTEAASDPQHPAASEGAAAAAASPPLLRCLVLTGFGGYDKVKLQSRPAAPPAPGPGQLTLRLRACGLNFADLMARQGLYDRLPPLPVTPGMEGAGVVIAVGEGVSDRKAGDRVMVLNRSGMWQEEVTVPSVQTFLIPEAMTFEEAAALLVNYITAYMVLFDFGNLQPGHSVLVHMAAGGVGMAAVQLCRTVENVTVFGTASASKHEALKENGVTHPIDYHTTDYVDEIKKISPKGVDIVMDPLGGSDTAKGYNLLKPMGKVVTYGMANLLTGPKRNLMALARTWWNQFSVTALQLLQANRAVCGFHLGYLDGEVELVSGVVARLLALYNQGHIKPHIDSVWPFEKVADAMKQMQEKKNVGKVLLVPGPEKEN"}};

FragmentIndex_test tolTest;


auto params = tolTest.getParameters();
params.setValue("fragment:min_mz", 0);
params.setValue("fragment:max_mz", 90000);
params.setValue("fragment:mass_tolerance", 0.05);
params.setValue("fragment:mass_tolerance_unit", "Da");
params.setValue("precursor:mass_tolerance", 2.0);
params.setValue("precursor:mass_tolerance_unit", "Da");
params.setValue("modifications:variable", std::vector<std::string>{});
params.setValue("modifications:fixed", std::vector<std::string>{});
tolTest.setParameters(params);

tolTest.build(entries);


TheoreticalSpectrumGenerator tsg;
PeakSpectrum b_y_ions;

AASequence peptide = AASequence::fromString("EVAEAATGEDASSPPPK");

tsg.getSpectrum(b_y_ions, peptide,1,1);

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
TEST_TRUE(found)

END_SECTION


END_TEST




