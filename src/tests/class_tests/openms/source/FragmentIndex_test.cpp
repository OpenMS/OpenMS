// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////////
#include <OpenMS/ANALYSIS/ID/FragmentIndex.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/DigestionEnzyme.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

class FragmentIndex_test: public FragmentIndex
{
public:
  bool testDigestion(const std::vector<FragmentIndex::Peptide>& fe)
  {
    bool testbuilda = true;
    for (auto pepa : fe)
    {
      bool found = false;
      for (auto pepb : fi_peptides_)
      {
        if ((pepa.sequence_ == pepb.sequence_) && (pepa.modification_idx_ == pepb.modification_idx_))
          found = true;
      }
      testbuilda = testbuilda && found;
    }
    return testbuilda && (fe.size() == fi_peptides_.size());
  }
  bool peptidesSorted()
  {
    float last_mz = -1;
    for (auto& pep: fi_peptides_)
    {
      if (pep.precursor_mz_ >= last_mz)
      {
        last_mz = pep.precursor_mz_;
      }
      else
      {
        return  false;
      }
    }
    return true;
  }

  bool fragmentsSorted()
  {

    for (size_t fi_idx = 0; fi_idx < fi_fragments_.size(); fi_idx += bucketsize_)
    {
      UInt32 last_idx = 0;
      for (size_t bucket_idx = fi_idx; bucket_idx < fi_fragments_.size(); bucket_idx++)
      {
        if(bucket_idx >= fi_idx + bucketsize_)
          break;
        if(fi_fragments_[bucket_idx].peptide_idx_ >= last_idx)
          last_idx = fi_fragments_[bucket_idx].peptide_idx_;
        else
          return false;
      }
    }
    return true;
  }

  bool testQuery (UInt32 charge, bool precursor_mz_known, std::vector<FASTAFile::FASTAEntry> entries)
  {
    //fetch the parameters for the modifications
    auto param = getParameters();
    StringList modifications_fixed_ = ListUtils::toStringList<std::string>(param.getValue("modifications_fixed"));
    StringList modifications_variable_ = ListUtils::toStringList<std::string>(param.getValue("modifications_variable"));
    ModifiedPeptideGenerator::MapToResidueType fixed_modifications = ModifiedPeptideGenerator::getModifications(modifications_fixed_);
    ModifiedPeptideGenerator::MapToResidueType variable_modifications = ModifiedPeptideGenerator::getModifications(modifications_variable_);

    //Create theoretical spectra for different charges
    TheoreticalSpectrumGenerator tsg;
    vector<AASequence> mod_peptides;
    PeakSpectrum b_y_ions;
    MSSpectrum spec_theo;
    Precursor prec_theo;

    vector<FragmentIndex::Peptide> peptides = getPeptides();
    bool test = true;

    // Create different ms/ms spectra with different charges

      int peptide_idx = 0;
      // For each peptide that was created, we now generate a theoretical spectra for the given charge
      // Each peptide should hit its own entry in the db. In this case the test returns true
      for (auto& pep : peptides)
      {
        FragmentIndex::SpectrumMatchesTopN sms;
        b_y_ions.clear(true);
        mod_peptides.clear();
        spec_theo.clear(true);

        prec_theo.clearMetaInfo();
        AASequence unmod_peptide = AASequence::fromString(entries[0].sequence.substr(pep.sequence_.first, pep.sequence_.second));
        AASequence mod_peptide = AASequence(unmod_peptide); // copy the peptide
        ModifiedPeptideGenerator::applyFixedModifications(fixed_modifications, mod_peptide);
        ModifiedPeptideGenerator::applyVariableModifications(variable_modifications, mod_peptide, param.getValue("max_variable_mods_per_peptide"), mod_peptides);
        mod_peptide = mod_peptides[pep.modification_idx_];
        tsg.getSpectrum(b_y_ions, mod_peptide, charge, charge);
        prec_theo.setMZ(mod_peptide.getMZ(charge));
        if(precursor_mz_known)
          prec_theo.setCharge(charge);
        spec_theo.setMSLevel(2);
        spec_theo.setPrecursors({prec_theo});
        for (auto ion : b_y_ions)
        {
          spec_theo.push_back(ion);
        }

        querySpectrum(spec_theo, sms);
        bool found = false;
        bool found_all_peaks = false;
        for (auto s : sms.hits_)
        {
          if (s.peptide_idx_ == peptide_idx)
          {
            // The peak was found: All peaks which were created were found, and the correct charge was identified
            found = (s.num_matched_ >= spec_theo.size()) && (s.precursor_charge_, charge);
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
std::vector<FASTAFile::FASTAEntry> entries0{{"t", "t", "ARGEPADSSRKDFDMDMDM"},
                                             {"t2", "t2", "HALLORTSCHSM"}};
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
std::vector<FragmentIndex::Peptide> peptides_unmod_no_minmax{{0,0,{0,2},5},
                                                              {0,0,{2,8},5},
                                                              {0,0,{10,1},5},
                                                              {0,0,{11,8},5},
                                                              {1,0,{0,6},5},
                                                              {1,0,{6,6},5}};

std::vector<FragmentIndex::Peptide> peptides_unmod_minmax{{0,0,{0,2},5},
                                                           {1,0,{0,6},5},
                                                           {1,0,{6,6},5}};
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
auto param = buildTest.getParameters();
param.setValue("digestor_enzyme", "Trypsin");
param.setValue("missed_cleavages", 0);
param.setValue("peptide_min_mass", 0);
param.setValue("peptide_min_length", 0);
param.setValue("peptide_max_mass", 5000);
param.setValue("modifications_variable", std::vector<std::string>{});
param.setValue("modifications_fixed", std::vector<std::string>{});
buildTest.setParameters(param);

buildTest.build(entries0);
TEST_TRUE(buildTest.testDigestion(peptides_unmod_no_minmax))
TEST_TRUE(buildTest.peptidesSorted())
TEST_TRUE(buildTest.fragmentsSorted())

buildTest.clear();
param.setValue("peptide_min_length", 2);
param.setValue("peptide_max_length", 6);
buildTest.setParameters(param);
buildTest.build(entries0);
TEST_TRUE(buildTest.testDigestion(peptides_unmod_minmax))
TEST_TRUE(buildTest.peptidesSorted())
TEST_TRUE(buildTest.fragmentsSorted())

buildTest.clear();
param.setValue("peptide_max_length", 100);
param.setValue("missed_cleavages", 1);
buildTest.setParameters(param);
buildTest.build(entries0);
TEST_TRUE(buildTest.testDigestion(peptides_unmod_minmax_missed_cleavage))
TEST_TRUE(buildTest.peptidesSorted())
TEST_TRUE(buildTest.fragmentsSorted())

buildTest.clear();
param.setValue("digestor_enzyme", "Trypsin");
param.setValue("missed_cleavages", 0);
param.setValue("peptide_min_mass", 0);
param.setValue("peptide_min_length", 6);
param.setValue("modifications_variable", std::vector<std::string>{"Oxidation (M)"});
param.setValue("modifications_fixed", std::vector<std::string>{"Carbamidomethyl (C)"});
buildTest.setParameters(param);
buildTest.build(entries0);
TEST_TRUE(buildTest.testDigestion(peptides_we_should_hit_mod))
TEST_TRUE(buildTest.peptidesSorted())
TEST_TRUE(buildTest.fragmentsSorted())

END_SECTION

START_SECTION(clear())
std::vector<FASTAFile::FASTAEntry> entries0{{"t", "t", "ARGEPADSSRKDFDMDMDM"},
                                             {"t2", "t2", "HALLORTSCHS"}};
FragmentIndex clearTest;
clearTest.build(entries0);
clearTest.clear();

TEST_TRUE(clearTest.getPeptides().empty())

END_SECTION

START_SECTION(setParameters())



END_SECTION




////TEST Different Charges of the query Spectrum ////
START_SECTION(void querySpectrum(const MSSpectrum& spectrum, SpectrumMatchesTopN& sms))



std::vector<FASTAFile::FASTAEntry> entries{{"test1", "test1","MSDEREVAEAATGEDASSPPPKTEAASDPQHPAASEGAAAAAASPPLLRCLVLTGFGGYDKVKLQSRPAAPPAPGPGQLTLRLRACGLNFADLMARQGLYDRLPPLPVTPGMEGAGVVIAVGEGVSDRKAGDRVMVLNRSGMWQEEVTVPSVQTFLIPEAMTFEEAAALLVNYITAYMVLFDFGNLQPGHSVLVHMAAGGVGMAAVQLCRTVENVTVFGTASASKHEALKENGVTHPIDYHTTDYVDEIKKISPKGVDIVMDPLGGSDTAKGYNLLKPMGKVVTYGMANLLTGPKRNLMALARTWWNQFSVTALQLLQANRAVCGFHLGYLDGEVELVSGVVARLLALYNQGHIKPHIDSVWPFEKVADAMKQMQEKKNVGKVLLVPGPEKEN"}};
AASequence protein = AASequence::fromString(entries[0].sequence);

FragmentIndex_test queryTest;
queryTest.build(entries);

auto param = queryTest.getParameters();
param.setValue("max_fragment_charge", 4);
param.setValue("min_precursor_charge", 1);
param.setValue("max_precursor_charge", 4);

param.setValue("fragment_max_mz", 5000000); // That we definitively create all peptides
queryTest.setParameters(param);


// Create different ms/ms spectra with different charges

for(uint16_t charge = 1; charge <= 4; charge++)
{
  TEST_TRUE(queryTest.testQuery(charge, false, entries))
  TEST_TRUE(queryTest.testQuery(charge, true, entries))
}

END_SECTION

START_SECTION(isotope_error)
std::vector<FASTAFile::FASTAEntry> entries{{"test1", "test1","MSDEREVAEAATGEDASSPPPKTEAASDPQHPAASEGAAAAAASPPLLRCLVLTGFGGYDKVKLQSRPAAPPAPGPGQLTLRLRACGLNFADLMARQGLYDRLPPLPVTPGMEGAGVVIAVGEGVSDRKAGDRVMVLNRSGMWQEEVTVPSVQTFLIPEAMTFEEAAALLVNYITAYMVLFDFGNLQPGHSVLVHMAAGGVGMAAVQLCRTVENVTVFGTASASKHEALKENGVTHPIDYHTTDYVDEIKKISPKGVDIVMDPLGGSDTAKGYNLLKPMGKVVTYGMANLLTGPKRNLMALARTWWNQFSVTALQLLQANRAVCGFHLGYLDGEVELVSGVVARLLALYNQGHIKPHIDSVWPFEKVADAMKQMQEKKNVGKVLLVPGPEKEN"}};

FragmentIndex_test isoTest;
isoTest.build(entries);

auto param = isoTest.getParameters();
param.setValue("min_isotope_error", -3);
param.setValue("max_isotope_error", 3);
isoTest.setParameters(param);

TheoreticalSpectrumGenerator tsg;
PeakSpectrum b_y_ions;
AASequence peptide = AASequence::fromString("EVAEAATGEDASSPPPK");
tsg.getSpectrum(b_y_ions, peptide,1,1);
MSSpectrum theo_spec;
Precursor theo_prec;
theo_prec.setCharge(1);
theo_spec.setMSLevel(2);
for (auto& peak : b_y_ions)
{

  theo_spec.push_back(peak);
}
for (int iso = -3; iso <= 3; iso++)
{
  theo_prec.setMZ(peptide.getMZ(1) + iso );
  theo_spec.setPrecursors({theo_prec});
  FragmentIndex::SpectrumMatchesTopN sms;
  isoTest.querySpectrum(theo_spec, sms);
  bool found = false;
  for (auto& hit : sms.hits_){
      auto result = isoTest.getPeptides()[hit.peptide_idx_];
      auto psize = peptide.size();
      TEST_EQUAL(result.sequence_.first, 5)
      TEST_EQUAL(result.sequence_.second, psize)
      found = true;
  }
  TEST_TRUE(found);
}

END_SECTION

START_SECTION(tolerance)
std::vector<FASTAFile::FASTAEntry> entries{{"test1", "test1","MSDEREVAEAATGEDASSPPPKTEAASDPQHPAASEGAAAAAASPPLLRCLVLTGFGGYDKVKLQSRPAAPPAPGPGQLTLRLRACGLNFADLMARQGLYDRLPPLPVTPGMEGAGVVIAVGEGVSDRKAGDRVMVLNRSGMWQEEVTVPSVQTFLIPEAMTFEEAAALLVNYITAYMVLFDFGNLQPGHSVLVHMAAGGVGMAAVQLCRTVENVTVFGTASASKHEALKENGVTHPIDYHTTDYVDEIKKISPKGVDIVMDPLGGSDTAKGYNLLKPMGKVVTYGMANLLTGPKRNLMALARTWWNQFSVTALQLLQANRAVCGFHLGYLDGEVELVSGVVARLLALYNQGHIKPHIDSVWPFEKVADAMKQMQEKKNVGKVLLVPGPEKEN"}};

FragmentIndex_test tolTest;
tolTest.build(entries);

auto param = tolTest.getParameters();
param.setValue("fragment_mz_tolerance", 0.05);
param.setValue("precursor_mz_tolerance", 2.0);
tolTest.setParameters(param);

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
for (auto& peak : b_y_ions)
{
  float factor = (rand() % 90 - 45) * 0.001;
  peak.setMZ(peak.getMZ()+factor); //
  theo_spec.push_back(peak);
}

FragmentIndex::SpectrumMatchesTopN sms;
tolTest.querySpectrum(theo_spec, sms);
bool found = false;
for (auto& hit : sms.hits_)
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




