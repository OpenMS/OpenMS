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
//////////////////////////////
START_TEST(FragmentIndex, "$Id")

//////////////////////////////
using namespace OpenMS;
using namespace std;



/// Test the build for peptides
START_SECTION(build)
std::vector<FASTAFile::FASTAEntry> entries0{{"t", "t", "ARGEPADSSRKDFDMDMDM"}};
std::vector<FragmentIndex::Peptide> peptides_we_should_hit{{0,0,{2,8},5},
                                                            {0,0,{11,8},5},
                                                            {0,1,{11,8},5},
                                                            {0,2,{11,8},5},
                                                            {0,3,{11,8},5},
                                                            {0,4,{11,8},5},
                                                            {0,5,{11,8},5},
                                                            {0,6,{11,8},5}

};
FragmentIndex buildTest;
buildTest.build(entries0);

vector<FragmentIndex::Peptide> dbpeptides = buildTest.getPeptides();
TEST_EQUAL(dbpeptides.size(), peptides_we_should_hit.size())

bool testbuilda = true;
for(auto pepa: peptides_we_should_hit){
  bool found = false;
  for(auto pepb : dbpeptides){
    if((pepa.sequence_ == pepb.sequence_) && (pepa.modification_idx_ == pepb.modification_idx_)) found = true;
  }
  testbuilda = testbuilda && found;
}

TEST_TRUE(testbuilda)
END_SECTION


////TEST Different Charges of the query Spectrum ////
START_SECTION(void querySpectrum(const MSSpectrum& spectrum, SpectrumMatchesTopN& sms))



std::vector<FASTAFile::FASTAEntry> entries{{"test1", "test1","MSDEREVAEAATGEDASSPPPKTEAASDPQHPAASEGAAAAAASPPLLRCLVLTGFGGYDKVKLQSRPAAPPAPGPGQLTLRLRACGLNFADLMARQGLYDRLPPLPVTPGMEGAGVVIAVGEGVSDRKAGDRVMVLNRSGMWQEEVTVPSVQTFLIPEAMTFEEAAALLVNYITAYMVLFDFGNLQPGHSVLVHMAAGGVGMAAVQLCRTVENVTVFGTASASKHEALKENGVTHPIDYHTTDYVDEIKKISPKGVDIVMDPLGGSDTAKGYNLLKPMGKVVTYGMANLLTGPKRNLMALARTWWNQFSVTALQLLQANRAVCGFHLGYLDGEVELVSGVVARLLALYNQGHIKPHIDSVWPFEKVADAMKQMQEKKNVGKVLLVPGPEKEN"}};
AASequence protein = AASequence::fromString(entries[0].sequence);

FragmentIndex queryTest;
queryTest.build(entries);

auto param = queryTest.getParameters();
param.setValue("max_fragment_charge", 4);
param.setValue("fragment_max_mz", 5000000); // That we definitively create all peptides
queryTest.setParameters(param);
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



vector<FragmentIndex::Peptide> peptides = queryTest.getPeptides();
bool test = true;
bool every_peak_found_its_counter_part = true;

// Create different ms/ms spectra with different charges

for(uint16_t charge = 1; charge <= 4; charge++)
{
  int peptide_idx = 0;
// For each peptide that was created, we now generate a theoretical spectra for the given charge
  //Each peptide should hit its own entry in the db. In this case the test returns true
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
    spec_theo.setMSLevel(2);
    spec_theo.setPrecursors({prec_theo});
    for (auto ion : b_y_ions)
    {
      spec_theo.push_back(ion);
    }

    queryTest.querySpectrum(spec_theo, sms);
    bool found = false;
    bool found_all_peaks = false;
    for (auto s : sms.hits_)
    {
      if (s.peptide_idx_ == peptide_idx)
      {
        TEST_EQUAL(s.num_matched_, spec_theo.size())
        found_all_peaks = (s.num_matched_ == spec_theo.size());
        found = true;
      }
    }
    test = test && found;
    every_peak_found_its_counter_part = every_peak_found_its_counter_part && found_all_peaks;
    peptide_idx++;
  }
  TEST_TRUE(test);

}

END_SECTION


END_TEST




