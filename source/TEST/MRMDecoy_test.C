// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger, Hannes Roest, Witold Wolski $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/TraMLFile.h>

#include <boost/assign/std/vector.hpp>
#include <boost/assign/list_of.hpp>

///////////////////////////
#define private public
#include <OpenMS/ANALYSIS/OPENSWATH/MRMDecoy.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MRMDecoy, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MRMDecoy* ptr = 0;
MRMDecoy* nullPointer = 0;

START_SECTION(MRMDecoy())
{
  ptr = new MRMDecoy();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~MRMDecoy())
{
  delete ptr;
}
END_SECTION

START_SECTION(pair getDecoyIon(String ionid, map<String, map<String, DoubleReal> > & decoy_ionseries))
{
	MRMDecoy gen;

  OpenMS::TargetedExperiment::Peptide peptide;
  peptide.sequence = "TESTPEPTIDE";

  OpenMS::AASequence aas = gen.getAASequence(peptide);

	double ProductMZ = 371.66692;
	double mz_threshold = 0.8;
	int precursor_charge = 2;

	MRMDecoy::IonSeries reference_ionseries = gen.getIonSeries(aas, precursor_charge);
	std::pair<String, double> targetion = gen.getTargetIon(ProductMZ,mz_threshold, reference_ionseries);
	std::pair<String, double> decoyion = gen.getDecoyIon("b7/2+", reference_ionseries);
	std::pair<String, double> decoyion_missing = gen.getDecoyIon("b17/2+", reference_ionseries);

	TEST_EQUAL(targetion.first,"b7/2+")
	TEST_REAL_SIMILAR(targetion.second,decoyion.second)
	TEST_REAL_SIMILAR(decoyion_missing.second,-1)
}
END_SECTION

START_SECTION(vector find_all_tryptic(string sequence))
{    
  MRMDecoy gen;

  String sequence = "TRESTPEPTIKDE";
  std::vector< pair<std::string::size_type,std::string> > tryptic_results = gen.find_all_tryptic(sequence);
  std::vector< pair<std::string::size_type,std::string> > tryptic_control = boost::assign::list_of(std::make_pair(1,"R"))(std::make_pair(5,"P"))(std::make_pair(7,"P"))(std::make_pair(10,"K"));

  for(Size i = 0; i < tryptic_results.size(); i++)
  {
    pair<std::string::size_type,std::string> result = tryptic_results[i];
    pair<std::string::size_type,std::string> control = tryptic_control[i];
    TEST_EQUAL(result.first,control.first)
    TEST_EQUAL(result.second,control.second)
  }
}
END_SECTION

START_SECTION(OpenMS::AASequence getAASequence(OpenMS::TargetedExperiment::Peptide peptide))
{
  MRMDecoy gen;

  OpenMS::TargetedExperiment::Peptide peptide;
  peptide.sequence = "TESTPEPTIDE";
  OpenMS::TargetedExperiment::Peptide::Modification modification;
  modification.avg_mass_delta = 79.9799;
  modification.location = 2;
  modification.mono_mass_delta = 79.966331;
  peptide.mods.push_back(modification);

  OpenMS::AASequence aas = gen.getAASequence(peptide);
  OpenMS::String modified_sequence = "TES(Phospho)TPEPTIDE";
  TEST_EQUAL(aas.toUnmodifiedString(),peptide.sequence)
  //TEST_EQUAL(aas.toString(),modified_sequence)

  OpenMS::TargetedExperiment::Peptide peptide2;
  peptide2.sequence = "TESTPEPTIDER";
  OpenMS::TargetedExperiment::Peptide::Modification modification2;
  modification2.avg_mass_delta = 9.9296;
  modification2.location = 11;
  modification2.mono_mass_delta = 10.008269;
  peptide2.mods.push_back(modification2);

  OpenMS::AASequence aas2 = gen.getAASequence(peptide2);
  OpenMS::String modified_sequence2 = "TESTPEPTIDER(Label:13C(6)15N(4))";
  TEST_EQUAL(aas2.toUnmodifiedString(),peptide2.sequence)
  //TEST_EQUAL(aas2.toString(),modified_sequence2)

  OpenMS::TargetedExperiment::Peptide peptide3;
  peptide3.sequence = "TESTMPEPTIDE";
  OpenMS::TargetedExperiment::Peptide::Modification modification3;
  modification3.avg_mass_delta = 15.9994;
  modification3.location = 4;
  modification3.mono_mass_delta = 15.994915;
  peptide3.mods.push_back(modification3);

  OpenMS::AASequence aas3 = gen.getAASequence(peptide3);
  OpenMS::String modified_sequence3 = "TESTM(Oxidation)PEPTIDE";
  TEST_EQUAL(aas3.toUnmodifiedString(),peptide3.sequence)
  //TEST_EQUAL(aas3.toString(),modified_sequence3)
}
END_SECTION

START_SECTION(OpenMS::TargetedExperiment::Peptide shufflePeptide(OpenMS::TargetedExperiment::Peptide peptide, double identity_threshold))
{
  MRMDecoy gen;
  OpenMS::TargetedExperiment::Peptide peptide;
  peptide.sequence = "TESTPEPTIDE";
  OpenMS::TargetedExperiment::Peptide::Modification modification;
  modification.avg_mass_delta = 79.9799;
  modification.location = 2;
  modification.mono_mass_delta = 79.966331;
  peptide.mods.push_back(modification);

  OpenMS::String expected_sequence = "TETTPEPESID";
  OpenMS::Size expected_location = 8;

  OpenMS::TargetedExperiment::Peptide shuffled = gen.shufflePeptide(peptide, 0.7 , 43);

  TEST_EQUAL(shuffled.sequence,expected_sequence)
  TEST_EQUAL(shuffled.mods[0].location,expected_location)

  OpenMS::TargetedExperiment::Peptide shuffleAASequence_target_sequence_11;
  shuffleAASequence_target_sequence_11.sequence = "TESTPEPTIDE";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_expected_11;
  shuffleAASequence_expected_11.sequence = "TESTPEPTIDE";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_result_11;
  shuffleAASequence_result_11 = gen.shufflePeptide(shuffleAASequence_target_sequence_11, 1.1 , 42);
  TEST_EQUAL(shuffleAASequence_result_11.sequence, shuffleAASequence_expected_11.sequence)

  OpenMS::TargetedExperiment::Peptide shuffleAASequence_target_sequence_07;
  shuffleAASequence_target_sequence_07.sequence = "TESTPEPTIDE";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_expected_07;
  shuffleAASequence_expected_07.sequence = "ETSTPDPEETI";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_result_07;
  shuffleAASequence_result_07 = gen.shufflePeptide(shuffleAASequence_target_sequence_07, 0.7, 42);
  TEST_EQUAL(shuffleAASequence_result_07.sequence, shuffleAASequence_expected_07.sequence)

  OpenMS::TargetedExperiment::Peptide shuffleAASequence_target_sequence_09;
  shuffleAASequence_target_sequence_09.sequence = "TESTPEPTIDE";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_expected_09;
  shuffleAASequence_expected_09.sequence = "ETSTPDPEETI";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_result_09;
  shuffleAASequence_result_09 = gen.shufflePeptide(shuffleAASequence_target_sequence_09, 0.9 , 42);
  TEST_EQUAL(shuffleAASequence_result_09.sequence, shuffleAASequence_expected_09.sequence)

  OpenMS::TargetedExperiment::Peptide shuffleAASequence_target_sequence_01;
  shuffleAASequence_target_sequence_01.sequence = "TESTPEPTIDE";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_expected_01;
  shuffleAASequence_expected_01.sequence = "ITEEPCPETGTDS";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_result_01;
  shuffleAASequence_result_01 = gen.shufflePeptide(shuffleAASequence_target_sequence_01, 0.2 , 42, 10000);
  TEST_EQUAL(shuffleAASequence_result_01.sequence, shuffleAASequence_expected_01.sequence)

  OpenMS::TargetedExperiment::Peptide shuffleAASequence_target_sequence_00;
  shuffleAASequence_target_sequence_00.sequence = "TESTPEPTIDE";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_expected_00;
  shuffleAASequence_expected_00.sequence = "EEDEPTPTGICST";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_result_00;
	TEST_EXCEPTION(Exception::IllegalArgument,gen.shufflePeptide(shuffleAASequence_target_sequence_00, 0.0 , 42))
}
END_SECTION

START_SECTION(shuffle_peptide_with_modifications_and2attempts)
{
  // Regression test for JIRA issue ABL-749
  // A peptide with modifications that was shuffled twice did not get its
  // modifications shuffled as well.
  MRMDecoy gen;
  OpenMS::TargetedExperiment::Peptide peptide;
  peptide.sequence = "GPPSEDGPGVPPPSPR";
  OpenMS::TargetedExperiment::Peptide::Modification modification;
  modification.avg_mass_delta = 79.9799;
  modification.location = 3;
  modification.mono_mass_delta = 79.966331;
  peptide.mods.push_back(modification);
  modification.avg_mass_delta = 79.9799;
  modification.location = 13;
  modification.mono_mass_delta = 79.966331;
  peptide.mods.push_back(modification);

  OpenMS::String expected_sequence = "GPPEVSGPGSPPPDPR";
  OpenMS::Size expected_location_1 = 5;
  OpenMS::Size expected_location_2 = 9;

  OpenMS::TargetedExperiment::Peptide shuffled = gen.shufflePeptide(peptide, 0.7 , 130);

  // the two modifcations get switched
  TEST_EQUAL(shuffled.sequence, expected_sequence)
  TEST_EQUAL(shuffled.mods[1].location,expected_location_1)
  TEST_EQUAL(shuffled.mods[0].location,expected_location_2)
}
END_SECTION

START_SECTION(shuffle_peptide_with_KPR)
{
  MRMDecoy gen;
  OpenMS::TargetedExperiment::Peptide peptide;
  peptide.sequence = "KPRKPRPK";
  OpenMS::String expected_sequence = "KPRKPRPKLN";
	TEST_EXCEPTION(Exception::IllegalArgument,gen.shufflePeptide(peptide, 0.7 , 130))
}
END_SECTION

START_SECTION(float AASequenceIdentity(const String & sequence, const String & decoy))
{
  MRMDecoy gen;

  String AASequenceIdentity_target_sequence = "TESTPEPTIDE";
  String AASequenceIdentity_decoy_sequence = "EDITPEPTSET";
  float AASequenceIdentity_result = gen.AASequenceIdentity(AASequenceIdentity_target_sequence,AASequenceIdentity_decoy_sequence);
  float AASequenceIdentity_expected = 0.454545;
  TEST_REAL_SIMILAR(AASequenceIdentity_result, AASequenceIdentity_expected)
}
END_SECTION

START_SECTION(OpenMS::TargetedExperiment::Peptide pseudoreversePeptide(OpenMS::TargetedExperiment::Peptide peptide))
{
  MRMDecoy gen;

  OpenMS::TargetedExperiment::Peptide peptide;
  peptide.sequence = "TESTPEPTIDE";
  OpenMS::TargetedExperiment::Peptide::Modification modification;
  modification.avg_mass_delta = 79.9799;
  modification.location = 2;
  modification.mono_mass_delta = 79.966331;
  peptide.mods.push_back(modification);

  OpenMS::String expected_sequence = "DITPEPTSETE";
  OpenMS::Size expected_location = 7;

  OpenMS::TargetedExperiment::Peptide pseudoreverse = gen.pseudoreversePeptide(peptide);
  TEST_EQUAL(pseudoreverse.sequence,expected_sequence)
  TEST_EQUAL(pseudoreverse.mods[0].location,expected_location)

  OpenMS::TargetedExperiment::Peptide pseudoreverseAASequence_target_sequence;
  pseudoreverseAASequence_target_sequence.sequence = "TESTPEPTIDE";
  OpenMS::TargetedExperiment::Peptide pseudoreverseAASequence_expected;
  pseudoreverseAASequence_expected.sequence = "DITPEPTSETE";
  OpenMS::TargetedExperiment::Peptide pseudoreverseAASequence_result;
  pseudoreverseAASequence_result = gen.pseudoreversePeptide(pseudoreverseAASequence_target_sequence);
  TEST_EQUAL(pseudoreverseAASequence_result.sequence,pseudoreverseAASequence_expected.sequence)
}
END_SECTION

START_SECTION(OpenMS::TargetedExperiment::Peptide reversePeptide(OpenMS::TargetedExperiment::Peptide peptide))
{
  MRMDecoy gen;

  OpenMS::TargetedExperiment::Peptide peptide;
  peptide.sequence = "TESTPEPTIDE";
  OpenMS::TargetedExperiment::Peptide::Modification modification;
  modification.avg_mass_delta = 79.9799;
  modification.location = 2;
  modification.mono_mass_delta = 79.966331;
  peptide.mods.push_back(modification);

  OpenMS::String expected_sequence = "EDITPEPTSET";
  OpenMS::Size expected_location = 8;

  OpenMS::TargetedExperiment::Peptide reverse = gen.reversePeptide(peptide);
  TEST_EQUAL(reverse.sequence,expected_sequence)
  TEST_EQUAL(reverse.mods[0].location,expected_location)

  OpenMS::TargetedExperiment::Peptide reverseAASequence_target_sequence;
  reverseAASequence_target_sequence.sequence = "TESTPEPTIDE";
  OpenMS::TargetedExperiment::Peptide reverseAASequence_expected;
  reverseAASequence_expected.sequence = "EDITPEPTSET";
  OpenMS::TargetedExperiment::Peptide reverseAASequence_result;
  reverseAASequence_result = gen.reversePeptide(reverseAASequence_target_sequence);
  TEST_EQUAL(reverseAASequence_result.sequence,reverseAASequence_expected.sequence)
}
END_SECTION

START_SECTION(void generateDecoys(OpenMS::TargetedExperiment& exp, OpenMS::TargetedExperiment& dec, String method, String decoy_tag, double identity_threshold, double mz_threshold, bool theoretical, double mz_shift))
{
  String method = "pseudo-reverse";
  DoubleReal identity_threshold = 0.7;
	Int max_attempts = 10;
  DoubleReal mz_threshold = 0.8;
  DoubleReal mz_shift = 20;
  String decoy_tag = "DECOY_";
  Int min_transitions = 2;
  Int max_transitions = 6;
  bool theoretical = 1;
  String in = "MRMDecoyGenerator_input.TraML";
  String out = "MRMDecoyGenerator_output.TraML";
  String test;
  NEW_TMP_FILE(test);

  TraMLFile traml;
  TargetedExperiment targeted_exp;
  TargetedExperiment targeted_decoy;
        
  traml.load(OPENMS_GET_TEST_DATA_PATH(in), targeted_exp);
    
  MRMDecoy decoys = MRMDecoy();
  decoys.restrictTransitions(targeted_exp, min_transitions, max_transitions);
  decoys.generateDecoys(targeted_exp, targeted_decoy, method, decoy_tag, identity_threshold, max_attempts, mz_threshold, theoretical, mz_shift);
  traml.store(OPENMS_GET_TEST_DATA_PATH(test), targeted_decoy);
  
  TEST_FILE_EQUAL(OPENMS_GET_TEST_DATA_PATH(test), OPENMS_GET_TEST_DATA_PATH(out))
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



