// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  Copyright The OpenMS team, Eberhard Karls University Tübingen,
//  ETH Zürich and FU Berlin 2001-2012.
//  This software is released under a BSD license. For a full list of
//  authors, refer to the file AUTHORS. For full licensing conditions
//  refer to the file LICENSE.
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
  OpenMS::String modified_sequence3 = "TESTM(Oxidation)PEPTIDER";
  TEST_EQUAL(peptide3.sequence, aas3.toUnmodifiedString())
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
  shuffleAASequence_expected_01.sequence = "EEDEPTPTGICST";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_result_01;
  shuffleAASequence_result_01 = gen.shufflePeptide(shuffleAASequence_target_sequence_01, 0.1 , 42);
  TEST_EQUAL(shuffleAASequence_result_01.sequence, shuffleAASequence_expected_01.sequence)

  OpenMS::TargetedExperiment::Peptide shuffleAASequence_target_sequence_00;
  shuffleAASequence_target_sequence_00.sequence = "TESTPEPTIDE";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_expected_00;
  shuffleAASequence_expected_00.sequence = "EEDEPTPTGICST";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_result_00;
  shuffleAASequence_result_00 = gen.shufflePeptide(shuffleAASequence_target_sequence_00, 0.0 , 42);
  TEST_EQUAL(shuffleAASequence_result_00.sequence, shuffleAASequence_expected_00.sequence)
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
  OpenMS::TargetedExperiment::Peptide shuffled = gen.shufflePeptide(peptide, 0.7 , 130);

  TEST_EQUAL(shuffled.sequence, expected_sequence)
  TEST_REAL_SIMILAR(gen.AASequenceIdentity(peptide.sequence, shuffled.sequence), 1.0)
  // TODO(georger) this is not a successfull execution of the shufflePeptide method! The
  // sequence identity is still 100% => how to indicate failure?
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

  OpenMS::String expected_sequence = "DITPEPTSETE";
  OpenMS::Size expected_location = 7;

  OpenMS::TargetedExperiment::Peptide reverse = gen.reversePeptide(peptide);
  TEST_EQUAL(reverse.sequence,expected_sequence)
  TEST_EQUAL(reverse.mods[0].location,expected_location)

  OpenMS::TargetedExperiment::Peptide reverseAASequence_target_sequence;
  reverseAASequence_target_sequence.sequence = "TESTPEPTIDE";
  OpenMS::TargetedExperiment::Peptide reverseAASequence_expected;
  reverseAASequence_expected.sequence = "DITPEPTSETE";
  OpenMS::TargetedExperiment::Peptide reverseAASequence_result;
  reverseAASequence_result = gen.reversePeptide(reverseAASequence_target_sequence);
  TEST_EQUAL(reverseAASequence_result.sequence,reverseAASequence_expected.sequence)
}
END_SECTION

START_SECTION(OpenMS::TargetedExperiment::Peptide trypticreversePeptide(OpenMS::TargetedExperiment::Peptide peptide))
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

  OpenMS::TargetedExperiment::Peptide trypticreverse = gen.trypticreversePeptide(peptide);
  TEST_EQUAL(trypticreverse.sequence,expected_sequence)
  TEST_EQUAL(trypticreverse.mods[0].location,expected_location)

  OpenMS::TargetedExperiment::Peptide trypticreverseAASequence_target_sequence;
  trypticreverseAASequence_target_sequence.sequence = "TESTPEPTIDE";
  OpenMS::TargetedExperiment::Peptide trypticreverseAASequence_expected;
  trypticreverseAASequence_expected.sequence = "EDITPEPTSET";
  OpenMS::TargetedExperiment::Peptide trypticreverseAASequence_result;
  trypticreverseAASequence_result = gen.trypticreversePeptide(trypticreverseAASequence_target_sequence);
  TEST_EQUAL(trypticreverseAASequence_result.sequence,trypticreverseAASequence_expected.sequence)
}
END_SECTION

START_SECTION(void generateDecoys(OpenMS::TargetedExperiment& exp, OpenMS::TargetedExperiment& dec, String method, String decoy_tag, double identity_threshold, double mz_threshold, bool theoretical))
{
  String method = "reverse";
  DoubleReal identity_threshold = 0.7;
  DoubleReal mz_threshold = 0.8;
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
  decoys.generateDecoys(targeted_exp, targeted_decoy, method, decoy_tag, identity_threshold, mz_threshold, theoretical);
  traml.store(OPENMS_GET_TEST_DATA_PATH(test), targeted_decoy);
  
  TEST_FILE_EQUAL(OPENMS_GET_TEST_DATA_PATH(test), OPENMS_GET_TEST_DATA_PATH(out))
}
END_SECTION

    
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



