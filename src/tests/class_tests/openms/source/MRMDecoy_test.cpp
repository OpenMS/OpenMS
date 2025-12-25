// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger, Hannes Roest, Witold Wolski $
// --------------------------------------------------------------------------

#include <OpenMS/test_config.h>
#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/MRMDecoy.h>
///////////////////////////

#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

#include <boost/assign/std/vector.hpp>

#ifdef __clang__
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wshadow"
#endif

using namespace OpenMS;
using namespace std;

class MRMDecoyHelper : public MRMDecoy
{
public:
  OpenMS::TargetedExperiment::Peptide pseudoreversePeptide_helper(
    const OpenMS::TargetedExperiment::Peptide& peptide) const
  {
    return pseudoreversePeptide_(peptide);
  }
  OpenMS::TargetedExperiment::Peptide reversePeptide_helper(
    const OpenMS::TargetedExperiment::Peptide& peptide) const
  {
    return reversePeptide_(peptide);
  }
  IndexType findFixedResidues_helper(const std::string& sequence) const {return findFixedResidues_(sequence);}
  IndexType findFixedAndTermResidues_helper(const std::string& sequence) const {return findFixedAndTermResidues_(sequence);}

  // Helper method to expose protected pseudoreversePeptideLight_
  std::pair<std::string, std::vector<OpenSwath::LightModification>> pseudoreversePeptideLight_helper(
    const std::string& sequence,
    const std::vector<OpenSwath::LightModification>& modifications) const
  {
    return pseudoreversePeptideLight_(sequence, modifications);
  }

  // Helper method to expose protected hasCNterminalModsLight_
  static bool hasCNterminalModsLight_helper(
    const std::vector<OpenSwath::LightModification>& modifications,
    size_t sequence_length,
    bool checkCterminalAA)
  {
    return hasCNterminalModsLight_(modifications, sequence_length, checkCterminalAA);
  }
};

START_TEST(MRMDecoy, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MRMDecoy * ptr = nullptr;
MRMDecoy* nullPointer = nullptr;

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

START_SECTION((std::vector<std::pair<std::string::size_type, std::string> > findFixedResidues(std::string sequence)))
{
  MRMDecoyHelper gen;

  String sequence = "TRESTPEPTIKDE";
  MRMDecoy::IndexType tryptic_results = gen.findFixedResidues_helper(sequence);
  MRMDecoy::IndexType tryptic_expect = {1, 5, 7, 10};
  TEST_TRUE(tryptic_results == tryptic_expect)
}

END_SECTION

START_SECTION((std::vector<std::pair<std::string::size_type, std::string> > findFixedAndTermResidues(std::string sequence)))
{
  MRMDecoyHelper gen;

  String sequence = "TRESTPEPTIKDE";
  MRMDecoy::IndexType tryptic_results = gen.findFixedAndTermResidues_helper(sequence);
  MRMDecoy::IndexType tryptic_expect = {0, 1, 5, 7, 10, 12};
  TEST_TRUE(tryptic_results == tryptic_expect)
}

END_SECTION

START_SECTION(OpenMS::TargetedExperiment::Peptide shufflePeptide(OpenMS::TargetedExperiment::Peptide peptide, double identity_threshold, int seed = -1, int max_attempts = 10))
{
  MRMDecoy gen;
  OpenMS::TargetedExperiment::Peptide peptide;
  peptide.sequence = "TESTPEPTIDE";
  OpenMS::TargetedExperiment::Peptide::Modification modification;
  modification.avg_mass_delta = 79.9799;
  modification.location = 2;
  modification.mono_mass_delta = 79.966331;
  peptide.mods.push_back(modification);

  OpenMS::String expected_sequence = "TIDEPEPSTTE";
  OpenMS::Size expected_location = 7;

  OpenMS::TargetedExperiment::Peptide shuffled = gen.shufflePeptide(peptide, 0.7, 43);

  TEST_EQUAL(shuffled.sequence, expected_sequence)
  TEST_EQUAL(shuffled.mods.size(), 1)
  TEST_EQUAL(shuffled.mods[0].location, expected_location)

  OpenMS::TargetedExperiment::Peptide shuffleAASequence_target_sequence_11;
  shuffleAASequence_target_sequence_11.sequence = "TESTPEPTIDE";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_expected_11;
  shuffleAASequence_expected_11.sequence = "TESTPEPTIDE";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_result_11;
  shuffleAASequence_result_11 = gen.shufflePeptide(shuffleAASequence_target_sequence_11, 1.1, 42);
  TEST_EQUAL(shuffleAASequence_result_11.sequence, shuffleAASequence_expected_11.sequence)

  OpenMS::TargetedExperiment::Peptide shuffleAASequence_target_sequence_07;
  shuffleAASequence_target_sequence_07.sequence = "TESTPEPTIDE";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_expected_07;
  shuffleAASequence_expected_07.sequence = "TTETPEPIDSE";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_result_07;
  shuffleAASequence_result_07 = gen.shufflePeptide(shuffleAASequence_target_sequence_07, 0.7, 42);
  TEST_EQUAL(shuffleAASequence_result_07.sequence, shuffleAASequence_expected_07.sequence)

  OpenMS::TargetedExperiment::Peptide shuffleAASequence_target_sequence_09;
  shuffleAASequence_target_sequence_09.sequence = "TESTPEPTIDE";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_expected_09;
  shuffleAASequence_expected_09.sequence = "TTETPEPIDSE";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_result_09;
  shuffleAASequence_result_09 = gen.shufflePeptide(shuffleAASequence_target_sequence_09, 0.9, 42);
  TEST_EQUAL(shuffleAASequence_result_09.sequence, shuffleAASequence_expected_09.sequence)

  OpenMS::TargetedExperiment::Peptide shuffleAASequence_target_sequence_01;
  shuffleAASequence_target_sequence_01.sequence = "TESTPEPTIDE";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_expected_01;
  shuffleAASequence_expected_01.sequence = "TNGCADQQEAE";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_result_01;
  shuffleAASequence_result_01 = gen.shufflePeptide(shuffleAASequence_target_sequence_01, 0.2, 42, 10000);
  TEST_EQUAL(shuffleAASequence_result_01.sequence, shuffleAASequence_expected_01.sequence)

  OpenMS::TargetedExperiment::Peptide shuffleAASequence_target_sequence_00;
  shuffleAASequence_target_sequence_00.sequence = "TESTPEPTIDE";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_expected_00;
  shuffleAASequence_expected_00.sequence = "TEIEPAPTQTE";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_result_00;
  shuffleAASequence_result_00 = gen.shufflePeptide(shuffleAASequence_target_sequence_00, 0.0, 42, 20);
  TEST_EQUAL(shuffleAASequence_result_00.sequence, shuffleAASequence_expected_00.sequence)

  OpenMS::TargetedExperiment::Peptide shuffleAASequence_target_sequence_01b;
  shuffleAASequence_target_sequence_01b.sequence = "TESTPEPTIDE";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_expected_01b;
  shuffleAASequence_expected_01b.sequence = "TNGCADQQEAE";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_result_01b;
  shuffleAASequence_result_01b = gen.shufflePeptide(shuffleAASequence_target_sequence_01b, 0.2, 42, 10000);
  TEST_EQUAL(shuffleAASequence_result_01b.sequence, shuffleAASequence_expected_01b.sequence)

  OpenMS::TargetedExperiment::Peptide shuffleAASequence_target_sequence_00b;
  shuffleAASequence_target_sequence_00b.sequence = "TESTPEPTIDE";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_expected_00b;
  shuffleAASequence_expected_00b.sequence = "TNDQIADNNEE";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_result_00b;
  shuffleAASequence_result_00b = gen.shufflePeptide(shuffleAASequence_target_sequence_00b, 0.0, 42, 2000);
  TEST_EQUAL(shuffleAASequence_result_00b.sequence, shuffleAASequence_expected_00b.sequence)
  // ensure that C-terminal K and R are preserved
  {
    OpenMS::TargetedExperiment::Peptide original_input;
    original_input.sequence = "TESTPEPTIDEK";
    expected_sequence = "TTETPEPEDSIK";
    OpenMS::TargetedExperiment::Peptide shuffleAASequence_result_00;
    shuffled = gen.shufflePeptide(original_input, 0.7, 42, 20);
    TEST_EQUAL(shuffled.sequence[shuffled.sequence.size() - 1], 'K')
    TEST_EQUAL(shuffled.sequence, expected_sequence)
  }

  // ensure that C-terminal K and R are preserved
  {
    OpenMS::TargetedExperiment::Peptide original_input;
    original_input.sequence = "TESTPEPTIDER";
    expected_sequence = "TTETPEPEDSIR";
    OpenMS::TargetedExperiment::Peptide shuffleAASequence_result_00;
    shuffled = gen.shufflePeptide(original_input, 0.7, 42, 20);
    TEST_EQUAL(shuffled.sequence[shuffled.sequence.size() - 1], 'R')
    TEST_EQUAL(shuffled.sequence, expected_sequence)
  }

  {
    OpenMS::TargetedExperiment::Peptide original_input;
    OpenMS::TargetedExperiment::Peptide::Modification mod;
    // std::vector<TargetedExperiment::Peptide::Modification> mods;
    std::vector<TargetedExperiment::Peptide::Modification> mods;
    // original_input.sequence = "EPAHLMSLFGGKPM(UniMod:35)";
    original_input.sequence = "EPAHLMSLFGGKPM";
    mod.location = 13; // non-C terminal
    mod.unimod_id = 35;
    mods.push_back(mod);
    original_input.mods = mods;
    expected_sequence = "EPSALMGGHLFKPM";
    OpenMS::TargetedExperiment::Peptide shuffleAASequence_result_00;
    shuffled = gen.shufflePeptide(original_input, 0.7, 42, 20);
    TEST_EQUAL(shuffled.sequence[shuffled.sequence.size() - 1], 'M')
    TEST_EQUAL(shuffled.sequence, expected_sequence)
    TEST_EQUAL(shuffled.mods.size(), 1)
    TEST_EQUAL(shuffled.mods[0].location, 13) // the second M remained at position 13
  }

  {
    OpenMS::TargetedExperiment::Peptide original_input;
    OpenMS::TargetedExperiment::Peptide::Modification mod;
    // std::vector<TargetedExperiment::Peptide::Modification> mods;
    std::vector<TargetedExperiment::Peptide::Modification> mods;
    // original_input.sequence = "EPAHLMSLFGGKPM(UniMod:35)";
    original_input.sequence = "EPAHLMSLFGGKPM";
    mod.location = 14; // C terminal
    mod.unimod_id = 35;
    mods.push_back(mod);
    original_input.mods = mods;
    expected_sequence = "EPSALMGGHLFKPM";
    OpenMS::TargetedExperiment::Peptide shuffleAASequence_result_00;
    shuffled = gen.shufflePeptide(original_input, 0.7, 42, 20);
    TEST_EQUAL(shuffled.sequence[shuffled.sequence.size() - 1], 'M')
    TEST_EQUAL(shuffled.sequence, expected_sequence)
    TEST_EQUAL(shuffled.mods.size(), 1)
    TEST_EQUAL(shuffled.mods[0].location, 14) // Problem: this modification cannot be C terminal any more for F!
    // TODO: report and fix this
  }

}

END_SECTION

START_SECTION([EXTRA] shuffle_peptide_with_modifications_and2attempts)
{
  // Regression test for JIRA issue ABL-749
  // A peptide with modifications that was shuffled twice did not get its
  // modifications shuffled as well.
  MRMDecoy gen;
  OpenMS::TargetedExperiment::Peptide peptide;
  peptide.sequence = "GPPSEDGPGVPPPSPR";
  OpenMS::TargetedExperiment::Peptide::Modification modification;

  // modification on the fourth S (counting starts at zero)
  modification.avg_mass_delta = 79.9799;
  modification.location = 3;
  modification.mono_mass_delta = 79.966331;
  peptide.mods.push_back(modification);

  // modification on the second to last S
  modification.avg_mass_delta = 79.9799;
  modification.location = 13;
  modification.mono_mass_delta = 79.966331;
  peptide.mods.push_back(modification);

  OpenMS::String expected_sequence = "GPPGDSEPGSPPPVPR";
  OpenMS::Size expected_location_1 = 9;
  OpenMS::Size expected_location_2 = 5;

  OpenMS::TargetedExperiment::Peptide shuffled = gen.shufflePeptide(peptide, 0.7, 130);

  // the two modifications get switched (the first S now comes after the second S)
  TEST_EQUAL(shuffled.sequence, expected_sequence)
  TEST_EQUAL(shuffled.mods.size(), 2)
  TEST_EQUAL(shuffled.mods[1].location, expected_location_1)
  TEST_EQUAL(shuffled.mods[0].location, expected_location_2)
}

END_SECTION

START_SECTION([EXTRA] shuffle_peptide_with_terminal_modifications)
{
  // Shuffle a peptide with C/N terminal modifications
  MRMDecoy gen;
  AASequence original_sequence = AASequence::fromString("(UniMod:272)TESTPEPTIDE(UniMod:193)");
  TEST_EQUAL(original_sequence.hasNTerminalModification(), true)
  TEST_EQUAL(original_sequence.hasCTerminalModification(), true)

  OpenMS::TargetedExperiment::Peptide peptide;
  OpenMS::TargetedExperiment::Peptide::Modification modification;
  peptide.sequence = original_sequence.toUnmodifiedString();

  // "sulfonation of N-terminus"
  modification.avg_mass_delta = 136.1265;
  modification.location = -1;
  modification.mono_mass_delta = 135.983029;
  peptide.mods.push_back(modification);

  //O18 label at both C-terminal oxygens
  modification.avg_mass_delta = 3.9995;
  modification.location = peptide.sequence.size();
  modification.mono_mass_delta = 4.008491;
  peptide.mods.push_back(modification);

  OpenMS::String expected_sequence = "TIDEPEPSTTE";

  OpenMS::TargetedExperiment::Peptide shuffled = gen.shufflePeptide(peptide, 0.7, 43);

  TEST_EQUAL(shuffled.sequence, expected_sequence)
  TEST_EQUAL(shuffled.mods.size(), 2)
  TEST_EQUAL(shuffled.mods[0].location, -1)
  TEST_EQUAL(shuffled.mods[1].location, shuffled.sequence.size())
}

END_SECTION

START_SECTION([EXTRA] shuffle_peptide_with_KPR)
{
  MRMDecoy gen;
  OpenMS::TargetedExperiment::Peptide peptide;
  peptide.sequence = "KPRKPRPK";
  OpenMS::String expected_sequence = "KNRKPRPK";
  OpenMS::TargetedExperiment::Peptide shuffled = gen.shufflePeptide(peptide, 0.7, 130, 17);
  TEST_EQUAL(shuffled.sequence, expected_sequence)
}

END_SECTION

START_SECTION(float AASequenceIdentity(const String& sequence, const String& decoy))
{
  MRMDecoy gen;

  String AASequenceIdentity_target_sequence = "TESTPEPTIDE";
  String AASequenceIdentity_decoy_sequence = "EDITPEPTSET";
  float AASequenceIdentity_result = gen.AASequenceIdentity(AASequenceIdentity_target_sequence, AASequenceIdentity_decoy_sequence);
  float AASequenceIdentity_expected = static_cast<float>(0.454545);
  TEST_REAL_SIMILAR(AASequenceIdentity_result, AASequenceIdentity_expected)
}

END_SECTION

START_SECTION((OpenMS::TargetedExperiment::Peptide MRMDecoy::reversePeptide(
      OpenMS::TargetedExperiment::Peptide peptide, const bool keepN, const bool keepC, 
      const String& const_pattern) const))
{
  MRMDecoy gen;

  OpenMS::TargetedExperiment::Peptide peptide;
  peptide.sequence = "TESTPEPTIDE";
  OpenMS::TargetedExperiment::Peptide::Modification modification;
  modification.avg_mass_delta = 79.9799;
  modification.location = 2;
  modification.mono_mass_delta = 79.966331;
  peptide.mods.push_back(modification);

  {
    OpenMS::String expected_sequence = "DITPEPTSETE";
    OpenMS::Size expected_location = 7;

    OpenMS::TargetedExperiment::Peptide pseudoreverse = MRMDecoy::reversePeptide(peptide, false, true);
    TEST_EQUAL(pseudoreverse.sequence, expected_sequence)
    TEST_EQUAL(pseudoreverse.mods.size(), 1)
    TEST_EQUAL(pseudoreverse.mods[0].location, expected_location)
  }

  {
    modification.avg_mass_delta = 49.9799;
    modification.mono_mass_delta = 49.966331;
    modification.location = 0;
    peptide.mods.push_back(modification);

    OpenMS::String expected_sequence = "TDITPEPTSEE";

    OpenMS::TargetedExperiment::Peptide pseudoreverse = MRMDecoy::reversePeptide(peptide, true, true);
    TEST_EQUAL(pseudoreverse.sequence, expected_sequence)
    TEST_EQUAL(pseudoreverse.mods.size(), 2)
    TEST_EQUAL(pseudoreverse.mods[0].location, 8)
    TEST_REAL_SIMILAR(pseudoreverse.mods[0].mono_mass_delta, 79.966331)

    TEST_EQUAL(pseudoreverse.mods[1].location, 0)
    TEST_REAL_SIMILAR(pseudoreverse.mods[1].mono_mass_delta, 49.966331)
  }

  {
    String const_pattern = "I";
    OpenMS::String expected_sequence = "TDTPEPTSIEE"; // "I" stays in place

    OpenMS::TargetedExperiment::Peptide pseudoreverse = MRMDecoy::reversePeptide(peptide, true, true, const_pattern);
    TEST_EQUAL(pseudoreverse.sequence, expected_sequence)
    TEST_EQUAL(pseudoreverse.mods.size(), 2)
    TEST_EQUAL(pseudoreverse.mods[0].location, 7)
    TEST_REAL_SIMILAR(pseudoreverse.mods[0].mono_mass_delta, 79.966331)

    TEST_EQUAL(pseudoreverse.mods[1].location, 0)
    TEST_REAL_SIMILAR(pseudoreverse.mods[1].mono_mass_delta, 49.966331)
  }

}
END_SECTION

START_SECTION(OpenMS::TargetedExperiment::Peptide pseudoreversePeptide(OpenMS::TargetedExperiment::Peptide peptide))
{
  MRMDecoyHelper gen;

  OpenMS::TargetedExperiment::Peptide peptide;
  peptide.sequence = "TESTPEPTIDE";
  OpenMS::TargetedExperiment::Peptide::Modification modification;
  modification.avg_mass_delta = 79.9799;
  modification.location = 2;
  modification.mono_mass_delta = 79.966331;
  peptide.mods.push_back(modification);

  OpenMS::String expected_sequence = "DITPEPTSETE";
  OpenMS::Size expected_location = 7;

  OpenMS::TargetedExperiment::Peptide pseudoreverse = gen.pseudoreversePeptide_helper(peptide);
  TEST_EQUAL(pseudoreverse.sequence, expected_sequence)
  TEST_EQUAL(pseudoreverse.mods.size(), 1)
  TEST_EQUAL(pseudoreverse.mods[0].location, expected_location)

  OpenMS::TargetedExperiment::Peptide pseudoreverseAASequence_target_sequence;
  pseudoreverseAASequence_target_sequence.sequence = "TESTPEPTIDE";
  OpenMS::TargetedExperiment::Peptide pseudoreverseAASequence_expected;
  pseudoreverseAASequence_expected.sequence = "DITPEPTSETE";
  OpenMS::TargetedExperiment::Peptide pseudoreverseAASequence_result;
  pseudoreverseAASequence_result = gen.pseudoreversePeptide_helper(pseudoreverseAASequence_target_sequence);
  TEST_EQUAL(pseudoreverseAASequence_result.sequence, pseudoreverseAASequence_expected.sequence)
}
END_SECTION

START_SECTION(OpenMS::TargetedExperiment::Peptide reversePeptide(OpenMS::TargetedExperiment::Peptide peptide))
{
  MRMDecoyHelper gen;

  OpenMS::TargetedExperiment::Peptide peptide;
  peptide.sequence = "TESTPEPTIDE";
  OpenMS::TargetedExperiment::Peptide::Modification modification;
  modification.avg_mass_delta = 79.9799;
  modification.location = 2;
  modification.mono_mass_delta = 79.966331;
  peptide.mods.push_back(modification);

  OpenMS::String expected_sequence = "EDITPEPTSET";
  OpenMS::Size expected_location = 8;

  OpenMS::TargetedExperiment::Peptide reverse = gen.reversePeptide_helper(peptide);
  TEST_EQUAL(reverse.sequence, expected_sequence)
  TEST_EQUAL(reverse.mods.size(), 1)
  TEST_EQUAL(reverse.mods[0].location, expected_location)

  OpenMS::TargetedExperiment::Peptide reverseAASequence_target_sequence;
  reverseAASequence_target_sequence.sequence = "TESTPEPTIDE";
  OpenMS::TargetedExperiment::Peptide reverseAASequence_expected;
  reverseAASequence_expected.sequence = "EDITPEPTSET";
  OpenMS::TargetedExperiment::Peptide reverseAASequence_result;
  reverseAASequence_result = gen.reversePeptide_helper(reverseAASequence_target_sequence);
  TEST_EQUAL(reverseAASequence_result.sequence, reverseAASequence_expected.sequence)
}

END_SECTION

/// Public methods

START_SECTION((void generateDecoys(const OpenMS::TargetedExperiment& exp,
                        OpenMS::TargetedExperiment& dec,
                        const String& method,
                        const double aim_decoy_fraction,
                        const bool switchKR,
                        const String& decoy_tag,
                        const int max_attempts,
                        const double identity_threshold,
                        const double precursor_mz_shift,
                        const double product_mz_shift,
                        const double product_mz_threshold,
                        const std::vector<String>& fragment_types,
                        const std::vector<size_t>& fragment_charges,
                        const bool enable_specific_losses,
                        const bool enable_unspecific_losses,
                        const int round_decPow = -4) const))
{
  String method = "pseudo-reverse";
  double identity_threshold = 0.7;
  Int max_attempts = 5;
  double product_mz_threshold = 0.8;
  double precursor_mz_shift = 0.1;
  double product_mz_shift = 20;
  String decoy_tag = "DECOY_";
  std::vector<String> fragment_types;
  fragment_types.push_back(String("b"));
  fragment_types.push_back(String("y"));
  fragment_types.push_back(String("a"));
  std::vector<size_t> fragment_charges;
  fragment_charges.push_back(1);
  fragment_charges.push_back(2);
  fragment_charges.push_back(3);
  fragment_charges.push_back(4);
  fragment_charges.push_back(5);
  bool enable_unspecific_losses = false;
  bool enable_specific_losses = true;

  String in = "MRMDecoyGenerator_input.TraML";
  String out = "MRMDecoyGenerator_output.TraML";
  String test;
  NEW_TMP_FILE(test);

  TraMLFile traml;
  TargetedExperiment targeted_exp;
  TargetedExperiment targeted_decoy;

  traml.load(OPENMS_GET_TEST_DATA_PATH(in), targeted_exp);

  MRMDecoy decoys = MRMDecoy();
  TEST_EQUAL(targeted_exp.getPeptides().size(), 13)
  TEST_EQUAL(targeted_exp.getTransitions().size(), 36)
  decoys.generateDecoys(targeted_exp, 
                        targeted_decoy,
                        method,
                        1.0,
                        false,
                        decoy_tag,
                        max_attempts,
                        identity_threshold,
                        precursor_mz_shift, 
                        product_mz_shift, 
                        product_mz_threshold,
                        fragment_types, 
                        fragment_charges,
                        enable_specific_losses, 
                        enable_unspecific_losses); 
  traml.store(test, targeted_decoy);

  TEST_FILE_SIMILAR(test.c_str(), OPENMS_GET_TEST_DATA_PATH(out))
}

END_SECTION

/////////////////////////////////////////////////////////////
// Tests comparing Light methods with Heavy methods
/////////////////////////////////////////////////////////////

START_SECTION([EXTRA] reversePeptideLight_matches_heavy)
{
  // Test that reversePeptideLight produces the same sequence and modification
  // positions as the heavy reversePeptide method

  // Test case 1: Simple peptide with one modification
  {
    OpenMS::TargetedExperiment::Peptide peptide;
    peptide.sequence = "TESTPEPTIDE";
    OpenMS::TargetedExperiment::Peptide::Modification modification;
    modification.avg_mass_delta = 79.9799;
    modification.location = 2;  // S at position 2
    modification.mono_mass_delta = 79.966331;
    modification.unimod_id = 21;  // Phospho
    peptide.mods.push_back(modification);

    // Create equivalent light modifications
    std::vector<OpenSwath::LightModification> light_mods;
    OpenSwath::LightModification light_mod;
    light_mod.location = 2;
    light_mod.unimod_id = 21;
    light_mods.push_back(light_mod);

    // Test keepN=false, keepC=true (pseudo-reverse style)
    auto heavy_result = MRMDecoy::reversePeptide(peptide, false, true);
    auto light_result = MRMDecoy::reversePeptideLight(peptide.sequence, light_mods, false, true);

    TEST_EQUAL(light_result.first, heavy_result.sequence)
    TEST_EQUAL(light_result.second.size(), heavy_result.mods.size())
    TEST_EQUAL(light_result.second[0].location, heavy_result.mods[0].location)
  }

  // Test case 2: Peptide with multiple modifications
  {
    OpenMS::TargetedExperiment::Peptide peptide;
    peptide.sequence = "TESTPEPTIDE";
    OpenMS::TargetedExperiment::Peptide::Modification mod1, mod2;
    mod1.location = 2;  // S
    mod1.unimod_id = 21;
    mod2.location = 0;  // N-term T (kept in place with keepN)
    mod2.unimod_id = 1;
    peptide.mods.push_back(mod1);
    peptide.mods.push_back(mod2);

    std::vector<OpenSwath::LightModification> light_mods;
    OpenSwath::LightModification lm1, lm2;
    lm1.location = 2;
    lm1.unimod_id = 21;
    lm2.location = 0;
    lm2.unimod_id = 1;
    light_mods.push_back(lm1);
    light_mods.push_back(lm2);

    // Test keepN=true, keepC=true
    auto heavy_result = MRMDecoy::reversePeptide(peptide, true, true);
    auto light_result = MRMDecoy::reversePeptideLight(peptide.sequence, light_mods, true, true);

    TEST_EQUAL(light_result.first, heavy_result.sequence)
    TEST_EQUAL(light_result.second.size(), heavy_result.mods.size())
    for (size_t i = 0; i < light_result.second.size(); ++i)
    {
      TEST_EQUAL(light_result.second[i].location, heavy_result.mods[i].location)
    }
  }

  // Test case 3: Full reverse (keepN=false, keepC=false)
  {
    OpenMS::TargetedExperiment::Peptide peptide;
    peptide.sequence = "TESTPEPTIDE";
    OpenMS::TargetedExperiment::Peptide::Modification mod;
    mod.location = 2;
    mod.unimod_id = 21;
    peptide.mods.push_back(mod);

    std::vector<OpenSwath::LightModification> light_mods;
    OpenSwath::LightModification lm;
    lm.location = 2;
    lm.unimod_id = 21;
    light_mods.push_back(lm);

    auto heavy_result = MRMDecoy::reversePeptide(peptide, false, false);
    auto light_result = MRMDecoy::reversePeptideLight(peptide.sequence, light_mods, false, false);

    TEST_EQUAL(light_result.first, heavy_result.sequence)
    TEST_EQUAL(light_result.second.size(), heavy_result.mods.size())
    TEST_EQUAL(light_result.second[0].location, heavy_result.mods[0].location)
  }

  // Test case 4: With const_pattern
  {
    OpenMS::TargetedExperiment::Peptide peptide;
    peptide.sequence = "TESTPEPTIDE";
    OpenMS::TargetedExperiment::Peptide::Modification mod;
    mod.location = 2;
    mod.unimod_id = 21;
    peptide.mods.push_back(mod);

    std::vector<OpenSwath::LightModification> light_mods;
    OpenSwath::LightModification lm;
    lm.location = 2;
    lm.unimod_id = 21;
    light_mods.push_back(lm);

    String const_pattern = "I";
    auto heavy_result = MRMDecoy::reversePeptide(peptide, true, true, const_pattern);
    auto light_result = MRMDecoy::reversePeptideLight(peptide.sequence, light_mods, true, true, const_pattern);

    TEST_EQUAL(light_result.first, heavy_result.sequence)
    TEST_EQUAL(light_result.second.size(), heavy_result.mods.size())
    TEST_EQUAL(light_result.second[0].location, heavy_result.mods[0].location)
  }

  // Test case 5: Empty modifications
  {
    OpenMS::TargetedExperiment::Peptide peptide;
    peptide.sequence = "TESTPEPTIDE";

    std::vector<OpenSwath::LightModification> light_mods;

    auto heavy_result = MRMDecoy::reversePeptide(peptide, false, true);
    auto light_result = MRMDecoy::reversePeptideLight(peptide.sequence, light_mods, false, true);

    TEST_EQUAL(light_result.first, heavy_result.sequence)
    TEST_EQUAL(light_result.second.size(), 0)
  }
}
END_SECTION

START_SECTION([EXTRA] pseudoreversePeptideLight_matches_heavy)
{
  MRMDecoyHelper gen;

  // Test case 1: Simple peptide with modification
  {
    OpenMS::TargetedExperiment::Peptide peptide;
    peptide.sequence = "TESTPEPTIDE";
    OpenMS::TargetedExperiment::Peptide::Modification modification;
    modification.avg_mass_delta = 79.9799;
    modification.location = 2;
    modification.mono_mass_delta = 79.966331;
    modification.unimod_id = 21;
    peptide.mods.push_back(modification);

    std::vector<OpenSwath::LightModification> light_mods;
    OpenSwath::LightModification lm;
    lm.location = 2;
    lm.unimod_id = 21;
    light_mods.push_back(lm);

    auto heavy_result = gen.pseudoreversePeptide_helper(peptide);
    auto light_result = gen.pseudoreversePeptideLight_helper(peptide.sequence, light_mods);

    TEST_EQUAL(light_result.first, heavy_result.sequence)
    TEST_EQUAL(light_result.second.size(), heavy_result.mods.size())
    TEST_EQUAL(light_result.second[0].location, heavy_result.mods[0].location)
  }

  // Test case 2: No modifications
  {
    OpenMS::TargetedExperiment::Peptide peptide;
    peptide.sequence = "TESTPEPTIDE";

    std::vector<OpenSwath::LightModification> light_mods;

    auto heavy_result = gen.pseudoreversePeptide_helper(peptide);
    auto light_result = gen.pseudoreversePeptideLight_helper(peptide.sequence, light_mods);

    TEST_EQUAL(light_result.first, heavy_result.sequence)
    TEST_EQUAL(light_result.second.size(), 0)
  }

  // Test case 3: Peptide ending with K (tryptic)
  {
    OpenMS::TargetedExperiment::Peptide peptide;
    peptide.sequence = "TESTPEPTIDEK";
    OpenMS::TargetedExperiment::Peptide::Modification modification;
    modification.location = 2;
    modification.unimod_id = 21;
    peptide.mods.push_back(modification);

    std::vector<OpenSwath::LightModification> light_mods;
    OpenSwath::LightModification lm;
    lm.location = 2;
    lm.unimod_id = 21;
    light_mods.push_back(lm);

    auto heavy_result = gen.pseudoreversePeptide_helper(peptide);
    auto light_result = gen.pseudoreversePeptideLight_helper(peptide.sequence, light_mods);

    TEST_EQUAL(light_result.first, heavy_result.sequence)
    // C-terminal AA should be preserved
    TEST_EQUAL(light_result.first[light_result.first.size()-1], 'K')
  }
}
END_SECTION

START_SECTION([EXTRA] shufflePeptideLight_matches_heavy)
{
  MRMDecoy gen;

  // Test case 1: Simple shuffle with modification
  {
    OpenMS::TargetedExperiment::Peptide peptide;
    peptide.sequence = "TESTPEPTIDE";
    OpenMS::TargetedExperiment::Peptide::Modification modification;
    modification.avg_mass_delta = 79.9799;
    modification.location = 2;
    modification.mono_mass_delta = 79.966331;
    modification.unimod_id = 21;
    peptide.mods.push_back(modification);

    std::vector<OpenSwath::LightModification> light_mods;
    OpenSwath::LightModification lm;
    lm.location = 2;
    lm.unimod_id = 21;
    light_mods.push_back(lm);

    // Use same seed for reproducibility
    int seed = 43;
    double identity_threshold = 0.7;

    auto heavy_result = gen.shufflePeptide(peptide, identity_threshold, seed);
    auto light_result = gen.shufflePeptideLight(peptide.sequence, light_mods, identity_threshold, seed);

    TEST_EQUAL(light_result.first, heavy_result.sequence)
    TEST_EQUAL(light_result.second.size(), heavy_result.mods.size())
    TEST_EQUAL(light_result.second[0].location, heavy_result.mods[0].location)
  }

  // Test case 2: Shuffle with multiple modifications
  {
    OpenMS::TargetedExperiment::Peptide peptide;
    peptide.sequence = "GPPSEDGPGVPPPSPR";
    OpenMS::TargetedExperiment::Peptide::Modification mod1, mod2;
    mod1.avg_mass_delta = 79.9799;
    mod1.location = 3;  // first S
    mod1.mono_mass_delta = 79.966331;
    mod1.unimod_id = 21;
    mod2.avg_mass_delta = 79.9799;
    mod2.location = 13;  // second S
    mod2.mono_mass_delta = 79.966331;
    mod2.unimod_id = 21;
    peptide.mods.push_back(mod1);
    peptide.mods.push_back(mod2);

    std::vector<OpenSwath::LightModification> light_mods;
    OpenSwath::LightModification lm1, lm2;
    lm1.location = 3;
    lm1.unimod_id = 21;
    lm2.location = 13;
    lm2.unimod_id = 21;
    light_mods.push_back(lm1);
    light_mods.push_back(lm2);

    int seed = 130;
    double identity_threshold = 0.7;

    auto heavy_result = gen.shufflePeptide(peptide, identity_threshold, seed);
    auto light_result = gen.shufflePeptideLight(peptide.sequence, light_mods, identity_threshold, seed);

    TEST_EQUAL(light_result.first, heavy_result.sequence)
    TEST_EQUAL(light_result.second.size(), heavy_result.mods.size())
    for (size_t i = 0; i < light_result.second.size(); ++i)
    {
      TEST_EQUAL(light_result.second[i].location, heavy_result.mods[i].location)
    }
  }

  // Test case 3: Shuffle with terminal K preserved
  {
    OpenMS::TargetedExperiment::Peptide peptide;
    peptide.sequence = "TESTPEPTIDEK";

    std::vector<OpenSwath::LightModification> light_mods;

    int seed = 42;
    double identity_threshold = 0.7;

    auto heavy_result = gen.shufflePeptide(peptide, identity_threshold, seed, 20);
    auto light_result = gen.shufflePeptideLight(peptide.sequence, light_mods, identity_threshold, seed, 20);

    TEST_EQUAL(light_result.first, heavy_result.sequence)
    // Terminal K should be preserved
    TEST_EQUAL(light_result.first[light_result.first.size()-1], 'K')
    TEST_EQUAL(heavy_result.sequence[heavy_result.sequence.size()-1], 'K')
  }

  // Test case 4: Shuffle with terminal R preserved
  {
    OpenMS::TargetedExperiment::Peptide peptide;
    peptide.sequence = "TESTPEPTIDER";

    std::vector<OpenSwath::LightModification> light_mods;

    int seed = 42;
    double identity_threshold = 0.7;

    auto heavy_result = gen.shufflePeptide(peptide, identity_threshold, seed, 20);
    auto light_result = gen.shufflePeptideLight(peptide.sequence, light_mods, identity_threshold, seed, 20);

    TEST_EQUAL(light_result.first, heavy_result.sequence)
    // Terminal R should be preserved
    TEST_EQUAL(light_result.first[light_result.first.size()-1], 'R')
  }

  // Test case 5: No modifications
  {
    OpenMS::TargetedExperiment::Peptide peptide;
    peptide.sequence = "TESTPEPTIDE";

    std::vector<OpenSwath::LightModification> light_mods;

    int seed = 42;
    double identity_threshold = 0.7;

    auto heavy_result = gen.shufflePeptide(peptide, identity_threshold, seed);
    auto light_result = gen.shufflePeptideLight(peptide.sequence, light_mods, identity_threshold, seed);

    TEST_EQUAL(light_result.first, heavy_result.sequence)
    TEST_EQUAL(light_result.second.size(), 0)
  }

  // Test case 6: K/R/P pattern preservation
  {
    OpenMS::TargetedExperiment::Peptide peptide;
    peptide.sequence = "KPRKPRPK";

    std::vector<OpenSwath::LightModification> light_mods;

    int seed = 130;
    double identity_threshold = 0.7;
    int max_attempts = 17;

    auto heavy_result = gen.shufflePeptide(peptide, identity_threshold, seed, max_attempts);
    auto light_result = gen.shufflePeptideLight(peptide.sequence, light_mods, identity_threshold, seed, max_attempts);

    TEST_EQUAL(light_result.first, heavy_result.sequence)
  }
}
END_SECTION

START_SECTION([EXTRA] switchKRLight_matches_heavy)
{
  MRMDecoy gen;

  // Test case 1: Switch K to R
  {
    OpenMS::TargetedExperiment::Peptide heavy_peptide;
    heavy_peptide.sequence = "TESTPEPTIDEK";
    gen.switchKR(heavy_peptide);

    std::string light_sequence = "TESTPEPTIDEK";
    MRMDecoy::switchKRLight(light_sequence);

    TEST_EQUAL(light_sequence, heavy_peptide.sequence)
    TEST_EQUAL(light_sequence[light_sequence.size()-1], 'R')
  }

  // Test case 2: Switch R to K
  {
    OpenMS::TargetedExperiment::Peptide heavy_peptide;
    heavy_peptide.sequence = "TESTPEPTIDER";
    gen.switchKR(heavy_peptide);

    std::string light_sequence = "TESTPEPTIDER";
    MRMDecoy::switchKRLight(light_sequence);

    TEST_EQUAL(light_sequence, heavy_peptide.sequence)
    TEST_EQUAL(light_sequence[light_sequence.size()-1], 'K')
  }

  // Test case 3: Non-K/R terminal - both should randomize (but may differ due to separate random calls)
  // We just check that the sequence length is preserved and terminal changed
  {
    std::string light_sequence = "TESTPEPTIDEX";
    MRMDecoy::switchKRLight(light_sequence);

    TEST_EQUAL(light_sequence.size(), 12)
    // Terminal should have changed from X to something else
    TEST_NOT_EQUAL(light_sequence[light_sequence.size()-1], 'X')
  }
}
END_SECTION

START_SECTION([EXTRA] hasCNterminalModsLight_matches_behavior)
{
  // Test that hasCNterminalModsLight_ correctly identifies terminal modifications

  // Test case 1: N-terminal modification (location = -1)
  {
    std::vector<OpenSwath::LightModification> mods;
    OpenSwath::LightModification mod;
    mod.location = -1;  // N-terminal
    mod.unimod_id = 272;  // sulfonation of N-terminus
    mods.push_back(mod);

    size_t seq_length = 11;
    bool has_terminal = MRMDecoyHelper::hasCNterminalModsLight_helper(mods, seq_length, false);
    TEST_EQUAL(has_terminal, true)
  }

  // Test case 2: C-terminal modification (location = sequence length)
  {
    std::vector<OpenSwath::LightModification> mods;
    OpenSwath::LightModification mod;
    mod.location = 11;  // C-terminal for length 11
    mod.unimod_id = 193;  // O18 label
    mods.push_back(mod);

    size_t seq_length = 11;
    bool has_terminal = MRMDecoyHelper::hasCNterminalModsLight_helper(mods, seq_length, false);
    TEST_EQUAL(has_terminal, true)
  }

  // Test case 3: Modification on C-terminal AA (location = sequence length - 1, checkCterminalAA=true)
  {
    std::vector<OpenSwath::LightModification> mods;
    OpenSwath::LightModification mod;
    mod.location = 10;  // C-terminal AA for length 11
    mod.unimod_id = 35;
    mods.push_back(mod);

    size_t seq_length = 11;
    bool has_terminal_with_check = MRMDecoyHelper::hasCNterminalModsLight_helper(mods, seq_length, true);
    bool has_terminal_without_check = MRMDecoyHelper::hasCNterminalModsLight_helper(mods, seq_length, false);

    TEST_EQUAL(has_terminal_with_check, true)
    TEST_EQUAL(has_terminal_without_check, false)
  }

  // Test case 4: Internal modification - should return false
  {
    std::vector<OpenSwath::LightModification> mods;
    OpenSwath::LightModification mod;
    mod.location = 5;  // Internal
    mod.unimod_id = 21;
    mods.push_back(mod);

    size_t seq_length = 11;
    bool has_terminal = MRMDecoyHelper::hasCNterminalModsLight_helper(mods, seq_length, false);
    TEST_EQUAL(has_terminal, false)
  }

  // Test case 5: No modifications - should return false
  {
    std::vector<OpenSwath::LightModification> mods;
    size_t seq_length = 11;
    bool has_terminal = MRMDecoyHelper::hasCNterminalModsLight_helper(mods, seq_length, false);
    TEST_EQUAL(has_terminal, false)
  }

  // Test case 6: Both N and C terminal modifications
  {
    std::vector<OpenSwath::LightModification> mods;
    OpenSwath::LightModification mod1, mod2;
    mod1.location = -1;  // N-terminal
    mod1.unimod_id = 272;
    mod2.location = 11;  // C-terminal
    mod2.unimod_id = 193;
    mods.push_back(mod1);
    mods.push_back(mod2);

    size_t seq_length = 11;
    bool has_terminal = MRMDecoyHelper::hasCNterminalModsLight_helper(mods, seq_length, false);
    TEST_EQUAL(has_terminal, true)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

#ifdef __clang__
  #pragma clang diagnostic pop
#endif