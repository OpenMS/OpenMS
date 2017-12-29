// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

#include <OpenMS/test_config.h>
#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/MRMDecoy.h>
///////////////////////////

#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

#include <boost/assign/std/vector.hpp>
#include <boost/assign/list_of.hpp>

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wshadow"

using namespace OpenMS;
using namespace std;

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
  MRMDecoy gen;

  String sequence = "TRESTPEPTIKDE";
  std::vector<pair<std::string::size_type, std::string> > tryptic_results = gen.findFixedResidues(sequence);
  std::vector<pair<std::string::size_type, std::string> > tryptic_control = boost::assign::list_of(std::make_pair(1, "R")) (std::make_pair(5, "P")) (std::make_pair(7, "P")) (std::make_pair(10, "K"));

  for (Size i = 0; i < tryptic_results.size(); i++)
  {
    pair<std::string::size_type, std::string> result = tryptic_results[i];
    pair<std::string::size_type, std::string> control = tryptic_control[i];
    TEST_EQUAL(result.first, control.first)
    TEST_EQUAL(result.second, control.second)
  }
}

END_SECTION

START_SECTION((std::vector<std::pair<std::string::size_type, std::string> > findFixedAndTermResidues(std::string sequence)))
{
  MRMDecoy gen;

  String sequence = "TRESTPEPTIKDE";
  std::vector<pair<std::string::size_type, std::string> > tryptic_results = gen.findFixedAndTermResidues(sequence);
  std::vector<pair<std::string::size_type, std::string> > tryptic_control = boost::assign::list_of(std::make_pair(0, "T")) (std::make_pair(1, "R")) (std::make_pair(5, "P")) (std::make_pair(7, "P")) (std::make_pair(10, "K")) (std::make_pair(12, "E"));

  for (Size i = 0; i < tryptic_results.size(); i++)
  {
    pair<std::string::size_type, std::string> result = tryptic_results[i];
    pair<std::string::size_type, std::string> control = tryptic_control[i];
    TEST_EQUAL(result.first, control.first)
    TEST_EQUAL(result.second, control.second)
  }
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
  TEST_EQUAL(pseudoreverse.sequence, expected_sequence)
  TEST_EQUAL(pseudoreverse.mods.size(), 1)
  TEST_EQUAL(pseudoreverse.mods[0].location, expected_location)

  OpenMS::TargetedExperiment::Peptide pseudoreverseAASequence_target_sequence;
  pseudoreverseAASequence_target_sequence.sequence = "TESTPEPTIDE";
  OpenMS::TargetedExperiment::Peptide pseudoreverseAASequence_expected;
  pseudoreverseAASequence_expected.sequence = "DITPEPTSETE";
  OpenMS::TargetedExperiment::Peptide pseudoreverseAASequence_result;
  pseudoreverseAASequence_result = gen.pseudoreversePeptide(pseudoreverseAASequence_target_sequence);
  TEST_EQUAL(pseudoreverseAASequence_result.sequence, pseudoreverseAASequence_expected.sequence)
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
  TEST_EQUAL(reverse.sequence, expected_sequence)
  TEST_EQUAL(reverse.mods.size(), 1)
  TEST_EQUAL(reverse.mods[0].location, expected_location)

  OpenMS::TargetedExperiment::Peptide reverseAASequence_target_sequence;
  reverseAASequence_target_sequence.sequence = "TESTPEPTIDE";
  OpenMS::TargetedExperiment::Peptide reverseAASequence_expected;
  reverseAASequence_expected.sequence = "EDITPEPTSET";
  OpenMS::TargetedExperiment::Peptide reverseAASequence_result;
  reverseAASequence_result = gen.reversePeptide(reverseAASequence_target_sequence);
  TEST_EQUAL(reverseAASequence_result.sequence, reverseAASequence_expected.sequence)
}

END_SECTION

/// Public methods



START_SECTION((void generateDecoys(OpenMS::TargetedExperiment& exp, OpenMS::TargetedExperiment& dec,
                                   String method, String decoy_tag, int max_attempts, double identity_threshold,
                                   double precursor_mz_shift, double product_mz_shift, double product_mz_threshold,
                                   std::vector<String> fragment_types, std::vector<size_t> fragment_charges,
                                   bool enable_specific_losses, bool enable_unspecific_losses, int round_decPow); ))
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

  TEST_FILE_EQUAL(test.c_str(), OPENMS_GET_TEST_DATA_PATH(out))
}

END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

#pragma clang diagnostic pop
