// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/TraMLFile.h>

#include <boost/assign/std/vector.hpp>
#include <boost/assign/list_of.hpp>

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/MRMDecoy.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
///////////////////////////

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wshadow"

using namespace OpenMS;
using namespace std;

START_TEST(MRMDecoy, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MRMDecoy * ptr = 0;
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

START_SECTION((std::pair<String, double> getDecoyIon(String ionid, std::map<String, std::map<String, double> >&decoy_ionseries)))
{
  MRMDecoy gen;

  OpenMS::TargetedExperiment::Peptide peptide;
  peptide.sequence = "TESTPEPTIDE";

  OpenMS::AASequence aas = TargetedExperimentHelper::getAASequence(peptide);

  double ProductMZ = 371.66692;
  double mz_threshold = 0.8;
  int precursor_charge = 2;

  MRMDecoy::IonSeries reference_ionseries = gen.getIonSeries(aas, precursor_charge);
  std::pair<String, double> targetion = gen.getTargetIon(ProductMZ, mz_threshold, reference_ionseries, 1);
  std::pair<String, double> decoyion = gen.getDecoyIon("b7^2", reference_ionseries);
  std::pair<String, double> decoyion_missing = gen.getDecoyIon("b17^2", reference_ionseries);

  TEST_EQUAL(targetion.first, "b7^2")
  TEST_REAL_SIMILAR(targetion.second, decoyion.second)
  TEST_REAL_SIMILAR(decoyion_missing.second, -1)
}

END_SECTION

START_SECTION((std::vector<std::pair<std::string::size_type, std::string> > find_all_tryptic(std::string sequence)))
{
  MRMDecoy gen;

  String sequence = "TRESTPEPTIKDE";
  std::vector<pair<std::string::size_type, std::string> > tryptic_results = gen.find_all_tryptic(sequence);
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

  OpenMS::String expected_sequence = "TETTPEPESID";
  OpenMS::Size expected_location = 8;

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
  shuffleAASequence_expected_07.sequence = "ETSTPDPEETI";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_result_07;
  shuffleAASequence_result_07 = gen.shufflePeptide(shuffleAASequence_target_sequence_07, 0.7, 42);
  TEST_EQUAL(shuffleAASequence_result_07.sequence, shuffleAASequence_expected_07.sequence)

  OpenMS::TargetedExperiment::Peptide shuffleAASequence_target_sequence_09;
  shuffleAASequence_target_sequence_09.sequence = "TESTPEPTIDE";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_expected_09;
  shuffleAASequence_expected_09.sequence = "ETSTPDPEETI";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_result_09;
  shuffleAASequence_result_09 = gen.shufflePeptide(shuffleAASequence_target_sequence_09, 0.9, 42);
  TEST_EQUAL(shuffleAASequence_result_09.sequence, shuffleAASequence_expected_09.sequence)

  OpenMS::TargetedExperiment::Peptide shuffleAASequence_target_sequence_01;
  shuffleAASequence_target_sequence_01.sequence = "TESTPEPTIDE";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_expected_01;
  shuffleAASequence_expected_01.sequence = "SIECPAPDEETTT";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_result_01;
  shuffleAASequence_result_01 = gen.shufflePeptide(shuffleAASequence_target_sequence_01, 0.2, 42, 10000);
  TEST_EQUAL(shuffleAASequence_result_01.sequence, shuffleAASequence_expected_01.sequence)

  OpenMS::TargetedExperiment::Peptide shuffleAASequence_target_sequence_00;
  shuffleAASequence_target_sequence_00.sequence = "TESTPEPTIDE";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_expected_00;
  shuffleAASequence_expected_00.sequence = "TEEDPTPDGATECIS";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_result_00;
  shuffleAASequence_result_00 = gen.shufflePeptide(shuffleAASequence_target_sequence_00, 0.0, 42, 20);
  TEST_EQUAL(shuffleAASequence_result_00.sequence, shuffleAASequence_expected_00.sequence)

  OpenMS::TargetedExperiment::Peptide shuffleAASequence_target_sequence_01b;
  shuffleAASequence_target_sequence_01b.sequence = "TESTPEPTIDE";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_expected_01b;
  shuffleAASequence_expected_01b.sequence = "TCTTPDPISEE";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_result_01b;
  shuffleAASequence_result_01b = gen.shufflePeptide(shuffleAASequence_target_sequence_01b, 0.2, 42, 10000, true);
  TEST_EQUAL(shuffleAASequence_result_01b.sequence, shuffleAASequence_expected_01b.sequence)

  OpenMS::TargetedExperiment::Peptide shuffleAASequence_target_sequence_00b;
  shuffleAASequence_target_sequence_00b.sequence = "TESTPEPTIDE";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_expected_00b;
  shuffleAASequence_expected_00b.sequence = "HGCGPDPCCHG";
  OpenMS::TargetedExperiment::Peptide shuffleAASequence_result_00b;
  shuffleAASequence_result_00b = gen.shufflePeptide(shuffleAASequence_target_sequence_00b, 0.0, 42, 2000, true);
  TEST_EQUAL(shuffleAASequence_result_00b.sequence, shuffleAASequence_expected_00b.sequence)
  // ensure that C-terminal K and R are preserved
  {
    OpenMS::TargetedExperiment::Peptide original_input;
    original_input.sequence = "TESTPEPTIDEK";
    expected_sequence = "ETSTPDPEETIK";
    OpenMS::TargetedExperiment::Peptide shuffleAASequence_result_00;
    shuffled = gen.shufflePeptide(original_input, 0.7, 42, 20);
    TEST_EQUAL(shuffled.sequence[shuffled.sequence.size() - 1], 'K')
    TEST_EQUAL(shuffled.sequence, expected_sequence)
  }

  // ensure that C-terminal K and R are preserved
  {
    OpenMS::TargetedExperiment::Peptide original_input;
    original_input.sequence = "TESTPEPTIDER";
    expected_sequence = "ETSTPDPEETIR";
    OpenMS::TargetedExperiment::Peptide shuffleAASequence_result_00;
    shuffled = gen.shufflePeptide(original_input, 0.7, 42, 20);
    TEST_EQUAL(shuffled.sequence[shuffled.sequence.size() - 1], 'R')
    TEST_EQUAL(shuffled.sequence, expected_sequence)
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

  OpenMS::String expected_sequence = "GPPEVSGPGSPPPDPR";
  OpenMS::Size expected_location_1 = 5;
  OpenMS::Size expected_location_2 = 9;

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

  OpenMS::String expected_sequence = "TETTPEPESID";

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
  OpenMS::String expected_sequence = "KPRKPRPKNL";
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
START_SECTION((void generateDecoys(OpenMS::TargetedExperiment & exp,
                                   OpenMS::TargetedExperiment & dec, String method, String decoy_tag,
                                   double identity_threshold, int max_attempts, double mz_threshold,
                                   bool theoretical, double mz_shift, bool exclude_similar,
                                   double similarity_threshold, bool remove_CNterm_mods, double precursor_mass_shift,
                                   bool enable_losses, bool skip_unannotated); ))
{
  String method = "pseudo-reverse";
  double identity_threshold = 1.0;
  Int max_attempts = 5;
  double mz_threshold = 0.8;
  double mz_shift = 20;
  String decoy_tag = "DECOY_";
  Int min_transitions = 2;
  Int max_transitions = 6;
  bool theoretical = 1;
  bool exclude_similar = 1;
  bool remove_CNterminal_mods = false;
  double similarity_threshold = 0.1;
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
  TEST_EQUAL(targeted_exp.getPeptides().size(), 13)
  TEST_EQUAL(targeted_exp.getTransitions().size(), 33)
  decoys.generateDecoys(targeted_exp, targeted_decoy, method, decoy_tag, identity_threshold, max_attempts, mz_threshold, theoretical, mz_shift, exclude_similar, similarity_threshold, remove_CNterminal_mods, 0.1, 1, 0);  traml.store(test, targeted_decoy);

  TEST_FILE_EQUAL(test.c_str(), OPENMS_GET_TEST_DATA_PATH(out))
}

END_SECTION

START_SECTION(void restrictTransitions(OpenMS::TargetedExperiment& exp, int min_transitions, int max_transitions))
{
  // see above
  NOT_TESTABLE
}

END_SECTION

START_SECTION((extra))
{
  MRMDecoy gen;

  AASequence target_sequence = AASequence::fromString("ADSTGTLVITDPTR(UniMod:267)");
  AASequence decoy_sequence = AASequence::fromString("ALDSTTGVDTTPIR(UniMod:267)");
  MRMDecoy::IonSeries target_ionseries = gen.getIonSeries(target_sequence, 2);
  MRMDecoy::IonSeries decoy_ionseries = gen.getIonSeries(decoy_sequence, 2);

  {
    double target_mz = 924.539;
    double mz_threshold = 0.1;

    std::pair<String, double> targetion = gen.getTargetIon(target_mz, mz_threshold, target_ionseries, 1);
    std::pair<String, double> decoyion = gen.getDecoyIon(targetion.first, decoy_ionseries);

    TEST_EQUAL(targetion.first, "y8^1")
    TEST_REAL_SIMILAR(targetion.second, 924.539)

    TEST_EQUAL(decoyion.first, "y8^1")
    TEST_REAL_SIMILAR(decoyion.second, 868.47682)
  }

  {
    double target_mz = 1082.608;
    double mz_threshold = 0.8;

    std::pair<String, double> targetion = gen.getTargetIon(target_mz, mz_threshold, target_ionseries, 1);
    std::pair<String, double> decoyion = gen.getDecoyIon(targetion.first, decoy_ionseries);

    TEST_EQUAL(targetion.first, "y10^1")
    TEST_REAL_SIMILAR(targetion.second, 1082.60856)

    TEST_EQUAL(decoyion.first, "y10^1")
    TEST_REAL_SIMILAR(decoyion.second, 1070.57217)
  }

  {
    AASequence peptide = AASequence::fromString("KVGLDPSQLPVGENGIV");
    MRMDecoy::IonSeries target_ionseries = gen.getIonSeries(peptide, 2);

    double target_mz_1 = 737.394;
    double target_mz_2 = 793.936;
    double mz_threshold = 0.8;

    std::pair<String, double> targetion1 = gen.getTargetIon(target_mz_1, mz_threshold, target_ionseries, 1);

    TEST_EQUAL(targetion1.first, "b15-18^2")

    std::pair<String, double> targetion2 = gen.getTargetIon(target_mz_2, mz_threshold, target_ionseries, 1);

    TEST_EQUAL(targetion2.first, "b16-18^2")
  }
  {
    // Neutral loss of H2ONH3, H2OH2O, CH3NO, HCOOH
    AASequence peptide = AASequence::fromString("AAAAAAAAAPAAAATAPTTAATTAATAAQ");
    MRMDecoy::IonSeries target_ionseries = gen.getIonSeries(peptide, 3);

    double mz_threshold = 0.05;

    std::pair<String, double> targetion3 = gen.getTargetIon(510.2660, mz_threshold, target_ionseries, 1);

    TEST_EQUAL(targetion3.first, "b20-35^3")

    std::pair<String, double> targetion4 = gen.getTargetIon(678.3557, mz_threshold, target_ionseries, 1);

    TEST_EQUAL(targetion4.first, "b18-36^2")
  }
  {
    // Neutral loss of H2O, NH3
    AASequence peptide = AASequence::fromString("AAAAAAALQAK");
    MRMDecoy::IonSeries target_ionseries = gen.getIonSeries(peptide, 2);

    double mz_threshold = 0.05;

    std::pair<String, double> targetion1 = gen.getTargetIon(328.1980, mz_threshold, target_ionseries, 1);

    TEST_EQUAL(targetion1.first, "y3-18^1")

    std::pair<String, double> targetion2 = gen.getTargetIon(329.1833, mz_threshold, target_ionseries, 1);

    TEST_EQUAL(targetion2.first, "y3-17^1")

    std::pair<String, double> targetion3 = gen.getTargetIon(885.5169, mz_threshold, target_ionseries, 1);

    TEST_EQUAL(targetion3.first, "y10^1")
  }
  {
    // Neutral loss of NH3NH3
    AASequence peptide = AASequence::fromString("AALQEELQLC(UniMod:4)K");
    MRMDecoy::IonSeries target_ionseries = gen.getIonSeries(peptide, 2);

    double mz_threshold = 0.05;

    std::pair<String, double> targetion1 = gen.getTargetIon(240.0975, mz_threshold, target_ionseries, 1);

    TEST_EQUAL(targetion1.first, "b5-34^2")
  }
  {
    // Neutral loss (oxidation) of methionine only
    AASequence peptide = AASequence::fromString("AFADALEVIPMALSENSGM(UniMod:35)NPIQTMTEVR");
    MRMDecoy::IonSeries target_ionseries = gen.getIonSeries(peptide, 3);

    double mz_threshold = 0.05;

    std::pair<String, double> targetion1 = gen.getTargetIon(933.4456, mz_threshold, target_ionseries, 1);

    TEST_EQUAL(targetion1.first, "y26-64^3")

    std::pair<String, double> targetion2 = gen.getTargetIon(1112.5415, mz_threshold, target_ionseries, 1);

    TEST_EQUAL(targetion2.first, "b22-64^2")
  }
  {
    // Neutral losses (phospho) of serine and threonine only
    AASequence peptide = AASequence::fromString("AFADALEVIPMALSENS(Phospho)GM(UniMod:35)NPIQTMTEVR");
    MRMDecoy::IonSeries target_ionseries = gen.getIonSeries(peptide, 3);

    double mz_threshold = 0.05;

    std::pair<String, double> targetion3 = gen.getTargetIon(1144.541, mz_threshold, target_ionseries, 1);

    TEST_EQUAL(targetion3.first, "b22-80^2")

    std::pair<String, double> targetion4 = gen.getTargetIon(1135.535, mz_threshold, target_ionseries, 1);
    TEST_EQUAL(targetion4.first, "b22-98^2")
  }
}

END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
  END_TEST



/*
Tests of unknown methods:
  - 'vector find_all_tryptic(string sequence)'
Missing tests:
  - 'std::pair<String, double> getTargetIon(double ProductMZ, double mz_threshold, std::map< String, std::map< String, double > > target_ionseries)'
  - 'std::map<String, std::map<String, double> > getIonSeries(AASequence sequence, int precursor_charge)'
  - 'std::vector<std::pair<std::string::size_type, std::string> > find_all_tryptic(std::string sequence)'
*/


#pragma clang diagnostic pop
