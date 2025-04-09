// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger $
// --------------------------------------------------------------------------
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/MRMAssay.h>
#include <OpenMS/FORMAT/TraMLFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MRMAssay, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MRMAssay * ptr = nullptr;
MRMAssay* nullPointer = nullptr;

class MRMAssay_test :
  public MRMAssay
{
public:
  std::vector<std::string> getMatchingPeptidoforms_test(const double fragment_ion, std::vector<std::pair<double, std::string> >& ions, const double mz_threshold)
  {
    return getMatchingPeptidoforms_(fragment_ion, ions, mz_threshold);
  }

  int getSwath_test(const std::vector<std::pair<double, double> >& swathes, const double precursor_mz)
  {
    return getSwath_(swathes, precursor_mz);
  }

  bool isInSwath_test(const std::vector<std::pair<double, double> >& swathes, const double precursor_mz, const double product_mz)
  {
    return isInSwath_(swathes, precursor_mz, product_mz);
  }

  std::string getRandomSequence_test(int sequence_size, boost::variate_generator<boost::mt19937&, boost::uniform_int<> > pseudoRNG)
  {
    return getRandomSequence_(sequence_size, pseudoRNG);
  }

  std::vector<std::vector<size_t> > nchoosekcombinations_test(std::vector<size_t> n, size_t k)
  {
    return nchoosekcombinations_(n, k);
  }

  std::vector<OpenMS::AASequence> addModificationsSequences_test(std::vector<OpenMS::AASequence> sequences, std::vector<std::vector<size_t> > mods_combs, OpenMS::String modification)
  {
    return addModificationsSequences_(sequences, mods_combs, modification);
  }

  std::vector<OpenMS::AASequence> generateTheoreticalPeptidoforms_test(OpenMS::AASequence sequence)
  {
    return generateTheoreticalPeptidoforms_(sequence);
  }

  std::vector<OpenMS::AASequence> generateTheoreticalPeptidoformsDecoy_test(OpenMS::AASequence sequence, OpenMS::AASequence decoy_sequence)
  {
    return generateTheoreticalPeptidoformsDecoy_(sequence, decoy_sequence);
  }

};

START_SECTION(MRMAssay())
{
  ptr = new MRMAssay();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~MRMAssay())
{
  delete ptr;
}

END_SECTION

START_SECTION(std::vector<std::string> MRMAssay::getMatchingPeptidoforms_(const double fragment_ion, std::vector<std::pair<double, std::string> >& ions, const double mz_threshold); )
{
  MRMAssay_test mrma;

  std::vector<std::pair<double, std::string> > ions;
  ions.push_back(std::make_pair(100.00, "PEPTIDEK"));
  ions.push_back(std::make_pair(100.01, "PEPTIDEK"));
  ions.push_back(std::make_pair(100.10, "PEPT(UniMod:21)IDEK"));
  ions.push_back(std::make_pair(100.12, "PEPTIDEK"));
  ions.push_back(std::make_pair(100.11, "PEPTIDEK"));

  std::vector<std::string> isoforms1 = mrma.getMatchingPeptidoforms_test(100.06, ions, 0.03);
  std::vector<std::string> isoforms2 = mrma.getMatchingPeptidoforms_test(100.06, ions, 0.06);

  TEST_EQUAL(isoforms1.size(), 0)

  TEST_EQUAL(isoforms2.size(), 2)
  TEST_EQUAL(isoforms2[0], "PEPT(UniMod:21)IDEK")
  TEST_EQUAL(isoforms2[1], "PEPTIDEK")
}

END_SECTION

START_SECTION(int MRMAssay::getSwath_(const std::vector<std::pair<double, double> > swathes, const double precursor_mz))
{
  MRMAssay_test mrma;
  std::vector<std::pair<double, double> > swathes;
  swathes.push_back(std::make_pair(400, 425));
  swathes.push_back(std::make_pair(424, 450));
  swathes.push_back(std::make_pair(449, 475));
  swathes.push_back(std::make_pair(474, 500));
  swathes.push_back(std::make_pair(499, 525));
  swathes.push_back(std::make_pair(524, 550));
  swathes.push_back(std::make_pair(549, 575));
  swathes.push_back(std::make_pair(574, 600));
  swathes.push_back(std::make_pair(599, 625));
  swathes.push_back(std::make_pair(624, 650));
  swathes.push_back(std::make_pair(649, 675));
  swathes.push_back(std::make_pair(674, 700));
  swathes.push_back(std::make_pair(699, 725));
  swathes.push_back(std::make_pair(724, 750));
  swathes.push_back(std::make_pair(749, 775));
  swathes.push_back(std::make_pair(774, 800));
  swathes.push_back(std::make_pair(799, 825));
  swathes.push_back(std::make_pair(824, 850));
  swathes.push_back(std::make_pair(849, 875));
  swathes.push_back(std::make_pair(874, 900));
  swathes.push_back(std::make_pair(899, 925));
  swathes.push_back(std::make_pair(924, 950));
  swathes.push_back(std::make_pair(949, 975));
  swathes.push_back(std::make_pair(974, 1000));
  swathes.push_back(std::make_pair(999, 1025));
  swathes.push_back(std::make_pair(1024, 1050));
  swathes.push_back(std::make_pair(1049, 1075));
  swathes.push_back(std::make_pair(1074, 1100));
  swathes.push_back(std::make_pair(1099, 1125));
  swathes.push_back(std::make_pair(1124, 1150));
  swathes.push_back(std::make_pair(1149, 1175));
  swathes.push_back(std::make_pair(1174, 1200));

  TEST_EQUAL(mrma.getSwath_test(swathes, 427.229959), 1)
  TEST_EQUAL(mrma.getSwath_test(swathes, 449.0), 2)
  TEST_EQUAL(mrma.getSwath_test(swathes, 449.229959), 2)
  TEST_EQUAL(mrma.getSwath_test(swathes, 685.8547721), 11)
  TEST_EQUAL(mrma.getSwath_test(swathes, 1685.8547721), -1)
  TEST_EQUAL(mrma.getSwath_test(swathes, -41.1), -1)
}

END_SECTION

START_SECTION(bool MRMAssay::isInSwath_(const std::vector<std::pair<double, double> > swathes, const double precursor_mz, const double product_mz))
{
  MRMAssay_test mrma;
  std::vector<std::pair<double, double> > swathes;
  swathes.push_back(std::make_pair(400, 425));
  swathes.push_back(std::make_pair(424, 450));
  swathes.push_back(std::make_pair(449, 475));
  swathes.push_back(std::make_pair(474, 500));
  swathes.push_back(std::make_pair(499, 525));
  swathes.push_back(std::make_pair(524, 550));
  swathes.push_back(std::make_pair(549, 575));
  swathes.push_back(std::make_pair(574, 600));
  swathes.push_back(std::make_pair(599, 625));
  swathes.push_back(std::make_pair(624, 650));
  swathes.push_back(std::make_pair(649, 675));
  swathes.push_back(std::make_pair(674, 700));
  swathes.push_back(std::make_pair(699, 725));
  swathes.push_back(std::make_pair(724, 750));
  swathes.push_back(std::make_pair(749, 775));
  swathes.push_back(std::make_pair(774, 800));
  swathes.push_back(std::make_pair(799, 825));
  swathes.push_back(std::make_pair(824, 850));
  swathes.push_back(std::make_pair(849, 875));
  swathes.push_back(std::make_pair(874, 900));
  swathes.push_back(std::make_pair(899, 925));
  swathes.push_back(std::make_pair(924, 950));
  swathes.push_back(std::make_pair(949, 975));
  swathes.push_back(std::make_pair(974, 1000));
  swathes.push_back(std::make_pair(999, 1025));
  swathes.push_back(std::make_pair(1024, 1050));
  swathes.push_back(std::make_pair(1049, 1075));
  swathes.push_back(std::make_pair(1074, 1100));
  swathes.push_back(std::make_pair(1099, 1125));
  swathes.push_back(std::make_pair(1124, 1150));
  swathes.push_back(std::make_pair(1149, 1175));
  swathes.push_back(std::make_pair(1174, 1200));

  TEST_EQUAL(mrma.isInSwath_test(swathes, 685.8547721, 427.229959), false)
  TEST_EQUAL(mrma.isInSwath_test(swathes, 685.8547721, 689), true)
}

END_SECTION

START_SECTION(std::string MRMAssay::getRandomSequence_(int sequence_size, boost::variate_generator<boost::mt19937&, boost::uniform
                                                                                                   _int<> > pseudoRNG))
{
  MRMAssay_test mrma;

  boost::mt19937 generator(42);
  boost::uniform_int<> uni_dist;
  boost::variate_generator<boost::mt19937&, boost::uniform_int<> > pseudoRNG(generator, uni_dist);

  std::string sequence1 = mrma.getRandomSequence_test(10, pseudoRNG);

  TEST_EQUAL(sequence1, "CHLNHHQQNE");
}

END_SECTION

START_SECTION(std::vector<std::vector<size_t> > MRMAssay::nchoosekcombinations_(std::vector<size_t> n, size_t k))
{
  MRMAssay_test mrma;
  std::vector<size_t> n;
  size_t k = 5;

  for (size_t i = 1; i <= 16; i++)
  {
    n.push_back(i);
  }

  TEST_EQUAL(mrma.nchoosekcombinations_test(n, k).size(), 4368)
  TEST_EQUAL(mrma.nchoosekcombinations_test(n, k)[0].size(), 5)
}

END_SECTION

START_SECTION(std::vector<std::vector<size_t> > MRMAssay::nchoosekcombinations_(std::vector<size_t> n, size_t k))
{
  MRMAssay_test mrma;
  std::vector<size_t> n;
  size_t k = 3;

  for (size_t i = 1; i <= 5; i++)
  {
    n.push_back(i);
  }

  std::vector<std::vector<size_t> > res = mrma.nchoosekcombinations_test(n, k);

  TEST_EQUAL(res.size(), 10)
  TEST_EQUAL(res[0].size(), 3)

  TEST_EQUAL(res[0][0], 1)
  TEST_EQUAL(res[0][1], 2)
  TEST_EQUAL(res[0][2], 3)

  TEST_EQUAL(res[1][0], 1)
  TEST_EQUAL(res[1][1], 2)
  TEST_EQUAL(res[1][2], 4)

  TEST_EQUAL(res[2][0], 1)
  TEST_EQUAL(res[2][1], 2)
  TEST_EQUAL(res[2][2], 5)

  TEST_EQUAL(res[3][0], 1)
  TEST_EQUAL(res[3][1], 3)
  TEST_EQUAL(res[3][2], 4)

  TEST_EQUAL(res[4][0], 1)
  TEST_EQUAL(res[4][1], 3)
  TEST_EQUAL(res[4][2], 5)

  TEST_EQUAL(res[5][0], 1)
  TEST_EQUAL(res[5][1], 4)
  TEST_EQUAL(res[5][2], 5)

  TEST_EQUAL(res[6][0], 2)
  TEST_EQUAL(res[6][1], 3)
  TEST_EQUAL(res[6][2], 4)

  TEST_EQUAL(res[7][0], 2)
  TEST_EQUAL(res[7][1], 3)
  TEST_EQUAL(res[7][2], 5)

  TEST_EQUAL(res[8][0], 2)
  TEST_EQUAL(res[8][1], 4)
  TEST_EQUAL(res[8][2], 5)

  TEST_EQUAL(res[9][0], 3)
  TEST_EQUAL(res[9][1], 4)
  TEST_EQUAL(res[9][2], 5)
}

END_SECTION

START_SECTION(std::vector<OpenMS::AASequence> MRMAssay::addModificationsSequences_(std::vector<OpenMS::AASequence> sequences, std::vector<std::vector<size_t> > mods_combs, OpenMS::String modification))
{
  MRMAssay_test mrma;

  AASequence sequence = AASequence::fromString("PEPTDIEK");
  std::vector<OpenMS::AASequence> sequences;
  sequences.push_back(sequence);

  std::vector<size_t> no;
  no.push_back(1);
  no.push_back(3);
  no.push_back(5);
  no.push_back(8);
  no.push_back(9);

  std::vector<std::vector<size_t> > mods_combs_o = mrma.nchoosekcombinations_test(no, 1);

  sequences = mrma.addModificationsSequences_test(sequences, mods_combs_o, String("Oxidation"));

  std::vector<size_t> np;
  np.push_back(4);
  np.push_back(5);
  np.push_back(8);

  std::vector<std::vector<size_t> > mods_combs_p = mrma.nchoosekcombinations_test(np, 1);

  sequences = mrma.addModificationsSequences_test(sequences, mods_combs_p, String("Phospho"));

  TEST_EQUAL(sequences.size(), 10)
  TEST_EQUAL(sequences[0].toString(), String("P(Oxidation)EPT(Phospho)DIEK"));
  TEST_EQUAL(sequences[1].toString(), String("P(Oxidation)EPTD(Phospho)IEK"));
  TEST_EQUAL(sequences[2].toString(), String("P(Oxidation)EPTDIEK(Phospho)"));
  TEST_EQUAL(sequences[3].toString(), String("PEP(Oxidation)T(Phospho)DIEK"));
  TEST_EQUAL(sequences[4].toString(), String("PEP(Oxidation)TD(Phospho)IEK"));
  TEST_EQUAL(sequences[5].toString(), String("PEP(Oxidation)TDIEK(Phospho)"));
  TEST_EQUAL(sequences[6].toString(), String("PEPT(Phospho)D(Oxidation)IEK"));
  TEST_EQUAL(sequences[7].toString(), String("PEPTD(Oxidation)IEK(Phospho)"));
  TEST_EQUAL(sequences[8].toString(), String("PEPT(Phospho)DIEK(Oxidation)"));
  TEST_EQUAL(sequences[9].toString(), String("PEPTD(Phospho)IEK(Oxidation)"));

  std::vector<std::string> sequence_list {};
  for (std::vector<OpenMS::AASequence>::const_iterator sq_it = sequences.begin(); sq_it != sequences.end(); ++sq_it)
  {
	  OpenMS::AASequence temp_sequence = *sq_it;
	  sequence_list.push_back(temp_sequence.toString());
  }

  bool check_terminal_mod_present = (std::find(sequence_list.begin(), sequence_list.end(), "PEPT(Phospho)DIEK.(Oxidation)") != sequence_list.end());
  TEST_EQUAL(check_terminal_mod_present, false)
}

END_SECTION

START_SECTION(std::vector<OpenMS::AASequence> MRMAssay::generateTheoreticalPeptidoforms_(OpenMS::AASequence sequence))
{
  MRMAssay_test mrma;

  std::vector<AASequence> sequences = mrma.generateTheoreticalPeptidoforms_test(AASequence::fromString(".(Acetyl)PEPT(Phospho)DIEK"));

  TEST_EQUAL(sequences.size(), 13)
  TEST_EQUAL(sequences[0], AASequence::fromString(".(Acetyl)PE(Phospho)PTDIEK"));
  TEST_EQUAL(sequences[1], AASequence::fromString(".(Acetyl)PEPT(Phospho)DIEK"));
  TEST_EQUAL(sequences[2], AASequence::fromString(".(Acetyl)PEPTD(Phospho)IEK"));
  TEST_EQUAL(sequences[3], AASequence::fromString(".(Acetyl)PEPTDIE(Phospho)K"));
  TEST_EQUAL(sequences[4], AASequence::fromString(".(Acetyl)PEPTDIEK(Phospho)"));
  TEST_EQUAL(sequences[5], AASequence::fromString("PE(Phospho)PT(Acetyl)DIEK"));
  TEST_EQUAL(sequences[6], AASequence::fromString("PEPT(Acetyl)D(Phospho)IEK"));
  TEST_EQUAL(sequences[7], AASequence::fromString("PEPT(Acetyl)DIE(Phospho)K"));
  TEST_EQUAL(sequences[8], AASequence::fromString("PEPT(Acetyl)DIEK(Phospho)"));
  TEST_EQUAL(sequences[9], AASequence::fromString("PE(Phospho)PTDIEK(Acetyl)"));
  TEST_EQUAL(sequences[10], AASequence::fromString("PEPT(Phospho)DIEK(Acetyl)"));
  TEST_EQUAL(sequences[11], AASequence::fromString("PEPTD(Phospho)IEK(Acetyl)"));
  TEST_EQUAL(sequences[12], AASequence::fromString("PEPTDIE(Phospho)K(Acetyl)"));
}

END_SECTION

START_SECTION(std::vector<OpenMS::AASequence> MRMAssay::generateTheoreticalPeptidoformsDecoy_(OpenMS::AASequence sequence, OpenMS::AASequence decoy_sequence))
{
  MRMAssay_test mrma;

  std::vector<AASequence> sequences = mrma.generateTheoreticalPeptidoformsDecoy_test(AASequence::fromString(".(Acetyl)PEPT(Phospho)DIEK"), AASequence::fromString("PESTDIEK"));

  TEST_EQUAL(sequences.size(), 13)
  TEST_EQUAL(sequences[0], AASequence::fromString(".(Acetyl)PE(Phospho)STDIEK"));
  TEST_EQUAL(sequences[1], AASequence::fromString(".(Acetyl)PEST(Phospho)DIEK"));
  TEST_EQUAL(sequences[2], AASequence::fromString(".(Acetyl)PESTD(Phospho)IEK"));
  TEST_EQUAL(sequences[3], AASequence::fromString(".(Acetyl)PESTDIE(Phospho)K"));
  TEST_EQUAL(sequences[4], AASequence::fromString(".(Acetyl)PESTDIEK(Phospho)"));
  TEST_EQUAL(sequences[5], AASequence::fromString("PE(Phospho)ST(Acetyl)DIEK"));
  TEST_EQUAL(sequences[6], AASequence::fromString("PEST(Acetyl)D(Phospho)IEK"));
  TEST_EQUAL(sequences[7], AASequence::fromString("PEST(Acetyl)DIE(Phospho)K"));
  TEST_EQUAL(sequences[8], AASequence::fromString("PEST(Acetyl)DIEK(Phospho)"));
  TEST_EQUAL(sequences[9], AASequence::fromString("PE(Phospho)STDIEK(Acetyl)"));
  TEST_EQUAL(sequences[10], AASequence::fromString("PEST(Phospho)DIEK(Acetyl)"));
  TEST_EQUAL(sequences[11], AASequence::fromString("PESTD(Phospho)IEK(Acetyl)"));
  TEST_EQUAL(sequences[12], AASequence::fromString("PESTDIE(Phospho)K(Acetyl)"));
}

END_SECTION

START_SECTION(void reannotateTransitions(OpenMS::TargetedExperiment& exp, double precursor_mz_threshold, double product_mz_threshold, std::vector<String> fragment_types, std::vector<size_t> fragment_charges, bool enable_reannotation, bool enable_specific_losses, bool enable_specific_losses))
{
  TraMLFile traml;
  TargetedExperiment targeted_exp;
  String in = "MRMAssay_reannotateTransitions_input.TraML";
  traml.load(OPENMS_GET_TEST_DATA_PATH(in), targeted_exp);
  MRMAssay mrma;

  double precursor_mz_threshold1 = 0.05;
  double product_mz_threshold1 = 0.05;
  std::vector<String> fragment_types1;
  fragment_types1.push_back(String("y"));
  std::vector<size_t> fragment_charges1;
  fragment_charges1.push_back(2);
  bool enable_losses1 = false;

  String out1 = "MRMAssay_reannotateTransitions_output_1.TraML";

  TargetedExperiment targeted_exp1 = targeted_exp;

  mrma.reannotateTransitions(targeted_exp1, precursor_mz_threshold1,
      product_mz_threshold1, fragment_types1, fragment_charges1,
      enable_losses1, enable_losses1);

  String test1;
  NEW_TMP_FILE(test1);
  traml.store(test1, targeted_exp1);

  TEST_FILE_SIMILAR(test1.c_str(), OPENMS_GET_TEST_DATA_PATH(out1))

  double precursor_mz_threshold2 = 0.05;
  double product_mz_threshold2 = 0.05;
  std::vector<String> fragment_types2;
  fragment_types2.push_back(String("y"));
  fragment_types2.push_back(String("b"));
  std::vector<size_t> fragment_charges2;
  fragment_charges2.push_back(2);
  fragment_charges2.push_back(3);
  bool enable_losses2 = true;

  String out2 = "MRMAssay_reannotateTransitions_output_2.TraML";

  TargetedExperiment targeted_exp2 = targeted_exp;

  mrma.reannotateTransitions(targeted_exp2, precursor_mz_threshold2, product_mz_threshold2, fragment_types2, fragment_charges2, enable_losses2, enable_losses2);

  String test2;
  NEW_TMP_FILE(test2);
  traml.store(test2, targeted_exp2);

  TEST_FILE_SIMILAR(test2.c_str(), OPENMS_GET_TEST_DATA_PATH(out2))
}

END_SECTION

START_SECTION(void restrictTransitions(OpenMS::TargetedExperiment& exp, double lower_mz_limit, double upper_mz_limit, std::vector<std::pair<double, double> > swathes))
{
  std::vector<std::pair<double, double> > swathes;
  swathes.push_back(std::make_pair(400, 425));
  swathes.push_back(std::make_pair(424, 450));
  swathes.push_back(std::make_pair(449, 475));
  swathes.push_back(std::make_pair(474, 500));
  swathes.push_back(std::make_pair(499, 525));
  swathes.push_back(std::make_pair(524, 550));
  swathes.push_back(std::make_pair(549, 575));
  swathes.push_back(std::make_pair(574, 600));
  swathes.push_back(std::make_pair(599, 625));
  swathes.push_back(std::make_pair(624, 650));
  swathes.push_back(std::make_pair(649, 675));
  swathes.push_back(std::make_pair(674, 700));
  swathes.push_back(std::make_pair(699, 725));
  swathes.push_back(std::make_pair(724, 750));
  swathes.push_back(std::make_pair(749, 775));
  swathes.push_back(std::make_pair(774, 800));
  swathes.push_back(std::make_pair(799, 825));
  swathes.push_back(std::make_pair(824, 850));
  swathes.push_back(std::make_pair(849, 875));
  swathes.push_back(std::make_pair(874, 900));
  swathes.push_back(std::make_pair(899, 925));
  swathes.push_back(std::make_pair(924, 950));
  swathes.push_back(std::make_pair(949, 975));
  swathes.push_back(std::make_pair(974, 1000));
  swathes.push_back(std::make_pair(999, 1025));
  swathes.push_back(std::make_pair(1024, 1050));
  swathes.push_back(std::make_pair(1049, 1075));
  swathes.push_back(std::make_pair(1074, 1100));
  swathes.push_back(std::make_pair(1099, 1125));
  swathes.push_back(std::make_pair(1124, 1150));
  swathes.push_back(std::make_pair(1149, 1175));
  swathes.push_back(std::make_pair(1174, 1200));

  TraMLFile traml;
  TargetedExperiment targeted_exp;
  String in = "MRMAssay_restrictTransitions_input.TraML";
  traml.load(OPENMS_GET_TEST_DATA_PATH(in), targeted_exp);
  MRMAssay mrma;

  double lower_mz_limit = 400;
  double upper_mz_limit = 2000;

  String out1 = "MRMAssay_restrictTransitions_output.TraML";

  TargetedExperiment targeted_exp1 = targeted_exp;

  mrma.restrictTransitions(targeted_exp1, lower_mz_limit, upper_mz_limit, swathes);

  String test1;
  NEW_TMP_FILE(test1);
  traml.store(test1, targeted_exp1);

  TEST_FILE_SIMILAR(test1.c_str(), OPENMS_GET_TEST_DATA_PATH(out1))

}

END_SECTION

START_SECTION(void detectingTransitions(OpenMS::TargetedExperiment& exp, int min_transitions, int max_transitions))
{
  TraMLFile traml;
  TargetedExperiment targeted_exp;
  String in = "MRMAssay_detectingTransitions_input.TraML";
  traml.load(OPENMS_GET_TEST_DATA_PATH(in), targeted_exp);
  MRMAssay mrma;

  int min_transitions = 4;
  int max_transitions = 6;

  String out1 = "MRMAssay_detectingTransitions_output.TraML";

  TargetedExperiment targeted_exp1 = targeted_exp;

  mrma.detectingTransitions(targeted_exp1, min_transitions, max_transitions);

  String test1;
  NEW_TMP_FILE(test1);
  traml.store(test1, targeted_exp1);

  TEST_FILE_SIMILAR(test1.c_str(), OPENMS_GET_TEST_DATA_PATH(out1))

}

END_SECTION

START_SECTION(void filterMinMaxTransitionsCompound(OpenMS::TargetedExperiment& exp, int min_transitions, int max_transitions))
{
  TraMLFile traml;
  TargetedExperiment targeted_exp;
  String in = "MRMAssay_detectingTransistionCompound_input.TraML";
  traml.load(OPENMS_GET_TEST_DATA_PATH(in), targeted_exp);
  MRMAssay mrma;

  int min_transitions = 3;
  int max_transitions = 6;

  String out1 = "MRMAssay_detectingTransitionCompound_output.TraML";

  TargetedExperiment targeted_exp1 = targeted_exp;

  mrma.filterMinMaxTransitionsCompound(targeted_exp1, min_transitions, max_transitions);

  String test1;
  NEW_TMP_FILE(test1);
  traml.store(test1, targeted_exp1);

  TEST_FILE_SIMILAR(test1.c_str(), OPENMS_GET_TEST_DATA_PATH(out1))

}

END_SECTION

START_SECTION(void uisTransitions(OpenMS::TargetedExperiment& exp, std::vector<String> fragment_types, std::vector<size_t> fragment_charges, bool enable_specific_losses, bool enable_unspecific_losses, double mz_threshold, std::vector<std::pair<double, double> > swathes, int round_decPow, size_t max_num_alternative_localizations, double shuffle_identity_threshold, int shuffle_max_attempts, int shuffle_seed))
{
  std::vector<std::pair<double, double> > swathes;
  swathes.push_back(std::make_pair(400, 425));
  swathes.push_back(std::make_pair(424, 450));
  swathes.push_back(std::make_pair(449, 475));
  swathes.push_back(std::make_pair(474, 500));
  swathes.push_back(std::make_pair(499, 525));
  swathes.push_back(std::make_pair(524, 550));
  swathes.push_back(std::make_pair(549, 575));
  swathes.push_back(std::make_pair(574, 600));
  swathes.push_back(std::make_pair(599, 625));
  swathes.push_back(std::make_pair(624, 650));
  swathes.push_back(std::make_pair(649, 675));
  swathes.push_back(std::make_pair(674, 700));
  swathes.push_back(std::make_pair(699, 725));
  swathes.push_back(std::make_pair(724, 750));
  swathes.push_back(std::make_pair(749, 775));
  swathes.push_back(std::make_pair(774, 800));
  swathes.push_back(std::make_pair(799, 825));
  swathes.push_back(std::make_pair(824, 850));
  swathes.push_back(std::make_pair(849, 875));
  swathes.push_back(std::make_pair(874, 900));
  swathes.push_back(std::make_pair(899, 925));
  swathes.push_back(std::make_pair(924, 950));
  swathes.push_back(std::make_pair(949, 975));
  swathes.push_back(std::make_pair(974, 1000));
  swathes.push_back(std::make_pair(999, 1025));
  swathes.push_back(std::make_pair(1024, 1050));
  swathes.push_back(std::make_pair(1049, 1075));
  swathes.push_back(std::make_pair(1074, 1100));
  swathes.push_back(std::make_pair(1099, 1125));
  swathes.push_back(std::make_pair(1124, 1150));
  swathes.push_back(std::make_pair(1149, 1175));
  swathes.push_back(std::make_pair(1174, 1200));

  TraMLFile traml;
  TargetedExperiment targeted_exp;
  String in = "MRMAssay_uisTransitions_input_1.TraML";
  traml.load(OPENMS_GET_TEST_DATA_PATH(in), targeted_exp);
  MRMAssay mrma;

  std::vector<String> fragment_types1;
  fragment_types1.push_back(String("y"));
  std::vector<size_t> fragment_charges1;
  fragment_charges1.push_back(2);
  bool enable_specific_losses1 = true;
  bool enable_unspecific_losses1 = false;
  bool enable_ms2_precursors1 = false;
  double product_mz_threshold1 = 0.05;

#if OPENMS_BOOST_VERSION_MINOR < 56
  String out1 = "MRMAssay_uisTransitions_output_1.TraML";
#else
  String out1 = "MRMAssay_uisTransitions_output_1_boost58.TraML";
#endif

  TargetedExperiment targeted_exp1 = targeted_exp;

  mrma.uisTransitions(targeted_exp1, fragment_types1, fragment_charges1, enable_specific_losses1, enable_unspecific_losses1, enable_ms2_precursors1, product_mz_threshold1, swathes, -4, 20, 42);

  String test1;
  NEW_TMP_FILE(test1);
  traml.store(test1, targeted_exp1);
#if !defined(__APPLE__) // currently fails on macOS likely due to different boost version and different random number generator
  TEST_FILE_SIMILAR(test1.c_str(), OPENMS_GET_TEST_DATA_PATH(out1))
#endif
  std::vector<String> fragment_types2;
  fragment_types2.push_back(String("y"));
  std::vector<size_t> fragment_charges2;
  fragment_charges2.push_back(2);
  bool enable_specific_losses2 = true;
  bool enable_unspecific_losses2 = true;
  bool enable_ms2_precursors2 = false;
  double product_mz_threshold2 = 0.05;

#if OPENMS_BOOST_VERSION_MINOR < 56
  String out2 = "MRMAssay_uisTransitions_output_2.TraML";
#else
  String out2 = "MRMAssay_uisTransitions_output_2_boost58.TraML";
#endif

  TargetedExperiment targeted_exp2 = targeted_exp;

  mrma.uisTransitions(targeted_exp2, fragment_types2, fragment_charges2, enable_specific_losses2, enable_unspecific_losses2, enable_ms2_precursors2, product_mz_threshold2, swathes, -4, 20, 42);

  String test2;
  NEW_TMP_FILE(test2);
  traml.store(test2, targeted_exp2);

#if !defined(__APPLE__) // currently fails on macOS likely due to different boost version and different random number generator
 TEST_FILE_SIMILAR(test2.c_str(), OPENMS_GET_TEST_DATA_PATH(out2))
#endif
}

END_SECTION

START_SECTION(void uisTransitions(OpenMS::TargetedExperiment& exp, std::vector<String> fragment_types, std::vector<size_t> fragment_charges, bool enable_specific_losses, bool enable_unspecific_losses, double mz_threshold, std::vector<std::pair<double, double> > swathes, int round_decPow, size_t max_num_alternative_localizations))
{
  std::vector<std::pair<double, double> > swathes;
  swathes.push_back(std::make_pair(400, 425));
  swathes.push_back(std::make_pair(424, 450));
  swathes.push_back(std::make_pair(449, 475));
  swathes.push_back(std::make_pair(474, 500));
  swathes.push_back(std::make_pair(499, 525));
  swathes.push_back(std::make_pair(524, 550));
  swathes.push_back(std::make_pair(549, 575));
  swathes.push_back(std::make_pair(574, 600));
  swathes.push_back(std::make_pair(599, 625));
  swathes.push_back(std::make_pair(624, 650));
  swathes.push_back(std::make_pair(649, 675));
  swathes.push_back(std::make_pair(674, 700));
  swathes.push_back(std::make_pair(699, 725));
  swathes.push_back(std::make_pair(724, 750));
  swathes.push_back(std::make_pair(749, 775));
  swathes.push_back(std::make_pair(774, 800));
  swathes.push_back(std::make_pair(799, 825));
  swathes.push_back(std::make_pair(824, 850));
  swathes.push_back(std::make_pair(849, 875));
  swathes.push_back(std::make_pair(874, 900));
  swathes.push_back(std::make_pair(899, 925));
  swathes.push_back(std::make_pair(924, 950));
  swathes.push_back(std::make_pair(949, 975));
  swathes.push_back(std::make_pair(974, 1000));
  swathes.push_back(std::make_pair(999, 1025));
  swathes.push_back(std::make_pair(1024, 1050));
  swathes.push_back(std::make_pair(1049, 1075));
  swathes.push_back(std::make_pair(1074, 1100));
  swathes.push_back(std::make_pair(1099, 1125));
  swathes.push_back(std::make_pair(1124, 1150));
  swathes.push_back(std::make_pair(1149, 1175));
  swathes.push_back(std::make_pair(1174, 1200));

  TraMLFile traml;
  TargetedExperiment targeted_exp;
  String in = "MRMAssay_uisTransitions_input_3.TraML";
  traml.load(OPENMS_GET_TEST_DATA_PATH(in), targeted_exp);
  MRMAssay mrma;

  std::vector<String> fragment_types1;
  fragment_types1.push_back(String("b"));
  std::vector<size_t> fragment_charges1;
  fragment_charges1.push_back(3);
  bool enable_losses1 = true;
  bool enable_ms2_precursors1 = false;
  double product_mz_threshold1 = 0.05;

#if OPENMS_BOOST_VERSION_MINOR < 56
  String out1 = "MRMAssay_uisTransitions_output_3.TraML";
#else
  String out1 = "MRMAssay_uisTransitions_output_3_boost58.TraML";
#endif

  TargetedExperiment targeted_exp1 = targeted_exp;

  mrma.uisTransitions(targeted_exp1, fragment_types1, fragment_charges1, enable_losses1, enable_losses1, enable_ms2_precursors1, product_mz_threshold1, swathes, -4, 20, 42);

  String test1;
  NEW_TMP_FILE(test1);
  traml.store(test1, targeted_exp1);
	   
#if !defined(__APPLE__) // currently fails on macOS likely due to different boost version and different random number generator
 TEST_FILE_SIMILAR(test1.c_str(), OPENMS_GET_TEST_DATA_PATH(out1)) 
#endif
	
  std::vector<String> fragment_types2;
  fragment_types2.push_back(String("y"));
  fragment_types2.push_back(String("b"));
  std::vector<size_t> fragment_charges2;
  fragment_charges2.push_back(2);
  fragment_charges2.push_back(3);
  bool enable_losses2 = true;
  bool enable_ms2_precursors2 = false;
  double product_mz_threshold2 = 0.05;

#if OPENMS_BOOST_VERSION_MINOR < 56
  String out2 = "MRMAssay_uisTransitions_output_4.TraML";
#else
  String out2 = "MRMAssay_uisTransitions_output_4_boost58.TraML";
#endif

  TargetedExperiment targeted_exp2 = targeted_exp;

  mrma.uisTransitions(targeted_exp2, fragment_types2, fragment_charges2, enable_losses2, enable_losses2, enable_ms2_precursors2, product_mz_threshold2, swathes, -4, 20, 42);

  String test2;
  NEW_TMP_FILE(test2);
  traml.store(test2, targeted_exp2);

#if !defined(__APPLE__) // currently fails on macOS likely due to different boost version and different random number generator
 TEST_FILE_SIMILAR(test2.c_str(), OPENMS_GET_TEST_DATA_PATH(out2)) 
#endif
	
  std::vector<String> fragment_types3;
  fragment_types3.push_back(String("y"));
  fragment_types3.push_back(String("b"));
  std::vector<size_t> fragment_charges3;
  fragment_charges3.push_back(2);
  fragment_charges3.push_back(3);
  bool enable_losses3 = true;
  bool enable_ms2_precursors3 = true;
  double product_mz_threshold3 = 0.05;

#if OPENMS_BOOST_VERSION_MINOR < 56
  String out3 = "MRMAssay_uisTransitions_output_5.TraML";
#else
  String out3 = "MRMAssay_uisTransitions_output_5_boost58.TraML";
#endif

  TargetedExperiment targeted_exp3 = targeted_exp;

  mrma.uisTransitions(targeted_exp3, fragment_types3, fragment_charges3, enable_losses3, enable_losses3, enable_ms2_precursors3, product_mz_threshold3, swathes, -4, 20, 42);

  String test3;
  NEW_TMP_FILE(test3);
  traml.store(test3, targeted_exp3);

#if !defined(__APPLE__) // currently fails on macOS likely due to different boost version and different random number generator
 TEST_FILE_SIMILAR(test3.c_str(), OPENMS_GET_TEST_DATA_PATH(out3)) 
#endif	
}

END_SECTION

END_TEST
