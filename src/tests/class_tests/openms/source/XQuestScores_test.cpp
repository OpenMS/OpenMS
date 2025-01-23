// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/ANALYSIS/XLMS/XQuestScores.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGeneratorXLMS.h>
#include <OpenMS/ANALYSIS/XLMS/OPXLSpectrumProcessingAlgorithms.h>

using namespace OpenMS;

START_TEST(XQuestScores, "$Id$")
TheoreticalSpectrumGeneratorXLMS specGen;
Param param = specGen.getParameters();
param.setValue("add_isotopes", "false");
param.setValue("add_metainfo", "true");
param.setValue("add_first_prefix_ion", "false");
specGen.setParameters(param);

PeakSpectrum theo_spec_1, theo_spec_2, theo_spec_3;
AASequence peptide1 = AASequence::fromString("PEPTIDEPEPTIDEPEPTIDE");
AASequence peptide2 = AASequence::fromString("PEPTIDEEDITPEPTIDE");
AASequence peptide3 = AASequence::fromString("EDITPEPTIDE");
specGen.getLinearIonSpectrum(theo_spec_1, peptide1, 3, true, 2);
specGen.getLinearIonSpectrum(theo_spec_2, peptide2, 3, true, 2);
specGen.getLinearIonSpectrum(theo_spec_3, peptide3, 5, true, 2);
// specGen.getLinearIonSpectrum(theo_spec_4, peptide3, 2, true);
std::vector <std::pair <Size, Size> > alignment1;
std::vector <std::pair <Size, Size> > alignment2;

DataArrays::FloatDataArray dummy_array1;
DataArrays::FloatDataArray dummy_array2;

OPXLSpectrumProcessingAlgorithms::getSpectrumAlignmentFastCharge(alignment1, 20, true, theo_spec_1, theo_spec_2, theo_spec_1.getIntegerDataArrays()[0], theo_spec_2.getIntegerDataArrays()[0], dummy_array1);
OPXLSpectrumProcessingAlgorithms::getSpectrumAlignmentFastCharge(alignment2, 20, true, theo_spec_1, theo_spec_3, theo_spec_1.getIntegerDataArrays()[0], theo_spec_3.getIntegerDataArrays()[0], dummy_array2);

START_SECTION(static float preScore(Size matched_alpha, Size ions_alpha, Size matched_beta, Size ions_beta))
	TEST_REAL_SIMILAR(XQuestScores::preScore(1, 1, 1, 1), 1.0)
	TEST_REAL_SIMILAR(XQuestScores::preScore(2, 4, 3, 6), 0.5)
	TEST_REAL_SIMILAR(XQuestScores::preScore(3, 2, 9, 6), 1.5) // more matched peaks, than theoretical peaks. practically impossible
	TEST_REAL_SIMILAR(XQuestScores::preScore(0, 5, 0, 5), 0.0)
	TEST_REAL_SIMILAR(XQuestScores::preScore(0, 5, 3, 5), 0.10954)
	TEST_REAL_SIMILAR(XQuestScores::preScore(2, 5, 0, 5), 0.08944)
	TEST_REAL_SIMILAR(XQuestScores::preScore(0, 50, 0, 50), 0.0)
	TEST_REAL_SIMILAR(XQuestScores::preScore(0, 50, 3, 50), 0.01095)
	TEST_REAL_SIMILAR(XQuestScores::preScore(2, 50, 0, 50), 0.00894)
	TEST_REAL_SIMILAR(XQuestScores::preScore(5, 50, 0, 50), 0.01414)
	TEST_REAL_SIMILAR(XQuestScores::preScore(45, 50, 0, 50), 0.04242)
	TEST_REAL_SIMILAR(XQuestScores::preScore(2, 50, 3, 50), 0.04898)
	TEST_REAL_SIMILAR(XQuestScores::preScore(1, 50, 1, 50), 0.02)
	TEST_REAL_SIMILAR(XQuestScores::preScore(2, 50, 2, 50), 0.04)
	TEST_REAL_SIMILAR(XQuestScores::preScore(45, 50, 5, 50), 0.3)
	TEST_REAL_SIMILAR(XQuestScores::preScore(25, 50, 25, 50), 0.5)
END_SECTION

START_SECTION(static float preScore(Size matched_alpha, Size ions_alpha))
	TEST_REAL_SIMILAR(XQuestScores::preScore(1, 1), 1.0)
	TEST_REAL_SIMILAR(XQuestScores::preScore(2, 1), 2.0)
	TEST_REAL_SIMILAR(XQuestScores::preScore(0, 2), 0.0)
	TEST_REAL_SIMILAR(XQuestScores::preScore(0, 50), 0.0)
	TEST_REAL_SIMILAR(XQuestScores::preScore(1, 50), 0.02)
	TEST_REAL_SIMILAR(XQuestScores::preScore(3, 50), 0.06)
	TEST_REAL_SIMILAR(XQuestScores::preScore(9, 18), 0.5)
END_SECTION

START_SECTION(static double matchOddsScore(const PeakSpectrum& theoretical_spec,  const Size matched_size, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, bool is_xlink_spectrum = false, Size n_charges = 1))
        TEST_EQUAL(theo_spec_1.size(), 46)
        TEST_EQUAL(alignment1.size(), 28)
        TEST_EQUAL(alignment2.size(), 10)
        TEST_REAL_SIMILAR(theo_spec_1.back().getMZ() - theo_spec_1[0].getMZ(), 1903.33405)
        TEST_REAL_SIMILAR(std::log(theo_spec_1.back().getMZ()) - std::log(theo_spec_1[0].getMZ()), 3.99930)

        TEST_REAL_SIMILAR(XQuestScores::matchOddsScore(theo_spec_1, alignment1.size(), 0.1, false), 106.63674);
        TEST_REAL_SIMILAR(XQuestScores::matchOddsScore(theo_spec_2, alignment1.size(), 0.1, false), 111.87796);
        TEST_REAL_SIMILAR(XQuestScores::matchOddsScore(theo_spec_1, alignment2.size(), 0.1, false), 28.07671);
        TEST_REAL_SIMILAR(XQuestScores::matchOddsScore(theo_spec_3, alignment2.size(), 0.1, false), 24.22081);
        TEST_REAL_SIMILAR(XQuestScores::matchOddsScore(theo_spec_1, alignment1.size(), 0.2, false, true, 2), 106.63373);
        TEST_REAL_SIMILAR(XQuestScores::matchOddsScore(theo_spec_2, alignment1.size(), 0.2, false, true, 2), 111.87432);
        TEST_REAL_SIMILAR(XQuestScores::matchOddsScore(theo_spec_3, alignment2.size(), 0.2, false), 17.11504);
        TEST_REAL_SIMILAR(XQuestScores::matchOddsScore(theo_spec_1, alignment1.size(), 10, true), 187.24386);
        TEST_REAL_SIMILAR(XQuestScores::matchOddsScore(theo_spec_2, alignment1.size(), 10, true), 198.42811);
        TEST_REAL_SIMILAR(XQuestScores::matchOddsScore(theo_spec_1, alignment2.size(), 10, true), 58.41773);
        TEST_REAL_SIMILAR(XQuestScores::matchOddsScore(theo_spec_3, alignment2.size(), 10, true), 63.85680);
        TEST_REAL_SIMILAR(XQuestScores::matchOddsScore(theo_spec_1, alignment1.size(), 20, true, true, 2), 187.24367);
        TEST_REAL_SIMILAR(XQuestScores::matchOddsScore(theo_spec_3, alignment2.size(), 20, true), 56.24576);
END_SECTION

START_SECTION(static double logOccupancyProb(const PeakSpectrum& theoretical_spec,  const Size matched_size, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm))
        TEST_REAL_SIMILAR(XQuestScores::logOccupancyProb(theo_spec_1, alignment1.size(), 0.1, false), 126.59011);
        TEST_REAL_SIMILAR(XQuestScores::logOccupancyProb(theo_spec_3, alignment2.size(), 0.1, false), 31.58523);
        TEST_REAL_SIMILAR(XQuestScores::logOccupancyProb(theo_spec_1, alignment2.size(), 0.1, false), 35.52062);
        TEST_REAL_SIMILAR(XQuestScores::logOccupancyProb(theo_spec_1, alignment1.size(), 0.2, false), 106.63674);
        TEST_REAL_SIMILAR(XQuestScores::logOccupancyProb(theo_spec_1, alignment2.size(), 0.2, false), 28.07671);
        TEST_REAL_SIMILAR(XQuestScores::logOccupancyProb(theo_spec_3, alignment2.size(), 0.2, false), 24.22081);
        TEST_REAL_SIMILAR(XQuestScores::logOccupancyProb(theo_spec_1, alignment1.size(), 10, true), 214.75707);
        TEST_REAL_SIMILAR(XQuestScores::logOccupancyProb(theo_spec_1, alignment2.size(), 10, true), 68.84436);
        TEST_REAL_SIMILAR(XQuestScores::logOccupancyProb(theo_spec_3, alignment2.size(), 10, true), 73.47408);
        TEST_REAL_SIMILAR(XQuestScores::logOccupancyProb(theo_spec_1, alignment1.size(), 20, true), 194.66285);
        TEST_REAL_SIMILAR(XQuestScores::logOccupancyProb(theo_spec_1, alignment2.size(), 20, true), 61.22836);
        TEST_REAL_SIMILAR(XQuestScores::logOccupancyProb(theo_spec_3, alignment2.size(), 20, true), 65.85512);
        TEST_REAL_SIMILAR(XQuestScores::logOccupancyProb(theo_spec_1, alignment1.size(), 200, true), 128.01463);
        TEST_REAL_SIMILAR(XQuestScores::logOccupancyProb(theo_spec_1, alignment2.size(), 200, true), 36.05495);
        TEST_REAL_SIMILAR(XQuestScores::logOccupancyProb(theo_spec_3, alignment2.size(), 200, true), 40.62847);
END_SECTION

START_SECTION(static double weightedTICScoreXQuest(Size alpha_size, Size beta_size, double intsum_alpha, double intsum_beta, double total_current, bool type_is_cross_link))
    TEST_REAL_SIMILAR(XQuestScores::weightedTICScoreXQuest(20, 10, 500.0, 500.0, 1500.0, true), 0.13636)
    TEST_REAL_SIMILAR(XQuestScores::weightedTICScoreXQuest(20, 10, 1000.0, 500.0, 1500.0, true), 0.18181)
    TEST_REAL_SIMILAR(XQuestScores::weightedTICScoreXQuest(20, 10, 500.0, 1000.0, 1500.0, true), 0.22727)
    TEST_REAL_SIMILAR(XQuestScores::weightedTICScoreXQuest(20, 10, 1450.0, 50.0, 1500.0, true), 0.14090)
    TEST_REAL_SIMILAR(XQuestScores::weightedTICScoreXQuest(20, 10, 50.0, 1450.0, 1500.0, true), 0.26818)
    TEST_REAL_SIMILAR(XQuestScores::weightedTICScoreXQuest(20, 0, 500.0, 0.0, 1500.0, false), 0.08333)
END_SECTION

START_SECTION(static double weightedTICScore(Size alpha_size, Size beta_size, double intsum_alpha, double intsum_beta, double total_current, bool type_is_cross_link))
    TEST_REAL_SIMILAR(XQuestScores::weightedTICScore(20, 10, 500.0, 500.0, 1500.0, true), 0.5)
    TEST_REAL_SIMILAR(XQuestScores::weightedTICScore(20, 10, 1000.0, 500.0, 1500.0, true), 0.66666)
    TEST_REAL_SIMILAR(XQuestScores::weightedTICScore(20, 10, 500.0, 1000.0, 1500.0, true), 0.83333)
    TEST_REAL_SIMILAR(XQuestScores::weightedTICScore(20, 10, 1450.0, 50.0, 1500.0, true), 0.51666)
    TEST_REAL_SIMILAR(XQuestScores::weightedTICScore(20, 10, 50.0, 1450.0, 1500.0, true), 0.98333)
    TEST_REAL_SIMILAR(XQuestScores::weightedTICScore(20, 0, 500.0, 0.0, 1500.0, false), 0.33333)
END_SECTION

START_SECTION(static double matchedCurrentChain(const std::vector< std::pair< Size, Size > >& matched_spec_linear, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks, const PeakSpectrum& spectrum_linear_peaks, const PeakSpectrum& spectrum_xlink_peaks))
    TEST_REAL_SIMILAR(XQuestScores::matchedCurrentChain(alignment1, alignment2, theo_spec_2, theo_spec_3), 38.0)
END_SECTION

START_SECTION(static double totalMatchedCurrent(const std::vector< std::pair< Size, Size > >& matched_spec_linear_alpha, const std::vector< std::pair< Size, Size > >& matched_spec_linear_beta, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks_alpha, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks_beta, const PeakSpectrum& spectrum_linear_peaks, const PeakSpectrum& spectrum_xlink_peaks))
    TEST_REAL_SIMILAR(XQuestScores::totalMatchedCurrent(alignment1, alignment1, alignment2, alignment2, theo_spec_2, theo_spec_3), 38.0)
END_SECTION

START_SECTION(static std::vector< double > xCorrelation(const PeakSpectrum & spec1, const PeakSpectrum & spec2, Int maxshift, double tolerance))
    std::vector <double> xcorr_scores = XQuestScores::xCorrelation(theo_spec_1, theo_spec_2, 2, 0.2);

    TEST_EQUAL(xcorr_scores [0] < 0.5, true)
    TEST_EQUAL(xcorr_scores [1] < 0, true)
    TEST_REAL_SIMILAR(xcorr_scores [2], 0.65121)
    TEST_EQUAL(xcorr_scores [3] < 0, true)
    TEST_EQUAL(xcorr_scores [4] < 0.5, true)
END_SECTION

START_SECTION(static double XQuestScores::xCorrelationPrescore(const PeakSpectrum & spec1, const PeakSpectrum & spec2, double tolerance))
    double xcorr_fast = XQuestScores::xCorrelationPrescore(theo_spec_1, theo_spec_2, 0.2);
    TEST_REAL_SIMILAR(xcorr_fast, 0.7)
END_SECTION

END_TEST
