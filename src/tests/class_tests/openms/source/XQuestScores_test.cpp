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

PeakSpectrum theo_spec_1, theo_spec_2, theo_spec_3, theo_spec_4;
specGen.getCommonIonSpectrum(theo_spec_1, AASequence::fromString("PEPTIDE"), 3, true);
specGen.getCommonIonSpectrum(theo_spec_2, AASequence::fromString("PEPTIDE"), 2, true);
specGen.getCommonIonSpectrum(theo_spec_3, AASequence::fromString("PEPTIDE"), 3, true);
specGen.getCommonIonSpectrum(theo_spec_4, AASequence::fromString("PEPTEDI"), 2, true);
std::vector <std::pair <Size, Size> > alignment1;
std::vector <std::pair <Size, Size> > alignment2;

OPXLSpectrumProcessingAlgorithms::getSpectrumAlignment(alignment1, theo_spec_1, theo_spec_2, 20, true);
OPXLSpectrumProcessingAlgorithms::getSpectrumAlignment(alignment2, theo_spec_3, theo_spec_4, 20, true);

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
        TEST_REAL_SIMILAR(XQuestScores::matchOddsScore(theo_spec_1, alignment1.size(), 0.1, false), 28.53678);
        TEST_REAL_SIMILAR(XQuestScores::matchOddsScore(theo_spec_3, alignment2.size(), 0.1, false), 6.82721);
        TEST_REAL_SIMILAR(XQuestScores::matchOddsScore(theo_spec_1, alignment1.size(), 0.2, false, true, 5), 34.02875);
        TEST_REAL_SIMILAR(XQuestScores::matchOddsScore(theo_spec_3, alignment2.size(), 0.2, false), 5.47099);
        TEST_REAL_SIMILAR(XQuestScores::matchOddsScore(theo_spec_1, alignment1.size(), 10, true), 708.39641);
        TEST_REAL_SIMILAR(XQuestScores::matchOddsScore(theo_spec_3, alignment2.size(), 10, true), 14.26185);
        TEST_REAL_SIMILAR(XQuestScores::matchOddsScore(theo_spec_1, alignment1.size(), 20, true, true, 5), 708.39641);
        TEST_REAL_SIMILAR(XQuestScores::matchOddsScore(theo_spec_3, alignment2.size(), 20, true), 12.87628);
END_SECTION

START_SECTION(static double logOccupancyProb(const PeakSpectrum& theoretical_spec,  const Size matched_size, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm))
        TEST_REAL_SIMILAR(XQuestScores::logOccupancyProb(theo_spec_1, alignment1.size(), 0.1, false), 32.67635);
        TEST_REAL_SIMILAR(XQuestScores::logOccupancyProb(theo_spec_3, alignment2.size(), 0.1, false), 8.19844);
        TEST_REAL_SIMILAR(XQuestScores::logOccupancyProb(theo_spec_1, alignment1.size(), 0.2, false), 28.53678);
        TEST_REAL_SIMILAR(XQuestScores::logOccupancyProb(theo_spec_3, alignment2.size(), 0.2, false), 6.82721);
        TEST_REAL_SIMILAR(XQuestScores::logOccupancyProb(theo_spec_1, alignment1.size(), 10, true), 708.39641);
        TEST_REAL_SIMILAR(XQuestScores::logOccupancyProb(theo_spec_3, alignment2.size(), 10, true), 15.94029);
        TEST_REAL_SIMILAR(XQuestScores::logOccupancyProb(theo_spec_1, alignment1.size(), 20, true), 708.39641);
        TEST_REAL_SIMILAR(XQuestScores::logOccupancyProb(theo_spec_3, alignment2.size(), 20, true), 14.55431);
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

START_SECTION(static double matchedCurrentChain(const std::vector< std::pair< Size, Size > >& matched_spec_common, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks, const PeakSpectrum& spectrum_common_peaks, const PeakSpectrum& spectrum_xlink_peaks))
    TEST_REAL_SIMILAR(XQuestScores::matchedCurrentChain(alignment1, alignment1, theo_spec_2, theo_spec_4), 10.0)
END_SECTION

START_SECTION(static double totalMatchedCurrent(const std::vector< std::pair< Size, Size > >& matched_spec_common_alpha, const std::vector< std::pair< Size, Size > >& matched_spec_common_beta, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks_alpha, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks_beta, const PeakSpectrum& spectrum_common_peaks, const PeakSpectrum& spectrum_xlink_peaks))
    TEST_REAL_SIMILAR(XQuestScores::totalMatchedCurrent(alignment1, alignment2, alignment1, alignment2, theo_spec_2, theo_spec_4), 10.0)
END_SECTION

START_SECTION(static std::vector< double > xCorrelation(const PeakSpectrum & spec1, const PeakSpectrum & spec2, Int maxshift, double tolerance))
    std::vector <double> xcorr_scores = XQuestScores::xCorrelation(theo_spec_1, theo_spec_2, 2, 0.2);

    TEST_EQUAL(xcorr_scores [0] < 0, true)
    TEST_EQUAL(xcorr_scores [1] < 0, true)
    TEST_REAL_SIMILAR(xcorr_scores [2], 0.83291)
    TEST_EQUAL(xcorr_scores [3] < 0, true)
    TEST_EQUAL(xcorr_scores [4] < 0, true)
END_SECTION

END_TEST
