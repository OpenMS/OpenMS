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
std::cout << "WTF0" << std::endl;
specGen.getCommonIonSpectrum(theo_spec_1, AASequence::fromString("PEPTIDE"), 2, true);
std::cout << "WTF1" << std::endl;
specGen.getCommonIonSpectrum(theo_spec_2, AASequence::fromString("PEPTEDI"), 4, true);
std::cout << "WTF2" << std::endl;
specGen.getCommonIonSpectrum(theo_spec_3, AASequence::fromString("PEPTIDE"), 3, true);
std::cout << "WTF3" << std::endl;
specGen.getCommonIonSpectrum(theo_spec_4, AASequence::fromString("PEPTEDI"), 1, true);
std::cout << "WTF4" << std::endl;
std::vector <std::pair <Size, Size> > alignment1;
std::vector <std::pair <Size, Size> > alignment2;

OPXLSpectrumProcessingAlgorithms::getSpectrumAlignment(alignment1, theo_spec_1, theo_spec_2, 50, true);
OPXLSpectrumProcessingAlgorithms::getSpectrumAlignment(alignment2, theo_spec_3, theo_spec_4, 50, true);

START_SECTION(preScore())
   /* @brief compute a simple and fast to compute pre-score for a cross-link spectrum match
    * @param number of experimental peaks matched to theoretical common ions from the alpha peptide
    * @param number of theoretical ions from the alpha peptide
    * @param number of experimental peaks matched to theoretical common ions from the beta peptide
    * @param number of theoretical ions from the beta peptide
    */
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

	TEST_REAL_SIMILAR(XQuestScores::preScore(1, 1), 1.0)
	TEST_REAL_SIMILAR(XQuestScores::preScore(2, 1), 2.0)
	TEST_REAL_SIMILAR(XQuestScores::preScore(0, 2), 0.0)
	TEST_REAL_SIMILAR(XQuestScores::preScore(0, 50), 0.0)
	TEST_REAL_SIMILAR(XQuestScores::preScore(1, 50), 0.02)
	TEST_REAL_SIMILAR(XQuestScores::preScore(3, 50), 0.06)
	TEST_REAL_SIMILAR(XQuestScores::preScore(9, 18), 0.5)
END_SECTION

START_SECTION(matchOddsScore())
   /* @brief compute the match-odds score, a score based on the probability of getting the given number of matched peaks by chance
    * @param theoretical spectrum, sorted by position
    * @param alignment between the theoretical and the experimental spectra
    * @param fragment mass tolerance of the alignment
    * @param fragment mass tolerance unit of the alignment, true = ppm, false = Da
    * @param type of cross-link, true = cross-link, false = mono-link
    * @param number of considered charges in the theoretical spectrum
    */
        TEST_REAL_SIMILAR(XQuestScores::matchOddsScore(theo_spec_1, alignment1, 0.2, false, true, 1), 6.07050);
        TEST_REAL_SIMILAR(XQuestScores::matchOddsScore(theo_spec_1, alignment1, 0.2, false, true, 2), 7.42116);
        TEST_REAL_SIMILAR(XQuestScores::matchOddsScore(theo_spec_1, alignment1, 0.2, false, true, 3), 8.20397);
        TEST_REAL_SIMILAR(XQuestScores::matchOddsScore(theo_spec_1, alignment1, 0.2, false, false, 3), 6.07050);
END_SECTION

START_SECTION(weightedTICScore())
   /* @brief compute the weighted total ion current score for a cross-link. Reimplementation from xQuest.
    * @param sequence length of alpha peptide
    * @param sequence length of beta peptide
    * @param intensity sum of matched peaks from alpha peptide
    * @param intensity sum of matched peaks from beta peptide
    * @param type of cross-link, true = cross-link, false = mono-link
    * @param sum of peak intensities of the experimental spectrum
    * @param true = cross-link, false = mono-link. in case of a mono-link, beta_size and intsum_beta should be 0
    */
    TEST_REAL_SIMILAR(XQuestScores::weightedTICScoreXQuest(20, 10, 500.0, 500.0, 1500.0, true), 0.13636)
    TEST_REAL_SIMILAR(XQuestScores::weightedTICScoreXQuest(20, 10, 1000.0, 500.0, 1500.0, true), 0.18181)
    TEST_REAL_SIMILAR(XQuestScores::weightedTICScoreXQuest(20, 10, 500.0, 1000.0, 1500.0, true), 0.22727)
    TEST_REAL_SIMILAR(XQuestScores::weightedTICScoreXQuest(20, 10, 1450.0, 50.0, 1500.0, true), 0.14090)
    TEST_REAL_SIMILAR(XQuestScores::weightedTICScoreXQuest(20, 10, 50.0, 1450.0, 1500.0, true), 0.26818)
    TEST_REAL_SIMILAR(XQuestScores::weightedTICScoreXQuest(20, 0, 500.0, 0.0, 1500.0, false), 0.08333)

    TEST_REAL_SIMILAR(XQuestScores::weightedTICScore(20, 10, 500.0, 500.0, 1500.0, true), 0.5)
    TEST_REAL_SIMILAR(XQuestScores::weightedTICScore(20, 10, 1000.0, 500.0, 1500.0, true), 0.66666)
    TEST_REAL_SIMILAR(XQuestScores::weightedTICScore(20, 10, 500.0, 1000.0, 1500.0, true), 0.83333)
    TEST_REAL_SIMILAR(XQuestScores::weightedTICScore(20, 10, 1450.0, 50.0, 1500.0, true), 0.51666)
    TEST_REAL_SIMILAR(XQuestScores::weightedTICScore(20, 10, 50.0, 1450.0, 1500.0, true), 0.98333)
    TEST_REAL_SIMILAR(XQuestScores::weightedTICScore(20, 0, 500.0, 0.0, 1500.0, false), 0.33333)

END_SECTION

START_SECTION(matchedCurrentChain())
   /* @brief computes sum of peak intensities of matched peaks for either the alpha or the beta peptide
    * @param alignment between common alpha or beta ions and common experimental peaks
    * @param alignment between xlink alpha or beta ions and xlink experimental peaks
    * @param experimental common ion spectrum
    * @param experimental xlink spectrum
    */
    TEST_REAL_SIMILAR(XQuestScores::matchedCurrentChain(alignment1, alignment1, theo_spec_2, theo_spec_4), 4.0)

END_SECTION

START_SECTION(totalMatchedCurrent())
   /* @brief computes sum of peak intensities of all matched peaks
    * @param alignment between common alpha ions and common experimental peaks
    * @param alignment between common beta ions and common experimental peaks
    * @param alignment between xlink alpha ions and xlink experimental peaks
    * @param alignment between xlink beta ions and xlink experimental peaks
    * @param experimental common ion spectrum
    * @param experimental xlink spectrum
    */
    TEST_REAL_SIMILAR(XQuestScores::totalMatchedCurrent(alignment1, alignment2, alignment1, alignment2, theo_spec_2, theo_spec_4), 6.0)
END_SECTION

START_SECTION(xCorrelation())
    // Cross-correlation, with shifting the second spectrum from -maxshift to +maxshift of tolerance bins (Tolerance in Da, a constant binsize)
   /* @brief computes a crude cross-correlation between two spectra. Crude, because it uses a static binsize based on a tolerance in Da and it uses equal intensities for all peaks
    * @param first spectrum
    * @param second spectrum
    * @param number of bins, that should be considered for shifting the second spectrum. the second spectrum is shifted from -maxshift to +maxshift of tolerance bins and a correlation is computed for each position.
    * @param tolerance or binsize in Da
    */
    std::vector <double> xcorr_scores = XQuestScores::xCorrelation(theo_spec_1, theo_spec_3, 2, 0.2);

    TEST_EQUAL(xcorr_scores [0] < 1.0, true)
    TEST_EQUAL(xcorr_scores [1] < 1.0, true)
    TEST_REAL_SIMILAR(xcorr_scores [2], 0.83291)
    TEST_EQUAL(xcorr_scores [3] < 1.0, true)
    TEST_EQUAL(xcorr_scores [4] < 1.0, true)
END_SECTION

END_TEST

