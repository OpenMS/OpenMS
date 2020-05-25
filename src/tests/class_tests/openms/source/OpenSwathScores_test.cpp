// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathScores.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(OpenSwathScores, "$Id$")
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

OpenSwath_Scores* nullPointer = nullptr;
OpenSwath_Scores* ptr = nullptr;

START_SECTION(OpenSwath_Scores())
{
  ptr = new OpenSwath_Scores();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~OpenSwath_Scores())
{
  delete ptr;
}
END_SECTION

START_SECTION(EXTRA OpenSwath_Scores())
{
  OpenSwath_Scores scores;

  TEST_REAL_SIMILAR( scores.ms1_xcorr_coelution_score, -1.0)
  TEST_REAL_SIMILAR( scores.ms1_xcorr_shape_score, -1.0)
  TEST_REAL_SIMILAR( scores.ms1_mi_score, -1.0)
}
END_SECTION

START_SECTION((double get_quick_lda_score(double library_corr_, double library_norm_manhattan_, double norm_rt_score_, double xcorr_coelution_score_,
                                          double xcorr_shape_score_, double log_sn_score_) const))
{
  OpenSwath_Scores scores;

  TEST_REAL_SIMILAR( scores.get_quick_lda_score(1.0, 1.0, 1.0, 1.0, 1.0, 1.0), 
                                                                                 -0.5319046 +
                                                                                  2.1643962 +
                                                                                  8.0353047 +
                                                                                  0.1458914 +
                                                                                 -1.6901925 +
                                                                                 -0.8002824)
}
END_SECTION

START_SECTION((double calculate_lda_prescore(const OpenSwath_Scores& scores) const))
{
  OpenSwath_Scores scores;
  TEST_REAL_SIMILAR( scores.calculate_lda_prescore(scores), 0.0)

  scores.library_corr                     = 1.0;
  scores.library_norm_manhattan           = 1.0;
  scores.norm_rt_score                    = 1.0;
  scores.xcorr_coelution_score            = 1.0;
  scores.xcorr_shape_score                = 1.0;
  scores.log_sn_score                     = 1.0;
  scores.elution_model_fit_score          = 1.0;

  TEST_REAL_SIMILAR( scores.calculate_lda_prescore(scores), 
                                                             -0.34664267 + 
                                                              2.98700722 + 
                                                              7.05496384 + 
                                                              0.09445371 + 
                                                             -5.71823862 + 
                                                             -0.72989582 + 
                                                              1.88443209)  


}
END_SECTION

START_SECTION((double calculate_lda_single_transition(const OpenSwath_Scores& scores) const))
{
  OpenSwath_Scores scores;
  TEST_REAL_SIMILAR( scores.calculate_lda_single_transition(scores), 0.0)

  scores.library_corr                     = 1.0;
  scores.library_norm_manhattan           = 1.0;
  scores.norm_rt_score                    = 1.0;
  scores.xcorr_coelution_score            = 1.0;
  scores.xcorr_shape_score                = 1.0;
  scores.log_sn_score                     = 1.0;
  scores.elution_model_fit_score          = 1.0;

  TEST_REAL_SIMILAR( scores.calculate_lda_single_transition(scores), 
                                                                     7.05496384 +
                                                                    -0.72989582 +
                                                                     -1.08443209)
}
END_SECTION
    
START_SECTION((double calculate_swath_lda_prescore(const OpenSwath_Scores& scores) const))
{
  OpenSwath_Scores scores;
  TEST_REAL_SIMILAR( scores.calculate_swath_lda_prescore(scores), 0.0)

  scores.library_corr              = 1.0;
  scores.library_norm_manhattan    = 1.0;
  scores.norm_rt_score             = 1.0;
  scores.isotope_correlation       = 1.0;
  scores.isotope_overlap           = 1.0;
  scores.massdev_score             = 1.0;
  scores.xcorr_coelution_score     = 1.0;
  scores.xcorr_shape_score         = 1.0;
  scores.yseries_score             = 1.0;
  scores.log_sn_score              = 1.0;


  TEST_REAL_SIMILAR( scores.calculate_swath_lda_prescore(scores), 
                                                                   -0.19011762 + 
                                                                    2.47298914 + 
                                                                    5.63906731 + 
                                                                   -0.62640133 + 
                                                                    0.36006925 + 
                                                                    0.08814003 + 
                                                                    0.13978311 + 
                                                                   -1.16475032 + 
                                                                   -0.19267813 + 
                                                                   -0.61712054)  
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

