// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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

#include <OpenMS/KERNEL/MRMFeature.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(MRMFeature, "$Id$")

/////////////////////////////////////////////////////////////

MRMFeature* ptr = nullptr;
MRMFeature* nullPointer = nullptr;

START_SECTION(MRMFeature())
{
	ptr = new MRMFeature();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~MRMFeature())
{
  delete ptr;
}
END_SECTION

START_SECTION(MRMFeature(const MRMFeature &rhs))
{
  MRMFeature tmp;
  tmp.setIntensity(100.0);
  tmp.setExpectedRT(400.0);

  MRMFeature tmp2 (tmp);
  TEST_REAL_SIMILAR(tmp2.getIntensity(), 100.0)
  TEST_REAL_SIMILAR(tmp2.getExpectedRT(), 400.0)
}
END_SECTION

START_SECTION(MRMFeature& operator=(const MRMFeature &rhs))
{
  MRMFeature tmp;
  tmp.setIntensity(100.0);
  tmp.setExpectedRT(400.0);

  MRMFeature tmp2;
  tmp2 = tmp;

  TEST_REAL_SIMILAR(tmp2.getIntensity(), 100.0)
  TEST_REAL_SIMILAR(tmp2.getExpectedRT(), 400.0)
}
END_SECTION

START_SECTION(const double & getExpectedRT() const)
{
  MRMFeature tmp;
  tmp.setExpectedRT(400.0);

  TEST_REAL_SIMILAR(tmp.getExpectedRT(), 400.0)
}
END_SECTION

START_SECTION(void setExpectedRT(double rt))
{
  // tested above
  NOT_TESTABLE
}
END_SECTION
    
START_SECTION (const OpenSwath_Scores & getScores() const)
{
  // tested with set/add score
  NOT_TESTABLE
  MRMFeature tmp;
  OpenSwath_Scores tmp2 = tmp.getScores();
  
  TEST_REAL_SIMILAR(tmp2.library_manhattan, 0)
}
END_SECTION

START_SECTION (Feature & getFeature(String key))
{
  MRMFeature mrmfeature;
  Feature f1;
  f1.setMetaValue("dummy", 1);
  Feature f2;
  mrmfeature.addFeature(f1, "chromatogram1");
  mrmfeature.addFeature(f1, "chromatogram2");
  TEST_EQUAL(mrmfeature.getFeature("chromatogram1").getMetaValue("dummy"), 1)
}
END_SECTION

START_SECTION (void addFeature(Feature & feature, const String & key))
{
  // tested in getFeature
  NOT_TESTABLE
}
END_SECTION

START_SECTION (const std::vector<Feature> & getFeatures() const)
{
  MRMFeature mrmfeature;
  Feature f1;
  f1.setMetaValue("dummy", 1);
  Feature f2;
  mrmfeature.addFeature(f1, "chromatogram1");
  mrmfeature.addFeature(f1, "chromatogram2");
  TEST_EQUAL(mrmfeature.getFeatures().size(), 2)
}
END_SECTION

START_SECTION (void getFeatureIDs(std::vector<String> & result) const)
{
  MRMFeature mrmfeature;
  Feature f1;
  f1.setMetaValue("dummy", 1);
  Feature f2;
  mrmfeature.addFeature(f1, "chromatogram1");
  mrmfeature.addFeature(f1, "chromatogram2");
  std::vector<String> result;
  mrmfeature.getFeatureIDs(result);
  TEST_EQUAL(result.size(), 2)
  TEST_EQUAL(result[0], "chromatogram1")
  TEST_EQUAL(result[1], "chromatogram2")
}
END_SECTION

START_SECTION (void addPrecursorFeature(Feature & feature, const String & key))
{
  // Initially, there should be no feature present
  MRMFeature mrmfeature;
  {
    std::vector<String> result;
    mrmfeature.getPrecursorFeatureIDs(result);
    TEST_EQUAL(result.size(), 0)
  }

  // After adding a feature, there should be one feature present
  Feature f1;
  mrmfeature.addPrecursorFeature(f1, "precursor_chromatogram1");
  {
    std::vector<String> result;
    mrmfeature.getPrecursorFeatureIDs(result);
    TEST_EQUAL(result.size(), 1)
  }
}
END_SECTION

START_SECTION (void getPrecursorFeatureIDs(std::vector<String> & result) const)
{
  MRMFeature mrmfeature;
  Feature f1;
  f1.setMetaValue("dummy", 1);
  Feature f2;
  mrmfeature.addPrecursorFeature(f1, "chromatogram1");
  mrmfeature.addPrecursorFeature(f1, "chromatogram2");
  std::vector<String> result;
  mrmfeature.getPrecursorFeatureIDs(result);
  TEST_EQUAL(result.size(), 2)
  TEST_EQUAL(result[0], "chromatogram1")
  TEST_EQUAL(result[1], "chromatogram2")
}
END_SECTION

START_SECTION (Feature & getPrecursorFeature(String key))
{
  MRMFeature mrmfeature;
  Feature f1;
  f1.setMetaValue("dummy", 1);
  Feature f2;
  mrmfeature.addPrecursorFeature(f1, "chromatogram1");
  mrmfeature.addPrecursorFeature(f1, "chromatogram2");
  TEST_EQUAL(mrmfeature.getPrecursorFeature("chromatogram1").getMetaValue("dummy"), 1)
}
END_SECTION

START_SECTION (void MRMFeature::IDScoresAsMetaValue(bool decoy, const OpenSwath_Ind_Scores& idscores))
{

  {
    MRMFeature mrmfeature;
    OpenSwath_Ind_Scores s;
    mrmfeature.IDScoresAsMetaValue(false, s);
    TEST_EQUAL(mrmfeature.metaValueExists("id_target_transition_names"), true)
    TEST_EQUAL(mrmfeature.metaValueExists("id_decoy_transition_names"), false)
    TEST_EQUAL(mrmfeature.metaValueExists("id_target_ind_mi_score"), true)
  }

  {
    MRMFeature mrmfeature;
    OpenSwath_Ind_Scores s;
    mrmfeature.IDScoresAsMetaValue(true, s);
    TEST_EQUAL(mrmfeature.metaValueExists("id_decoy_transition_names"), true)
    TEST_EQUAL(mrmfeature.metaValueExists("id_target_transition_names"), false)
    TEST_EQUAL(mrmfeature.metaValueExists("id_decoy_ind_mi_score"), true)
  }

}
END_SECTION

START_SECTION ( void MRMFeature::scoresAsMetaValue(bool ms1only, const OpenSwath_Scores_Usage& su_) )
{

  MRMFeature mrmfeature;
  mrmfeature.getScores().xcorr_coelution_score = 5.0;
  mrmfeature.getScores().norm_rt_score = 8.0;
  {
    OpenSwath_Scores_Usage s;
    mrmfeature.scoresAsMetaValue(false, s);
    TEST_EQUAL(mrmfeature.metaValueExists("var_xcorr_coelution"), true)
    TEST_REAL_SIMILAR(mrmfeature.getMetaValue("var_xcorr_coelution"), 5.0)

    TEST_EQUAL(mrmfeature.metaValueExists("var_library_corr"), true)
    TEST_EQUAL(mrmfeature.metaValueExists("var_norm_rt_score"), true)
    TEST_REAL_SIMILAR(mrmfeature.getMetaValue("var_norm_rt_score"), 8.0)
  }

}
END_SECTION

/////////////////////////////////////////////////////////////
END_TEST
