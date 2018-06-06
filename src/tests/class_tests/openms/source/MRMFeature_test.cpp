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
  tmp.addScore("testscore", 200);

  MRMFeature tmp2 (tmp);

  TEST_EQUAL(tmp2.getScore("testscore"), 200)
  TEST_REAL_SIMILAR(tmp2.getIntensity(), 100.0)
}
END_SECTION

START_SECTION(MRMFeature& operator=(const MRMFeature &rhs))
{
  MRMFeature tmp;
  tmp.setIntensity(100.0);
  tmp.addScore("testscore", 200);

  MRMFeature tmp2;
  tmp2 = tmp;

  TEST_EQUAL(tmp2.getScore("testscore"), 200)
  TEST_REAL_SIMILAR(tmp2.getIntensity(), 100.0)
}
END_SECTION

START_SECTION (const PGScoresType & getScores() const)
{
  // tested with set/add score
  NOT_TESTABLE
}
END_SECTION

START_SECTION (double getScore(const String & score_name))
{
  // tested with set/add score
  NOT_TESTABLE
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

START_SECTION (void setScores(const PGScoresType & scores))
{
  MRMFeature mrmfeature;
  MRMFeature::PGScoresType scores;
  scores["score1"] = 1;
  scores["score2"] = 2;
  mrmfeature.setScores(scores);
  scores = mrmfeature.getScores();
  TEST_EQUAL(mrmfeature.getScore("score1"), 1)
  TEST_EQUAL(mrmfeature.getScore("score2"), 2)
  TEST_EQUAL(scores[String("score1")], 1)
  TEST_EQUAL(scores[String("score2")], 2)
}
END_SECTION

START_SECTION (void addScore(const String & score_name, double score))
{
  MRMFeature mrmfeature;
  mrmfeature.addScore("score1",1);
  mrmfeature.addScore("score2",2);
  MRMFeature::PGScoresType scores = mrmfeature.getScores();
  TEST_EQUAL(mrmfeature.getScore("score1"), 1)
  TEST_EQUAL(mrmfeature.getScore("score2"), 2)
  TEST_EQUAL(scores[String("score1")], 1)
  TEST_EQUAL(scores[String("score2")], 2)
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

/////////////////////////////////////////////////////////////
END_TEST
