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
// $Maintainer: Hendrik Weisser $
// $Authors: Marc Sturm, Andreas Bertsch, Sven Nahnsen, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmWorst.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(ConsensusIDAlgorithmWorst, "$Id$")

/////////////////////////////////////////////////////////////

ConsensusIDAlgorithm* ptr = nullptr;
ConsensusIDAlgorithm* null_pointer = nullptr;
START_SECTION(ConsensusIDAlgorithmWorst())
{
  ptr = new ConsensusIDAlgorithmWorst();
  TEST_NOT_EQUAL(ptr, null_pointer);
}
END_SECTION

START_SECTION(~ConsensusIDAlgorithmWorst())
{
  delete(ptr);
}
END_SECTION

// create 3 ID runs:
PeptideIdentification temp;
temp.setScoreType("Posterior Error Probability");
temp.setHigherScoreBetter(false);
vector<PeptideIdentification> ids(3, temp);
vector<PeptideHit> hits;
// the first ID has 5 hits
hits.resize(5);
hits[0].setSequence(AASequence::fromString("A"));
hits[0].setScore(0.1);
hits[1].setSequence(AASequence::fromString("B"));
hits[1].setScore(0.2);
hits[2].setSequence(AASequence::fromString("C"));
hits[2].setScore(0.3);
hits[3].setSequence(AASequence::fromString("D"));
hits[3].setScore(0.4);
hits[4].setSequence(AASequence::fromString("E"));
hits[4].setScore(0.5);
ids[0].setHits(hits);
// the second ID has 3 hits
hits.resize(3);
hits[0].setSequence(AASequence::fromString("C"));
hits[0].setScore(0.2);
hits[1].setSequence(AASequence::fromString("A"));
hits[1].setScore(0.4);
hits[2].setSequence(AASequence::fromString("B"));
hits[2].setScore(0.6);
ids[1].setHits(hits);
// the third ID has 10 hits
hits.resize(10);
hits[0].setSequence(AASequence::fromString("F"));
hits[0].setScore(0.0);
hits[1].setSequence(AASequence::fromString("C"));
hits[1].setScore(0.1);
hits[2].setSequence(AASequence::fromString("G"));
hits[2].setScore(0.2);
hits[3].setSequence(AASequence::fromString("D"));
hits[3].setScore(0.3);
hits[4].setSequence(AASequence::fromString("B"));
hits[4].setScore(0.4);
hits[5].setSequence(AASequence::fromString("E"));
hits[5].setScore(0.5);
hits[6].setSequence(AASequence::fromString("H"));
hits[6].setScore(0.6);
hits[7].setSequence(AASequence::fromString("I"));
hits[7].setScore(0.7);
hits[8].setSequence(AASequence::fromString("J"));
hits[8].setScore(0.8);
hits[9].setSequence(AASequence::fromString("K"));
hits[9].setScore(0.9);
ids[2].setHits(hits);

START_SECTION(void apply(std::vector<PeptideIdentification>& ids))
{
  TOLERANCE_ABSOLUTE(0.01)

  ConsensusIDAlgorithmWorst consensus;
  // define parameters:
  Param param;
  param.setValue("filter:considered_hits", 0);
  consensus.setParameters(param);
  // apply:
  vector<PeptideIdentification> f = ids;
  consensus.apply(f);

  TEST_EQUAL(f.size(), 1);
  hits = f[0].getHits();
  TEST_EQUAL(hits.size(), 11);

  TEST_EQUAL(hits[0].getRank(), 1);
  TEST_EQUAL(hits[0].getSequence(), AASequence::fromString("F"));
  TEST_REAL_SIMILAR(hits[0].getScore(), 0.0);

  TEST_EQUAL(hits[1].getRank(), 2);
  TEST_EQUAL(hits[1].getSequence(), AASequence::fromString("G"));
  TEST_REAL_SIMILAR(hits[1].getScore(), 0.2);

  TEST_EQUAL(hits[2].getRank(), 3);
  TEST_EQUAL(hits[2].getSequence(), AASequence::fromString("C"));
  TEST_REAL_SIMILAR(hits[2].getScore(), 0.3);

  // hits with the same score get assigned the same rank:
  TEST_EQUAL(hits[3].getRank(), 4);
  TEST_EQUAL(hits[3].getSequence(), AASequence::fromString("A"));
  TEST_REAL_SIMILAR(hits[3].getScore(), 0.4);

  TEST_EQUAL(hits[4].getRank(), 4);
  TEST_EQUAL(hits[4].getSequence(), AASequence::fromString("D"));
  TEST_REAL_SIMILAR(hits[4].getScore(), 0.4);

  TEST_EQUAL(hits[5].getRank(), 5);
  TEST_EQUAL(hits[5].getSequence(), AASequence::fromString("E"));
  TEST_REAL_SIMILAR(hits[5].getScore(), 0.5);

  TEST_EQUAL(hits[6].getRank(), 6);
  TEST_EQUAL(hits[6].getSequence(), AASequence::fromString("B"));
  TEST_REAL_SIMILAR(hits[6].getScore(), 0.6);

  TEST_EQUAL(hits[7].getRank(), 6);
  TEST_EQUAL(hits[7].getSequence(), AASequence::fromString("H"));
  TEST_REAL_SIMILAR(hits[7].getScore(), 0.6);

  TEST_EQUAL(hits[8].getRank(), 7);
  TEST_EQUAL(hits[8].getSequence(), AASequence::fromString("I"));
  TEST_REAL_SIMILAR(hits[8].getScore(), 0.7);

  TEST_EQUAL(hits[9].getRank(), 8);
  TEST_EQUAL(hits[9].getSequence(), AASequence::fromString("J"));
  TEST_REAL_SIMILAR(hits[9].getScore(), 0.8);

  TEST_EQUAL(hits[10].getRank(), 9);
  TEST_EQUAL(hits[10].getSequence(), AASequence::fromString("K"));
  TEST_REAL_SIMILAR(hits[10].getScore(), 0.9);


  ids[2].setHigherScoreBetter(true);
  TEST_EXCEPTION(Exception::InvalidValue, consensus.apply(ids));
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
