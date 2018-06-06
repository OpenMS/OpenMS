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
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/ANALYSIS/PIP/PeakIntensityPredictor.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(PeakIntensityPredictor, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TOLERANCE_ABSOLUTE(0.001)

AASequence seq1 = AASequence::fromString("LTSEAR");
AASequence seq2 = AASequence::fromString("AEAQIR");
AASequence seq3 = AASequence::fromString("TLEDAR");

vector<AASequence> vec;
vec.push_back(seq1);
vec.push_back(seq2);
vec.push_back(seq3);

PeakIntensityPredictor* ptr;
PeakIntensityPredictor* nullPointer = nullptr;

START_SECTION(PeakIntensityPredictor())
  ptr = new PeakIntensityPredictor();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~PeakIntensityPredictor()))
  delete ptr;
END_SECTION

START_SECTION(double predict(const AASequence& sequence))
  PeakIntensityPredictor pip;
  TEST_REAL_SIMILAR(pip.predict(seq1), -0.531675)
  TEST_REAL_SIMILAR(pip.predict(seq2), 0.0171194)
  TEST_REAL_SIMILAR(pip.predict(seq3), -0.595362)
END_SECTION


START_SECTION(double predict(const AASequence& sequence, std::vector<double>& add_info))
  PeakIntensityPredictor pip;
  std::vector<double> add_info;
  pip.predict(seq1,add_info);
  TEST_EQUAL(add_info.size(),3)
  TEST_REAL_SIMILAR(add_info[0],0.0)
  TEST_REAL_SIMILAR(add_info[1],1.0)
  TEST_REAL_SIMILAR(add_info[2],2.04653)
END_SECTION

START_SECTION(std::vector<double> predict(const std::vector<AASequence>& sequences))
  PeakIntensityPredictor pip;
  vector<double> ref = pip.predict(vec);
  TEST_REAL_SIMILAR(ref[0], -0.531675)
  TEST_REAL_SIMILAR(ref[1], 0.0171194)
  TEST_REAL_SIMILAR(ref[2], -0.595362)
END_SECTION

START_SECTION(std::vector<double> predict(const std::vector<AASequence>& sequences, std::vector<std::vector<double> >& add_info))
  PeakIntensityPredictor pip;
  vector<vector<double> > add_info;
  pip.predict(vec,add_info);
  TEST_EQUAL(add_info.size(),3)
  TEST_EQUAL(add_info[0].size(),3)
  TEST_EQUAL(add_info[1].size(),3)
  TEST_EQUAL(add_info[2].size(),3)
  TEST_REAL_SIMILAR(add_info[0][0],0.0)
  TEST_REAL_SIMILAR(add_info[0][1],1.0)
  TEST_REAL_SIMILAR(add_info[0][2],2.04653)
  TEST_REAL_SIMILAR(add_info[1][0],0.0)
  TEST_REAL_SIMILAR(add_info[1][1],1.0)
  TEST_REAL_SIMILAR(add_info[1][2],2.30648)
  TEST_REAL_SIMILAR(add_info[2][0],0.0)
  TEST_REAL_SIMILAR(add_info[2][1],1.0)
  TEST_REAL_SIMILAR(add_info[2][2],2.24984)
END_SECTION


END_TEST


