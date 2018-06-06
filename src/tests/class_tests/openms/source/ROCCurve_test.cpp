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
// $Authors: Volker Mosthaf, Andreas Bertsch $
// --------------------------------------------------------------------------
//



#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/MATH/STATISTICS/ROCCurve.h>

///////////////////////////

#include <cmath>
#include <ctime>
#include <vector>

///////////////////////////
START_TEST(ROCCurve, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;
using namespace OpenMS::Math;

ROCCurve* rcp = nullptr;
ROCCurve* rcp_nullPointer = nullptr;

START_SECTION((ROCCurve()))
  rcp = new ROCCurve();
  TEST_NOT_EQUAL(rcp, rcp_nullPointer)
END_SECTION

START_SECTION((void insertPair(double score, bool clas)))
  srand((unsigned)time(nullptr));
  for (Size i = 0; i < 1000; ++i)
  {
    double score = (double)rand()/RAND_MAX;
    bool clas = (rand() > RAND_MAX/2);
    rcp->insertPair(score, clas);
  }
	NOT_TESTABLE
END_SECTION

START_SECTION((double AUC()))
  // random test:
  double auc = rcp->AUC();
  bool inBounds = ( auc >= 0 && auc <= 1 );
  TEST_EQUAL(inBounds,true)

  // some real data:
  ROCCurve rc;
  TEST_EQUAL(rc.AUC(),0.5)
END_SECTION

START_SECTION((std::vector<std::pair<double, double> > curve(UInt resolution = 10)))
  vector<pair<double,double> > curvePoints = rcp->curve(100);
  TEST_EQUAL(curvePoints.size(),100)
END_SECTION

START_SECTION((double cutoffPos(double fraction=0.95)))
  double cop = rcp->cutoffPos();
  bool inBounds( cop >=0 && cop <= 1 );
  TEST_EQUAL(inBounds,true)
END_SECTION

START_SECTION((double cutoffNeg(double fraction=0.95)))
  double con = rcp->cutoffNeg();
  bool inBounds( con >=0 && con <= 1 );
  TEST_EQUAL(inBounds,true)
END_SECTION

START_SECTION((ROCCurve(const ROCCurve& source)))
  ROCCurve crc(*rcp);
  double ccop = crc.cutoffPos();
  double cop = rcp->cutoffPos();
  TEST_REAL_SIMILAR(ccop,cop)
END_SECTION

START_SECTION((ROCCurve& operator = (const ROCCurve& source)))
  ROCCurve crc = *rcp;
  double ccop = crc.cutoffPos();
  double cop = rcp->cutoffPos();
  TEST_REAL_SIMILAR(cop,ccop)
END_SECTION

START_SECTION((virtual ~ROCCurve()))
  delete rcp;
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
