// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Volker Mosthaf, Andreas Bertsch $
// --------------------------------------------------------------------------
//



#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/ML/ROCCURVE/ROCCurve.h>

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

  // deterministic check that cutoffNeg picks a NEGATIVE-class score (issue #9488, ML-21)
  ROCCurve rc2;
  // negatives at low scores, positives at high scores
  rc2.insertPair(0.10, false);
  rc2.insertPair(0.20, false);
  rc2.insertPair(0.30, false);
  rc2.insertPair(0.40, false);
  rc2.insertPair(0.70, true);
  rc2.insertPair(0.80, true);
  rc2.insertPair(0.90, true);
  rc2.insertPair(0.95, true);
  // with the fix, cutoffNeg(0.95) walks negatives only and returns a
  // negative-class score (<= 0.40), never a positive-class score (>= 0.70)
  double con2 = rc2.cutoffNeg(0.95);
  bool fromNeg = (con2 <= 0.40 + 1e-9 && con2 >= 0.0);
  TEST_EQUAL(fromNeg, true)
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
