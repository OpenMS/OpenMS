// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//



#include <OpenMS/CONCEPT/ClassTest.h>

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

ROCCurve* rcp;

CHECK(ROCCurve())
  rcp = new ROCCurve();
  TEST_NOT_EQUAL(rcp, 0)
RESULT

CHECK(void insertPair(double score, bool clas))
  srand((unsigned)time(NULL));
  for (uint i = 0; i < 1000; ++i)
  {
    double score = (double)rand()/RAND_MAX;
    bool clas = (rand() > RAND_MAX/2);
    rcp->insertPair(score, clas);
  }
RESULT

#ifdef OPENMS_HAS_CGAL 

CHECK(double AUC())
  double auc = rcp->AUC();
  bool inBounds = ( auc >= 0 && auc <= 1 );
  TEST_EQUAL(inBounds,1)
RESULT

#endif

CHECK((std::vector<std::pair<double, double> curve(uint resolution = 10)))
  vector<pair<double,double> > curvePoints = rcp->curve(100);
  TEST_EQUAL(curvePoints.size(),100)
RESULT

CHECK(double cutoffPos(double fraction = 0.95))
  double cop = rcp->cutoffPos();
  bool inBounds( cop >=0 && cop <= 1 );
  TEST_EQUAL(inBounds,1)
RESULT

CHECK(double cutoffNeg(double fraction = 0.95))
  double con = rcp->cutoffNeg();
  bool inBounds( con >=0 && con <= 1 );
  TEST_EQUAL(inBounds,1)
RESULT

CHECK(ROCCurve(const ROCCurve& source))
  ROCCurve crc(*rcp);
  double ccop = crc.cutoffPos();
  double cop = rcp->cutoffPos();
  TEST_REAL_EQUAL(ccop,cop)
RESULT

CHECK(ROCCurve& operator = (const ROCCurve& source))
  ROCCurve crc = *rcp;
  double ccop = crc.cutoffPos();
  double cop = rcp->cutoffPos();
  TEST_REAL_EQUAL(cop,ccop)
RESULT

CHECK((bool operator () (const std::pair<double, bool>& a, const std::pair<double, bool>& b)))
	// this is a private predicate; this test is just to avoid the checker warnings
RESULT


CHECK(~ROCCurve())
  delete rcp;
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
