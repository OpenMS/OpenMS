// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
START_TEST(ROCCurve, "$Id: ROCCurve_test.C,v 1.4 2006/03/28 12:53:13 marc_sturm Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

ROCCurve* rcp;

CHECK(ROCCurve::ROCCurve())
  rcp = new ROCCurve();
  TEST_NOT_EQUAL(rcp, 0)
RESULT

CHECK(ROCCurve::insertPair())
  srand( (unsigned)time( NULL ) );
  for ( uint i = 0; i < 1000; ++i )
  {
    double score = (double)rand()/RAND_MAX;
    bool clas = ( rand() > RAND_MAX/2 );
    rcp->insertPair(score,clas);
  }
RESULT

#ifdef OPENMS_HAS_CGAL 

CHECK(ROCCurve::AUC())
  double auc = rcp->AUC();
  bool inBounds = ( auc >= 0 && auc <= 1 );
  TEST_EQUAL(inBounds,1)
RESULT

#endif

CHECK(ROCCurve::curve())
  vector<pair<double,double> > curvePoints = rcp->curve(100);
  TEST_EQUAL(curvePoints.size(),100)
RESULT

CHECK(ROCCurve::cutoffPos())
  double cop = rcp->cutoffPos();
  bool inBounds( cop >=0 && cop <= 1 );
  TEST_EQUAL(inBounds,1)
RESULT

CHECK(ROCCurve::cutoffNeg())
  double con = rcp->cutoffNeg();
  bool inBounds( con >=0 && con <= 1 );
  TEST_EQUAL(inBounds,1)
RESULT

CHECK(ROCCurve::ROCCurve(const ROCCurve&))
  ROCCurve crc(*rcp);
  double ccop = crc.cutoffPos();
  double cop = rcp->cutoffPos();
  TEST_REAL_EQUAL(ccop,cop)
RESULT

CHECK(ROCCurve::operator=(const ROCCurve&))
  ROCCurve crc = *rcp;
  double ccop = crc.cutoffPos();
  double cop = rcp->cutoffPos();
  TEST_REAL_EQUAL(cop,ccop)
RESULT

CHECK(ROCCurve::~ROCCurve())
  delete rcp;
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
