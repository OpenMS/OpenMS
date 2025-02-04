// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Witold Wolski $
// --------------------------------------------------------------------------

#include "OpenMS/OPENSWATHALGO/OpenSwathAlgoConfig.h"

#include "OpenMS/OPENSWATHALGO/ALGO/StatsHelpers.h"
#include "OpenMS/OPENSWATHALGO/DATAACCESS/SpectrumHelpers.h"

#include <OpenMS/CONCEPT/ClassTest.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(DIAHelpers, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(testDotProdScore)
{
	double arr1[] = { 100., 200., 4., 30., 20. };
	double arr2[] = { 100., 100., 4., 100., 200. };
	std::vector<double> vec1;
	std::vector<double> vec2;
	vec1.assign(arr1, arr1 + sizeof(arr1) / sizeof(double));
	vec2.assign(arr2, arr2 + sizeof(arr2) / sizeof(double));
	/*
	x<-c(100., 200., 4., 30., 20.)
	y<-c(100., 100., 4., 100., 200.)
	xs<-sqrt(x)
	ys<-sqrt(y)
	xsn<-xs/sqrt(sum(xs*xs))
	ysn<-ys/sqrt(sum(ys*ys))
	sum(xsn*ysn)
	*/
	//0.8604286
	double scor = OpenSwath::dotprodScoring(vec1,vec2);

	TEST_REAL_SIMILAR (scor, 0.8604286);

	//xsm <- xs/sum(xs)
	//ysm <-ys/sum(ys)
	//sum(fabs(ysm-xsm))
	scor = OpenSwath::manhattanScoring(vec1,vec2);
	TEST_REAL_SIMILAR (scor, 0.4950837);
	//0.4950837
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
