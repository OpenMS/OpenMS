// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ML/NNLS/NonNegativeLeastSquaresSolver.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(NonNegativeLeastSquaresSolver, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

NonNegativeLeastSquaresSolver* ptr = nullptr;
NonNegativeLeastSquaresSolver* nullPointer = nullptr;
START_SECTION(NonNegativeLeastSquaresSolver())
{
	ptr = new NonNegativeLeastSquaresSolver();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~NonNegativeLeastSquaresSolver())
{
	delete ptr;
}
END_SECTION

START_SECTION((static Int solve(const Matrix< double > &A, const Matrix< double > &b, Matrix< double > &x)))
{
	
	// CASE 1	
	double A_1[3][4] = 
		{
			{1, 10, 4, 10},
			{4, 5 , 1, 12},
			{5, 1 , 9, 20},
		};
	double b_1[3][1] = {{4},{7},{4}};
	double x_1[4][1] = {{0.931153},{0.36833},{0},{0}};

	Matrix<double> A,b,x;
	A.setMatrix<double,3,4>(A_1);
	b.setMatrix<double,3,1>(b_1);
	x.getEigenMatrix().resize(4,1);
	
	TOLERANCE_ABSOLUTE(0.0005);
	
	NonNegativeLeastSquaresSolver::solve(A,b,x);
	for (size_t i = 0;i < x.rows(); ++i)
	{
		TEST_REAL_SIMILAR(x(i,0), x_1[i][0]);
	}	
	
	
		
	// CASE 2
	double A_2[4][4] = 
		{
			{0.9290,    0.0200,         0,         0},
			{0.0590,    0.9230,    0.0300,    0.0010},
			{0.0020,    0.0560,    0.9240,    0.0400},
			{		  0,    0.0010,    0.0450,    0.9240}
		};
	double b_2[4][1] = {{5},{45},{4},{31}};
	double x_2[4][1] = {{4.3395},{48.4364},{0},{33.4945}};	
	
	A.setMatrix<double,4,4>(A_2);
	b.setMatrix<double,4,1>(b_2);
	x.getEigenMatrix().resize(4,1);
	
	NonNegativeLeastSquaresSolver::solve(A,b,x);
	for (size_t i=0;i<x.rows();++i)
	{
		TEST_REAL_SIMILAR(x(i,0), x_2[i][0]);
	}	
	
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



