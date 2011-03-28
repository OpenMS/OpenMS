// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/MATH/MISC/NonNegativeLeastSquaresSolver.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(NonNegativeLeastSquaresSolver, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

NonNegativeLeastSquaresSolver* ptr = 0;
NonNegativeLeastSquaresSolver* nullPointer = 0;
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
	A.setMatrix<3,4>(A_1);
	b.setMatrix<3,1>(b_1);
	x.resize(4,1);
	
	TOLERANCE_ABSOLUTE(0.0005);
	
	NonNegativeLeastSquaresSolver::solve(A,b,x);
	for (size_t i=0;i<x.rows();++i)
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
	
	A.setMatrix<4,4>(A_2);
	b.setMatrix<4,1>(b_2);
	x.resize(4,1);
	
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



