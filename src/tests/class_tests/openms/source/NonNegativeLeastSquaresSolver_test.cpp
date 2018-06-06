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
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/MATH/MISC/NonNegativeLeastSquaresSolver.h>
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



