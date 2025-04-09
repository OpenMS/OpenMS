// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ML/NNLS/NNLS.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(NNLS, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

/*NNLS* ptr = 0;
NNLS* null_ptr = 0;
START_SECTION(NNLS())
{
	ptr = new NNLS();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~NNLS())
{
	delete ptr;
}
END_SECTION
*/

START_SECTION([EXTRA]int nnls_(double *a, integer *mda, integer *m, integer *n, double *b, double *x, double *rnorm, double *w, double *zz, integer *index, integer *mode))

	// translate A to array a (column major order)
	double a_vec[4]= {1, 0, 0, 1};
	int a_rows = 2;
	int a_cols = 2;
	
	// translate b
	double b_vec[2] = {2, 3};
	
	// prepare solution array (directly copied from example)
	double *x_vec = new double[2+1];
	double rnorm;
	double *w = new double[2+1];
	double *zz = new double[2+1];
	int *indx = new int[2+1];
	int mode;
	
	NNLS::nnls_(a_vec, &a_rows, &a_rows, &a_cols, b_vec, x_vec, &rnorm, w, zz, indx, &mode);

	TEST_EQUAL(mode, 1)
	
	double x_solution[2] = {2, 3};
	for (Size i=0; i<2; ++i)
	{
		TEST_EQUAL(x_vec[i], x_solution[i])
	}
	delete[] x_vec;
	delete[] w;
	delete[] zz;
	delete[] indx;
			
END_SECTION			

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



