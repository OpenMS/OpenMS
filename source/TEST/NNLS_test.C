// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/MATH/MISC/NNLS/NNLS.h>
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

START_SECTION([EXTRA]int nnls_(doublereal *a, integer *mda, integer *m, integer *n, doublereal *b, doublereal *x, doublereal *rnorm, doublereal *w, doublereal *zz, integer *index, integer *mode))

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
			
END_SECTION			

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



