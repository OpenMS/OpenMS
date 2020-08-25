// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

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
			
END_SECTION			

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



