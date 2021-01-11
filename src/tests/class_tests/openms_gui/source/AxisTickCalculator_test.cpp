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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/VISUAL/AxisTickCalculator.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
///////////////////////////

using namespace OpenMS;

START_TEST(AxisTickCalculator, "$Id$")

/////////////////////////////////////////////////////////////


START_SECTION((static void calcGridLines(double x1, double x2, GridVector& grid)))
	std::vector<std::vector<double> > vector1;
  AxisTickCalculator::calcGridLines(1.0,4.0,vector1);

  TEST_EQUAL(2,vector1.size());
	TEST_EQUAL(4,vector1[0].size());

	TEST_REAL_SIMILAR(1.0,vector1[0][0]);
	TEST_REAL_SIMILAR(2.0,vector1[0][1]);
	TEST_REAL_SIMILAR(3.0,vector1[0][2]);
	TEST_REAL_SIMILAR(4.0,vector1[0][3]);

END_SECTION

START_SECTION((static void calcLogGridLines(double x1, double x2, GridVector& grid)))
	std::vector<std::vector<double> > vector1;
	AxisTickCalculator::calcLogGridLines(log10(10.0),log10(100.0),vector1);

	TEST_EQUAL(1,vector1[0].size());
	TEST_EQUAL(8,vector1[1].size());

	TEST_REAL_SIMILAR(1.0,vector1[0][0]);

	TEST_REAL_SIMILAR( 1.30103 ,vector1[1][0]);
	TEST_REAL_SIMILAR( 1.47712 ,vector1[1][1]);
	TEST_REAL_SIMILAR( 1.60206 ,vector1[1][2]);
	TEST_REAL_SIMILAR( 1.69897 ,vector1[1][3]);
	TEST_REAL_SIMILAR( 1.77815 ,vector1[1][4]);
	TEST_REAL_SIMILAR( 1.8451 ,vector1[1][5]);
	TEST_REAL_SIMILAR( 1.90309 ,vector1[1][6]);
	TEST_REAL_SIMILAR( 1.95425 ,vector1[1][7]);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



