// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

// /////////////////////////

#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <iostream>
#include <vector>

using namespace OpenMS;
using namespace OpenMS::Math;
using namespace std;

// /////////////////////////

START_TEST(Distribution, "$Id$")

// ///////////////////////////////////////////////////////////

START_SECTION((ceilDecimal))
	TEST_REAL_SIMILAR(ceilDecimal(12345.671,-2),12345.68)
	TEST_REAL_SIMILAR(ceilDecimal(12345.67,-1),12345.7)
	TEST_REAL_SIMILAR(ceilDecimal(12345.67,0),12346.0)
	TEST_REAL_SIMILAR(ceilDecimal(12345.67,1),12350.0)
	TEST_REAL_SIMILAR(ceilDecimal(12345.67,2),12400.0)
END_SECTION

START_SECTION((roundDecimal))
	TEST_REAL_SIMILAR(roundDecimal(12345.671,-2),12345.67)
	TEST_REAL_SIMILAR(roundDecimal(12345.67,-1),12345.7)
	TEST_REAL_SIMILAR(roundDecimal(12345.67,0),12346.0)
	TEST_REAL_SIMILAR(roundDecimal(12345.67,1),12350.0)
	TEST_REAL_SIMILAR(roundDecimal(12345.67,2),12300.0)
END_SECTION

START_SECTION((intervalTransformation))
	TEST_REAL_SIMILAR(intervalTransformation(0.5,0.0,1.0,0.0,600.0),300.0)
	TEST_REAL_SIMILAR(intervalTransformation(0.5,0.25,1.0,0.0,600.0),200.0)
	TEST_REAL_SIMILAR(intervalTransformation(0.5,0.0,0.75,0.0,600.0),400.0)
	TEST_REAL_SIMILAR(intervalTransformation(0.5,0.0,1.0,150.0,600.0),375.0)
	TEST_REAL_SIMILAR(intervalTransformation(0.5,0.0,1.0,0.0,450.0),225.0)
END_SECTION 

START_SECTION((linear2log))
	TEST_REAL_SIMILAR(linear2log(0.0),0.0)
	TEST_REAL_SIMILAR(linear2log(9.0),1.0)
	TEST_REAL_SIMILAR(linear2log(99.0),2.0)
	TEST_REAL_SIMILAR(linear2log(999.0),3.0)
END_SECTION

START_SECTION((log2linear))
	TEST_REAL_SIMILAR(log2linear(0.0),0.0)
	TEST_REAL_SIMILAR(log2linear(1.0),9.0)
	TEST_REAL_SIMILAR(log2linear(2.0),99.0)
	TEST_REAL_SIMILAR(log2linear(3.0),999.0)
END_SECTION

START_SECTION((isOdd))
	TEST_EQUAL(isOdd(0),false)
	TEST_EQUAL(isOdd(1),true)
	TEST_EQUAL(isOdd(2),false)
	TEST_EQUAL(isOdd(3),true)
END_SECTION

START_SECTION((template <typename T> T round (T x)))
	float f_down=14.49f;		 // expected 14
	float f_up = 14.50f;		 // expected 15
	double d_up = -999.49;   // expected -999
	double d_down = -675.77; // expected -676
	TEST_REAL_SIMILAR(round(f_down), 14.0)
	TEST_REAL_SIMILAR(round(f_up), 15.0)
	TEST_REAL_SIMILAR(round(d_up), -999)
	TEST_REAL_SIMILAR(round(d_down), -676)
END_SECTION


START_SECTION((bool approximatelyEqual(DoubleReal a, DoubleReal b, DoubleReal tol)))
	TEST_EQUAL(approximatelyEqual(1.1, 1.1002, 0.1), true)
	TEST_EQUAL(approximatelyEqual(1.1, 1.1002, 0.01), true)
	TEST_EQUAL(approximatelyEqual(1.1, 1.1002, 0.001), true)
	TEST_EQUAL(approximatelyEqual(1.1, 1.1002, 0.0001), false)
END_SECTION

/////////////////////////////////////////////////////////////);
/////////////////////////////////////////////////////////////
END_TEST
