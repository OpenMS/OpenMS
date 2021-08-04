// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
