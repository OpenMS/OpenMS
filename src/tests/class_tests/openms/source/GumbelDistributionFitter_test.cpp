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
// $Maintainer: Timo Sachsenberg $
// $Authors: David Wojnar$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/MATH/STATISTICS/GumbelDistributionFitter.h>
///////////////////////////

using namespace OpenMS;
using namespace Math;
using namespace std;

START_TEST(GumbelDistributionFitter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

GumbelDistributionFitter* ptr = nullptr;
GumbelDistributionFitter* nullPointer = nullptr;
START_SECTION(GumbelDistributionFitter())
{
	ptr = new GumbelDistributionFitter();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((virtual ~GumbelDistributionFitter()))
{
  delete ptr;
	NOT_TESTABLE
}
END_SECTION

START_SECTION((GumbelDistributionFitResult fit(std::vector<DPosition<2> >& points)))
{

	DPosition<2> pos;
  vector<DPosition<2> > points;

	pos.setX(-2.7); pos.setY(0.017); points.push_back(pos);
	pos.setX(-2.5); pos.setY(0.025); points.push_back(pos);
	pos.setX(-2); pos.setY(0.052); points.push_back(pos);
	pos.setX(-1); pos.setY(0.127); points.push_back(pos);
	pos.setX(-0.7); pos.setY(0.147); points.push_back(pos);
	pos.setX(-0.01); pos.setY(0.178); points.push_back(pos);
	pos.setX(0); pos.setY(0.178); points.push_back(pos);
	pos.setX(0.2); pos.setY(0.182); points.push_back(pos);
	pos.setX(0.5); pos.setY(0.184); points.push_back(pos);
	pos.setX(1); pos.setY(0.179); points.push_back(pos);
	pos.setX(1.3); pos.setY(0.171); points.push_back(pos);
	pos.setX(1.9); pos.setY(0.151); points.push_back(pos);
	pos.setX(2.5); pos.setY(0.127); points.push_back(pos);
	pos.setX(2.6); pos.setY(0.123); points.push_back(pos);
	pos.setX(2.7); pos.setY(0.119); points.push_back(pos);
	pos.setX(2.8); pos.setY(0.115); points.push_back(pos);
	pos.setX(2.9); pos.setY(0.111); points.push_back(pos);
	pos.setX(3); pos.setY(0.108); points.push_back(pos);
	pos.setX(3.5); pos.setY(0.089); points.push_back(pos);
	pos.setX(3.9); pos.setY(0.076); points.push_back(pos);
	pos.setX(4.01); pos.setY(0.073); points.push_back(pos);
	pos.setX(4.22); pos.setY(0.067); points.push_back(pos);
	pos.setX(4.7); pos.setY(0.054); points.push_back(pos);
	pos.setX(4.9); pos.setY(0.05); points.push_back(pos);
	pos.setX(5); pos.setY(0.047); points.push_back(pos);
	pos.setX(6); pos.setY(0.03); points.push_back(pos);
	pos.setX(7); pos.setY(0.017); points.push_back(pos);
	pos.setX(7.5); pos.setY(0.015); points.push_back(pos);
	pos.setX(7.9); pos.setY(0.012); points.push_back(pos);
	pos.setX(8.03); pos.setY(0.011); points.push_back(pos);
	//a= 0.5, b = 2
	

	ptr = new GumbelDistributionFitter;
		GumbelDistributionFitter::GumbelDistributionFitResult init_param;
	init_param.a = 1.0;
	init_param.b = 3.0;
	ptr->setInitialParameters(init_param);
	GumbelDistributionFitter::GumbelDistributionFitResult result = ptr->fit(points);

	TOLERANCE_ABSOLUTE(0.1)
	TEST_REAL_SIMILAR(result.a, 0.5)
	TEST_REAL_SIMILAR(result.b, 2.0)	
	
	vector<DPosition<2> > points2;
	pos.setX(0); pos.setY(0.18); points2.push_back(pos);
	pos.setX(0.2); pos.setY(0.24); points2.push_back(pos);
	pos.setX(0.5); pos.setY(0.32); points2.push_back(pos);
	pos.setX(1); pos.setY(0.37); points2.push_back(pos);
	pos.setX(1.3); pos.setY(0.35); points2.push_back(pos);
	pos.setX(1.9); pos.setY(0.27); points2.push_back(pos);
	pos.setX(2.5); pos.setY(0.18); points2.push_back(pos);
	pos.setX(2.6); pos.setY(0.16); points2.push_back(pos);
	pos.setX(3); pos.setY(0.12); points2.push_back(pos);
	pos.setX(5); pos.setY(0.02); points2.push_back(pos);
	//a = 1, b = 1
	
	init_param.a = 3.0;
	init_param.b = 3.0;
	ptr->setInitialParameters(init_param);
	GumbelDistributionFitter::GumbelDistributionFitResult result2 = ptr->fit(points2);

	TOLERANCE_ABSOLUTE(0.1)
	TEST_REAL_SIMILAR(result2.a, 1.0)
	TEST_REAL_SIMILAR(result2.b, 1.0)	
	
}
END_SECTION

START_SECTION((void setInitialParameters(const GumbelDistributionFitResult& result)))
{
  GumbelDistributionFitter f1;
  GumbelDistributionFitter::GumbelDistributionFitResult result;
  f1.setInitialParameters(result);
	
	NOT_TESTABLE //implicitly tested in fit method
}
END_SECTION

START_SECTION((GumbelDistributionFitter(const GumbelDistributionFitter& rhs)))
NOT_TESTABLE
END_SECTION

START_SECTION((GumbelDistributionFitter& operator = (const GumbelDistributionFitter& rhs)))
NOT_TESTABLE
END_SECTION
GumbelDistributionFitter::GumbelDistributionFitResult* p = nullptr;
START_SECTION((GumbelDistributionFitter::GumbelDistributionFitResult()))
p =  new GumbelDistributionFitter::GumbelDistributionFitResult;
	TEST_NOT_EQUAL(ptr, nullPointer)
TEST_REAL_SIMILAR(p->a, 1.0)
TEST_REAL_SIMILAR(p->b, 2.0)
END_SECTION

START_SECTION((GumbelDistributionFitter::GumbelDistributionFitResult(const GumbelDistributionFitter::GumbelDistributionFitResult& rhs)))
p-> a = 5.0;
p->b = 4.0;
GumbelDistributionFitter::GumbelDistributionFitResult obj(*p);
TEST_REAL_SIMILAR(obj.a, 5.0)
TEST_REAL_SIMILAR(obj.b, 4.0)
END_SECTION

START_SECTION((GumbelDistributionFitter::GumbelDistributionFitResult& operator = (const GumbelDistributionFitter::GumbelDistributionFitResult& rhs)))
p-> a = 3.0;
p->b = 2.2;
GumbelDistributionFitter::GumbelDistributionFitResult obj = *p;
TEST_REAL_SIMILAR(obj.a, 3.0)
TEST_REAL_SIMILAR(obj.b, 2.2)
END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



