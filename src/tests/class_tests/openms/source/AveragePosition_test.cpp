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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/MATH/STATISTICS/AveragePosition.h>
///////////////////////////

using namespace OpenMS;
using namespace std;
using OpenMS::Math::AveragePosition;

START_TEST(AveragePosition, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

AveragePosition<3>* ptr = nullptr;
AveragePosition<3>* nullPointer = nullptr;
START_SECTION(AveragePosition())
{
	ptr = new AveragePosition<3>();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~AveragePosition())
{
	delete ptr;
}
END_SECTION

START_SECTION((AveragePosition(AveragePosition const &rhs)))
{
	DPosition<4> pos1;
	pos1[0] = 1.0;
	pos1[1] = 2.0;
	pos1[2] = 3.0;
	pos1[3] = 4.0;

	DPosition<4> pos2;
	pos2[0] = 13.0;
	pos2[1] = 10.0;
	pos2[2] = 7.0;
	pos2[3] = 4.0;

	AveragePosition<4> avg;
	avg.add(pos1,6);
	avg.add(pos2);

	AveragePosition<4> avg_cpy = avg;

	TEST_REAL_SIMILAR(avg.getPosition()[0],avg_cpy.getPosition()[0]);
	TEST_REAL_SIMILAR(avg.getPosition()[1],avg_cpy.getPosition()[1]);
	TEST_REAL_SIMILAR(avg.getPosition()[2],avg_cpy.getPosition()[2]);
	TEST_REAL_SIMILAR(avg.getPosition()[3],avg_cpy.getPosition()[3]);
	TEST_REAL_SIMILAR(avg.getWeight(),avg_cpy.getWeight());
}
END_SECTION

START_SECTION((PositionType const& getPosition() const))
{
	DPosition<4> pos1;
	pos1[0] = 1.0;
	pos1[1] = 2.0;
	pos1[2] = 3.0;
	pos1[3] = 4.0;

	DPosition<4> pos2;
	pos2[0] = 13.0;
	pos2[1] = 10.0;
	pos2[2] = 7.0;
	pos2[3] = 4.0;

	AveragePosition<4> avg;
	avg.add(pos1,-1);
	avg.add(pos2);

	TEST_REAL_SIMILAR(avg.getPosition()[0],0);
	TEST_REAL_SIMILAR(avg.getPosition()[1],0);
	TEST_REAL_SIMILAR(avg.getPosition()[2],0);
	TEST_REAL_SIMILAR(avg.getPosition()[3],0);
	TEST_REAL_SIMILAR(avg.getWeight(),0);

	avg.add(pos1,4);
	avg.add(pos2);

	TEST_REAL_SIMILAR(avg.getPosition()[0],5.8);
	TEST_REAL_SIMILAR(avg.getPosition()[1],5.2);
	TEST_REAL_SIMILAR(avg.getPosition()[2],4.6);
	TEST_REAL_SIMILAR(avg.getPosition()[3],4);
	TEST_REAL_SIMILAR(avg.getWeight(),5);
}
END_SECTION

START_SECTION((CoordinateType const& getWeight() const))
{
	AveragePosition<1> avg;
	avg.add(DPosition<1>(9),2);
	TEST_REAL_SIMILAR(avg.getWeight(),2);
	TEST_REAL_SIMILAR(avg.getPosition()[0],9);
	avg.add(DPosition<1>(9),3);
	TEST_REAL_SIMILAR(avg.getWeight(),5);
	TEST_REAL_SIMILAR(avg.getPosition()[0],9);
	avg.add(DPosition<1>(6),10);
	TEST_REAL_SIMILAR(avg.getWeight(),15);
	TEST_REAL_SIMILAR(avg.getPosition()[0],7);
}
END_SECTION

START_SECTION((void clear()))
{
	DPosition<4> pos1;
	pos1[0] = 1.0;
	pos1[1] = 2.0;
	pos1[2] = 3.0;
	pos1[3] = 4.0;
	AveragePosition<4> avg;
	avg.add(pos1,2);
	TEST_EQUAL(avg.getPosition(),pos1);
	TEST_REAL_SIMILAR(avg.getWeight(),2);
	avg.clear();
	TEST_EQUAL(avg.getPosition(),DPosition<4>::zero());
	TEST_EQUAL(avg.getWeight(),0);
}
END_SECTION

START_SECTION((void add(PositionType position, CoordinateType const weight=1)))
{
	// already tested above
	NOT_TESTABLE;
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



