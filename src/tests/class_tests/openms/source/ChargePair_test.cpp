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
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/ChargePair.h>
///////////////////////////

#include <OpenMS/DATASTRUCTURES/Compomer.h>
#include <OpenMS/DATASTRUCTURES/Adduct.h>

using namespace OpenMS;
using namespace std;

START_TEST(ChargePair, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ChargePair* ptr = nullptr;
ChargePair* nullPointer = nullptr;
START_SECTION(ChargePair())
{
	ptr = new ChargePair();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~ChargePair())
{
	delete ptr;
}
END_SECTION

Compomer cmp;
cmp.setID(99);

START_SECTION((ChargePair(const Size &index0, const Size &index1, const Int &charge0, const Int &charge1, const Compomer &compomer, const double &mass_diff, const bool active)))
{
	ChargePair cp(34,45, 4,5, cmp, 12.34, false);
	TEST_EQUAL(cp.getElementIndex(0), 34);
	TEST_EQUAL(cp.getElementIndex(1), 45);	
	TEST_EQUAL(cp.getCharge(0), 4);
	TEST_EQUAL(cp.getCharge(1), 5);	
	TEST_EQUAL(cp.getCompomer(), cmp);
	TEST_REAL_SIMILAR(cp.getMassDiff(), 12.34);	
	TEST_EQUAL(cp.isActive(), false);	
}
END_SECTION

START_SECTION((ChargePair(const ChargePair &rhs)))
{
	ChargePair cp2(34,45, 4,5, cmp, 12.34, false);
	ChargePair cp (cp2);
	TEST_EQUAL(cp.getElementIndex(0), 34);
	TEST_EQUAL(cp.getElementIndex(1), 45);	
	TEST_EQUAL(cp.getCharge(0), 4);
	TEST_EQUAL(cp.getCharge(1), 5);	
	TEST_EQUAL(cp.getCompomer(), cmp);
	TEST_REAL_SIMILAR(cp.getMassDiff(), 12.34);
	TEST_EQUAL(cp.getEdgeScore(), 1);	
	TEST_EQUAL(cp.isActive(), false);	
}
END_SECTION

START_SECTION((ChargePair& operator=(const ChargePair &rhs)))
{
	ChargePair cp2(34,45, 4,5, cmp, 12.34, false);
	ChargePair cp = cp2;
	TEST_EQUAL(cp.getElementIndex(0), 34);
	TEST_EQUAL(cp.getElementIndex(1), 45);	
	TEST_EQUAL(cp.getCharge(0), 4);
	TEST_EQUAL(cp.getCharge(1), 5);	
	TEST_EQUAL(cp.getCompomer(), cmp);
	TEST_REAL_SIMILAR(cp.getMassDiff(), 12.34);
	TEST_EQUAL(cp.getEdgeScore(), 1);
	TEST_EQUAL(cp.isActive(), false);	
}
END_SECTION


START_SECTION((Int getCharge(UInt pairID) const ))
{
	NOT_TESTABLE //well.. tested below...	
}
END_SECTION

START_SECTION((void setCharge(UInt pairID, Int e)))
{
  ChargePair cp;
	cp.setCharge(0,123);
	cp.setCharge(1,321);
	TEST_EQUAL(cp.getCharge(0), 123)
	TEST_EQUAL(cp.getCharge(1), 321)
}
END_SECTION

START_SECTION((Size getElementIndex(UInt pairID) const ))
{
	NOT_TESTABLE //well.. tested below...
}
END_SECTION

START_SECTION((void setElementIndex(UInt pairID, Size e)))
{
  ChargePair cp;
	cp.setElementIndex(0,123);
	cp.setElementIndex(1,321);
	TEST_EQUAL(cp.getElementIndex(0), 123)
	TEST_EQUAL(cp.getElementIndex(1), 321)
}
END_SECTION

START_SECTION((const Compomer& getCompomer() const))
{
	NOT_TESTABLE //well.. tested below...
}
END_SECTION

START_SECTION((void setCompomer(const Compomer &compomer)))
{
  ChargePair cp;
	cp.setCompomer(cmp);
	TEST_EQUAL(cp.getCompomer(), cmp)
}
END_SECTION

START_SECTION((double getMassDiff() const))
{
	NOT_TESTABLE //well.. tested below...
}
END_SECTION

START_SECTION((void setMassDiff(double mass_diff)))
{
  ChargePair cp;
	cp.setMassDiff(123.432);
	TEST_REAL_SIMILAR(cp.getMassDiff(), 123.432)
}
END_SECTION

START_SECTION((double getEdgeScore() const))
{
	NOT_TESTABLE //well.. tested below...
}
END_SECTION

START_SECTION((void setEdgeScore(double score)))
{
  ChargePair cp;
	cp.setEdgeScore(1123.432f);
	TEST_REAL_SIMILAR(cp.getEdgeScore(), 1123.432)
}
END_SECTION
		

START_SECTION((bool isActive() const))
{
	NOT_TESTABLE //well.. tested below...
}
END_SECTION

START_SECTION((void setActive(const bool active)))
{
  ChargePair cp;
	cp.setActive(true);
	TEST_EQUAL(cp.isActive(), true)
	cp.setActive(false);
	TEST_EQUAL(cp.isActive(), false)
}
END_SECTION

START_SECTION((virtual bool operator==(const ChargePair &i) const))
{
	ChargePair cp1(34,45, 4,5, cmp, 12.34, false);
	ChargePair cp2(34,15, 4,5, cmp, 12.34, false);
	TEST_EQUAL(cp1==cp2, false);
	ChargePair cp3(34,15, 4,5, cmp, 12.34, true);
	ChargePair cp4(34,15, 4,5, cmp, 12.34, false);
	TEST_EQUAL(cp3==cp4, false);
	ChargePair cp5(34,15, 4,5, cmp, 12.34, false);
	ChargePair cp6(34,15, 4,5, cmp, 12.34, false);
	TEST_EQUAL(cp5==cp6, true);
	
}
END_SECTION

START_SECTION((virtual bool operator!=(const ChargePair &i) const))
{
	ChargePair cp1(34,45, 4,5, cmp, 12.34, false);
	ChargePair cp2(34,15, 4,5, cmp, 12.34, false);
	TEST_EQUAL(cp1!=cp2, true);
	ChargePair cp3(34,15, 4,5, cmp, 12.34, true);
	ChargePair cp4(34,15, 4,5, cmp, 12.34, false);
	TEST_EQUAL(cp3!=cp4, true);
	ChargePair cp5(34,15, 4,5, cmp, 12.34, false);
	ChargePair cp6(34,15, 4,5, cmp, 12.34, false);
	TEST_EQUAL(cp5!=cp6, false);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



