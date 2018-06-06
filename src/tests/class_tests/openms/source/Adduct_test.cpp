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
#include <OpenMS/DATASTRUCTURES/Adduct.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(Adduct, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Adduct* ptr = nullptr;
Adduct* nullPointer = nullptr;
START_SECTION(Adduct())
{
	ptr = new Adduct();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~Adduct())
{
	delete ptr;
}
END_SECTION

START_SECTION((Adduct(Int charge)))
{
	Adduct a(123);
	TEST_EQUAL(a.getCharge(), 123);
}
END_SECTION

START_SECTION((Adduct(Int charge, Int amount, double singleMass, String formula, double log_prob, double rt_shift, const String label="")))
{
	Adduct a(123, 43, 123.456f, "S", -0.3453, -10);
	TEST_EQUAL(a.getCharge(), 123);
	TEST_EQUAL(a.getAmount(), 43);
	TEST_REAL_SIMILAR(a.getSingleMass(), 123.456);
	TEST_EQUAL(a.getFormula()=="S1", true);
	TEST_REAL_SIMILAR(a.getLogProb(), -0.3453);
  TEST_REAL_SIMILAR(a.getRTShift(), -10);
  TEST_EQUAL(a.getLabel(), "");

	Adduct a2(123, 43, 123.456f, "S", -0.3453, -10, "testlabel");
  TEST_EQUAL(a2.getLabel(), "testlabel");	
}
END_SECTION

START_SECTION([EXTRA] friend OPENMS_DLLAPI bool operator==(const Adduct& a, const Adduct& b))
{
	
	Adduct a(123,  3, 123.456f, "S", -0.3453f, 0);
	Adduct b(a);

	TEST_EQUAL(a==b, true);
	a.setAmount(22);
	TEST_EQUAL(a==b, false);
	
}
END_SECTION

START_SECTION((const Int& getCharge() const))
{
	NOT_TESTABLE //well.. tested below...
}
END_SECTION

START_SECTION((void setCharge(const Int &charge)))
{
	Adduct a;
	a.setCharge(123);
	TEST_EQUAL(a.getCharge(), 123);
}
END_SECTION

START_SECTION((const Int& getAmount() const))
{
	NOT_TESTABLE //well.. tested below...
}
END_SECTION

START_SECTION((void setAmount(const Int &amount)))
{
	Adduct a;
	a.setAmount(43);
  TEST_EQUAL(a.getAmount(), 43);
}
END_SECTION

START_SECTION((const double& getSingleMass() const))
{
	NOT_TESTABLE //well.. tested below...
}
END_SECTION

START_SECTION((void setSingleMass(const double &singleMass)))
{
	Adduct a;
	a.setSingleMass(43.21);
  TEST_REAL_SIMILAR(a.getSingleMass(), 43.21);
}
END_SECTION

START_SECTION((const double& getLogProb() const))
{
	NOT_TESTABLE //well.. tested below...
}
END_SECTION

START_SECTION((void setLogProb(const double &log_prob)))
{
	Adduct a;
	a.setLogProb(43.21f);
  TEST_REAL_SIMILAR(a.getLogProb(), 43.21);
}
END_SECTION

START_SECTION((const String& getFormula() const))
{
	NOT_TESTABLE //well.. tested below...
}
END_SECTION

START_SECTION((void setFormula(const String &formula)))
	Adduct a;
	a.setFormula("S");
  TEST_EQUAL(a.getFormula()=="S1", true);
END_SECTION

START_SECTION((const double& getRTShift() const))
	Adduct a(123, 43, 123.456f, "S", -0.3453, -10);
  TEST_REAL_SIMILAR(a.getRTShift(), -10);
	Adduct a1(123, 43, 123.456f, "S", -0.3453, 11);
  TEST_REAL_SIMILAR(a1.getRTShift(), 11);
END_SECTION

START_SECTION((const String& getLabel() const ))
	Adduct a(123, 43, 123.456f, "S", -0.3453, -10);
  TEST_EQUAL(a.getLabel(), "");
	Adduct a1(123, 43, 123.456f, "S", -0.3453, 11, "mylabel");
  TEST_EQUAL(a1.getLabel(), "mylabel");
END_SECTION



START_SECTION((Adduct operator *(const Int m) const))
{
	Adduct a_p(123, 43, 123.456, "S", -0.3453, 0);
	Adduct a = a_p*4;
	TEST_EQUAL(a.getCharge(), 123);
	TEST_EQUAL(a.getAmount(), 43*4);
	TEST_REAL_SIMILAR(a.getSingleMass(), 123.456f);
	TEST_EQUAL(a.getFormula()=="S1", true);
	TEST_REAL_SIMILAR(a.getLogProb(), -0.3453);
}
END_SECTION

START_SECTION((Adduct operator+(const Adduct &rhs)))
{
	Adduct a_p(123, 43, 123.456f, "S", -0.3453f, 0);
	Adduct a_p2(123, 40, 123.456f, "S", -0.3453f, 0);
	Adduct a = a_p + a_p2;
	TEST_EQUAL(a.getCharge(), 123);
	TEST_EQUAL(a.getAmount(), 43+40);
	TEST_REAL_SIMILAR(a.getSingleMass(), 123.456);
	TEST_EQUAL(a.getFormula()=="S1", true);
	TEST_REAL_SIMILAR(a.getLogProb(), -0.3453);
}
END_SECTION

START_SECTION((void operator+=(const Adduct &rhs)))
{
	Adduct a_p(123, 43, 123.456f, "S", -0.3453f, 0);
	Adduct a(a_p);
	a.setAmount(10);
	a	+= a_p;
	TEST_EQUAL(a.getAmount(), 43+10);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



