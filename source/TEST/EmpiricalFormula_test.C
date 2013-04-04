// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(ElementDB, "$Id$")

/////////////////////////////////////////////////////////////

EmpiricalFormula* e_ptr = 0;
EmpiricalFormula* e_nullPointer = 0;
START_SECTION(EmpiricalFormula())
	e_ptr = new EmpiricalFormula;
	TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION(~EmpiricalFormula())
	delete e_ptr;
END_SECTION

START_SECTION(EmpiricalFormula(const String& rhs))
	e_ptr = new EmpiricalFormula("C4");
	TEST_NOT_EQUAL(e_ptr, e_nullPointer)
        EmpiricalFormula e0("C5(13)C4H2");
        EmpiricalFormula e1("C5(13)C4");
        EmpiricalFormula e2("(12)C5(13)C4");
        EmpiricalFormula e3("C9");
	TEST_REAL_SIMILAR(e1.getMonoWeight(), e2.getMonoWeight())
	TEST_REAL_SIMILAR(e1.getMonoWeight(), 112.013419)
	TEST_REAL_SIMILAR(e2.getMonoWeight(), 112.013419)
END_SECTION

START_SECTION(EmpiricalFormula(const EmpiricalFormula& rhs))
	EmpiricalFormula ef(*e_ptr);
	TEST_EQUAL(ef == *e_ptr, true)
END_SECTION

START_SECTION((EmpiricalFormula(SignedSize number, const Element* element, SignedSize charge=0)))
	EmpiricalFormula ef(4, e_ptr->getElement("C"));
	TEST_EQUAL(ef == *e_ptr, true)
	TEST_EQUAL(ef.getCharge(), 0)
END_SECTION

START_SECTION(const Element* getElement(UInt atomic_number) const)
	const Element* e = e_ptr->getElement(6);
	TEST_EQUAL(e->getSymbol(), "C")
END_SECTION

START_SECTION(const Element* getElement(const String& name) const)
	const Element* e = e_ptr->getElement("C");
	TEST_EQUAL(e->getSymbol(), "C")
END_SECTION

START_SECTION(Size getNumberOf(UInt atomic_number) const)
	Size num1 = e_ptr->getNumberOf(6);
	TEST_EQUAL(num1, 4);
END_SECTION

START_SECTION(Size getNumberOf(const String& name) const)
	Size num2 = e_ptr->getNumberOf("C");
	TEST_EQUAL(num2, 4);
END_SECTION

START_SECTION(Size getNumberOf(const Element* element) const)
	const Element* e = e_ptr->getElement(6);
	Size num3 = e_ptr->getNumberOf(e);
	TEST_EQUAL(num3, 4);
END_SECTION

START_SECTION(Size getNumberOfAtoms() const)
	Size num4 = e_ptr->getNumberOfAtoms();
	TEST_EQUAL(num4, 4);
END_SECTION

START_SECTION(EmpiricalFormula& operator = (const EmpiricalFormula& rhs))
	EmpiricalFormula ef;
	ef = *e_ptr;
	TEST_EQUAL(*e_ptr == ef, true)
END_SECTION

START_SECTION(EmpiricalFormula& operator = (const String& rhs))
	EmpiricalFormula ef;
	ef = "C4";
	TEST_EQUAL(*e_ptr == ef, true)
	TEST_EXCEPTION(Exception::ParseError, ef = "2C4")
END_SECTION


START_SECTION(EmpiricalFormula operator * (const SignedSize& times) const)
	EmpiricalFormula ef("C3H8");
	ef = ef * 3;
	TEST_EQUAL(ef, "C9H24")
END_SECTION


START_SECTION(EmpiricalFormula& operator += (const EmpiricalFormula& rhs))
	EmpiricalFormula ef("C3");
	ef += ef;
	TEST_EQUAL(ef, "C6")
	EmpiricalFormula ef2("C-6H2");
	ef += ef2;
	TEST_EQUAL(ef, "H2");
END_SECTION

START_SECTION(EmpiricalFormula& operator += (const String& rhs))
	EmpiricalFormula ef;
	ef += "C";
	TEST_EQUAL(ef, "C")
	ef += "C5";
	TEST_EQUAL(ef, "C6")
	ef += "C-5";
	TEST_EQUAL(ef, "C")
	ef += "C-1H2";
	TEST_EQUAL(ef, "H2")
END_SECTION

START_SECTION(EmpiricalFormula operator + (const EmpiricalFormula& rhs) const)
	EmpiricalFormula ef("C2");
	EmpiricalFormula ef2;
	ef2 = ef + ef;
	TEST_EQUAL(ef2, "C4")
	ef2 = ef2 + EmpiricalFormula("C-4H2");
	TEST_EQUAL(ef2, "H2")
END_SECTION

START_SECTION(EmpiricalFormula operator + (const String& rhs) const)
	EmpiricalFormula ef1("C2");
	EmpiricalFormula ef2;
	ef2 = ef1 + "C2";
	TEST_EQUAL(ef2, "C4")
	ef2 = ef2 + "C-4H2";
	TEST_EQUAL(ef2, "H2")
END_SECTION

START_SECTION(EmpiricalFormula& operator -= (const EmpiricalFormula& rhs))
	EmpiricalFormula ef1("C5H12"), ef2("CH12");
	ef1 -= ef2;
	TEST_EQUAL(*e_ptr == ef1, true)
	ef1 -= EmpiricalFormula("C4H-2");
	TEST_EQUAL(ef1, "H2");
END_SECTION

START_SECTION(EmpiricalFormula& operator -= (const String& rhs))
	EmpiricalFormula ef1("C5H12");
	ef1 -= "CH12";
	TEST_EQUAL(*e_ptr == ef1, true)
	ef1 -= "C4H-2";
	TEST_EQUAL(ef1, "H2");
END_SECTION

START_SECTION(EmpiricalFormula operator - (const EmpiricalFormula& rhs) const)
	EmpiricalFormula ef1("C5H12"), ef2("CH12");
	EmpiricalFormula ef3, ef4;
	ef3 = ef1 - ef2;
	cerr << *e_ptr << " " << ef3 << endl;
	TEST_EQUAL(*e_ptr == ef3, true)
	ef3 = ef3 - EmpiricalFormula("C4H-2");
	TEST_EQUAL(ef3, "H2");
END_SECTION

START_SECTION(EmpiricalFormula operator - (const String& rhs) const)
	EmpiricalFormula ef1("C5H12"), ef2("CH12"), ef4;
	ef4 = ef1 - "CH12";
	TEST_EQUAL(*e_ptr == ef4, true)
	TEST_EXCEPTION(Exception::ParseError, ef1-"BLUBB")
	ef4 = ef4 - "C4H-2";
	TEST_EQUAL(ef4, "H2");
END_SECTION

START_SECTION(bool isEmpty() const)
	EmpiricalFormula ef;
	TEST_EQUAL(ef.isEmpty(), true)
	TEST_EQUAL(e_ptr->isEmpty(), false)
END_SECTION

START_SECTION(bool hasElement(const String& name) const)
	TEST_EQUAL(e_ptr->hasElement("C"), true)
	TEST_EQUAL(e_ptr->hasElement("N"), false)
END_SECTION

START_SECTION(bool hasElement(UInt atomic_number) const)
	TEST_EQUAL(e_ptr->hasElement(6), true)
	TEST_EQUAL(e_ptr->hasElement(7), false)
END_SECTION

START_SECTION(bool hasElement(const Element* element) const)
	const Element* e = e_ptr->getElement(6);
	TEST_EQUAL(e_ptr->hasElement(e), true)
	e = e_ptr->getElement(1);
	TEST_EQUAL(e_ptr->hasElement(e), false)
END_SECTION

START_SECTION(void setCharge(SignedSize charge))
	e_ptr->setCharge(1);
	NOT_TESTABLE // will be tested in next check
END_SECTION

START_SECTION(SignedSize getCharge() const)
	TEST_EQUAL(e_ptr->getCharge(), 1)
	EmpiricalFormula ef1("C2+");
	TEST_EQUAL(ef1.getCharge(), 1)
	EmpiricalFormula ef2("C2+3");
	TEST_EQUAL(ef2.getCharge(), 3)
END_SECTION

START_SECTION(bool isCharged() const)
	TEST_EQUAL(e_ptr->isCharged(), true)
	e_ptr->setCharge(0);
	TEST_EQUAL(e_ptr->isCharged(), false)
END_SECTION

START_SECTION(DoubleReal getAverageWeight() const)
	EmpiricalFormula ef("C2");
	const Element* e = e_ptr->getElement("C");
	TEST_REAL_SIMILAR(ef.getAverageWeight(), e->getAverageWeight() * 2)
END_SECTION

START_SECTION(DoubleReal getMonoWeight() const)
	EmpiricalFormula ef("C2");
	const Element* e = e_ptr->getElement("C");
	TEST_REAL_SIMILAR(ef.getMonoWeight(), e->getMonoWeight() * 2)
END_SECTION

START_SECTION(String getString() const)
	EmpiricalFormula ef("C2H5");
	String str = ef.getString();
	TEST_EQUAL(String(str).hasSubstring("H5"), true)
	TEST_EQUAL(String(str).hasSubstring("C2"), true)
END_SECTION

START_SECTION(const ElementDB* getElementDB() const)
	const ElementDB* db = 0;
  const ElementDB* db_nullPointer = 0;
	db = e_ptr->getElementDB();
  TEST_NOT_EQUAL(db, db_nullPointer);
	TEST_EQUAL(db->getElement("C")->getSymbol(), "C")
END_SECTION

START_SECTION([EXTRA](friend std::ostream& operator << (std::ostream&, const EmpiricalFormula&)))
	stringstream ss;
	EmpiricalFormula ef("C2H5");
	ss << ef;
	TEST_EQUAL(String(ss.str()).hasSubstring("H5"), true);
	TEST_EQUAL(String(ss.str()).hasSubstring("C2"), true);
END_SECTION

START_SECTION(bool operator != (const EmpiricalFormula& rhs) const)
	EmpiricalFormula ef1("C2H5"), ef2(*e_ptr);
	TEST_EQUAL(ef1 != ef2, true)
	TEST_EQUAL(ef1 != ef1, false)
	ef2.setCharge(1);
	TEST_EQUAL(ef2 != *e_ptr, true)
END_SECTION

START_SECTION(bool operator != (const String& rhs) const)
	EmpiricalFormula ef1("C2H5");
	TEST_EQUAL(ef1 != "C2", true)
	TEST_EQUAL(ef1 != "C2H5", false)
END_SECTION

START_SECTION(bool operator == (const EmpiricalFormula& rhs) const)
	EmpiricalFormula ef1("C2H5"), ef2(*e_ptr);
	TEST_EQUAL(ef1 == ef2, false)
	TEST_EQUAL(ef1 == ef1, true)
	ef2.setCharge(1);
	TEST_EQUAL(ef2 == *e_ptr, false)
END_SECTION

START_SECTION(bool operator == (const String& rhs) const)
	EmpiricalFormula ef1("C2H5");
	TEST_EQUAL(ef1 == "C2", false)
	TEST_EQUAL(ef1 == "C2H5", true)
END_SECTION

START_SECTION(ConstIterator begin() const)
	EmpiricalFormula ef("C6H12O6");
	Map<String, SignedSize> formula;
	formula["C"] = 6;
	formula["H"] = 12;
	formula["O"] = 6;
	for (EmpiricalFormula::ConstIterator it = ef.begin(); it != ef.end(); ++it)
	{
		TEST_EQUAL(it->second, formula[it->first->getSymbol()])
	}

END_SECTION

START_SECTION(ConstIterator end() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(IsotopeDistribution getIsotopeDistribution(UInt max_depth) const)
	EmpiricalFormula ef("C");
	IsotopeDistribution iso = ef.getIsotopeDistribution(20);
	double result[] = { 0.9893, 0.0107};
	Size i = 0;
	for (IsotopeDistribution::ConstIterator it = iso.begin(); it != iso.end(); ++it, ++i)
	{
		TEST_REAL_SIMILAR(it->second, result[i])
	}
END_SECTION

START_SECTION(([EXTRA] Check correct charge semantics))
	EmpiricalFormula ef1("H4C+"); // CH4 +1 charge
	TEST_EQUAL(ef1.getNumberOf("H"), 4)
	TEST_EQUAL(ef1.getNumberOf("C"), 1)
	TEST_EQUAL(ef1.getCharge(), 1)
	EmpiricalFormula ef2("H4C1+"); // ""
  TEST_EQUAL(ef2.getNumberOf("H"), 4)
  TEST_EQUAL(ef2.getNumberOf("C"), 1)
  TEST_EQUAL(ef2.getCharge(), 1)
	EmpiricalFormula ef3("H4C-1+"); // C-1 H4 +1 charge
  TEST_EQUAL(ef3.getNumberOf("H"), 4)
  TEST_EQUAL(ef3.getNumberOf("C"), -1)
  TEST_EQUAL(ef3.getCharge(), 1)
	EmpiricalFormula ef4("H4C-1"); // C-1 H4 0 charge
  TEST_EQUAL(ef4.getNumberOf("H"), 4)
  TEST_EQUAL(ef4.getNumberOf("C"), -1)
  TEST_EQUAL(ef4.getCharge(), 0)
	EmpiricalFormula ef5("H4C1-1"); // C1 H4 -1 charge
  TEST_EQUAL(ef5.getNumberOf("H"), 4)
  TEST_EQUAL(ef5.getNumberOf("C"), 1)
  TEST_EQUAL(ef5.getCharge(), -1)
	EmpiricalFormula ef6("H4C-1-1"); // C-1 H4 -1 charge
  TEST_EQUAL(ef6.getNumberOf("H"), 4)
  TEST_EQUAL(ef6.getNumberOf("C"), -1)
  TEST_EQUAL(ef6.getCharge(), -1)
	EmpiricalFormula ef7("H4C-1-"); // C-1 H4 -1 charge
  TEST_EQUAL(ef7.getNumberOf("H"), 4)
  TEST_EQUAL(ef7.getNumberOf("C"), -1)
  TEST_EQUAL(ef7.getCharge(), -1)
	EmpiricalFormula ef8("-"); // -1 Charge
  TEST_EQUAL(ef8.getNumberOf("H"), 0)
  TEST_EQUAL(ef8.getNumberOf("C"), 0)
  TEST_EQUAL(ef8.getCharge(), -1)
	EmpiricalFormula ef9("+"); // +1 Charge
  TEST_EQUAL(ef9.getNumberOf("H"), 0)
  TEST_EQUAL(ef9.getNumberOf("C"), 0)
  TEST_EQUAL(ef9.getCharge(), 1)
	EmpiricalFormula ef10("-3"); // -3 Charge
  TEST_EQUAL(ef10.getNumberOf("H"), 0)
  TEST_EQUAL(ef10.getNumberOf("C"), 0)
  TEST_EQUAL(ef10.getCharge(), -3)
	EmpiricalFormula ef11("+3"); // +3 Charge
  TEST_EQUAL(ef11.getNumberOf("H"), 0)
  TEST_EQUAL(ef11.getNumberOf("C"), 0)
  TEST_EQUAL(ef11.getCharge(), 3)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

