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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/Container.h>
#include <OpenMS/DATASTRUCTURES/String.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(Element, "$Id$")

/////////////////////////////////////////////////////////////

Element* e_ptr = nullptr;
Element* e_nullPointer = nullptr;
START_SECTION(Element())
	e_ptr = new Element;
	TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION(~Element())
	delete e_ptr;
END_SECTION

IsotopeDistribution dist;
String name("Name"), symbol("Symbol");
UInt atomic_number(43);
double average_weight(0.12345);
double mono_weight(0.123456789);

e_ptr = nullptr;
START_SECTION((Element(const String& name, const String& symbol, UInt atomic_number, double average_weight, double mono_weight, const IsotopeDistribution& isotopes)))
	e_ptr = new Element(name, symbol, atomic_number, average_weight, mono_weight, dist);	
	TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION(Element(const Element& element))
	Element copy(*e_ptr);
	TEST_EQUAL(*e_ptr == copy, true)
END_SECTION

delete e_ptr;
e_ptr = new Element;

START_SECTION(void setAtomicNumber(UInt atomic_number))
	e_ptr->setAtomicNumber(atomic_number);
	NOT_TESTABLE
END_SECTION

START_SECTION(UInt getAtomicNumber() const)
	TEST_EQUAL(e_ptr->getAtomicNumber(), atomic_number)
END_SECTION

START_SECTION(void setName(const String& name))
	e_ptr->setName(name);
	NOT_TESTABLE
END_SECTION

START_SECTION(const String& getName() const)
	TEST_EQUAL(e_ptr->getName(), name)
END_SECTION

START_SECTION(void setSymbol(const String& symbol))
	e_ptr->setSymbol(symbol);
	NOT_TESTABLE
END_SECTION

START_SECTION(const String& getSymbol() const)
	TEST_EQUAL(e_ptr->getSymbol(), symbol)
END_SECTION

START_SECTION(void setIsotopeDistribution(const IsotopeDistribution& isotopes))
	e_ptr->setIsotopeDistribution(dist);
	NOT_TESTABLE
END_SECTION

START_SECTION((const IsotopeDistribution& getIsotopeDistribution() const))
	TEST_EQUAL(e_ptr->getIsotopeDistribution() == dist, true)
END_SECTION

START_SECTION(void setAverageWeight(double weight))
	e_ptr->setAverageWeight(average_weight);
	NOT_TESTABLE
END_SECTION

START_SECTION(double getAverageWeight() const)
	TEST_REAL_SIMILAR(e_ptr->getAverageWeight(), average_weight)
END_SECTION

START_SECTION(void setMonoWeight(double weight))
	e_ptr->setMonoWeight(2.333);
	NOT_TESTABLE
END_SECTION

START_SECTION(double getMonoWeight() const)
	TEST_REAL_SIMILAR(e_ptr->getMonoWeight(), 2.333)
END_SECTION

START_SECTION(Element& operator = (const Element& element))
	Element e = *e_ptr;
	TEST_EQUAL(e == *e_ptr, true)
END_SECTION

START_SECTION(bool operator != (const Element& element) const)
	Element e(*e_ptr);
	TEST_EQUAL(e != *e_ptr, false)
	e.setAverageWeight(0.54321);
	TEST_EQUAL(e != *e_ptr, true)
END_SECTION

START_SECTION(bool operator == (const Element& element) const)
	Element e(*e_ptr);
	TEST_EQUAL(e == *e_ptr, true)
	e.setAverageWeight(0.54321);
	TEST_EQUAL(e == *e_ptr, false)
END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
