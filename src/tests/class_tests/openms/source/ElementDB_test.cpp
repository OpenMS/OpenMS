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

#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/DATASTRUCTURES/Map.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(ElementDB, "$Id$")

/////////////////////////////////////////////////////////////

const ElementDB* e_ptr = nullptr;
const ElementDB* e_nullPointer = nullptr;
const Element * elem_nullPointer = nullptr;

START_SECTION(static const ElementDB* getInstance())
	e_ptr = ElementDB::getInstance();
	TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION((const Map<String, const Element*>& getNames() const))
	Map<String, const Element*> names = e_ptr->getNames();
	const Element * e = e_ptr->getElement("Carbon");
	TEST_EQUAL(e, names["Carbon"])
  TEST_NOT_EQUAL(e, elem_nullPointer)
END_SECTION


START_SECTION((const Map<String, const Element*>& getSymbols() const))
	Map<String, const Element*> symbols = e_ptr->getSymbols();
	const Element * e = e_ptr->getElement("Carbon");
	TEST_EQUAL(e, symbols["C"])
  TEST_NOT_EQUAL(e, elem_nullPointer)
END_SECTION

START_SECTION((const Map<UInt, const Element*>& getAtomicNumbers() const))
	Map<UInt, const Element*> atomic_numbers = e_ptr->getAtomicNumbers();
	const Element * e = e_ptr->getElement("Carbon");
	TEST_EQUAL(e, atomic_numbers[6])
  TEST_NOT_EQUAL(e, elem_nullPointer)
END_SECTION

START_SECTION(const Element* getElement(const String& name) const)
	const Element * e1 = e_ptr->getElement("Hydrogen");
	const Element * e2 = e_ptr->getElement("H");
	TEST_EQUAL(e1, e2);
  TEST_NOT_EQUAL(e1, elem_nullPointer);
END_SECTION

START_SECTION(const Element* getElement(UInt atomic_number) const)
	const Element * e1 = e_ptr->getElement("Carbon");
	const Element * e2 = e_ptr->getElement(6);
	TEST_EQUAL(e1, e2)
  TEST_NOT_EQUAL(e1, elem_nullPointer)
END_SECTION

START_SECTION(bool hasElement(const String& name) const)
	TEST_EQUAL(e_ptr->hasElement("Carbon"), true)
END_SECTION

START_SECTION(bool hasElement(UInt atomic_number) const)
	TEST_EQUAL(e_ptr->hasElement(6), true)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
