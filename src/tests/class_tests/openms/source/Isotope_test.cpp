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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/Isotope.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(Isotope, "$Id$")

/////////////////////////////////////////////////////////////

Isotope* e_ptr = nullptr;
Isotope* e_nullPointer = nullptr;
START_SECTION(Isotope())
	e_ptr = new Isotope;
	TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION(~Isotope())
	delete e_ptr;
END_SECTION

IsotopeDistribution dist;
string name("Name"), symbol("Symbol");
unsigned int atomic_number(15);
double mono_weight(0.123456789);

e_ptr = nullptr;
START_SECTION((Isotope(const std::string & name,
            const std::string & symbol,
            unsigned int atomic_number,
            unsigned int neutrons,
            double mono_weight,
            double abundance,
            double half_life,
            Isotope::DecayMode dm)))
	e_ptr = new Isotope(name, symbol, atomic_number, 10u, mono_weight, 0.6, 42, Isotope::DecayMode::UNKNOWN);	
	TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION( const Element* getElement() const)
  const Element* el = e_ptr->getElement();
  TEST_NOT_EQUAL(el, e_nullPointer)
  TEST_EQUAL(el->getSymbol(), "P")
END_SECTION

START_SECTION(Isotope(const Isotope& Isotope))
	Isotope copy(*e_ptr);
	TEST_EQUAL(*e_ptr == copy, true)
END_SECTION

delete e_ptr;
e_ptr = new Isotope;

START_SECTION(void setHalfLife(double hl))
	e_ptr->setHalfLife(8.6);
	NOT_TESTABLE
END_SECTION

START_SECTION(double getHalfLife() const)
	TEST_REAL_SIMILAR(e_ptr->getHalfLife(), 8.6)
END_SECTION

START_SECTION(void setAbundance(double hl))
	e_ptr->setAbundance(0.6);
	NOT_TESTABLE
END_SECTION

START_SECTION(double getAbundance() const)
	TEST_REAL_SIMILAR(e_ptr->getAbundance(), 0.6)
END_SECTION

START_SECTION(void setNeutrons(int hl))
	e_ptr->setNeutrons(10);
	NOT_TESTABLE
END_SECTION

START_SECTION(int getNeutrons() const)
	TEST_EQUAL(e_ptr->getNeutrons(), 10)
END_SECTION

START_SECTION(void setDecayMode(int hl))
	e_ptr->setDecayMode(Isotope::DecayMode::ALPHA);
	NOT_TESTABLE
END_SECTION

START_SECTION(int getDecayMode() const)
	TEST_EQUAL(e_ptr->getDecayMode(), Isotope::DecayMode::ALPHA)
END_SECTION

START_SECTION( virtual bool isIsotope() )
  TEST_EQUAL(e_ptr->isIsotope(), true)
END_SECTION

START_SECTION(bool isStable() const)
  TEST_EQUAL(e_ptr->isStable(), false)
END_SECTION

START_SECTION(Isotope& operator = (const Isotope& Isotope))
	Isotope e = *e_ptr;
	TEST_EQUAL(e == *e_ptr, true)
END_SECTION

START_SECTION(bool operator != (const Isotope& Isotope) const)
	Isotope e(*e_ptr);
	TEST_EQUAL(e != *e_ptr, false)
	e.setAverageWeight(0.54321);
	TEST_EQUAL(e != *e_ptr, true)
END_SECTION

START_SECTION(bool operator == (const Isotope& Isotope) const)
	Isotope e(*e_ptr);
	TEST_EQUAL(e == *e_ptr, true)
	e.setAverageWeight(0.54321);
	TEST_EQUAL(e == *e_ptr, false)
END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


