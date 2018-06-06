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
// $Authors: Clemens Groepl, Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------
//

///////////////////////////

// This one is going to be tested.
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
///////////////////////////

// More headers

#include <iostream>
#include <iterator>
#include <utility>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

/////////////////////////////////////////////////////////////

START_TEST(IsotopeDistribution, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wshadow"

using namespace OpenMS;
using namespace std;

IsotopeDistribution* nullPointer = nullptr;

START_SECTION(CoarseIsotopePatternGenerator())
	IsotopeDistribution* ptr = nullptr;
	ptr = new IsotopeDistribution();
	Size container_size = ptr->size();
  TEST_EQUAL(container_size, 1)
	TEST_NOT_EQUAL(ptr, nullPointer)
	delete ptr;
END_SECTION



IsotopeDistribution* iso = new IsotopeDistribution();

START_SECTION(IsotopeDistribution(const IsotopeDistribution& isotope_distribution))
	IsotopeDistribution copy;
	copy = *iso;
  for (Size i = 0; i != copy.getContainer().size(); ++i)
  {
    TEST_EQUAL(copy.getContainer()[i].getMZ(), iso->getContainer()[i].getMZ())
    TEST_EQUAL(copy.getContainer()[i].getIntensity(), iso->getContainer()[i].getIntensity())
  }
	TEST_EQUAL(copy.getMin(), iso->getMin())
	TEST_EQUAL(copy.getMax(), iso->getMax())
	TEST_EQUAL(copy.size(), iso->size())
END_SECTION

START_SECTION(~IsotopeDistribution())
	IsotopeDistribution* ptr = new IsotopeDistribution();
	delete ptr;
END_SECTION

START_SECTION(IsotopeDistribution& operator = (const CoarseIsotopePatternGenerator& isotope_distribution))
	IsotopeDistribution copy;
	copy = *iso;
	for (Size i = 0; i != copy.getContainer().size(); ++i)
	{
		TEST_EQUAL(copy.getContainer()[i].getMZ(), iso->getContainer()[i].getMZ())
		TEST_EQUAL(copy.getContainer()[i].getIntensity(), iso->getContainer()[i].getIntensity())
	}
	TEST_EQUAL(copy.getMin(), iso->getMin())
	TEST_EQUAL(copy.getMax(), iso->getMax())
	TEST_EQUAL(copy.size(), iso->size())
END_SECTION


START_SECTION(bool operator==(const IsotopeDistribution &isotope_distribution) const)
	IsotopeDistribution iso1, iso2;
	TEST_EQUAL(iso1 == iso2, true)
	IsotopeDistribution iso3(EmpiricalFormula("C4").getIsotopeDistribution(CoarseIsotopePatternGenerator(11))),
    iso4(EmpiricalFormula("C4").getIsotopeDistribution(CoarseIsotopePatternGenerator(11)));
	TEST_EQUAL(iso3 == iso4, true)
END_SECTION

START_SECTION(void set(const ContainerType &distribution))
	IsotopeDistribution iso1(EmpiricalFormula("C4").getIsotopeDistribution(CoarseIsotopePatternGenerator(11))), iso2;
	TEST_EQUAL(iso1 == iso2, false)
	IsotopeDistribution::ContainerType container = iso1.getContainer();
	iso2.set(container);
	TEST_EQUAL(iso1.getContainer() == iso2.getContainer(), true)
	TEST_EQUAL(iso1 == iso2, true)
END_SECTION

START_SECTION(const ContainerType& getContainer() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(Size getMax() const)
	IsotopeDistribution iso(EmpiricalFormula("H2").getIsotopeDistribution(CoarseIsotopePatternGenerator(11)));
	TEST_EQUAL(iso.getMax(), 6)
END_SECTION

START_SECTION(Size getMin() const)
	IsotopeDistribution iso(EmpiricalFormula("H2").getIsotopeDistribution(CoarseIsotopePatternGenerator(11)));
	TEST_EQUAL(iso.getMin(), 2)
	IsotopeDistribution iso2(EmpiricalFormula("C4").getIsotopeDistribution(CoarseIsotopePatternGenerator(11)));
	TEST_EQUAL(iso2.getMin(), 48)
END_SECTION

START_SECTION(Size size() const)
	IsotopeDistribution iso1, iso2(EmpiricalFormula("C4").getIsotopeDistribution(CoarseIsotopePatternGenerator(11)));
	TEST_EQUAL(iso1.size(), 1)
	TEST_EQUAL(iso2.size(), 5)
END_SECTION

START_SECTION(void clear())
	IsotopeDistribution iso2(EmpiricalFormula("C4").getIsotopeDistribution(CoarseIsotopePatternGenerator(11)));
	TEST_EQUAL(iso2.size(), 5)
	iso2.clear();
	TEST_EQUAL(iso2.size(), 0)
END_SECTION


START_SECTION(void trimRight(double cutoff))
	IsotopeDistribution iso(EmpiricalFormula("C160").getIsotopeDistribution(CoarseIsotopePatternGenerator(10)));
	TEST_NOT_EQUAL(iso.size(),3)
	iso.trimRight(0.2);
	TEST_EQUAL(iso.size(),3)
END_SECTION

START_SECTION(void trimLeft(double cutoff))
  IsotopeDistribution iso(EmpiricalFormula("C160").getIsotopeDistribution(CoarseIsotopePatternGenerator(10)));
	iso.trimRight(0.2);
	iso.trimLeft(0.2);
	TEST_EQUAL(iso.size(),2)
END_SECTION

START_SECTION(void renormalize())
	IsotopeDistribution iso(EmpiricalFormula("C160").getIsotopeDistribution(CoarseIsotopePatternGenerator(10)));
	iso.trimRight(0.2);
	iso.trimLeft(0.2);
	iso.renormalize();
	double sum = 0;
	for (IsotopeDistribution::ConstIterator it = iso.begin(); it != iso.end(); ++it)
	{
		sum += it->getIntensity();
	}

	TEST_REAL_SIMILAR(sum, 1.0)
END_SECTION



START_SECTION(bool operator!=(const IsotopeDistribution &isotope_distribution) const)
  IsotopeDistribution iso1;
  IsotopeDistribution iso2;
  TEST_EQUAL(iso1 != iso2, false)
  IsotopeDistribution iso3(EmpiricalFormula("C4").getIsotopeDistribution(CoarseIsotopePatternGenerator(11))),
                      iso4(EmpiricalFormula("C4").getIsotopeDistribution(CoarseIsotopePatternGenerator(11)));
  TEST_EQUAL(iso3 != iso4, false)
  TEST_EQUAL(iso2 != iso3, true)
END_SECTION

START_SECTION(Iterator begin())
	NOT_TESTABLE
END_SECTION

START_SECTION(Iterator end())
	NOT_TESTABLE
END_SECTION

START_SECTION(ConstIterator begin() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(ConstIterator end() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(ReverseIterator rbegin())
	NOT_TESTABLE
END_SECTION

START_SECTION(ReverseIterator rend())
	NOT_TESTABLE
END_SECTION

START_SECTION(ConstReverseIterator rbegin() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(ConstReverseIterator rend() const)
	NOT_TESTABLE
END_SECTION

delete iso;


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

#pragma clang diagnostic pop

