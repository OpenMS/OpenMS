// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>

///////////////////////////

// More headers

#include <iostream>
#include <iterator>
#include <utility>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

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

IsotopeDistribution* nullPointer = 0;

START_SECTION(IsotopeDistribution())
	IsotopeDistribution* ptr = 0;
	ptr = new IsotopeDistribution();
	Size max_isotope = ptr->getMaxIsotope();
  TEST_EQUAL(max_isotope, 0)
	TEST_NOT_EQUAL(ptr, nullPointer)
	delete ptr;
END_SECTION

START_SECTION(IsotopeDistribution(Size max_isotope))
	IsotopeDistribution* ptr = new IsotopeDistribution(117);
	Size max_isotope = ptr->getMaxIsotope();
  TEST_EQUAL(max_isotope, 117)
	TEST_NOT_EQUAL(ptr, nullPointer)
	delete ptr;
END_SECTION

IsotopeDistribution* iso = new IsotopeDistribution();

START_SECTION(IsotopeDistribution(const IsotopeDistribution& isotope_distribution))
	IsotopeDistribution copy;
	copy = *iso;
  for (Size i = 0; i != copy.getContainer().size(); ++i)
  {
    TEST_EQUAL(copy.getContainer()[i].first, iso->getContainer()[i].first)
    TEST_EQUAL(copy.getContainer()[i].second, iso->getContainer()[i].second)
  }
	TEST_EQUAL(copy.getMin(), iso->getMin())
	TEST_EQUAL(copy.getMax(), iso->getMax())
	TEST_EQUAL(copy.size(), iso->size())
	TEST_EQUAL(copy.getMaxIsotope(), iso->getMaxIsotope())
END_SECTION

START_SECTION(~IsotopeDistribution())
	IsotopeDistribution* ptr = new IsotopeDistribution(117);
	delete ptr;
END_SECTION

START_SECTION(IsotopeDistribution& operator = (const IsotopeDistribution& isotope_distribution))
	IsotopeDistribution copy;
	copy = *iso;
	for (Size i = 0; i != copy.getContainer().size(); ++i)
	{
		TEST_EQUAL(copy.getContainer()[i].first, iso->getContainer()[i].first)
		TEST_EQUAL(copy.getContainer()[i].second, iso->getContainer()[i].second)
	}
	TEST_EQUAL(copy.getMin(), iso->getMin())
	TEST_EQUAL(copy.getMax(), iso->getMax())
	TEST_EQUAL(copy.size(), iso->size())
	TEST_EQUAL(copy.getMaxIsotope(), iso->getMaxIsotope())
END_SECTION

START_SECTION(void setMaxIsotope(Size max_isotope))
	IsotopeDistribution iso2;
	iso2.estimateFromPeptideWeight(1234.2);
	TEST_EQUAL(iso->getMaxIsotope(), 0)
	TEST_EQUAL(iso2.getContainer().size(), 275)
	iso->setMaxIsotope(117);
	TEST_EQUAL(iso->getMaxIsotope(), 117)
END_SECTION

START_SECTION(Size getMaxIsotope() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(IsotopeDistribution operator + (const IsotopeDistribution& isotope_distribution) const)
	IsotopeDistribution iso1(1), iso2(1);
	IsotopeDistribution result = iso1 + iso2;
	TEST_EQUAL(result.size(), 1)
	IsotopeDistribution::ContainerType container = result.getContainer();
	TEST_EQUAL(container[0].first, 0)
	TEST_EQUAL(container[0].second, 1)
END_SECTION

START_SECTION(IsotopeDistribution& operator *= (Size factor))
	EmpiricalFormula ef("C222N190O110");
	IsotopeDistribution id = ef.getIsotopeDistribution(11);
	IsotopeDistribution::ContainerType container;
	container.push_back(make_pair<Size, double>(7084, 0.0349429));
	container.push_back(make_pair<Size, double>(7085, 0.109888));
	container.push_back(make_pair<Size, double>(7086, 0.180185));
	container.push_back(make_pair<Size, double>(7087, 0.204395));
	container.push_back(make_pair<Size, double>(7088, 0.179765));
	container.push_back(make_pair<Size, double>(7089, 0.130358));
	container.push_back(make_pair<Size, double>(7090, 0.0809864));
	container.push_back(make_pair<Size, double>(7091, 0.0442441));
	container.push_back(make_pair<Size, double>(7092, 0.0216593));
	container.push_back(make_pair<Size, double>(7093, 0.00963707));
	container.push_back(make_pair<Size, double>(7094, 0.0039406));

	for (Size i = 0; i != id.size(); ++i)
	{
		TEST_EQUAL(id.getContainer()[i].first, container[i].first)
		TEST_REAL_SIMILAR(id.getContainer()[i].second, container[i].second)
	}

  // test gapped isotope distributions, e.g. bromide 79,81 (missing 80)
  {
    EmpiricalFormula ef("Br2");
    IsotopeDistribution id = ef.getIsotopeDistribution(5);
    container.clear();
    // the expected results as pairs of
    // [nominal mass, probability]
    // derived via convolution of elemental probabilities; the sum of all probabilities is 1
    // For Br2, this is simply the product of Bromine x Bromine, which
    // has a light isotope (79 Da, ~50% probability) and a heavy isotope (81 Da, ~50% probability)
    container.push_back(make_pair<Size, double>(158, 0.2569476));  // 79+79, ~ 0.5 * 0.5
    container.push_back(make_pair<Size, double>(159, 0.0));        // this mass cannot be explained by two Br atoms
    container.push_back(make_pair<Size, double>(160, 0.49990478)); // 79+81 (or 81+79), ~ 0.5 * 0.5 + 0.5 * 0.5
    container.push_back(make_pair<Size, double>(161, 0.0));        // same as mass 159
    container.push_back(make_pair<Size, double>(162, 0.24314761)); // 81+81, ~ 0.5 * 0.5
    for (Size i = 0; i != id.size(); ++i)
    {
      TEST_EQUAL(id.getContainer()[i].first, container[i].first)
      TEST_REAL_SIMILAR(id.getContainer()[i].second, container[i].second)
    }
  }
  {
    // testing a formula which has more than one element (here: C and Br), since the internal computation is different
    // The convolution is similar to the one above, but add another convolution step with Carbon (hence the lightest mass is 12 Da heavier)
    EmpiricalFormula ef("CBr2");
    IsotopeDistribution id = ef.getIsotopeDistribution(7);
    container.clear();
    container.push_back(make_pair<Size, double>(170, 0.254198270573));
    container.push_back(make_pair<Size, double>(171, 0.002749339427));
    container.push_back(make_pair<Size, double>(172, 0.494555798854));
    container.push_back(make_pair<Size, double>(173, 0.005348981146));
    container.push_back(make_pair<Size, double>(174, 0.240545930573));
    container.push_back(make_pair<Size, double>(175, 0.002601679427));
    for (Size i = 0; i != id.size(); ++i)
    {
      TEST_EQUAL(id.getContainer()[i].first, container[i].first)
      TEST_REAL_SIMILAR(id.getContainer()[i].second, container[i].second)
    }
  }

END_SECTION

START_SECTION(bool operator==(const IsotopeDistribution &isotope_distribution) const)
	IsotopeDistribution iso1(1);
	IsotopeDistribution iso2(2);
	TEST_EQUAL(iso1 == iso2, false)
	iso2.setMaxIsotope(1);
	TEST_EQUAL(iso1 == iso2, true)
	IsotopeDistribution iso3(EmpiricalFormula("C4").getIsotopeDistribution(11)),
											iso4(EmpiricalFormula("C4").getIsotopeDistribution(11));
	TEST_EQUAL(iso3 == iso4, true)
END_SECTION

START_SECTION(void set(const ContainerType &distribution))
	IsotopeDistribution iso1(EmpiricalFormula("C4").getIsotopeDistribution(11)), iso2;
	TEST_EQUAL(iso1 == iso2, false)
	IsotopeDistribution::ContainerType container = iso1.getContainer();
	iso2.set(container);
	TEST_EQUAL(iso1.getContainer() == iso2.getContainer(), true)
	iso2.setMaxIsotope(iso1.getMaxIsotope());
	TEST_EQUAL(iso1 == iso2, true)
END_SECTION

START_SECTION(const ContainerType& getContainer() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(Size getMax() const)
	IsotopeDistribution iso(EmpiricalFormula("H2").getIsotopeDistribution(11));
	TEST_EQUAL(iso.getMax(), 6)
END_SECTION

START_SECTION(Size getMin() const)
	IsotopeDistribution iso(EmpiricalFormula("H2").getIsotopeDistribution(11));
	TEST_EQUAL(iso.getMin(), 2)
	IsotopeDistribution iso2(EmpiricalFormula("C4").getIsotopeDistribution(11));
	TEST_EQUAL(iso2.getMin(), 48)
END_SECTION

START_SECTION(Size size() const)
	IsotopeDistribution iso1, iso2(EmpiricalFormula("C4").getIsotopeDistribution(11));
	TEST_EQUAL(iso1.size(), 1)
	TEST_EQUAL(iso2.size(), 5)
END_SECTION

START_SECTION(void clear())
	IsotopeDistribution iso2(EmpiricalFormula("C4").getIsotopeDistribution(11));
	TEST_EQUAL(iso2.size(), 5)
	iso2.clear();
	TEST_EQUAL(iso2.size(), 0)
END_SECTION

START_SECTION(void estimateFromPeptideWeight(double average_weight))
	// hard to test as this is an rough estimate
	IsotopeDistribution iso(3);
	iso.estimateFromPeptideWeight(100.0);
	iso.renormalize();
	TEST_REAL_SIMILAR(iso.begin()->second, 0.95137)

	iso.estimateFromPeptideWeight(1000.0);
	TEST_REAL_SIMILAR(iso.begin()->second, 0.572779)

	iso.estimateFromPeptideWeight(10000.0);
	TEST_REAL_SIMILAR(iso.begin()->second, 0.00291426)
END_SECTION

START_SECTION(void trimRight(double cutoff))
	IsotopeDistribution iso(EmpiricalFormula("C160").getIsotopeDistribution(10));
	TEST_NOT_EQUAL(iso.size(),3)
	iso.trimRight(0.2);
	TEST_EQUAL(iso.size(),3)
END_SECTION

START_SECTION(void trimLeft(double cutoff))
	IsotopeDistribution iso(EmpiricalFormula("C160").getIsotopeDistribution(10));
	iso.trimRight(0.2);
	iso.trimLeft(0.2);
	TEST_EQUAL(iso.size(),2)
END_SECTION

START_SECTION(void renormalize())
	IsotopeDistribution iso(EmpiricalFormula("C160").getIsotopeDistribution(10));
	iso.trimRight(0.2);
	iso.trimLeft(0.2);
	iso.renormalize();
	double sum = 0;
	for (IsotopeDistribution::ConstIterator it = iso.begin(); it != iso.end(); ++it)
	{
		sum += it->second;
	}

	TEST_REAL_SIMILAR(sum, 1.0)
END_SECTION

START_SECTION(IsotopeDistribution& operator+=(const IsotopeDistribution &isotope_distribution))
	IsotopeDistribution iso1(EmpiricalFormula("H1").getIsotopeDistribution(11)),
											iso2(EmpiricalFormula("H2").getIsotopeDistribution(11));
	TEST_EQUAL(iso1 == iso2, false)
	iso1 += IsotopeDistribution(EmpiricalFormula("H1").getIsotopeDistribution(11));
	TEST_EQUAL(iso1.size() == iso2.size(), true)
	IsotopeDistribution::ConstIterator it1(iso1.begin()), it2(iso2.begin());

	for (; it1 != iso1.end(); ++it1, ++it2)
	{
		TEST_EQUAL(it1->first, it2->first)
		TEST_REAL_SIMILAR(it2->second, it2->second)
	}

END_SECTION

START_SECTION(IsotopeDistribution operator *(Size factor) const)
	IsotopeDistribution iso1(EmpiricalFormula("H1").getIsotopeDistribution(11)),
											iso2(EmpiricalFormula("H5").getIsotopeDistribution(11));
	TEST_EQUAL(iso1 == iso2, false)
	IsotopeDistribution iso3 = iso1 * 5;
	iso3.renormalize();
	iso2.renormalize();

	TEST_EQUAL(iso2.size(), iso3.size())
	IsotopeDistribution::ConstIterator it1(iso2.begin()), it2(iso3.begin());

	for (; it1 != iso2.end(); ++it1, ++it2)
	{
		TEST_EQUAL(it1->first, it2->first)
		TEST_REAL_SIMILAR(it1->second, it2->second)
	}
END_SECTION

START_SECTION(bool operator!=(const IsotopeDistribution &isotope_distribution) const)
  IsotopeDistribution iso1(1);
  IsotopeDistribution iso2(2);
  TEST_EQUAL(iso1 != iso2, true)
  iso2.setMaxIsotope(1);
  TEST_EQUAL(iso1 != iso2, false)
  IsotopeDistribution iso3(EmpiricalFormula("C4").getIsotopeDistribution(11)),
                      iso4(EmpiricalFormula("C4").getIsotopeDistribution(11));
  TEST_EQUAL(iso3 != iso4, false)
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

delete iso;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

#pragma clang diagnostic pop

