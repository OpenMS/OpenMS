// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Clemens Groepl, Andreas Bertsch $
// $Authors: Clemens Groepl, Andreas Bertsch $
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

/////////////////////////////////////////////////////////////

START_TEST(IsotopeDistribution, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

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

START_SECTION(void trimRight(DoubleReal cutoff))
	IsotopeDistribution iso(EmpiricalFormula("C160").getIsotopeDistribution(10));
	TEST_NOT_EQUAL(iso.size(),3)
	iso.trimRight(0.2);
	TEST_EQUAL(iso.size(),3)
END_SECTION

START_SECTION(void trimLeft(DoubleReal cutoff))
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
