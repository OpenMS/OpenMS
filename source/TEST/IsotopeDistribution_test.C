// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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

CHECK(IsotopeDistribution())
	IsotopeDistribution* ptr = 0;
	ptr = new IsotopeDistribution();
	Size max_isotope = ptr->getMaxIsotope();
  TEST_EQUAL(max_isotope, 0)
	TEST_NOT_EQUAL(ptr, 0)
	delete ptr;
RESULT

CHECK(IsotopeDistribution(Size max_isotope))
	IsotopeDistribution* ptr = new IsotopeDistribution(117);
	Size max_isotope = ptr->getMaxIsotope();
  TEST_EQUAL(max_isotope, 117)
	TEST_NOT_EQUAL(ptr, 0)
	delete ptr;
RESULT

IsotopeDistribution* iso = new IsotopeDistribution();

CHECK(IsotopeDistribution(const IsotopeDistribution& isotope_distribution))
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
RESULT

CHECK(~IsotopeDistribution())
	IsotopeDistribution* ptr = new IsotopeDistribution(117);
	delete ptr;
RESULT

CHECK(IsotopeDistribution& operator = (const IsotopeDistribution& isotope_distribution))
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
RESULT

CHECK(void setMaxIsotope(Size max_isotope))
	iso->setMaxIsotope(117);
	TEST_EQUAL(iso->getMaxIsotope(), 117)
RESULT

CHECK(IsotopeDistribution operator + (const IsotopeDistribution& isotope_distribution) const)
	IsotopeDistribution iso1(1), iso2(1);
	IsotopeDistribution result = iso1 + iso2;
	TEST_EQUAL(result.size(), 1)
	IsotopeDistribution::ContainerType container = result.getContainer();
	TEST_EQUAL(container[0].first, 0)
	TEST_EQUAL(container[0].second, 1)
RESULT

CHECK(IsotopeDistribution& operator *= (UnsignedInt factor))
	EmpiricalFormula ef("C222N190O110");
	IsotopeDistribution id = ef.getIsotopeDistribution(11);
	IsotopeDistribution::ContainerType container;
	container.push_back(make_pair<UnsignedInt, double>(7084, 0.0349429));
	container.push_back(make_pair<UnsignedInt, double>(7085, 0.109888));
	container.push_back(make_pair<UnsignedInt, double>(7086, 0.180185));
	container.push_back(make_pair<UnsignedInt, double>(7087, 0.204395));
	container.push_back(make_pair<UnsignedInt, double>(7088, 0.179765));
	container.push_back(make_pair<UnsignedInt, double>(7089, 0.130358));
	container.push_back(make_pair<UnsignedInt, double>(7090, 0.0809864));
	container.push_back(make_pair<UnsignedInt, double>(7091, 0.0442441));
	container.push_back(make_pair<UnsignedInt, double>(7092, 0.0216593));
	container.push_back(make_pair<UnsignedInt, double>(7093, 0.00963707));
	container.push_back(make_pair<UnsignedInt, double>(7094, 0.0039406));

	for (Size i = 0; i != id.size(); ++i)
	{
		TEST_EQUAL(id.getContainer()[i].first, container[i].first)
		TEST_REAL_EQUAL(id.getContainer()[i].second, container[i].second)
	}
	
RESULT


CHECK(ConstIterator begin() const)

RESULT

CHECK(ConstIterator end() const)

RESULT

CHECK(IsotopeDistribution operator * (UnsignedInt factor) const)

RESULT

CHECK(IsotopeDistribution& operator += (const IsotopeDistribution& isotope_distribution))

RESULT

CHECK(Iterator begin())

RESULT

CHECK(Iterator end())

RESULT

CHECK(Size getMaxIsotope() const)

RESULT

CHECK(Size size() const)

RESULT

CHECK(UnsignedInt getMax() const)

RESULT

CHECK(UnsignedInt getMin() const)

RESULT

CHECK(bool operator != (const IsotopeDistribution& isotope_distribution) const)

RESULT

CHECK(bool operator == (const IsotopeDistribution& isotope_distribution) const)

RESULT

CHECK(const ContainerType& getContainer() const)

RESULT

CHECK(double getTrimRightCutoff())

RESULT

CHECK(void clear())

RESULT

CHECK(void estimateFromPeptideWeight(double weight))

RESULT

CHECK(void renormalize())

RESULT

CHECK(void set(const ContainerType& distribution))
	// TODO
RESULT

CHECK(void setTrimRightCutoff(double const cutoff))

RESULT

CHECK(void trimRight())

RESULT

delete iso;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
