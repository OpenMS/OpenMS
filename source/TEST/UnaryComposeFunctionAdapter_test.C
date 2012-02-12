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
// $Maintainer: Stephan Aiche$
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/CONCEPT/UnaryComposeFunctionAdapter.h>
///////////////////////////

#include <string>
#include <vector>
#include <utility>
#include <algorithm>

using namespace OpenMS;
using namespace std;

// test class
class Element {
 public:
  Element(const string& a) : _a(a) {}
  const string& getA() const { return _a; }
 private:
  string _a;
};

START_TEST(UnaryComposeFunctionAdapter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

typedef equal_to<string> T_EQUAL_TO;
typedef binder2nd<T_EQUAL_TO> T_BINDER_2ND;
typedef std::const_mem_fun_ref_t<const std::string&, Element> T_MEMBER_FUNCTIOM;

typedef UnaryComposeFunctionAdapter<T_BINDER_2ND, T_MEMBER_FUNCTIOM> UCFA;

UCFA* ptr = 0;
UCFA* null_ptr = 0;
START_SECTION((UnaryComposeFunctionAdapter(const OP1 &o1, const OP2 &o2)))
{
  ptr = new UCFA(bind2nd(equal_to<string>(), "3"), mem_fun_ref(&Element::getA));
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~UnaryComposeFunctionAdapter())
{
	delete ptr;
}
END_SECTION

START_SECTION((OP1::result_type operator()(const typename OP2::argument_type &x) const ))
{
  Element a("4"), b("3"), c("2"), d("1");

  vector<Element> elements;
  elements.push_back(a);
  elements.push_back(b);
  elements.push_back(c);
  elements.push_back(d);

  vector<Element>::const_iterator found = find_if(elements.begin(), elements.end(),
    unaryCompose(bind2nd(equal_to<string>(), "3"), mem_fun_ref(&Element::getA)));

  TEST_EQUAL(found != elements.end(), true)
  TEST_EQUAL(found - elements.begin(), 1)

  vector<Element>::const_iterator not_found = find_if(elements.begin(), elements.end(),
    unaryCompose(bind2nd(equal_to<string>(), "10"), mem_fun_ref(&Element::getA)));

  TEST_EQUAL(not_found == elements.end(), true)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



