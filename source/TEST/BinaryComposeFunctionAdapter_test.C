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
#include <OpenMS/CONCEPT/BinaryComposeFunctionAdapter.h>
///////////////////////////

#include <string>
#include <vector>
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

START_TEST(BinaryComposeFunctionAdapter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
typedef std::less<std::basic_string<char, std::char_traits<char>, std::allocator<char> > > T_STD_LESS;
typedef std::const_mem_fun_ref_t<const std::string&, Element> T_MEMBER_FUNCTIOM;

typedef BinaryComposeFunctionAdapter<T_STD_LESS,T_MEMBER_FUNCTIOM,T_MEMBER_FUNCTIOM > BCFA;

BCFA* ptr = 0;
BCFA* null_ptr = 0;

START_SECTION((BinaryComposeFunctionAdapter(const OP1 &o1, const OP2 &o2, const OP3 &o3)))
{
  ptr = new BCFA(less<string>(), mem_fun_ref(&Element::getA), mem_fun_ref(&Element::getA));
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~BinaryComposeFunctionAdapter())
{
	delete ptr;
}
END_SECTION

START_SECTION((OP1::result_type operator()(const typename OP2::argument_type &x, const typename OP3::argument_type &y) const ))
{
  Element a("Matthias"), b("Marcel"), c("Anton"), d("Henner");

  vector<Element> elements;
  elements.push_back(a);
  elements.push_back(b);
  elements.push_back(c);
  elements.push_back(d);

  // the below function sorts elements based on the &Element::getA result value
  sort(elements.begin(), elements.end(),
              binaryCompose(
                      less<string>(),
                      mem_fun_ref(&Element::getA),
                      mem_fun_ref(&Element::getA)));

  TEST_EQUAL(elements.size(), 4)
  TEST_EQUAL(elements[0].getA(), "Anton")
  TEST_EQUAL(elements[1].getA(), "Henner")
  TEST_EQUAL(elements[2].getA(), "Marcel")
  TEST_EQUAL(elements[3].getA(), "Matthias")
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



