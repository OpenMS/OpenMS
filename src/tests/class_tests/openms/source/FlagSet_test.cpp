// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/FlagSet.h>

#include <iostream>
//////////////////////////

using namespace OpenMS;

enum Enum_test
{
  e_red,
  e_green,
  e_blue
};

enum class Enum_broken
{
  eb_red = -1, // too small
  eb_green,
  eb_blue = 64 // too large
};

using FST = FlagSet<Enum_test>;

using FSB = FlagSet<Enum_broken>;

std::ostream& operator<<(std::ostream& str, FST s)
{
  return str << s.value();
}

START_TEST(FlagSet, "$Id$")
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////
//start Section
/////////////////////////////////////////////////////////////////



FlagSet<Enum_test>* ptr = nullptr;
FlagSet<Enum_test>* nulpt = nullptr;
START_SECTION(FlagSet())
{
  ptr = new FlagSet<Enum_test>();
  TEST_NOT_EQUAL(ptr, nulpt)
}
END_SECTION


START_SECTION(explicit FlagSet(const ENUM& en))
{
  TEST_NOT_EQUAL(FST(e_red), FST())
  TEST_EQUAL(FST(e_red).value(), 1)
  TEST_PRECONDITION_VIOLATED(FSB{Enum_broken::eb_red}) // negative value
  TEST_PRECONDITION_VIOLATED(FSB{Enum_broken::eb_blue}) // too large for uint64
}
END_SECTION

START_SECTION(FlagSet(const FlagSet& stat))
  TEST_EQUAL(FST(FST(e_green)).value(), 2)
END_SECTION

START_SECTION(FlagSet& operator=(const FlagSet & stat))
  FST gg(e_green);
  FST target = gg;
  TEST_EQUAL(target, gg)
END_SECTION


START_SECTION(~FlagSet())
  NOT_TESTABLE
END_SECTION

/// Equality
START_SECTION(bool operator==(const FlagSet & stat) const)
{
    FST gg(e_green);
    FST target = gg;
    TEST_EQUAL(target == gg, true)
    TEST_NOT_EQUAL(gg, FST())
}
END_SECTION

/// bitwise AND
START_SECTION(FlagSet operator&(const ENUM & en) const)
{
  FST gg(e_green);
  FST empty = gg & e_red;
  FST just_green = gg & e_green;
  TEST_EQUAL(empty, FST())
  TEST_EQUAL(just_green, gg)
}
END_SECTION

/// bitwise AND
START_SECTION(FlagSet operator&(const FlagSet & rhs) const)
{
  FST gg(e_green);
  FST rr(e_red);
  FST empty = gg & rr;
  FST just_green = gg & gg;
  TEST_EQUAL(empty, FST())
  TEST_EQUAL(just_green, gg)
}
END_SECTION

START_SECTION(FlagSet& operator&=(const ENUM & en))
{
  FST gg(e_green);
  FST empty = gg;
  empty &= e_red;
  FST just_green = gg;
  just_green &= e_green;
  TEST_EQUAL(empty, FST())
  TEST_EQUAL(just_green, gg)
}
END_SECTION


START_SECTION(FlagSet& operator&=(const FlagSet & rhs))
{
  FST gg(e_green);
  FST rr(e_red);
  FST empty = gg;
  empty &= rr;
  FST just_green = gg;
  just_green &= gg;
  TEST_EQUAL(empty, FST())
  TEST_EQUAL(just_green, gg)
}
END_SECTION

START_SECTION(FlagSet operator|(const ENUM & en) const)
{
  FST gg(e_green);
  FST green_or_red = gg | e_red;
  FST green_or_green = gg | e_green;
  TEST_EQUAL(green_or_red.value(), 3)
  TEST_EQUAL(green_or_green, gg)
}
END_SECTION


START_SECTION(FlagSet operator|(const FlagSet & rhs) const)
{
  FST gg;
  FST empty_or_red = gg | e_red;
  FST red_or_green = empty_or_red | e_green;
  TEST_EQUAL(empty_or_red, FST(e_red))
  TEST_EQUAL(red_or_green.value(), 3)
}
END_SECTION

START_SECTION(FlagSet& operator|=(const ENUM & en))
{
  FST gg(e_green);
  FST green_or_red = gg;
  green_or_red |= e_red;
  FST green_or_green = gg;
  green_or_green |= e_green;
  TEST_EQUAL(green_or_red.value(), 3)
  TEST_EQUAL(green_or_green, gg)
}
END_SECTION

START_SECTION(FlagSet& operator|=(const FlagSet & rhs))
{
  FST gg;
  FST empty_or_red = gg;
  empty_or_red |= e_red;
  FST red_or_green = empty_or_red;
  red_or_green |= e_green;
  TEST_EQUAL(empty_or_red, FST(e_red))
  TEST_EQUAL(red_or_green.value(), 3)
}
END_SECTION

START_SECTION(FlagSet operator+(const ENUM & en) const)
{
  FST gg(e_green);
  FST green_or_red = gg + e_red;
  FST green_or_green = gg + e_green;
  TEST_EQUAL(green_or_red.value(), 3)
  TEST_EQUAL(green_or_green, gg)
}
END_SECTION

START_SECTION(FlagSet operator+(const FlagSet & en) const)
{
  FST gg;
  FST empty_or_red = gg + e_red;
  FST red_or_green = empty_or_red + e_green;
  TEST_EQUAL(empty_or_red, FST(e_red))
  TEST_EQUAL(red_or_green.value(), 3)
}
END_SECTION

START_SECTION(FlagSet& operator+=(const ENUM & rhs))
{
  FST gg(e_green);
  FST green_or_red = gg;
  green_or_red += e_red;
  FST green_or_green = gg;
  green_or_green += e_green;
  TEST_EQUAL(green_or_red.value(), 3)
  TEST_EQUAL(green_or_green, gg)
}
END_SECTION

START_SECTION(FlagSet& operator+=(const FlagSet & rhs))
{
  FST gg;
  FST empty_or_red = gg;
  empty_or_red += e_red;
  FST red_or_green = empty_or_red;
  red_or_green += e_green;
  TEST_EQUAL(empty_or_red, FST(e_red))
  TEST_EQUAL(red_or_green.value(), 3)
}
END_SECTION

START_SECTION(FlagSet operator-(const FlagSet & rhs))
{
  FST gg;
  FST empty = gg - FST(e_red);
  TEST_EQUAL(FST(), empty)
  FST red_or_green = (FST(e_red) + e_green);
  FST red_or_green_no_blue = red_or_green - FST(e_blue);
  TEST_EQUAL(red_or_green, red_or_green_no_blue)
  FST red_only = red_or_green - e_green;
  TEST_EQUAL(red_only, FST(e_red))
}
END_SECTION

START_SECTION(FlagSet& operator-=(const FlagSet & rhs))
{
  FST gg;
  FST empty = gg;
  empty -= FST(e_red);
  TEST_EQUAL(FST(), empty)
  FST red_or_green = (FST(e_red) + e_green);
  FST red_or_green_no_blue = red_or_green;
  red_or_green_no_blue -= FST(e_blue);
  TEST_EQUAL(red_or_green, red_or_green_no_blue)
  FST red_only = red_or_green;
  red_only -= FST(e_green);
  TEST_EQUAL(red_only, FST(e_red))
}
END_SECTION

START_SECTION(FlagSet operator-(const ENUM & rhs))
{
  FST gg;
  FST empty = gg - e_red;
  TEST_EQUAL(FST(), empty)
  FST red_or_green = (FST(e_red) + e_green);
  FST red_or_green_no_blue = red_or_green - e_blue;
  TEST_EQUAL(red_or_green, red_or_green_no_blue)
  FST red_only = red_or_green - e_green;
  TEST_EQUAL(red_only, FST(e_red))
}
END_SECTION

START_SECTION(FlagSet& operator-=(const ENUM & rhs))
{
  FST gg;
  FST empty = gg;
  empty -= e_red;
  TEST_EQUAL(FST(), empty)
  FST red_or_green = (FST(e_red) + e_green);
  FST red_or_green_no_blue = red_or_green;
  red_or_green_no_blue -= e_blue;
  TEST_EQUAL(red_or_green, red_or_green_no_blue)
  FST red_only = red_or_green;
  red_only -= e_green;
  TEST_EQUAL(red_only, FST(e_red))
}
END_SECTION

START_SECTION(bool isSuperSetOf(const FlagSet & required) const)
{
  FST gg;
  FST empty = gg - e_red;
  TEST_EQUAL(gg.isSuperSetOf(empty), true)
  TEST_EQUAL(empty.isSuperSetOf(gg), true)
  FST red_or_green = (FST(e_red) + e_green);
  FST red_or_green_or_blue = red_or_green + e_blue;
  TEST_EQUAL(red_or_green_or_blue.isSuperSetOf(red_or_green), true)
  TEST_EQUAL(red_or_green_or_blue.isSuperSetOf(red_or_green_or_blue), true)
  TEST_EQUAL(red_or_green_or_blue.isSuperSetOf(empty), true)

  TEST_EQUAL(red_or_green.isSuperSetOf(red_or_green_or_blue), false)
  TEST_EQUAL(empty.isSuperSetOf(red_or_green_or_blue), false)
  TEST_EQUAL(FST(e_red).isSuperSetOf(red_or_green_or_blue), false)
}
END_SECTION

START_SECTION(bool isSuperSetOf(const ENUM & required) const)
{
  FST empty;
  TEST_EQUAL(empty.isSuperSetOf(e_red), false)
  FST red_or_green = (FST(e_red) + e_green);
  FST red_or_green_or_blue = red_or_green + e_blue;
  TEST_EQUAL(red_or_green_or_blue.isSuperSetOf(e_red), true)
  TEST_EQUAL(red_or_green_or_blue.isSuperSetOf(e_blue), true)
  TEST_EQUAL(red_or_green_or_blue.isSuperSetOf(e_green), true)
}
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
