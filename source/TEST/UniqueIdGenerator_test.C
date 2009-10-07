// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(UniqueIdGenerator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((UniqueIdGenerator()))
{
  // singleton has private ctor
  NOT_TESTABLE;
}
END_SECTION

START_SECTION((~UniqueIdGenerator()))
{
  // singleton has private dtor
  NOT_TESTABLE;
}
END_SECTION

START_SECTION((static UInt64 getUniqueId()))
{
  STATUS("OpenMS::UniqueIdGenerator::getUniqueId(): " << OpenMS::UniqueIdGenerator::getUniqueId());
  // the actual values are unpredictable, but see setSeed() below
  NOT_TESTABLE;
}
END_SECTION

START_SECTION((static void setSeed(const DateTime &)))
{
  OpenMS::DateTime one_moment_in_time;
  one_moment_in_time.set(5,4,666,3,2,1);
  OpenMS::UniqueIdGenerator::setSeed(one_moment_in_time);

  // hoping that your compiler already supports the ull suffix for unsigned long long (aka Int64Type) integer literals
  OpenMS::UInt64 unique_ids[] =
                  { 1663028827116059880ull, 16317018546938646277ull, 4499304917001700489ull, 14007980103328265649ull, 8982057078544736839ull };

  const int num_num = sizeof(unique_ids)/sizeof(*unique_ids);

  for ( int i = 0; i < num_num; ++i )
  {
    OpenMS::UInt64 uid = OpenMS::UniqueIdGenerator::getUniqueId();
    TEST_EQUAL(uid,unique_ids[i]);
  }

  OpenMS::UniqueIdGenerator::setSeed(one_moment_in_time);

  for ( int i = 0; i < num_num; ++i )
  {
    OpenMS::UInt64 uid = OpenMS::UniqueIdGenerator::getUniqueId();
    TEST_EQUAL(uid,unique_ids[i]);
  }

}
END_SECTION

START_SECTION((static Param const& getInfo()))
{
  STATUS(std::endl << OpenMS::UniqueIdGenerator::getInfo());
  Param const & param(OpenMS::UniqueIdGenerator::getInfo());
  TEST_STRING_EQUAL(param.getValue("generator_type"),"mt19937");
  TEST_STRING_EQUAL(param.getValue("generator_min"),"0");
  TEST_STRING_EQUAL(param.getValue("generator_max"),"4294967295");
  TEST_STRING_EQUAL(param.getValue("initialization_date_time_as_string"),"0666-05-04 03:02:01");
  TEST_STRING_EQUAL(param.getValue("initialization_date_time_as_longlong"),"6660504030201000");
  TEST_STRING_EQUAL(param.getValue("actually_used_seed"),"1890151705");
  STATUS(OpenMS::UniqueIdGenerator::getInfo());
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
