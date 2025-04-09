// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <ctime>
#include <algorithm> // for std::sort and std::adjacent_find
// array_wrapper needs to be included before it is used
// only in boost1.64+. See issue #2790
#if OPENMS_BOOST_VERSION_MINOR >= 64
#include <boost/serialization/array_wrapper.hpp>
#endif
#include <boost/accumulators/statistics/covariance.hpp>
#include <boost/typeof/incr_registration_group.hpp>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(UniqueIdGenerator, "$Id$")

unsigned nofIdsToGenerate = 100000;
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
  /* test for collisions, test will be different for every test execution */
  OpenMS::UniqueIdGenerator::setSeed(std::time(nullptr));
  std::vector<OpenMS::UInt64> ids;
  ids.reserve(nofIdsToGenerate);
  for (unsigned i=0; i<nofIdsToGenerate; ++i)
  {
    ids.push_back(OpenMS::UniqueIdGenerator::getUniqueId());
  }
  std::sort(ids.begin(), ids.end());
  // check if the generated ids contain (at least) two equal ones
  std::vector<OpenMS::UInt64>::iterator iter = std::adjacent_find(ids.begin(), ids.end());
  TEST_EQUAL(iter == ids.end(), true);
}
END_SECTION

START_SECTION((static void setSeed(UInt seed)))
{
  UInt one_moment_in_time = 546666321;
//  OpenMS::DateTime one_moment_in_time;
//  one_moment_in_time.set(5,4,6666,3,2,1);
//  OpenMS::UniqueIdGenerator::setSeed(one_moment_in_time);

  /* check if the generator changed */
  UInt64 large_int = 0;
  std::vector<UInt64> unique_ids;
  large_int = 4039984684862977299U;
  unique_ids.push_back(large_int);
  large_int = 11561668883169444769U;
  unique_ids.push_back(large_int);
  large_int = 8153960635892418594U;
  unique_ids.push_back(large_int);
  large_int = 12940485248168291983U;
  unique_ids.push_back(large_int);
  large_int = 11522917731873626020U;
  unique_ids.push_back(large_int);
  large_int = 4387255872055054320U;
  unique_ids.push_back(large_int);

  OpenMS::UniqueIdGenerator::setSeed(one_moment_in_time);
  for (Size i=0; i<unique_ids.size(); ++i)
  {
    OpenMS::UInt64 uid = OpenMS::UniqueIdGenerator::getUniqueId();
    TEST_EQUAL(uid,unique_ids[i]);
  }

  /* check if the same sequence is generated form the same seed */
  std::vector<OpenMS::UInt64> ids;
  ids.reserve(nofIdsToGenerate);
  OpenMS::UniqueIdGenerator::setSeed(one_moment_in_time);
  for (unsigned i=0; i<nofIdsToGenerate; ++i)
  {
    ids.push_back(OpenMS::UniqueIdGenerator::getUniqueId());
  }
  std::vector<OpenMS::UInt64> ids2;
  ids2.reserve(nofIdsToGenerate);
  OpenMS::UniqueIdGenerator::setSeed(one_moment_in_time);
  for (unsigned i=0; i<nofIdsToGenerate; ++i)
  {
    ids2.push_back(OpenMS::UniqueIdGenerator::getUniqueId());
  }

  for ( unsigned i = 0; i < nofIdsToGenerate; ++i )
  {
    if(ids[i] != ids2[i])
    {
      TEST_EQUAL(ids[i],  ids2[i]);
    }
  }
}
END_SECTION

START_SECTION([EXTRA] multithreaded example)
{

   /* test for collisions, test will be different for every test execution */
  OpenMS::UniqueIdGenerator::setSeed(std::time(nullptr));
  std::vector<OpenMS::UInt64> ids;
  ids.reserve(nofIdsToGenerate);
#pragma omp parallel for
  for (int i = 0; i < static_cast<int>(nofIdsToGenerate); ++i)
  {
    OpenMS::UInt64 tmp = OpenMS::UniqueIdGenerator::getUniqueId();
#pragma omp critical (add_test)
    {
      ids.push_back(tmp);
    }
  }
  std::sort(ids.begin(), ids.end());
  // check if the generated ids contain (at least) two equal ones
  std::vector<OpenMS::UInt64>::iterator iter = std::adjacent_find(ids.begin(), ids.end());
  TEST_EQUAL(iter == ids.end(), true);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
