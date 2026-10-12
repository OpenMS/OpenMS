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
#include <algorithm> // for std::sort and std::adjacent_find
#include <array>
#include <barrier>
#include <ctime>
#include <limits>
#include <random>
#include <set>
#include <string>
#include <thread>
#include <vector>
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
  for (unsigned i = 0; i < nofIdsToGenerate; ++i)
  {
    ids.push_back(OpenMS::UniqueIdGenerator::getUniqueId());
  }
  std::sort(ids.begin(), ids.end());
  // check if the generated ids contain (at least) two equal ones
  std::vector<OpenMS::UInt64>::iterator iter = std::adjacent_find(ids.begin(), ids.end());
  TEST_EQUAL(iter == ids.end(), true);
}
END_SECTION

START_SECTION((static std::string getUUID()))
{
  std::string u = OpenMS::UniqueIdGenerator::getUUID();
  STATUS("OpenMS::UniqueIdGenerator::getUUID(): " << u);
  TEST_EQUAL(u.size(), 36)
  TEST_TRUE(u[8] == '-' && u[13] == '-' && u[18] == '-' && u[23] == '-')

  // version nibble must be '4' (UUIDv4)
  TEST_EQUAL(u[14], '4')

  // variant nibble encodes the 10xx bits -> one of 8, 9, a, b
  const char var = u[19];
  TEST_TRUE(var == '8' || var == '9' || var == 'a' || var == 'b' || var == 'A' || var == 'B')

  // every non-hyphen character is a hex digit
  bool all_hex = true;
  for (size_t i = 0; i < u.size(); ++i)
  {
    if (i == 8 || i == 13 || i == 18 || i == 23) continue;
    const char c = u[i];
    const bool hex = (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f') || (c >= 'A' && c <= 'F');
    if (! hex)
    {
      all_hex = false;
      break;
    }
  }
  TEST_TRUE(all_hex)

  // uniqueness across many draws (no reseeding per call -> collisions negligible)
  std::set<std::string> seen;
  for (int i = 0; i < 1000; ++i)
    seen.insert(OpenMS::UniqueIdGenerator::getUUID());
  TEST_EQUAL(seen.size(), 1000)
}
END_SECTION

START_SECTION((static void setSeed(UInt64 seed)))
{
  constexpr UInt64 one_moment_in_time = 546666321;
  const std::array<UInt64, 6> unique_ids {4039984684862977299ULL,  11561668883169444769ULL, 8153960635892418594ULL,
                                          12940485248168291983ULL, 11522917731873626020ULL, 4387255872055054320ULL};

  OpenMS::UniqueIdGenerator::setSeed(one_moment_in_time);
  for (Size i = 0; i < unique_ids.size(); ++i)
  {
    OpenMS::UInt64 uid = OpenMS::UniqueIdGenerator::getUniqueId();
    TEST_EQUAL(uid, unique_ids[i]);
  }

  // Reseeding must reproduce the sequence, including a partially consumed cache.
  std::vector<OpenMS::UInt64> ids;
  ids.reserve(nofIdsToGenerate);
  OpenMS::UniqueIdGenerator::setSeed(one_moment_in_time);
  for (unsigned i = 0; i < nofIdsToGenerate; ++i)
  {
    ids.push_back(OpenMS::UniqueIdGenerator::getUniqueId());
  }
  OpenMS::UniqueIdGenerator::setSeed(one_moment_in_time);
  for (UInt64 expected : ids)
  {
    const UInt64 actual = OpenMS::UniqueIdGenerator::getUniqueId();
    if (actual != expected) { TEST_EQUAL(actual, expected); }
  }
}
END_SECTION

START_SECTION([EXTRA] multithreaded example)
{
  // Cross many cache boundaries, with a partially consumed final batch per worker.
  constexpr Size count = 1048589;
  UniqueIdGenerator::setSeed(12345);
  std::vector<UInt64> ids(count);
#pragma omp parallel for num_threads(8)
  for (SignedSize i = 0; i < static_cast<SignedSize>(count); ++i)
  {
    ids[i] = UniqueIdGenerator::getUniqueId();
  }
  std::sort(ids.begin(), ids.end());
  TEST_NOT_EQUAL(ids.front(), 0)
  TEST_TRUE(std::adjacent_find(ids.begin(), ids.end()) == ids.end())
}
END_SECTION

START_SECTION((static UInt64 getSeed()))
{
  for (UInt64 seed : {UInt64(0), UInt64(546666321), std::numeric_limits<UInt64>::max()})
  {
    UniqueIdGenerator::setSeed(seed);
    TEST_EQUAL(UniqueIdGenerator::getSeed(), seed)
    UniqueIdGenerator::getUniqueId();
    TEST_EQUAL(UniqueIdGenerator::getSeed(), seed)
  }
}
END_SECTION

START_SECTION([EXTRA] seeded sequence over more than one million ids)
{
  for (UInt64 seed : {UInt64(0), UInt64(546666321), std::numeric_limits<UInt64>::max()})
  {
    UniqueIdGenerator::setSeed(seed);
    std::mt19937_64 reference(seed);
    Size mismatches = 0;
    for (Size i = 0; i < 1000017; ++i)
    {
      UInt64 expected;
      do
      {
        expected = reference();
      } while (expected == 0);
      mismatches += UniqueIdGenerator::getUniqueId() != expected;
    }
    TEST_EQUAL(mismatches, 0)
  }
}
END_SECTION

START_SECTION([EXTRA] UUID and integer ids share the seeded sequence)
{
  constexpr UInt64 seed = 546666321;
  UniqueIdGenerator::setSeed(seed);
  std::mt19937_64 reference(seed);
  // A UUID spans two batches when requested after these integer ids.
  for (Size i = 0; i < 1023; ++i)
  {
    UniqueIdGenerator::getUniqueId();
    reference();
  }
  const std::string uuid = UniqueIdGenerator::getUUID();
  reference();
  reference();
  TEST_EQUAL(UniqueIdGenerator::getUniqueId(), reference())

  UniqueIdGenerator::setSeed(seed);
  for (Size i = 0; i < 1023; ++i)
  {
    UniqueIdGenerator::getUniqueId();
  }
  TEST_EQUAL(UniqueIdGenerator::getUUID(), uuid)
}
END_SECTION

START_SECTION([EXTRA] standard threads generate more than one million ids)
{
  constexpr Size worker_count = 8;
  constexpr Size per_worker = 131073;
  UniqueIdGenerator::setSeed(654321);
  std::vector<UInt64> ids(worker_count * per_worker);
  std::vector<std::thread> workers;
  for (Size worker = 0; worker < worker_count; ++worker)
  {
    workers.emplace_back([&, worker]() {
      for (Size i = 0; i < per_worker; ++i)
      {
        ids[worker * per_worker + i] = UniqueIdGenerator::getUniqueId();
      }
    });
  }
  for (auto& worker : workers)
  {
    worker.join();
  }
  std::sort(ids.begin(), ids.end());
  TEST_NOT_EQUAL(ids.front(), 0)
  TEST_TRUE(std::adjacent_find(ids.begin(), ids.end()) == ids.end())
}
END_SECTION

START_SECTION([EXTRA] reseeding invalidates caches of persistent standard threads)
{
  constexpr Size worker_count = 8;
  constexpr Size per_worker = 3;
  const std::array<UInt64, 4> seeds {12345, 12345, 98765, 98765};
  std::array<std::vector<UInt64>, seeds.size()> ids;
  for (auto& epoch_ids : ids)
  {
    epoch_ids.resize(worker_count * per_worker);
  }
  std::barrier phase(static_cast<std::ptrdiff_t>(worker_count + 1));
  std::vector<std::thread> workers;
  for (Size worker = 0; worker < worker_count; ++worker)
  {
    workers.emplace_back([&, worker]() {
      for (Size epoch = 0; epoch < seeds.size(); ++epoch)
      {
        phase.arrive_and_wait();
        for (Size i = 0; i < per_worker; ++i)
        {
          ids[epoch][worker * per_worker + i] = UniqueIdGenerator::getUniqueId();
        }
        phase.arrive_and_wait();
      }
    });
  }
  for (UInt64 seed : seeds)
  {
    // Both barriers ensure reseeding never overlaps worker calls. The worker
    // threads survive each reset and retain any partially consumed TLS cache.
    UniqueIdGenerator::setSeed(seed);
    phase.arrive_and_wait();
    phase.arrive_and_wait();
  }
  for (auto& worker : workers)
  {
    worker.join();
  }
  for (auto& epoch_ids : ids)
  {
    std::sort(epoch_ids.begin(), epoch_ids.end());
    TEST_NOT_EQUAL(epoch_ids.front(), 0)
    TEST_TRUE(std::adjacent_find(epoch_ids.begin(), epoch_ids.end()) == epoch_ids.end())
  }
  TEST_TRUE(ids[0] == ids[1])
  TEST_TRUE(ids[2] == ids[3])
  TEST_TRUE(ids[0] != ids[2])
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
