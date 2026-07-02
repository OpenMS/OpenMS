// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/UniqueIdGenerator.h>

#include <OpenMS/SYSTEM/SysInfo.h>

#include <boost/date_time/posix_time/posix_time_types.hpp> //no i/o just types

#include <cstdint>
#include <cstdio>
#include <cstring>

namespace OpenMS
{
  UInt64 UniqueIdGenerator::seed_ = 0;
  UniqueIdGenerator* UniqueIdGenerator::instance_ = nullptr;
  boost::mt19937_64* UniqueIdGenerator::rng_ = nullptr;
  boost::uniform_int<UInt64>* UniqueIdGenerator::dist_ = nullptr;

  UInt64 UniqueIdGenerator::getUniqueId()
  {
    UniqueIdGenerator& instance = getInstance_();
#ifdef _OPENMP
    UInt64 val;
#pragma omp critical (OPENMS_UniqueIdGenerator_getUniqueId)
    {
      val = (*instance.dist_)(*instance.rng_);
    }
    // note: OpenMP can only work on a structured block, return needs to be outside that block
    return val; 
#else
    return (*instance.dist_)(*instance.rng_);
#endif
  }

  std::string UniqueIdGenerator::getUUID()
  {
    // Two independent 64-bit draws from the shared, once-seeded generator give us
    // the 128 bits of a UUID without reseeding a fresh RNG per call (which would
    // depend on std::random_device's entropy quality on every single call).
    UInt64 hi = getUniqueId();
    UInt64 lo = getUniqueId();
    uint8_t bytes[16];
    std::memcpy(bytes, &hi, 8);
    std::memcpy(bytes + 8, &lo, 8);
    bytes[6] = (bytes[6] & 0x0F) | 0x40; // version 4
    bytes[8] = (bytes[8] & 0x3F) | 0x80; // variant 1
    char buf[37];
    std::snprintf(buf, sizeof(buf),
      "%02x%02x%02x%02x-%02x%02x-%02x%02x-%02x%02x-%02x%02x%02x%02x%02x%02x",
      bytes[0], bytes[1], bytes[2], bytes[3],
      bytes[4], bytes[5], bytes[6], bytes[7],
      bytes[8], bytes[9], bytes[10], bytes[11],
      bytes[12], bytes[13], bytes[14], bytes[15]);
    return std::string(buf);
  }

  UInt64 UniqueIdGenerator::getSeed()
  {
    return getInstance_().seed_;
  }

  void UniqueIdGenerator::setSeed(UInt64 seed)
  {
  // modifies static members
#ifdef _OPENMP
#pragma omp critical (OPENMS_UniqueIdGenerator_setSeed)
#endif
    {
      UniqueIdGenerator& instance = getInstance_();
      instance.seed_ = seed;
      instance.rng_->seed( instance.seed_ );
      instance.dist_->reset();
    }
  }

  UniqueIdGenerator::UniqueIdGenerator() = default;

  UniqueIdGenerator & UniqueIdGenerator::getInstance_()
  {
  // modifies static members
#ifdef _OPENMP
#pragma omp critical (OPENMS_UniqueIdGenerator_getInstance_)
#endif
    {
      if (!instance_)
      {
        instance_ = new UniqueIdGenerator();
        instance_->init_();
      }
    }
    return *instance_;
  }

  void UniqueIdGenerator::init_()
  {
  // modifies static members
#ifdef _OPENMP
#pragma omp critical (OPENMS_UniqueIdGenerator_init_)
#endif
    {
      // find a seed:
      // get something with high resolution (around microseconds) -- its hard to do better on Windows --
      // which has absolute system time (there is higher resolution available for the time since program startup, but
      // we do not want this here since this seed usually gets initialized at the same program uptime).
      // Reason for high-res: in pipelines, instances of TOPP tools can get initialized almost simultaneously (i.e., resolution in seconds is not enough),
      // leading to identical random numbers (e.g. feature-IDs) in two or more distinct files.
      // C++11 note: C++ build-in alternative once C++11 can be presumed: 'std::chrono::high_resolution_clock'
      boost::posix_time::ptime t(boost::posix_time::microsec_clock::local_time() );
      seed_ = t.time_of_day().ticks();  // independent of implementation; as opposed to nanoseconds(), which need not be available on every platform

      // Mix in process ID to ensure different seeds for parallel processes started at the same microsecond.
      // This prevents duplicate UniqueIds when multiple TOPP tools run simultaneously in a pipeline.
      seed_ ^= static_cast<UInt64>(SysInfo::getProcessId()) << 32;

      rng_ = new boost::mt19937_64 (seed_);
      dist_ = new boost::uniform_int<UInt64> (0, std::numeric_limits<UInt64>::max());
    }
  }

  UniqueIdGenerator::~UniqueIdGenerator()
  {
    delete rng_;
    delete dist_;
  }

}
