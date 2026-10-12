// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/SYSTEM/SysInfo.h>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <mutex>
#include <random>

namespace OpenMS
{
namespace
{
  UInt64 initialSeed()
  {
    // Microsecond time and the process ID distinguish tools started close together.
    const auto now = std::chrono::system_clock::now().time_since_epoch();
    return static_cast<UInt64>(std::chrono::duration_cast<std::chrono::microseconds>(now).count())
           ^ (static_cast<UInt64>(SysInfo::getProcessId()) << 32);
  }

  struct GeneratorState
  {
    std::mutex mutex;
    UInt64 seed {initialSeed()};
    std::mt19937_64 rng {seed};
  };

  GeneratorState& generatorState()
  {
    // Preserve the original singleton's process lifetime for shutdown callers.
    static auto* const state = new GeneratorState;
    return *state;
  }
} // namespace

UInt64 UniqueIdGenerator::getUniqueId()
{
  auto& state = generatorState();
  // Every draw advances the one shared engine, for both OpenMP and standard
  // threads, so the ids drawn for a seed do not depend on thread count or scheduling.
  const std::lock_guard lock(state.mutex);
  UInt64 value;
  do
  {
    value = state.rng();
  } while (value == 0); // UniqueIdInterface reserves zero for an invalid ID.
  return value;
}

std::string UniqueIdGenerator::getUUID()
{
  // Two independent 64-bit draws from the shared, once-seeded generator give us
  // the 128 bits of a UUID without reseeding a fresh RNG per call (which would
  // depend on std::random_device's entropy quality on every single call).
  UInt64 hi = getUniqueId();
  UInt64 lo = getUniqueId();
  Byte bytes[16];
  std::memcpy(bytes, &hi, 8);
  std::memcpy(bytes + 8, &lo, 8);
  bytes[6] = (bytes[6] & 0x0F) | 0x40; // version 4
  bytes[8] = (bytes[8] & 0x3F) | 0x80; // variant 1
  char buf[37];
  std::snprintf(buf, sizeof(buf), "%02x%02x%02x%02x-%02x%02x-%02x%02x-%02x%02x-%02x%02x%02x%02x%02x%02x", bytes[0], bytes[1], bytes[2], bytes[3],
                bytes[4], bytes[5], bytes[6], bytes[7], bytes[8], bytes[9], bytes[10], bytes[11], bytes[12], bytes[13], bytes[14], bytes[15]);
  return std::string(buf);
}

UInt64 UniqueIdGenerator::getSeed()
{
  auto& state = generatorState();
  const std::lock_guard lock(state.mutex);
  return state.seed;
}

void UniqueIdGenerator::setSeed(UInt64 seed)
{
  auto& state = generatorState();
  const std::lock_guard lock(state.mutex);
  state.seed = seed;
  state.rng.seed(seed);
}

UniqueIdGenerator::UniqueIdGenerator() = default;
UniqueIdGenerator::~UniqueIdGenerator() = default;

} // namespace OpenMS
