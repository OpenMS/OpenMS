// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <string>


namespace OpenMS
{

/**
  @brief  A generator for unique ids.

  The unique ids are nonzero 64-bit random unsigned integers.
  All callers share one generator with process lifetime.
  The random generator is std::mt19937_64 (MT19937-64), whose output sequence is
  fully specified by the algorithm and therefore identical on every platform and
  standard library. Zero engine outputs are skipped because they represent an
  invalid unique id. Random values can still collide; map-level uniqueness
  must be checked separately.

  Each thread caches 1024 engine outputs, acquiring a shared mutex only when
  refilling its cache. Generation is thread-safe with or without OpenMP.
  Calls made by one thread preserve the seeded sequence, except for skipped
  zeros. With multiple generating threads, values are reserved in batches;
  unused values are discarded when a thread terminates or setSeed() is called.
  Parallel output therefore depends on scheduling and thread lifetime.

  @ingroup Concept
*/
class OPENMS_DLLAPI UniqueIdGenerator
{

public:
  /// Returns a new, nonzero unique id.
  static UInt64 getUniqueId();

  /**
    @brief Returns a new random UUID (version 4)

    Formatted as the standard 36-character string
    "xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx".
  */
  static std::string getUUID();

  /**
    @brief Reseed the shared generator and invalidate every thread's cached ids.

    @param[in] seed Seed for the Mersenne Twister engine.

    @pre No other thread may be generating ids or UUIDs during reseeding.
  */
  static void setSeed(UInt64 seed);

  /// Get the seed
  static UInt64 getSeed();

protected:
  UniqueIdGenerator();
  ~UniqueIdGenerator();

private:
  UniqueIdGenerator(const UniqueIdGenerator&); // protect from c++ auto-generation
};

} // namespace OpenMS
