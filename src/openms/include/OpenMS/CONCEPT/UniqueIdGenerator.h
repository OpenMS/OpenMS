// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_int.hpp>

#include <string>


namespace OpenMS
{

  class DateTime;

  /**
    @brief  A generator for unique ids.

    The unique ids are 64-bit random unsigned random integers.
    The class is implemented as a singleton.
    The random generator is implemented using boost::random.

    @ingroup Concept
  */
  class OPENMS_DLLAPI UniqueIdGenerator
  {

public:

    /// Returns a new unique id
    static UInt64 getUniqueId();

    /**
      @brief Returns a new random UUID (version 4)

      Formatted as the standard 36-character string
      "xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx".
    */
    static std::string getUUID();

    /// Initializes random generator using the given value.
    static void setSeed(const UInt64);

    /// Get the seed
    static UInt64 getSeed();

protected:
    UniqueIdGenerator();
    ~UniqueIdGenerator();

private:
    static UInt64 seed_;
    static UniqueIdGenerator* instance_;
    static boost::mt19937_64* rng_;
    static boost::uniform_int<UInt64>* dist_;

    static UniqueIdGenerator& getInstance_();
    void init_();
    UniqueIdGenerator(const UniqueIdGenerator& );//protect from c++ auto-generation
  };

} // namespace OpenMS

