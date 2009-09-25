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

#include <OpenMS/CONCEPT/UniqueIdGenerator.h>

// debugging
#define V_UniqueIdGenerator(a) ;
// #define V_UniqueIdGenerator(a) std::cout << "line " << __LINE__ << ": " << a << std::endl;

namespace OpenMS
{

  UniqueIdGenerator::UniqueId
  UniqueIdGenerator::getUniqueId()
  {
    V_UniqueIdGenerator("UniqueIdGenerator::getUID()");
    gsl_rng * const rng = getInstance_().random_number_generator_;
    return (UniqueIdGenerator::UniqueId(gsl_rng_get(rng)) << 32)
      + UniqueIdGenerator::UniqueId(gsl_rng_get(rng));
  }

  const Param&
  UniqueIdGenerator::getInfo()
  {
    V_UniqueIdGenerator("UniqueIdGenerator::getInfo()");
    return getInstance_().info_;
  }

  UniqueIdGenerator&
  UniqueIdGenerator::getInstance_()
  {
    V_UniqueIdGenerator("UniqueIdGenerator::getInstance_()");
    static UniqueIdGenerator* instance_ = 0;
    if ( !instance_ )
    {
      instance_ = new UniqueIdGenerator();
    }
    return *instance_;
  }

  UniqueIdGenerator::UniqueIdGenerator()
  {
    V_UniqueIdGenerator("UniqueIdGenerator::UniqueIdGenerator()");
    random_number_generator_ = gsl_rng_alloc(gsl_rng_mt19937);
    init_(OpenMS::DateTime::now());
    return;
  }

  void
  UniqueIdGenerator::setSeed(const DateTime& date_time)
  {
    V_UniqueIdGenerator("UniqueIdGenerator::setSeed()");
    getInstance_().init_(date_time);
    return;
  }

  void
  UniqueIdGenerator::init_(const DateTime& date_time)
  {
    V_UniqueIdGenerator("UniqueIdGenerator::init_()");

    initialization_date_time_ = date_time;
    const UInt64 seed_64 = initialization_date_time_.toString("yyyyMMddhhmmsszzz").toLongLong();
    const unsigned long int actually_used_seed = ((UInt64(1)<<32)-1)&((seed_64>>32)^seed_64); // just to mix the bits a bit
    gsl_rng_set(random_number_generator_, actually_used_seed);

    info_.setValue("generator_type", gsl_rng_name(random_number_generator_));
    info_.setValue("generator_min", String(gsl_rng_min(random_number_generator_)));
    info_.setValue("generator_max", String(gsl_rng_max(random_number_generator_)));
    info_.setValue("initialization_date_time_as_string", initialization_date_time_.get());
    info_.setValue("initialization_date_time_as_longlong", String(seed_64));
    info_.setValue("actually_used_seed", String(actually_used_seed));
    V_UniqueIdGenerator("info:\n" << info_);

    return;
  }

  UniqueIdGenerator::~UniqueIdGenerator()
  {
    V_UniqueIdGenerator("UniqueIdGenerator::~UniqueIdGenerator()");
    gsl_rng_free( random_number_generator_);
    return;
  }

}
