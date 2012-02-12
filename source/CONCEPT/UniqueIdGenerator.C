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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/UniqueIdGenerator.h>

#include <gsl/gsl_rng.h>


// debugging
#define V_UniqueIdGenerator(a) ;
// #define V_UniqueIdGenerator(a) std::cout << "line " << __LINE__ << ": " << a << std::endl;

namespace OpenMS
{

namespace
{

// Defined outside so that we can initialize from different places.
// Note the "=0" and read the following comment.
static UniqueIdGenerator* instance_ = 0;

// Do not add a "=0" here:  The initialization is taken care of by init_() or setSeed().
// Remember: The order of execution of static initializers is unspecified.
// Well, actually that paragraph does not apply here (instance_ is just a pointer, not an object)...
// But don't blame me ... just in case ...
static gsl_rng * rng_;

}


UInt64
UniqueIdGenerator::getUniqueId()
{
  V_UniqueIdGenerator("UniqueIdGenerator::getUID()");
  getInstance_();
  UInt64 r;
#pragma omp critical  
  {r = (UInt64(gsl_rng_get(rng_)) << 32) + UInt64(gsl_rng_get(rng_));}
  return r;
}

const Param&
UniqueIdGenerator::getInfo()
{
  V_UniqueIdGenerator("UniqueIdGenerator::getInfo()");
  return getInstance_().info_;
}

// We can't add an underscore to a private constructor, sorry!  Even templates won't help.  ;-)
UniqueIdGenerator::UniqueIdGenerator()
{
  V_UniqueIdGenerator("UniqueIdGenerator::UniqueIdGenerator()");
  rng_ = gsl_rng_alloc(gsl_rng_mt19937);
  // The random seed is set by a call to init_()
  // from within either getInstance_() or setSeed(),
  // depending upon what is called first.
  return;
}

UniqueIdGenerator&
UniqueIdGenerator::getInstance_()
{
  V_UniqueIdGenerator("UniqueIdGenerator::getInstance_()");
  if ( !instance_ )
  {
    instance_ = new UniqueIdGenerator();
    instance_->init_(OpenMS::DateTime::now());
  }
  return *instance_;
}

void
UniqueIdGenerator::setSeed(const DateTime& date_time)
{
  V_UniqueIdGenerator("UniqueIdGenerator::setSeed()");
  if ( !instance_ )
  {
    instance_ = new UniqueIdGenerator();
  }
  instance_->init_(date_time);
  return;
}

void
UniqueIdGenerator::init_(const DateTime& date_time)
{
  V_UniqueIdGenerator("UniqueIdGenerator::init_()");

  const UInt64 seed_64 = date_time.toString("yyyyMMddhhmmsszzz").toLongLong();
  const unsigned long int actually_used_seed = ((UInt64(1) << 32) - 1) & ((seed_64 >> 32) ^ seed_64); // just to mix the bits a bit
  gsl_rng_set(rng_, actually_used_seed);

  info_.setValue("generator_type", gsl_rng_name(rng_));
  info_.setValue("generator_min", String(gsl_rng_min(rng_)));
  info_.setValue("generator_max", String(gsl_rng_max(rng_)));
  info_.setValue("initialization_date_time_as_string", date_time.get());
  info_.setValue("initialization_date_time_as_longlong", String(seed_64));
  info_.setValue("actually_used_seed", String(actually_used_seed));
  V_UniqueIdGenerator("info:\n" << info_);

  return;
}

UniqueIdGenerator::~UniqueIdGenerator()
{
  V_UniqueIdGenerator("UniqueIdGenerator::~UniqueIdGenerator()");
  gsl_rng_free(rng_);
  return;
}

}
