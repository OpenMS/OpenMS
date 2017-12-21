// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/UniqueIdGenerator.h>

#include <boost/date_time/posix_time/posix_time_types.hpp> //no i/o just types
#include <limits>
#include <iostream>

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

  UniqueIdGenerator::UniqueIdGenerator()
  {
  }

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
