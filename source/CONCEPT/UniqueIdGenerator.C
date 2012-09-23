// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
    static UniqueIdGenerator * instance_ = 0;

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
    {
      r = (UInt64(gsl_rng_get(rng_)) << 32) + UInt64(gsl_rng_get(rng_));
    }
    return r;
  }

  const Param &
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

  UniqueIdGenerator &
  UniqueIdGenerator::getInstance_()
  {
    V_UniqueIdGenerator("UniqueIdGenerator::getInstance_()");
    if (!instance_)
    {
      instance_ = new UniqueIdGenerator();
      instance_->init_(OpenMS::DateTime::now());
    }
    return *instance_;
  }

  void
  UniqueIdGenerator::setSeed(const DateTime & date_time)
  {
    V_UniqueIdGenerator("UniqueIdGenerator::setSeed()");
    if (!instance_)
    {
      instance_ = new UniqueIdGenerator();
    }
    instance_->init_(date_time);
    return;
  }

  void
  UniqueIdGenerator::init_(const DateTime & date_time)
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
