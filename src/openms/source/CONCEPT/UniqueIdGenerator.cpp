// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
#include <ctime>
#include <limits>
#include <iostream>

namespace OpenMS
{
  UInt64 UniqueIdGenerator::seed_ = 0;
  UniqueIdGenerator* UniqueIdGenerator::instance_ = NULL;
  boost::random::mt19937_64* UniqueIdGenerator::rng_ = NULL;
  boost::random::uniform_int_distribution<UInt64>* UniqueIdGenerator::dist_ = NULL;

  UInt64
  UniqueIdGenerator::getUniqueId()
  {
    UniqueIdGenerator& instance = getInstance_();
    return (*instance.dist_)(*instance.rng_);
  }

  UInt
  UniqueIdGenerator::getSeed()
  {
    return getInstance_().seed_;
  }

  void
  UniqueIdGenerator::setSeed(UInt seed)
  {
    UniqueIdGenerator& instance = getInstance_();
    instance.seed_ = seed;
    instance.rng_->seed( instance.seed_ );
    instance.dist_->reset();
  }

  UniqueIdGenerator::UniqueIdGenerator()
  {

  }

  UniqueIdGenerator &
  UniqueIdGenerator::getInstance_()
  {
    if (!instance_)
    {
      instance_ = new UniqueIdGenerator();
      instance_->init_();
    }
    return *instance_;
  }

  void
  UniqueIdGenerator::init_()
  {
    seed_ = std::time(0);
    rng_ = new boost::random::mt19937_64 (seed_);
    dist_ = new boost::random::uniform_int_distribution<UInt64> (0,std::numeric_limits<UInt64>::max());
  }

  UniqueIdGenerator::~UniqueIdGenerator()
  {
    delete rng_;
    delete dist_;
  }

}
