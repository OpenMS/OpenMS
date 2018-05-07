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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

/**
    @defgroup Multithreading Multithreading macros

    @brief Macros used for locking in multithreading environment

    These provide a simple interface to acquire locks of different kinds.
    Specifically of interest are read-write locks which allow multiple readers
    but only one writer. A sample implementation will look like this:

    @code

    STATIC_LOCK(SomeSingleton_mutex) 
    class SomeSingleton
    {
      int readData()
      {
        OPENMS_NONUNIQUELOCK(SomeSingleton_mutex, lock)
        // do some work that only reads data
      }
      int writeData()
      {
        OPENMS_UNIQUELOCK(SomeSingleton_mutex, lock)
        // do some work that requires unique access
      }
      int writeDataConditional()
      {
        OPENMS_UPGRADEABLE_UNIQUELOCK(SomeSingleton_mutex, lock)
        // do some work that only reads data
        OPENMS_UPGRADE_UNIQUELOCK(lock, uniqueLock)
        // do some work that requires unique access
      }
    }
    @endcode

    This will allow multiple threads to read simultaneously while only a single
    thread can write at one time. More importantly, while a single thread is
    writing, all other threads are blocked from reading (and the writing thread
    has to wait for all readers to finish before it can start writing). The
    example above thus implements the "single writer - multiple reader"
    paradigm for a singleton resource.

    @note The scope of the lock is controlled by RAII, therefore the lock is
    released as soon as it is out of scope.

    The current implementation uses boost, this can be replaced with C++17
    mutexes at some point.

    @ingroup Concept

    @{
*/


#ifdef _OPENMP
#define OPENMS_MULTITHREADING_ON
#include <boost/thread.hpp>
#include <mutex>
#include <omp.h>
#endif


/// Parameter: mutex vs openmp critical sections
///   replace all mutexes with OpenMP critical sections:
// #define USE_OPENMP_CRITICAL

/// Parameter: std::mutex vs boost::shared_mutex 
///   use std::mutex instead of boostd::shared_mutex (which allows more fine grained locking) 
///   note: if USE_OPENMP_CRITICAL is defined, neither will be selected
// #define USE_STD_MUTEX

#ifdef OPENMS_MULTITHREADING_ON

#ifndef USE_OPENMP_CRITICAL

#ifdef USE_STD_MUTEX

#define STATIC_LOCK(name) \
  static std::mutex name;

#define OPENMS_UNIQUELOCK(name, lockname) \
  std::lock_guard<std::mutex> lockname(name);

#define OPENMS_UPGRADEABLE_UNIQUELOCK(name, lockname) \
  std::lock_guard<std::mutex> lockname(name);

#define OPENMS_UPGRADE_UNIQUELOCK(name, lockname)

#define OPENMS_NONUNIQUELOCK(name, lockname) \
  std::lock_guard<std::mutex> lockname(name);

#else // USE_STD_MUTEX

/**
    @brief Initialize a lock (static).
*/
#define STATIC_LOCK(name) \
  static boost::shared_mutex name;

/**
    @brief Acquire a unique lock.
*/
#define OPENMS_UNIQUELOCK(name, lockname) \
    boost::unique_lock<boost::shared_mutex> lockname{name};

/**
    @brief Acquire a shared lock that can be upgraded to a unique lock.
*/
#define OPENMS_UPGRADEABLE_UNIQUELOCK(name, lockname) \
    boost::upgrade_lock<boost::shared_mutex> lockname{name};

/**
    @brief Upgrade a shared lock to a unique lock.
*/
#define OPENMS_UPGRADE_UNIQUELOCK(name, lockname) \
    boost::upgrade_to_unique_lock<boost::shared_mutex> lockname(name);

/**
    @brief Acquire a shared lock.
*/
#define OPENMS_NONUNIQUELOCK(name, lockname) \
    boost::shared_lock<boost::shared_mutex> lockname{name};
#endif // USE_STD_MUTEX

#else // USE_OPENMP_CRITICAL

#define STATIC_LOCK(name) 

#define OPENMS_UNIQUELOCK(name, lockname) \
  OPENMS_THREAD_CRITICAL(name)

#define OPENMS_UPGRADEABLE_UNIQUELOCK(name, lockname) \
  OPENMS_THREAD_CRITICAL(name)

#define OPENMS_UPGRADE_UNIQUELOCK(name, lockname) \
  OPENMS_THREAD_CRITICAL(name)

#define OPENMS_NONUNIQUELOCK(name, lockname) \
  OPENMS_THREAD_CRITICAL(name)

#endif // USE_OPENMP_CRITICAL


#define STRINGIFY(a) #a

#define OPENMS_THREAD_CRITICAL(name) \
  _Pragma( STRINGIFY( omp critical (name) ) )

#else

// NOP for a single threaded environment

#define STATIC_LOCK(name) 

#define OPENMS_UNIQUELOCK(name) 

#define OPENMS_UPGRADEABLE_UNIQUELOCK(name, lockname) 

#define OPENMS_UPGRADE_UNIQUELOCK(name, lockname) 

#define OPENMS_NONUNIQUELOCK(name) 

#define OPENMS_THREAD_CRITICAL(name)

#endif

/** @} */ // end of Multithreading



