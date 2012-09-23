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

#ifndef OPENMS_CONCEPT_UNIQUEIDGENERATOR_H
#define OPENMS_CONCEPT_UNIQUEIDGENERATOR_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/DATASTRUCTURES/Param.h>

#include <iostream>

namespace OpenMS
{

/**
 @brief  A generator for unique ids.

 The unique ids are 64-bit random unsigned random integers.
 The class is implemented as a singleton.
 The random generator is initialized upon startup using the current system time and date.
 Collisions are not excluded by design, but extremely unlikely.
 (To estimate the probability of collisions,
 note that \f$ 10^9*60*60*24*365*100 / 2^{64} \doteq 0.17 \f$,
 so it is unlikely you will see one in your lifetime.)

 @ingroup Concept
 */
  class OPENMS_DLLAPI UniqueIdGenerator
  {

public:

    /// Returns a new unique id
    static UInt64
    getUniqueId();

    /// Initializes random generator using the given DateTime instead of DateTime::now().  This is intended for debugging and testing.
    static void
    setSeed(const DateTime &);

    /// Returns a summary of internal status
    static Param const &
    getInfo();

private:

    UniqueIdGenerator();

    ~UniqueIdGenerator();

    static UniqueIdGenerator &
    getInstance_();

    void
    init_(const DateTime & date_time);

    Param info_;

  };

} // namespace OpenMS

#endif  // OPENMS_CONCEPT_UNIQUEIDGENERATOR_H
