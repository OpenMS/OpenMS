// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
    setSeed(const DateTime&);

    /// Returns a summary of internal status
    static Param const&
    getInfo();

  private:

    UniqueIdGenerator();

    ~UniqueIdGenerator();

    static UniqueIdGenerator&
    getInstance_();

    void
    init_(const DateTime& date_time);

    Param info_;

};

} // namespace OpenMS

#endif  // OPENMS_CONCEPT_UNIQUEIDGENERATOR_H
