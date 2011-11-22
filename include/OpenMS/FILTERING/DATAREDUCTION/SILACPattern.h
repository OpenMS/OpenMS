// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:expandtab
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2011 -- Bastian Blank
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
// $Maintainer: Lars Nilse $
// $Authors: Bastian Blank $
// --------------------------------------------------------------------------

#include <OpenMS/FILTERING/DATAREDUCTION/SILACPoint.h>

#ifndef OPENMS_FILTERING_DATAREDUCTION_SILACPATTERN_H
#define OPENMS_FILTERING_DATAREDUCTION_SILACPATTERN_H

namespace OpenMS
{
  /**
   * @brief A single SILAC pattern containing multiple found points
   * @see HashGrid
   * @ingroup Datastructures
   */
  class OPENMS_DLLAPI SILACPattern
    : public SILACPoint
  {
    public:
      /**
       * Points checked and found in the raw data.
       */
      std::vector<SILACPoint> points;
  };
}

#endif /* OPENMS_FILTERING_DATAREDUCTION_SILACPATTERN_H */
