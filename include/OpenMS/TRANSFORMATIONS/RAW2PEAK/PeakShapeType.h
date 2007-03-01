// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKSHAPETYPE_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKSHAPETYPE_H

namespace OpenMS
{
  /** 
    @brief Peak shape type (asymmetric lorentzian or asymmetric hyperbolic secans squared).

    The peak shape can represent an asymmetric lorentzian function, given by 
                  
    l(x) = height/(1.+pow(left_width*(x - mz_position), 2)) (x<=mz_position) 
                  
    l(x) = height/(1.+pow(right_width*(x - mz_position), 2)) (x>mz_position)
                  
    or an asymmetric hyperbolic secans squared function 
                  
    s(x) = height/pow(cosh(left_width*(x-mz_position)), 2) (x<=mz_position)
                  
    s(x) = height/pow(cosh(right_width*(x-mz_position)), 2) (x>mz_position)
  */
  struct PeakShapeType
  {
    enum Enum
    {
      LORENTZ_PEAK,
      SECH_PEAK,
      UNDEFINED
    };
  };
} // namespace OpenMS

#endif
