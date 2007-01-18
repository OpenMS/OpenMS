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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------


#ifndef OPENMS_KERNEL_DIMENSIONDESCRIPTION_H
#define OPENMS_KERNEL_DIMENSIONDESCRIPTION_H

namespace OpenMS
{

  //------------------------------------------------------------
  
  // You can't instantiate this in general. Only some template specializations are provided for tag classes.
  template < typename /* DimensionDescriptionTagClass */ > struct DimensionDescription;

  //------------------------------------------------------------

  //Tag class for LCMS experiments (it is never instanciated)
  struct LCMS_Tag;

  /**
  	@brief Dimension description for LCMS experiments
  	
  	@ingroup Kernel
  */
  template <>
  struct DimensionDescription < LCMS_Tag >
  {
    /// This maps symbolic names of the dimensions to numbers.
    enum DimensionId
      {
				
				RT, ///< Mass-to-charge dimension id (0 if used as a const int)
				MZ, ///< Retention time dimension id (1 if used as a const int)
				
				
				/** 
					This is the last value in the @p enum.  It is used as size
					parameter of the arrays, so that we will get a compile time error
					if a new id is added and the corresponding names are not provided
					in the static intitializers.
				*/
				DIMENSION 
      };
		
    /**
    	@brief Short name of the dimension (abbreviated form).
    	
    	By convention, it should be the same as the identifier in the @p enum,
			e.g. <code>dimension_name_short[MZ] == "MZ"</code>, etc.
		*/
    static char const * const dimension_name_short [DIMENSION];
		
    /// Long name of the dimension (self-explanatory form)
    static char const * const dimension_name_full  [DIMENSION];

    /// Unit of measurement (abbreviated form)
    static char const * const dimension_unit_short [DIMENSION];

    /// Unit of measurement (self-explanatory form)
    static char const * const dimension_unit_full  [DIMENSION];
    
  };
  
} // namespace OpenMS

#endif // OPENMS_KERNEL_DIMENSIONDESCRIPTION_H
