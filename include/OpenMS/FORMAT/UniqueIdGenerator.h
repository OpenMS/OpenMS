// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Knut Reinert $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_UNIQUEIDGENERATOR_H
#define OPENMS_FORMAT_UNIQUEIDGENERATOR_H

#include <OpenMS/CONCEPT/Types.h>

namespace OpenMS
{
  /**
    @brief Generator for unique IDs.
    
    This class is a singleton, so you have to access it via the method @p instance().
    
  	@ingroup Format
  */
	class UniqueIdGenerator
	{
		public:
			/// returns a regerence to the instance
			static UniqueIdGenerator& instance();
			/// returns a unique ID
			UID getUID();
			
		private:
			/// the actual unique ID
			static UID id_;
			/// Default constructor
			UniqueIdGenerator();
			
	};
} // namespace OpenMS

#endif // OPENMS_FOMAT_UNIQUEIDGENERATOR_H
