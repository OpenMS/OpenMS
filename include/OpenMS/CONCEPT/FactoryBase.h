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
// $Maintainer: Clemens Groepl, Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_CONCEPT_FACTORYBASE_H
#define OPENMS_CONCEPT_FACTORYBASE_H

#include <OpenMS/config.h>

namespace OpenMS
{
  /** 
  	@brief Base class for Factory<T>
		
		Just be able to use dynamic_cast on a pointer
		
		@ingroup Concept
	*/
  class OPENMS_DLLAPI FactoryBase
  {
		public:
			/// destructor 
			virtual ~FactoryBase(){}
		
  };

}
#endif // OPENMS_CONCEPT_FACTORYBASE_H
