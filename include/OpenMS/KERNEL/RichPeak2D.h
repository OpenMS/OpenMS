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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_RICHPEAK2D_H
#define OPENMS_KERNEL_RICHPEAK2D_H

#include <OpenMS/KERNEL/Peak2D.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

namespace OpenMS
{

	/**	
		@brief A 2-dimensional raw data point or peak mith meta information.
	 
		This datastructure is intended for continuous data or peak data.
		If wou do not need to annotated single peaks with meta data, use Peak2D instead.
	
		@ingroup Kernel
	*/
	class RichPeak2D 
		: public Peak2D, 
			public MetaInfoInterface
	{
		public:
		
		/// Default constructor
		RichPeak2D() 
			: Peak2D(),
				MetaInfoInterface()
		{
			
		}
		
		/// Copy constructor
		inline RichPeak2D(const RichPeak2D& p) 
			: Peak2D(p),
				MetaInfoInterface(p)
		{
		}

		/// Destructor
		~RichPeak2D() 
		{
		}
    
		/// Assignment operator
		inline RichPeak2D& operator = (const RichPeak2D& rhs)
		{
			if (this==&rhs) return *this;
			
			Peak2D::operator = (rhs);
			MetaInfoInterface::operator = (rhs);
			
			return *this;
		}

		/// Equality operator
		inline bool operator == (const RichPeak2D& rhs) const
		{
			return 
				Peak2D::operator == (rhs) &&
				MetaInfoInterface::operator == (rhs)
			;
		}

		/// Equality operator
		inline bool operator != (const RichPeak2D& rhs) const
		{
			return !(operator == (rhs));
		}
	
	};

} // namespace OpenMS

#endif // OPENMS_KERNEL_RICHPEAK2D_H
