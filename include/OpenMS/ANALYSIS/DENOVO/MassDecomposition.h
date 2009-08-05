// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYIS_DENOVO_DECOMPOSITION_H
#define OPENMS_ANALYIS_DENOVO_DECOMPOSITION_H

#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS
{
	/** @brief Class represents a decomposition of a mass into amino acids
			
			This class represents a mass decomposition into amino acids. A
			decomposition are amino acids given with frequencies which add
			up to a specific mass.
	*/
	class OPENMS_DLLAPI MassDecomposition
	{
		public:

			/** @name Constructors and destructors
			*/
			//@{
			/// default constructor
			MassDecomposition();

			/// copy constructor
			MassDecomposition(const MassDecomposition& deco);
		
			/// constructor with String as parameter
			MassDecomposition(const String& deco);
			//@}

			/** @name Operators
			*/
			//@{
			/// assignment operator
			MassDecomposition& operator = (const MassDecomposition& rhs);

			MassDecomposition& operator += (const MassDecomposition& d);

			String toString() const;

			String toExpandedString() const;
			
			MassDecomposition operator + (const MassDecomposition& rhs) const;
			
			UInt getNumberOfMaxAA() const;
			//@}

			/** @name Predicates
			*/
			//@{
			/// less than predicate
			bool operator < (const MassDecomposition& rhs) const;

			/// equality operator 
			bool operator == (const String& deco) const;
	
			/// returns true if tag is contained in the mass decomposition
			bool containsTag(const String& tag) const;

			/// returns true if the mass decomposition if contained in this instance
			bool compatible(const MassDecomposition& deco) const;
			//@}

		protected:

			Map<char, UInt> decomp;
			
			UInt number_of;
			
			UInt number_of_max_aa;

			UInt min_number;
			
			UInt max_number;

	};

} // namespace OpenMS

#endif
