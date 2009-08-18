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

#ifndef OPENMS_ANALYSIS_DENOVO_MASSDECOMPOSITION_H
#define OPENMS_ANALYSIS_DENOVO_MASSDECOMPOSITION_H

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

			/** @name Operators and accessors
			*/
			//@{
			/// assignment operator
			MassDecomposition& operator = (const MassDecomposition& rhs);

			/// adds the mass decomposition d to this object
			MassDecomposition& operator += (const MassDecomposition& d);

			/// returns the decomposition as a string
			String toString() const;

			/// returns the decomposition as a string; instead of frequencies the amino acids are repeated
			String toExpandedString() const;
			
			/// adds this decomposition and the decomposition given and returns a new composition
			MassDecomposition operator + (const MassDecomposition& rhs) const;
			
			/// returns the max frequency of this composition
			Size getNumberOfMaxAA() const;
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

			Map<char, Size> decomp_;
			
			Size number_of_max_aa_;

	};

} // namespace OpenMS

#endif
