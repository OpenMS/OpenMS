// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_CHEMISTRY_MODIFICATIONDEFINITION_H
#define OPENMS_CHEMISTRY_MODIFICATIONDEFINITION_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>

namespace OpenMS
{
	/** @ingroup Chemistry
	
			@brief Representation of modification definition

			This class defines a modification type e.g. a input parameter of a search engine.
			The modification is defined using an unique name of the modification present 
			in the modifications DB instance.
	*/
	class OPENMS_DLLAPI ModificationDefinition
	{
		public:

			/** @name Constructor and Destructors
			*/
			//@{
			/// default constructor
			ModificationDefinition();

			/// copy constructor
			ModificationDefinition(const ModificationDefinition& rhs);

			/// detailed constructor specifying the modifications name
			ModificationDefinition(const String& mod);
			
			/// destructor
			virtual ~ModificationDefinition();
			//@}

			/** @name Accessors
			*/
			//@{
			/// sets the allowed position of the modification
			void setTermSpecificity(ResidueModification::Term_Specificity pos);

			/// returns the allowed position of the modification
			ResidueModification::Term_Specificity getTermSpecificity() const;

			/// sets whether this modification definition is fixed or variable (modification must occur vs. can occur)
			void setFixedModification(bool fixed);

			/// returns if the modification if fixed true, else false
			bool isFixedModification() const;

			/// set the maximal number of occurences per peptide, unbound if 0
			void setMaxOccurences(UInt num);

			/// returns the maximal number of occurences per peptide
			UInt getMaxOccurences() const;

			/// returns the modification set
			String getModification() const;

			/// sets the modification, allowed are unique names provided by ModificationsDB
			void setModification(const String& modification);
			//@}

			/** @name Assignment
			*/
			//@{
			/// assignment operator
			ModificationDefinition& operator = (const ModificationDefinition& element);
			//@}

			/** @name Predicates
			*/
			//@{
			/// equality operator
			bool operator == (const ModificationDefinition& rhs) const;

			/// inequality operator
			bool operator != (const ModificationDefinition& rhs) const;

			/// less than operator for e.g. usage in maps; only mod FullIds are compared!
			bool operator < (const OpenMS::ModificationDefinition&) const;
			//@}

		protected:

			/// allowed position
			ResidueModification::Term_Specificity term_spec_;

			/// the modification
			const ResidueModification* mod_;

			/// fixed (true) or variable (false)
			bool fixed_modification_;

			/// maximal number of occurences per peptide
			UInt max_occurences_;
	};

} // namespace OpenMS

#endif

