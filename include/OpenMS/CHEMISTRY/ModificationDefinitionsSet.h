// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_CHEMISTRY_RESIDUEMODIFICATIONSSET_H
#define OPENMS_CHEMISTRY_RESIDUEMODIFICATIONSSET_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/ModificationDefinition.h>

namespace OpenMS
{
	/** @ingroup Chemistry
	
			@brief Representation of an element
	*/
	class ModificationDefinitionsSet
	{
		public:

			/** @name Constructor and Destructors
			*/
			//@{
			/// default constructor
			ModificationDefinitionsSet();

			/// copy constructor
			ModificationDefinitionsSet(const ModificationDefinitionsSet& rhs);

			/// destructor
			virtual ~ModificationDefinitionsSet();
			//@}

			/** @name Accessors
			*/
			//@{
			/// sets the maximal number of modifications allowed per peptide
			void setMaxModifications(UInt max_mod);

			/// return the maximal number of modifications allowed per peptide
			UInt getMaxModifications() const;

			/// returns the number of modifications stored in this set
			UInt getNumberOfModifications() const;

			UInt getNumberOfFixedModifications() const;

			UInt getNumberOfVariableModifications() const;

			void addModificationDefinition(const ModificationDefinition& mod_def);

			void setModificationDefinitions(const std::vector<ModificationDefinition>& mod_defs);
			//@}

			/** @name Assignment
			*/
			//@{
			/// assignment operator
			ModificationDefinitionsSet& operator = (const ModificationDefinitionsSet& element);
			//@}

			/** @name Predicates
			*/
			//@{
			/// equality operator
			bool operator == (const ModificationDefinitionsSet& element) const;

			/// inequality operator
			bool operator != (const ModificationDefinitionsSet& element) const;
			//@}


		protected:

			std::vector<ModificationDefinition> variable_mods_;

			std::vector<ModificationDefinition> fixed_mods_;

			UInt max_mods_per_peptide_;
	};


} // namespace OpenMS

#endif

