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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_CHEMISTRY_MODIFICATIONDEFINITIONSSET_H
#define OPENMS_CHEMISTRY_MODIFICATIONDEFINITIONSSET_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/ModificationDefinition.h>

#include <set>

namespace OpenMS
{
	/** @ingroup Chemistry
	
			@brief Representation of a set of modification definitions

			This class enhances the modification definitions as defined in the 
			class ModificationDefinition into a set of definitions. This is also
			e.g. used as input parameters  in search engines.
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

			/// detailed constructor with comma separated list of modifications
			ModificationDefinitionsSet(const String& fixed_modifications, const String& variable_modifications = "");
			
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

			/// returns the number of fixed modifications stored in this set
			UInt getNumberOfFixedModifications() const;

			/// returns the number of variable modifications stored in this set
			UInt getNumberOfVariableModifications() const;

			/// adds a modification definition to the set
			void addModification(const ModificationDefinition& mod_def);

			/// sets the modification definitions 
			void setModifications(const std::set<ModificationDefinition>& mod_defs);

			/** @brief set the modification definitions from a string
			 
					The strings should contain a comma separated list of modifications. The names
					can be PSI-MOD identifier or any other unique name supported by PSI-MOD. TermSpec
					definitions and other specific definitions are given by the modifications themselves.
			*/
			void setModifications(const String& fixed_modifications, const String& variable_modifications);
			
			/// returns the stored modification definitions
			std::set<ModificationDefinition> getModifications() const;

			/// returns the stored fixed modification definitions
			const std::set<ModificationDefinition>& getFixedModifications() const;

			/// returns the stored variable modification definitions
			const std::set<ModificationDefinition>& getVariableModifications() const;
			
			/// return only the names of the modifications stored in the set
			std::set<String> getModificationNames() const;

			/// return only the names of the fixed modifications
			std::set<String> getFixedModificationNames() const;
			
			/// return only the names of the variable modifications
			std::set<String> getVariableModificationNames() const;
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
			bool operator == (const ModificationDefinitionsSet& rhs) const;

			/// inequality operator
			bool operator != (const ModificationDefinitionsSet& rhs) const;
			//@}


		protected:

			std::set<ModificationDefinition> variable_mods_;

			std::set<ModificationDefinition> fixed_mods_;

			UInt max_mods_per_peptide_;
	};


} // namespace OpenMS

#endif

