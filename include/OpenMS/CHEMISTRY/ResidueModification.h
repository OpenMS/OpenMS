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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_CHEMISTRY_RESIDUEMODIFICATION_H
#define OPENMS_CHEMISTRY_RESIDUEMODIFICATION_H

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <set>

namespace OpenMS
{
	// forward declaration
	class Residue;

	/** @brief Representation of a modification

	*/
	class ResidueModification
	{
		public:
		
			/** @name Constructors and Destructors
			*/
			//@{
			/// default constructor
			ResidueModification();
			
			/// copy constructor
			ResidueModification(const ResidueModification& modification);

			/// destructor				
			virtual ~ResidueModification();
			//@}
			/** @name Assignment
			*/
			//@{
			/// assignment operator
			ResidueModification& operator = (const ResidueModification& modification);
			//@}

			/** @name Accessors
			*/
			//@{
			/// sets the name of the modification
			void setName(const String& name);

			/// returns the name of the modification
			const String& getName() const;

			/// sets the short name of the modification, this name is used in PeptideSequence as output
			void setShortName(const String& name);

			/// returns the short name of the modification
			const String& getShortName() const;

			/// sets the naming prefix of modified residues
			void setNamePrefix(const String& name_prefix);

			/// returns the naming prefix of modified residues
			const String& getNamePrefix() const;

			/// sets the synonyms of the modification
			void setSynonyms(const std::set<String>& synonyms);

			/// adds a synonym to the modification
			void addSynonym(const String& synonym);

			/// returns the synonym names of the modifications
			const std::set<String>& getSynonyms() const;
				
			/// sets the formula, which is added to the original residue
			void setAddFormula(const EmpiricalFormula& formula);
				
			/// returns the formula, which is added to the original residue
			const EmpiricalFormula& getAddFormula() const;

			/// sets the average weight of the added formula
			void setAddAverageWeight(Real weight);
				
			/// returns the weight of the added formula
			Real getAddAverageWeight() const;

			/// sets the mono isotopic weight of the added formula
			void setAddMonoWeight(Real weight);

			/// returns the mono isotopic weight of the added formula
			Real getAddMonoWeight() const;

			/// sets the formula which is deleted from the residue
			void setDelFormula(const EmpiricalFormula& formula);
				
			/// returns the formula which is deleted from the residue
			const EmpiricalFormula& getDelFormula() const;

			/// sets the average weight of the deletion 
			void setDelAverageWeight(Real weight);

			/// returns the average weight of the deletion
			Real getDelAverageWeight() const;

			/// sets the mono isotopic weight of the deletion
			void setDelMonoWeight(Real weight);

			/// returns the mono isotopic weight of the deletion
			Real getDelMonoWeight() const;

			/// sets the residues where the modification can be applied to
			void setValidResidues(const std::set<Residue*>& valid_residues);

			/// adds a valid residue
			void addValidResidue(Residue* valid_residue);
				
			/// returns the residues where the modifications can be applied to
			const std::set<Residue*>& getValidResidues() const;
			//@}

			/** @name Predicates
			*/
			//@{
			/// equality operator
			bool operator == (const ResidueModification& modification) const;

			/// inequality operator
			bool operator != (const ResidueModification& modification) const;
			//@}
				
		protected:

			// basic
			String name_;
			
			String short_name_;
			
			String name_prefix_;

			std::set<String> synonyms_;

			// additions
			EmpiricalFormula add_formula_;

			Real add_average_weight_;

			Real add_mono_weight_;
				
			// deletions
			EmpiricalFormula del_formula_;

			Real del_average_weight_;

			Real del_mono_weight_;

			// residues 
			std::set<Residue*> valid_residues_;
	};
}

#endif
