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

#ifndef OPENMS_CHEMISTRY_RESIDUEMODIFICATION2_H
#define OPENMS_CHEMISTRY_RESIDUEMODIFICATION2_H

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
	class ResidueModification2
	{
		public:
	
			/** Enums
			*/
			//@{
			/** @brief Position where the modification is allowed to occur
			
					The allowed sites are
					Any C-term
					Any N-term
					Anywhere
					Protein C-term
					Protein N-term

					This does not describe the amino acids which are valid for a 
					specific amino acid!

			*/
			enum AllowedPosition
			{
				ANY_C_TERM = 0,
				ANY_N_TERM = 1,
				ANYWHERE = 2,
				PROTEIN_C_TERM = 3,
				PROTEIN_N_TERM
			};
			//@}
						
			/** @name Constructors and Destructors
			*/
			//@{
			/// default constructor
			ResidueModification2();
			
			/// copy constructor
			ResidueModification2(const ResidueModification2& modification);

			/// destructor				
			virtual ~ResidueModification2();
			//@}

			/** @name Assignment operator
			*/
			//@{
			/// assignment operator
			ResidueModification2& operator = (const ResidueModification2& modification);
			//@}

			/** @name Accessors
			*/
			//@{

			void setTitle(const String& title);

			const String& getTitle() const;

			void setFullName(const String& full_name);

			const String& getFullName() const;

			void setAllowedPosition(AllowedPosition position);

			AllowedPosition getAllowedPosition() const;
			
			String getAllowedPositionName() const;
			
			/// site where the modification is allowed to occur; amino acids C-term and N-term are valid
			void setSite(const String& site);

			const String& getSite() const;

			/// classification as defined by the PSI-MOD
			void setClassification(const String& classification);

			/// returns the classification
			const String& getClassification() const;
			
			void setAverageMass(double mass);

			double getAverageMass() const;

			void setMonoMass(double mass);

			double getMonoMass() const;

			void setComposition(const String& composition);

			const String& getComposition() const;

			
			/*
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
			void setAddAverageWeight(DoubleReal weight);
				
			/// returns the weight of the added formula
			DoubleReal getAddAverageWeight() const;

			/// sets the mono isotopic weight of the added formula
			void setAddMonoWeight(DoubleReal weight);

			/// returns the mono isotopic weight of the added formula
			DoubleReal getAddMonoWeight() const;

			/// sets the formula which is deleted from the residue
			void setDelFormula(const EmpiricalFormula& formula);
				
			/// returns the formula which is deleted from the residue
			const EmpiricalFormula& getDelFormula() const;

			/// sets the average weight of the deletion 
			void setDelAverageWeight(DoubleReal weight);

			/// returns the average weight of the deletion
			DoubleReal getDelAverageWeight() const;

			/// sets the mono isotopic weight of the deletion
			void setDelMonoWeight(DoubleReal weight);

			/// returns the mono isotopic weight of the deletion
			DoubleReal getDelMonoWeight() const;

			/// sets the residues where the modification can be applied to
			void setValidResidues(const std::set<Residue*>& valid_residues);

			/// adds a valid residue
			void addValidResidue(Residue* valid_residue);
				
			/// returns the residues where the modifications can be applied to
			const std::set<Residue*>& getValidResidues() const;
			*/
			//@}

			/** @name Predicates
			*/
			//@{
			/// equality operator
			bool operator == (const ResidueModification2& modification) const;

			/// inequality operator
			bool operator != (const ResidueModification2& modification) const;
			//@}
				
		protected:

			String title_;

			String full_name_;

			AllowedPosition allowed_position_;
			
			String site_;

			String classification_;

			double average_mass_;
			
			double mono_mass_;
			
			String composition_;

			/*
			// basic
			String name_;
			
			String short_name_;
			
			String name_prefix_;

			std::set<String> synonyms_;

			// additions
			EmpiricalFormula add_formula_;

			DoubleReal add_average_weight_;

			DoubleReal add_mono_weight_;
				
			// deletions
			EmpiricalFormula del_formula_;

			DoubleReal del_average_weight_;

			DoubleReal del_mono_weight_;

			// residues 
			std::set<Residue*> valid_residues_;
			*/
	};
}

#endif
