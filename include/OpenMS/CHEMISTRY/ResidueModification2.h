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

			/** @brief Classification of the modification

					the PSI-MOD defined the following classes of modifications
						AA substitution
						Artefact
						Chemical derivative
						Co-translational
						Isotopic label
						Multiple
						N-linked glycosylation
						Non-standard residue
						O-linked glycosylation
						Other
						Other glycosylation
						Post-translational
						Pre-translational
						Synth. pep. protect. gp.
													 
				
			*/
			enum Classification
			{
				AA_SUBSTITUTION = 0,
				ARTEFACT,
				CHEMICAL_DERIVATIVE,
				CO_TRANSLATIONAL,
				ISOTOPIC_LABEL,
				MULTIPLE,
				N_LINKED_GLYCOSYLATION,
				NON_STANDARD_RESIDUE,
				OTHER,
				OTHER_GLYCOSYLATION,
				POST_TRANSLATIONAL,
				PRE_TRANSLATIONAL,
				SYNTH_PEP_PROTECT_GP				
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

			/// sets the valid residue(s) (in PSI-MOD this is referred as origin)
			void setValidResidues(const std::vector<String>& residue_names);

			/// returns the valid residues (in PSI-MOD terminology this is called origin)
			const std::vector<String>& getValidResidues() const;
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

			std::vector<String> valid_residues_;
	};
}

#endif
