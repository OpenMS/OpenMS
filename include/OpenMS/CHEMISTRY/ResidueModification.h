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
	class ResidueModification
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
			enum Term_Specificity
			{
				ANYWHERE = 0,
				C_TERM = 1,
				N_TERM =2,
				NUMBER_OF_TERM_SPECIFICITY
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
			enum Source_Classification
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
				SYNTH_PEP_PROTECT_GP,
				NUMBER_OF_SOURCE_CLASSIFICATION
			};
			//@}
			


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

			/** @name Assignment operator
			*/
			//@{
			/// assignment operator
			ResidueModification& operator = (const ResidueModification& modification);
			//@}

			/** @name Accessors
			*/
			//@{

			void setId(const String& id);

			const String& getId() const;

			void setFullName(const String& full_name);

			const String& getFullName() const;

			///
			void setName(const String& name);

			/// returns the PSI-MS-label if available; e.g. Mascot uses this name
			const String& getName() const;
			
			void setTermSpecificity(Term_Specificity term_spec);

			void setTermSpecificity(const String& name);
			
			Term_Specificity getTermSpecificity() const;
			
			String getTermSpecitificityName(Term_Specificity = NUMBER_OF_TERM_SPECIFICITY) const;
	
			void setOrigin(const String& origin);

			const String& getOrigin() const;

			/// classification as defined by the PSI-MOD
			void setSourceClassification(const String& classification);

			void setSourceClassification(Source_Classification classification);

			Source_Classification getSourceClassification() const;
			
			/// returns the classification
			String getSourceClassificationName(Source_Classification classification = NUMBER_OF_SOURCE_CLASSIFICATION) const;
			
			void setAverageMass(double mass);

			double getAverageMass() const;

			void setMonoMass(double mass);

			double getMonoMass() const;

			void setDiffAverageMass(double mass);

			double getDiffAverageMass() const;

			void setDiffMonoMass(double mass);

			double getDiffMonoMass() const;
			
			void setFormula(const String& composition);

			const String& getFormula() const;

			void setDiffFormula(const String& diff_formula);

			const String& getDiffFormula() const;
			
			/// sets the valid residue(s)
			//void setValidResidues(const std::vector<String>& residue_names);

			/// returns the valid residues (in PSI-MOD terminology this is called origin)
			//const std::vector<String>& getValidResidues() const;

			void setSynonyms(const std::set<String>& synonyms);

			void addSynonym(const String& synonym);

			const std::set<String>& getSynonyms() const;
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

			String id_;

			String full_name_;

			String name_;
			
			Term_Specificity term_spec_;
			
			String origin_;

			Source_Classification classification_;

			double average_mass_;
			
			double mono_mass_;
			
			double diff_average_mass_;

			double diff_mono_mass_;
			
			String formula_;

			String diff_formula_;

			//std::vector<String> valid_residues_;

			std::set<String> synonyms_;
	};
}

#endif
