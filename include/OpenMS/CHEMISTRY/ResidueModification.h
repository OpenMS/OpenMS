// -*- mode: C++; tab-width: 2; -*-
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
//

#ifndef OPENMS_CHEMISTRY_RESIDUEMODIFICATION_H
#define OPENMS_CHEMISTRY_RESIDUEMODIFICATION_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

#include <set>

namespace OpenMS
{
	// forward declaration
	class Residue;

	/** @brief Representation of a modification

			This class represents a modification of a residue. A residue modification
			has several attributes like the diff formula, a terminal specificity 
			a mass and maybe an origin which means a specific residue which it can
			be applied to. A residue modification can be represented by its PSI-MOD
			identifier, e.g. MOD:01214. This is a unique key which only occurs ones in 
			an OpenMS instance stored in the ModificationsDB. Some residue modifications
			have also unique synonyms which can be used instead.
	*/
	class OPENMS_DLLAPI ResidueModification
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
				PROTEIN_C_TERM = 3,
				PROTEIN_N_TERM = 4,
				NUMBER_OF_TERM_SPECIFICITY
			};

			/** @brief Classification of the modification
				
			*/
			enum Source_Classification
			{
				ARTIFACT = 0,
				HYPOTHETICAL, 
				NATURAL,
				NUMBER_OF_SOURCE_CLASSIFICATIONS
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
			/// set the identifier of the modification
			void setId(const String& id);

			/// returns the identifier of the modification
			const String& getId() const;

			/// sets the unimod accession
			void setUniModAccession(const String& id);

			/// returns the unimod accession if available
			const String& getUniModAccession() const;

			/// set the MOD:XXXXX accession of PSI-MOD
			void setPSIMODAccession(const String& id);

			/// returns the PSI-MOD accession if available
			const String& getPSIMODAccession() const;

			/// sets the full name of the modification
			void setFullName(const String& full_name);

			/// returns the full name of the modification
			const String& getFullName() const;

			/// sets the name of modification
			void setName(const String& name);

			/// returns the PSI-MS-label if available; e.g. Mascot uses this name
			const String& getName() const;
			
			/// sets the term specificity 
			void setTermSpecificity(Term_Specificity term_spec);

			/// sets the term specificity specified using a name
			void setTermSpecificity(const String& name);
			
			/// returns terminal specificity
			Term_Specificity getTermSpecificity() const;
			
			/// returns the terminal specificity name which is set or given as parameter
			String getTermSpecificityName(Term_Specificity = NUMBER_OF_TERM_SPECIFICITY) const;
	
			///	sets the origin 
			void setOrigin(const String& origin);

			/// returns the origin if set
			const String& getOrigin() const;

			/// classification as defined by the PSI-MOD
			void setSourceClassification(const String& classification);

			/// sets the source classification
			void setSourceClassification(Source_Classification classification);

			/// returns the source classification, if none was set, it is unspecific
			Source_Classification getSourceClassification() const;
			
			/// returns the classification
			String getSourceClassificationName(Source_Classification classification = NUMBER_OF_SOURCE_CLASSIFICATIONS) const;
			
			/// sets the average mass
			void setAverageMass(DoubleReal mass);

			/// returns the average mass if set 
			DoubleReal getAverageMass() const;

			/// sets the monoisotopic mass
			void setMonoMass(DoubleReal mass);

			/// return the monoisotopic mass, if set
			DoubleReal getMonoMass() const;

			/// set the difference average mass
			void setDiffAverageMass(DoubleReal mass);

			/// returns the difference average mass if set
			DoubleReal getDiffAverageMass() const;

			/// sets the difference monoisotopic mass 
			void setDiffMonoMass(DoubleReal mass);

			/// returns the diff monoisotopic mass if set
			DoubleReal getDiffMonoMass() const;
			
			/// set the formula 
			void setFormula(const String& composition);

			/// returns the chemical formula if set
			const String& getFormula() const;

			/// sets diff formula
			void setDiffFormula(const EmpiricalFormula& diff_formula);

			/// returns the diff formula if one was set
			const EmpiricalFormula& getDiffFormula() const;
			
			/// sets the synonyms of that modification
			void setSynonyms(const std::set<String>& synonyms);

			/// adds a synonym to the unique list
			void addSynonym(const String& synonym);

			/// returns the set of synonyms
			const std::set<String>& getSynonyms() const;

			/// sets the neutral loss formula
			void setNeutralLossDiffFormula(const EmpiricalFormula& loss);

			/// returns the neutral loss diff formula (if available)
			const EmpiricalFormula& getNeutralLossDiffFormula() const;

			/// set the neutral loss mono weight 
			void setNeutralLossMonoMass(DoubleReal mono_mass);

			/// returns the neutral loss mono weight
			DoubleReal getNeutralLossMonoMass() const;

			/// set the neutral loss average weight
			void setNeutralLossAverageMass(DoubleReal average_mass);

			/// returns the neutral loss average weight
			DoubleReal getNeutralLossAverageMass() const;
			//@}

			/** @name Predicates
			*/
			//@{
			/// returns true if a neutral loss formula is set
			bool hasNeutralLoss() const;
			
			/// equality operator
			bool operator == (const ResidueModification& modification) const;

			/// inequality operator
			bool operator != (const ResidueModification& modification) const;
			//@}
				
		protected:

			String id_;

			String psi_mod_accession_;

			String unimod_accession_;

			String full_name_;

			String name_;
			
			Term_Specificity term_spec_;
			
			String origin_;

			Source_Classification classification_;

			DoubleReal average_mass_;
			
			DoubleReal mono_mass_;
			
			DoubleReal diff_average_mass_;

			DoubleReal diff_mono_mass_;
			
			String formula_;

			EmpiricalFormula diff_formula_;

			std::set<String> synonyms_;

			EmpiricalFormula neutral_loss_diff_formula_;

			DoubleReal neutral_loss_mono_mass_;

			DoubleReal neutral_loss_average_mass_;
	};
}

#endif
