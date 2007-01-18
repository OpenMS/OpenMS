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

#ifndef OPENMS_CHEMISTRY_RESIDUE_H
#define OPENMS_CHEMISTRY_RESIDUE_H

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <iostream>
#include <set>
#include <vector>

namespace OpenMS
{
	/** 
		@ingroup Chemistry
		
		@brief Representation of a residue
	*/
	class Residue
	{
		public:
		
		/** @brief Representation of a modification
		*/
		class Modification
		{
			public:
			
				/** @name Constructors and Destructors
				*/
				//@{
				/// default constructor
				Modification();
				
				/// copy constructor
				Modification(const Modification& modification);

				/// destructor				
				virtual ~Modification();
				//@}

				/** @name Assignment
				*/
				//@{
				/// assignment operator
				Modification& operator = (const Modification& modification);
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
				bool operator == (const Modification& modification) const;

				/// inequality operator
				bool operator != (const Modification& modification) const;
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
		
		
			/** @name Typedefs and Constants
			*/
			//@{
			typedef std::vector<EmpiricalFormula>::const_iterator LowMassConstIterator;
			typedef std::vector<EmpiricalFormula>::iterator LowMassIterator;

			
			inline static const EmpiricalFormula& getInternalToFull()
			{
				static const EmpiricalFormula internal_to_full = EmpiricalFormula("H2O");
				return internal_to_full;
			}

			
			inline static Real getInternalToFullAverageWeight()
			{
				static const Real internal_to_full_average_weight = getInternalToFull().getAverageWeight();
				return internal_to_full_average_weight;
			}

			inline static Real getInternalToFullMonoWeight()
			{
				static const Real internal_to_full_mono_weight = getInternalToFull().getMonoWeight();
				return internal_to_full_mono_weight;
			}

			inline static const EmpiricalFormula& getNTerminalToFull()
			{
				static const EmpiricalFormula Nterminal_to_full = EmpiricalFormula("HO");
				return Nterminal_to_full;
			}

			inline static Real getNTerminalToFullAverageWeight()
			{
				static const Real Nterminal_to_full_average_weight = getNTerminalToFull().getAverageWeight();
				return Nterminal_to_full_average_weight;
			}

			inline static Real getNTerminalToFullMonoWeight()
			{
				static const Real Nterminal_to_full_mono_weight = getNTerminalToFull().getMonoWeight();
				return Nterminal_to_full_mono_weight;
			}
			
			inline static const EmpiricalFormula& getCTerminalToFull()
			{
				static const EmpiricalFormula Cterminal_to_full = EmpiricalFormula("H");
				return Cterminal_to_full;
			}

			inline static Real getCTerminalToFullAverageWeight()
			{
				static const Real Cterminal_to_full_average_weight = getCTerminalToFull().getAverageWeight();
				return Cterminal_to_full_average_weight;
			}
		
			inline static Real getCTerminalToFullMonoWeight()
			{
				static const Real Cterminal_to_full_mono_weight = getCTerminalToFull().getMonoWeight();
				return Cterminal_to_full_mono_weight;
			}
		
			inline static const EmpiricalFormula& getBIonToFull()
			{
				static const EmpiricalFormula b_ion_to_full = EmpiricalFormula("HO");
				return b_ion_to_full;
			}

			inline static Real getBIonToFullAverageWeight()
			{
				static const Real b_ion_to_full_average_weight = getBIonToFull().getAverageWeight();
				return b_ion_to_full_average_weight;
			}
			
			inline static Real getBIonToFullMonoWeight()
			{
				static const Real b_ion_to_full_mono_weight = getBIonToFull().getMonoWeight();
				return b_ion_to_full_mono_weight;
			}
			
			inline static const EmpiricalFormula& getAIonToFull()
			{
				static const EmpiricalFormula a_ion_to_full = EmpiricalFormula("HCO2");
				return a_ion_to_full;
			}
	
			inline static Real getAIonToFullAverageWeight()
			{
				static const Real a_ion_to_full_average_weight = getAIonToFull().getAverageWeight();
				return a_ion_to_full_average_weight;
			}

			inline static Real getAIonToFullMonoWeight()
			{
				static const Real a_ion_to_full_mono_weight = getAIonToFull().getMonoWeight();
				return a_ion_to_full_mono_weight;
			}

			inline static const EmpiricalFormula& getYIonToFull()
			{
				static const EmpiricalFormula y_ion_to_full = EmpiricalFormula("");
				return y_ion_to_full;
			}

			inline static Real getYIonToFullAverageWeight()
			{
				static const Real y_ion_to_full_average_weight = getYIonToFull().getAverageWeight();
				return y_ion_to_full_average_weight;
			}

			inline static Real getYIonToFullMonoWeight()
			{
				static const Real y_ion_to_full_mono_weight = getYIonToFull().getMonoWeight();
				return y_ion_to_full_mono_weight;
			}

			inline static const EmpiricalFormula& getCIonToFull()
			{
				static const EmpiricalFormula c_ion_to_full = EmpiricalFormula("");
				return c_ion_to_full;
			}

			inline static Real getCIonToFullAverageWeight()
			{
				static const Real c_ion_to_full_average_weight = getCIonToFull().getAverageWeight();
				return c_ion_to_full_average_weight;
			}

			inline static Real getCIonToFullMonoWeight()
			{
				static const Real c_ion_to_full_mono_weight = getCIonToFull().getMonoWeight();
				return c_ion_to_full_mono_weight;
			}
			
			inline static const EmpiricalFormula& getXIonToFull()
			{
				static const EmpiricalFormula x_ion_to_full = EmpiricalFormula("HCO");
				return x_ion_to_full;
			}

			inline static Real getXIonToFullAverageWeight()
			{
				static const Real x_ion_to_full_average_weight = getXIonToFull().getAverageWeight();
				return x_ion_to_full_average_weight;
			}

			inline static Real getXIonToFullMonoWeight()
			{
				static const Real x_ion_to_full_mono_weight = getXIonToFull().getMonoWeight();
				return x_ion_to_full_mono_weight;
			}

			inline static const EmpiricalFormula& getZIonToFull()
			{
				static const EmpiricalFormula z_ion_to_full = EmpiricalFormula("N");
				return z_ion_to_full;
			}

			inline static Real getZIonToFullAverageWeight()
			{
				static const Real z_ion_to_full_average_weight = getZIonToFull().getAverageWeight();
				return z_ion_to_full_average_weight;
			}

			inline static Real getZIonToFullMonoWeight()
			{
				static const Real z_ion_to_full_mono_weight = getZIonToFull().getMonoWeight();
				return z_ion_to_full_mono_weight;
			}
			//@}
			
			/** @name Enums
			*/
			//@{
			enum ResidueType
			{
				Full = 0,
				Internal,
				NTerminal,
				CTerminal,
				AIon,
				BIon,
				CIon,
				XIon,
				YIon,
				ZIon
			};
			//@}
			
			/** @name Constructors
			*/
			//@{
			/// default contructor
			Residue();
	
			/// copy constructor
			Residue(const Residue& residue);
	
			/// detailed constructor
			Residue(const String& name,
							const String& three_letter_code,
							const String& one_letter_code,
							const EmpiricalFormula& formula,
							const EmpiricalFormula& neutral_loss
							);
	
			/// destructor
			virtual ~Residue();
			//@}
			
			/** @name Assignment
			 */
			//@{
			/// assignment operator
			Residue& operator = (const Residue& residue);
			//@}
	
			/** Accessors
			*/
			//@{
			/// sets the name of the residue
			void setName(const String& name);
			
			/// returns the name of the residue
			const String& getName() const;

			/// sets the short name of the residue, this name is used in the PeptideSequence for output
			void setShortName(const String& short_name);

			/// returns the short name of the residue
			const String& getShortName() const;
	
			/// sets the synonyms
			void setSynonyms(const std::set<String>& synonyms);

			/// adds a synonym
			void addSynonym(const String& synonym);
			
			/// returns the sysnonyms
			const std::set<String>& getSynonyms() const;
			
			/// sets the name of the residue as three letter code
			void setThreeLetterCode(const String& three_letter_code);
	
			/// returns the name of the residue as three letter code
			const String& getThreeLetterCode() const;
	
			/// sets the name as one letter code
			void setOneLetterCode(const String& one_letter_code);
	
			/// returns the name as one letter code
			const String& getOneLetterCode() const;
	
			/// sets the neutral loss formula (if there is one)
			void setLossFormula(const EmpiricalFormula&);
	
			/// returns the neutral loss formula (or just an empty string if there is none)
			const EmpiricalFormula& getLossFormula() const;
	
			/// sets the average weight of the neutral loss molecule
			void setLossAverageWeight(Real weight);
	
			/// return the average weight of the neutral loss molecule
			Real getLossAverageWeight() const;
	
			/// sets the mono isotopic weight of the neutral loss molecule
			void setLossMonoWeight(Real weight);
	
			/// returns the mono isotopic weight of the neutral loss molecule
			Real getLossMonoWeight() const;
	
			/// set the neutral loss molecule weight (if there is one)
			void setLossName(const String& name);
			
			/// gets neutral loss name (if there is one, else returns an empty string)
			const String& getLossName() const;
			
			/// set empirical formula of the residue
			void setFormula(const EmpiricalFormula& formula, ResidueType res_type = Full);
		
			/// returns the empirical formula of thre residue
			EmpiricalFormula getFormula(ResidueType res_type = Full) const;
	
			/// sets average weight of the residue
			void setAverageWeight(Real weight, ResidueType res_type = Full);
	
			/// returns average weight of the residue
			Real getAverageWeight(ResidueType res_type = Full) const;
	
			/// sets mono weight of the residue
			void setMonoWeight(Real weight, ResidueType res_type = Full);
	
			/// returns mono weight of the residue
			Real getMonoWeight(ResidueType res_type = Full) const;

			/// sets the modification pointer
			void setModification(Modification* modification);
			
			/// returns a pointer to the modification, if the residue is unmodified, 0 is returned
			const Modification * getModification() const;
			
			/// sets the name of the unmodified residue
			void setUnmodifiedName(const String& name);

			/// returns the name of the unmodified residue
			const String& getUnmodifiedName() const;

			/// sets the low mass marker ions as a vector of formulas
			void setLowMassIons(const std::vector<EmpiricalFormula>& low_mass_ions);

			/// returns a vector of formulas with the low mass markers of the residue
			const std::vector<EmpiricalFormula>& getLowMassIons() const;
			//@}
			
			/** @name Predicates
			*/
			//@{
			/// true if the residue has neutral loss
			bool hasNeutralLoss() const;
			
			/// equality operator
			bool operator == (const Residue& residue) const;

			/// inequality operator
			bool operator != (const Residue& residue) const;

			/// equality operator for one letter code
			bool operator == (char one_letter_code) const;

			/// equality operator for one letter code
			bool operator != (char one_letter_code) const;

			/// returns the pka of the residue 
			DoubleReal getPka() const;
		
			/// returns the pkb of the residue 
			DoubleReal getPkb() const;
		
			/// returns the pkc of the residue if it exists otherwise -1
			DoubleReal getPkc() const;
		
			/// calculates the isoelectric point using the pk* values
			DoubleReal getPiValue() const;

			/// sets the pka of the residue
			void setPka(DoubleReal value);
		
			/// sets the pkb of the residue
			void setPkb(DoubleReal value);		
		
			/// sets the pkc of the residue
			void setPkc(DoubleReal value);		
		
			/// true if the residue is a modified one
			bool isModified() const;
			//@}
		
			/// ostream iterator to write the residue to a stream
			friend std::ostream& operator << (std::ostream& os, const Residue& residue);

		protected:

			// basic 
			String name_;

			String short_name_;
			
			std::set<String> synonyms_;
			
			String three_letter_code_;

			String one_letter_code_;

			EmpiricalFormula formula_;

			EmpiricalFormula internal_formula_;

			Real average_weight_;

			Real mono_weight_;

			// modification
			bool is_modified_;

			String pre_mod_name_;

			Modification * modification_;
			
			// loss
			String loss_name_;

			EmpiricalFormula loss_formula_;

			Real loss_average_weight_;

			Real loss_mono_weight_;

			// low mass markers like immonium ions
			std::vector<EmpiricalFormula> low_mass_ions_;
				
			// pka values
			DoubleReal pka_;

			// pkb values
			DoubleReal pkb_;
			
			// pkc values
			DoubleReal pkc_;
			
	};
	
	std::ostream& operator << (std::ostream& os, const Residue& residue);

	/// returns the ion name given as a residue type
	String getResidueTypeName(Residue::ResidueType res_type);
}

#endif
