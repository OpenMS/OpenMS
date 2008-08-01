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
	// forward declarations
	class ResidueModification;

	/** 
		@ingroup Chemistry
		
		@brief Representation of a residue

		This class represents residues. Residues can have many different attributes, like
		the formula physico-chemical values of properties and so on. 

		A very important property of residues are their modifications. By default no 
		modification is present. Any modification which is present in the ModificationsDB can 
		by applied, if appropriate. 
	*/
	class Residue
	{
		public:
		
			/** @name Typedefs and Constants
			*/
			//@{
			//typedef std::vector<EmpiricalFormula>::const_iterator LowMassConstIterator;
			//typedef std::vector<EmpiricalFormula>::iterator LowMassIterator;

			
			inline static const EmpiricalFormula& getInternalToFull()
			{
				static const EmpiricalFormula internal_to_full = EmpiricalFormula("H2O");
				return internal_to_full;
			}

			
			inline static DoubleReal getInternalToFullAverageWeight()
			{
				static const DoubleReal internal_to_full_average_weight = getInternalToFull().getAverageWeight();
				return internal_to_full_average_weight;
			}

			inline static DoubleReal getInternalToFullMonoWeight()
			{
				static const DoubleReal internal_to_full_mono_weight = getInternalToFull().getMonoWeight();
				return internal_to_full_mono_weight;
			}

			inline static const EmpiricalFormula& getNTerminalToFull()
			{
				static const EmpiricalFormula Nterminal_to_full = EmpiricalFormula("HO");
				return Nterminal_to_full;
			}

			inline static DoubleReal getNTerminalToFullAverageWeight()
			{
				static const DoubleReal Nterminal_to_full_average_weight = getNTerminalToFull().getAverageWeight();
				return Nterminal_to_full_average_weight;
			}

			inline static DoubleReal getNTerminalToFullMonoWeight()
			{
				static const DoubleReal Nterminal_to_full_mono_weight = getNTerminalToFull().getMonoWeight();
				return Nterminal_to_full_mono_weight;
			}
			
			inline static const EmpiricalFormula& getCTerminalToFull()
			{
				static const EmpiricalFormula Cterminal_to_full = EmpiricalFormula("H");
				return Cterminal_to_full;
			}

			inline static DoubleReal getCTerminalToFullAverageWeight()
			{
				static const DoubleReal Cterminal_to_full_average_weight = getCTerminalToFull().getAverageWeight();
				return Cterminal_to_full_average_weight;
			}
		
			inline static DoubleReal getCTerminalToFullMonoWeight()
			{
				static const DoubleReal Cterminal_to_full_mono_weight = getCTerminalToFull().getMonoWeight();
				return Cterminal_to_full_mono_weight;
			}
		
			inline static const EmpiricalFormula& getBIonToFull()
			{
				static const EmpiricalFormula b_ion_to_full = EmpiricalFormula("HO");
				return b_ion_to_full;
			}

			inline static DoubleReal getBIonToFullAverageWeight()
			{
				static const DoubleReal b_ion_to_full_average_weight = getBIonToFull().getAverageWeight();
				return b_ion_to_full_average_weight;
			}
			
			inline static DoubleReal getBIonToFullMonoWeight()
			{
				static const DoubleReal b_ion_to_full_mono_weight = getBIonToFull().getMonoWeight();
				return b_ion_to_full_mono_weight;
			}
			
			inline static const EmpiricalFormula& getAIonToFull()
			{
				static const EmpiricalFormula a_ion_to_full = EmpiricalFormula("HCO2");
				return a_ion_to_full;
			}
	
			inline static DoubleReal getAIonToFullAverageWeight()
			{
				static const DoubleReal a_ion_to_full_average_weight = getAIonToFull().getAverageWeight();
				return a_ion_to_full_average_weight;
			}

			inline static DoubleReal getAIonToFullMonoWeight()
			{
				static const DoubleReal a_ion_to_full_mono_weight = getAIonToFull().getMonoWeight();
				return a_ion_to_full_mono_weight;
			}

			inline static const EmpiricalFormula& getYIonToFull()
			{
				static const EmpiricalFormula y_ion_to_full = EmpiricalFormula("");
				return y_ion_to_full;
			}

			inline static DoubleReal getYIonToFullAverageWeight()
			{
				static const DoubleReal y_ion_to_full_average_weight = getYIonToFull().getAverageWeight();
				return y_ion_to_full_average_weight;
			}

			inline static DoubleReal getYIonToFullMonoWeight()
			{
				static const DoubleReal y_ion_to_full_mono_weight = getYIonToFull().getMonoWeight();
				return y_ion_to_full_mono_weight;
			}

			inline static const EmpiricalFormula& getCIonToFull()
			{
				static const EmpiricalFormula c_ion_to_full = EmpiricalFormula("");
				return c_ion_to_full;
			}

			inline static DoubleReal getCIonToFullAverageWeight()
			{
				static const DoubleReal c_ion_to_full_average_weight = getCIonToFull().getAverageWeight();
				return c_ion_to_full_average_weight;
			}

			inline static DoubleReal getCIonToFullMonoWeight()
			{
				static const DoubleReal c_ion_to_full_mono_weight = getCIonToFull().getMonoWeight();
				return c_ion_to_full_mono_weight;
			}
			
			inline static const EmpiricalFormula& getXIonToFull()
			{
				static const EmpiricalFormula x_ion_to_full = EmpiricalFormula("HCO");
				return x_ion_to_full;
			}

			inline static DoubleReal getXIonToFullAverageWeight()
			{
				static const DoubleReal x_ion_to_full_average_weight = getXIonToFull().getAverageWeight();
				return x_ion_to_full_average_weight;
			}

			inline static DoubleReal getXIonToFullMonoWeight()
			{
				static const DoubleReal x_ion_to_full_mono_weight = getXIonToFull().getMonoWeight();
				return x_ion_to_full_mono_weight;
			}

			inline static const EmpiricalFormula& getZIonToFull()
			{
				static const EmpiricalFormula z_ion_to_full = EmpiricalFormula("NH2");
				return z_ion_to_full;
			}

			inline static DoubleReal getZIonToFullAverageWeight()
			{
				static const DoubleReal z_ion_to_full_average_weight = getZIonToFull().getAverageWeight();
				return z_ion_to_full_average_weight;
			}

			inline static DoubleReal getZIonToFullMonoWeight()
			{
				static const DoubleReal z_ion_to_full_mono_weight = getZIonToFull().getMonoWeight();
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
							const EmpiricalFormula& neutral_loss);
	
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
			void setLossAverageWeight(DoubleReal weight);
	
			/// return the average weight of the neutral loss molecule
			DoubleReal getLossAverageWeight() const;
	
			/// sets the mono isotopic weight of the neutral loss molecule
			void setLossMonoWeight(DoubleReal weight);
	
			/// returns the mono isotopic weight of the neutral loss molecule
			DoubleReal getLossMonoWeight() const;
	
			/// set the neutral loss molecule weight (if there is one)
			void setLossName(const String& name);
			
			/// gets neutral loss name (if there is one, else returns an empty string)
			const String& getLossName() const;
			
			/// set empirical formula of the residue
			void setFormula(const EmpiricalFormula& formula, ResidueType res_type = Full);
		
			/// returns the empirical formula of thre residue
			EmpiricalFormula getFormula(ResidueType res_type = Full) const;
	
			/// sets average weight of the residue
			void setAverageWeight(DoubleReal weight, ResidueType res_type = Full);
	
			/// returns average weight of the residue
			DoubleReal getAverageWeight(ResidueType res_type = Full) const;
	
			/// sets mono weight of the residue
			void setMonoWeight(DoubleReal weight, ResidueType res_type = Full);
	
			/// returns mono weight of the residue
			DoubleReal getMonoWeight(ResidueType res_type = Full) const;

			/// sets by the name, this mod should be present in ModificationsDB
			void setModification(const String& name);
			
			/// returns the name of the modification to the modification
			const String& getModification() const;
			
			/// sets the name of the unmodified residue
			//void setUnmodifiedName(const String& name);

			/// returns the name of the unmodified residue
			//const String& getUnmodifiedName() const;

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
		
      /// returns the side chain basicity
      DoubleReal getSideChainBasicity() const;

      /// sets the side chain basicity
      void setSideChainBasicity(DoubleReal gb_sc);

      /// returns the backbone basicitiy if located in N-terminal direction
      DoubleReal getBackboneBasicityLeft() const;

      /// sets the N-terminal direction backbone basicitiy
      void setBackboneBasicityLeft(DoubleReal gb_bb_l);

      /// returns the C-terminal direction backbone basicitiy
      DoubleReal getBackboneBasicityRight() const;

      /// sets the C-terminal direction backbone basicity
      void setBackboneBasicityRight(DoubleReal gb_bb_r);

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

			DoubleReal average_weight_;

			DoubleReal mono_weight_;

			// modification
			bool is_modified_;

			String pre_mod_name_;

			String modification_;
			
			// loss
			String loss_name_;

			EmpiricalFormula loss_formula_;

			DoubleReal loss_average_weight_;

			DoubleReal loss_mono_weight_;

			// low mass markers like immonium ions
			std::vector<EmpiricalFormula> low_mass_ions_;
				
			// pka values
			DoubleReal pka_;

			// pkb values
			DoubleReal pkb_;
			
			// pkc values
			DoubleReal pkc_;
		
      DoubleReal gb_sc_;

      DoubleReal gb_bb_l_;

      DoubleReal gb_bb_r_;
	
	};
	
	std::ostream& operator << (std::ostream& os, const Residue& residue);

	/// returns the ion name given as a residue type
	String getResidueTypeName(Residue::ResidueType res_type);
}

#endif
