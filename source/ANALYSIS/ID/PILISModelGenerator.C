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


#include <OpenMS/ANALYSIS/ID/PILISModelGenerator.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

using namespace std;

#define PRECURSOR_MODEL_DEBUG
#undef  PRECURSOR_MODEL_DEBUG

namespace OpenMS 
{
	PILISModelGenerator::PILISModelGenerator()
		: DefaultParamHandler("PILISModelGenerator")
	{
		defaults_.setValue("model_depth", 4, "The number of explicitly modeled backbone cleavages from N-terminus and C-terminus, would be 9 for the default value");
		defaults_.setValue("visible_model_depth", 30, "The maximal possible size of a peptide to be modeled");
		defaults_.setValue("variable_modifications", StringList::create("MOD:00719,MOD:01214"), "Modifications which should be included in the model, represented by PSI-MOD accessions.");
		defaults_.setValue("fixed_modifications", StringList::create(""), "Modifications which should replace the unmodified amino acid, represented by PSI-MOD accessions.");
		defaultsToParam_();
	}

	PILISModelGenerator::PILISModelGenerator(const PILISModelGenerator& rhs)
		:	DefaultParamHandler(rhs)
	{
	}

	PILISModelGenerator& PILISModelGenerator::operator = (const PILISModelGenerator& rhs)
	{
		if (this != &rhs)
		{
			DefaultParamHandler::operator = (rhs);
		}
		return *this;
	}

	PILISModelGenerator::~PILISModelGenerator()
	{
	}

	void PILISModelGenerator::getModel(HiddenMarkovModel& model)
	{
		UInt visible_model_depth = (UInt)param_.getValue("visible_model_depth");
		UInt model_depth = (UInt)param_.getValue("model_depth");

		model.addNewState(new HMMState("endcenter", false));
		model.addNewState(new HMMState("end", false));
		
		model.addNewState("BBcenter");
		model.addNewState("AAcenter");
		model.addNewState("CRcenter");
		model.addNewState("Acenter");
		model.addNewState("SCcenter");
		model.addNewState("ASCcenter");
		
		model.addNewState("bxyz");
		model.addNewState("axyz");
		model.addNewState("D");
		model.addNewState("E");
		
		model.addNewState("AABase1");
		model.addNewState("AABase2");
		
		model.addNewState("K");
		model.addNewState("H");
		model.addNewState("R");

		// 
		set<const Residue*> residues(ResidueDB::getInstance()->getResidues(ResidueDB::NATURAL_20));
		StringList variable_modifications = param_.getValue("variable_modifications");
		for (StringList::const_iterator it = variable_modifications.begin(); it != variable_modifications.end(); ++it)
		{
			residues.insert(ResidueDB::getInstance()->getModifiedResidue(*it));
		}

		cerr << "Using " << residues.size() << " residues" << endl;
		
		//const String residues("ACDEFGHIKLMNPQRSTVWY");
	
		// create the residue states states
		for (set<const Residue*>::const_iterator it = residues.begin(); it != residues.end(); ++it)
    {
			AASequence first_aa;
			first_aa += *it;
      String first(first_aa.toString());
      
			model.addNewState(first + "_D");
			model.addNewState(first + "_E");

			model.addNewState(first + "_K");
			model.addNewState(first + "_H");
			model.addNewState(first + "_R");
			
			for (set<const Residue*>::const_iterator jt = residues.begin(); jt != residues.end(); ++jt)
			{
				AASequence second_aa;
				second_aa += *jt;
				String second(second_aa.toString());
      	
				model.addNewState(first+second +"_bxyz");
				model.addNewState(first+second +"_axyz");
			}
			model.addNewState(first + "_bk-1");
			model.addNewState(first + "_bk-2");
    }

		model.addNewState("bk-1");
		model.addNewState("bk-2");
		model.addNewState(new HMMState("bk-1_-ions", false));
    model.addNewState(new HMMState("bk-2_-ions", false));

		for (Size i = 1; i <= visible_model_depth; ++i)
		{
			// these states are really created
			// charge states
			String num(i);
			model.addNewState("BB" + num);
			model.addNewState("BBk-" + num);
			model.addNewState("CR" + num); 
			model.addNewState("CRk-" + num);

			model.addNewState("SC" + num);
			model.addNewState("SCk-" + num);

			// states for trans mapping 
			model.addNewState("AA"+num);
			model.addNewState("AAk-"+num);

			model.addNewState("A"+num);
			model.addNewState("Ak-"+num);

			model.addNewState("ASC"+num);
			model.addNewState("ASCk-"+num);

			// emitting ion states
			model.addNewState(new HMMState("axyz_" + num + "-ions", false));
			model.addNewState(new HMMState("axyz_k-" + num + "-ions", false));
			model.addNewState(new HMMState("bxyz_" + num + "-ions", false));
			model.addNewState(new HMMState("bxyz_k-" + num + "-ions", false));
			model.addNewState(new HMMState("D_" + num + "-ions", false));
			model.addNewState(new HMMState("D_k-" + num + "-ions", false));
			model.addNewState(new HMMState("E_" + num + "-ions", false));
			model.addNewState(new HMMState("E_k-" + num + "-ions", false));
			model.addNewState(new HMMState("K_" + num + "-ions", false));
			model.addNewState(new HMMState("K_k-" + num + "-ions", false));
			model.addNewState(new HMMState("R_" + num + "-ions", false));
			model.addNewState(new HMMState("R_k-" + num + "-ions", false));
			model.addNewState(new HMMState("H_" + num + "-ions", false));
			model.addNewState(new HMMState("H_k-" + num + "-ions", false));

			//model.addNewState(new HMMState("b" + num + "+", false));
			//model.addNewState(new HMMState("bk-" + num + "+", false));
			//model.addNewState(new HMMState("y" + num + "+", false));
			//model.addNewState(new HMMState("yk-" + num + "+", false));
			//model.addNewState(new HMMState("a" + num + "+", false));
			//model.addNewState(new HMMState("ak-" + num + "+", false));

			//model.addNewState(new HMMState("b"+num+"++", false));
      //model.addNewState(new HMMState("bk-" + num + "++", false));
      //model.addNewState(new HMMState("y" + num + "++", false));
      //model.addNewState(new HMMState("yk-" + num + "++", false));
      //model.addNewState(new HMMState("a" + num + "++", false));
      //model.addNewState(new HMMState("ak-" + num + "++", false));

			model.addNewState(new HMMState("end"+num, false));
			model.addNewState(new HMMState("endk-"+num, false));

			// post AA collector states
			model.addNewState("bxyz"+num);
			model.addNewState("bxyzk-"+num);

			model.addNewState("axyz"+num);
			model.addNewState("axyzk-"+num);
			
			model.addNewState("D"+num);
			model.addNewState("Dk-"+num);
			model.addNewState("E"+num);
			model.addNewState("Ek-"+num);

			model.addNewState("K"+num);
			model.addNewState("Kk-"+num);
			model.addNewState("H"+num);
			model.addNewState("Hk-"+num);
			model.addNewState("R"+num);
			model.addNewState("Rk-"+num);

			// map the residue states
			for (set<const Residue*>::const_iterator it = residues.begin(); it != residues.end(); ++it)
			{
				AASequence first_aa;
				first_aa += *it;
				String first(first_aa.toString());
				
				model.addNewState(first + "_D" + num);
				model.addNewState(first + "_Dk-" + num);

        model.addNewState(first + "_E" + num);
        model.addNewState(first + "_Ek-" + num);

				model.addNewState(first + "_K" + num);
				model.addNewState(first + "_Kk-" + num);

				model.addNewState(first + "_H" + num);
        model.addNewState(first + "_Hk-" + num);

				model.addNewState(first + "_R" + num);
        model.addNewState(first + "_Rk-" + num);
				
				for (set<const Residue*>::const_iterator jt = residues.begin(); jt != residues.end(); ++jt)
				{
					AASequence second_aa;
					second_aa += *jt;
					String second(second_aa.toString());
					
					model.addNewState(first + second + "_bxyz" + num);
					model.addNewState(first + second + "_bxyzk-" + num);

					model.addNewState(first + second + "_axyz" + num);
					model.addNewState(first + second + "_axyzk-" + num);
				}
			}
		}

		model.setTransitionProbability("AABase1", "AABase2", 1);

    // CR bk-1, bk-2
		for (set<const Residue*>::const_iterator it = residues.begin(); it != residues.end(); ++it)
		{
			AASequence first_aa;
			first_aa += *it;
			String first(first_aa.toString());

    	model.addSynonymTransition("AABase1", "AABase2", "Ak-1", first+"_bk-1");
    	model.setTransitionProbability(first+"_bk-1", "bk-1", 0.5);
    	model.setTransitionProbability(first+"_bk-1", "endk-1", 0.5);

			model.addSynonymTransition("AABase1", "AABase2", "Ak-2", first+"_bk-2");
			model.setTransitionProbability(first+"_bk-2", "bk-2", 0.5);
			model.setTransitionProbability(first+"_bk-2", "endk-2", 0.5);
			
		}

		// set the initial transitions
		for (Size i = 1; i <= visible_model_depth; ++i)
		{
			String num(i);
			if (i <= model_depth)
			{
		
				model.setTransitionProbability("BB"+num, "end"+num, 0.5);
				model.setTransitionProbability("BBk-"+num, "endk-"+num, 0.5);
			
				model.setTransitionProbability("SC"+num, "end"+num, 0.5);
				model.setTransitionProbability("SCk-"+num, "endk-"+num, 0.5);

				model.setTransitionProbability("CR"+num, "end"+num, 0.5);
				model.setTransitionProbability("CRk-"+num, "endk-"+num, 0.5);
	
				model.setTransitionProbability("BB"+num, "AA"+num, 0.5);
				model.setTransitionProbability("BBk-"+num, "AAk-"+num, 0.5);

				model.setTransitionProbability("CR"+num, "A"+num, 0.5);
				model.setTransitionProbability("CRk-"+num, "Ak-"+num, 0.5);

        model.setTransitionProbability("SC"+num, "ASC"+num, 0.5);
        model.setTransitionProbability("SCk-"+num, "ASCk-"+num, 0.5);
			}
			else
			{
				model.addSynonymTransition("BBcenter", "endcenter", "BB"+num, "end"+num);
				model.addSynonymTransition("BBcenter", "endcenter", "BBk-"+num, "endk-"+num);
				model.setTransitionProbability("BBcenter", "endcenter", 0.5);
			
				model.addSynonymTransition("CRcenter", "endcenter", "CR"+num, "end"+num);
				model.addSynonymTransition("CRcenter", "endcenter", "CRk-"+num, "endk-"+num);
				model.setTransitionProbability("CRcenter", "endcenter", 0.5);

				model.addSynonymTransition("SCcenter", "endcenter", "SC"+num, "end"+num);
				model.addSynonymTransition("SCcenter", "endcenter", "SCk-"+num, "endk-"+num);
				model.setTransitionProbability("SCcenter", "endcenter", 0.5);

				model.addSynonymTransition("BBcenter", "AAcenter", "BB"+num, "AA"+num);
				model.addSynonymTransition("BBcenter", "AAcenter", "BBk-"+num, "AAk-"+num);
				model.setTransitionProbability("BBcenter", "AAcenter", 0.5);

				model.addSynonymTransition("CRcenter", "Acenter", "CR"+num, "A"+num);
				model.addSynonymTransition("CRcenter", "Acenter", "CRk-"+num, "Ak-"+num);
				model.setTransitionProbability("CRcenter", "Acenter", 0.5);

        model.addSynonymTransition("SCcenter", "ASCcenter", "SC"+num, "ASC"+num);
        model.addSynonymTransition("SCcenter", "ASCcenter", "SCk-"+num, "ASCk-"+num);
        model.setTransitionProbability("SCcenter", "ASCcenter", 0.5);
			}
			
			for (set<const Residue*>::const_iterator it = residues.begin(); it != residues.end(); ++it)
			{
				AASequence first_aa;
				first_aa += *it;
				String first(first_aa.toString());

				// CR D
				model.addSynonymTransition("AABase1", "AABase2", "A"+num, first+"_D"+num);
				model.addSynonymTransition("AABase1", "AABase2", "Ak-"+num, first+"_Dk-"+num);
			
				model.addSynonymTransition(first+"_D", "D", first+"_D"+num, "D"+num);
				model.addSynonymTransition(first+"_D", "end", first+"_D"+num, "end"+num);
				model.addSynonymTransition(first+"_D", "D", first+"_Dk-"+num, "Dk-"+num);
				model.addSynonymTransition(first+"_D", "end", first+"_Dk-"+num, "endk-"+num);

				model.setTransitionProbability(first+"_D", "D", 0.5);
				model.setTransitionProbability(first+"_D", "end", 0.5);
		
				// CR E
				model.addSynonymTransition("AABase1", "AABase2", "A"+num, first+"_E"+num);
        model.addSynonymTransition("AABase1", "AABase2", "Ak-"+num, first+"_Ek-"+num);
        
        model.addSynonymTransition(first+"_E", "E", first+"_E"+num, "E"+num);
				model.addSynonymTransition(first+"_E", "end", first+"_E"+num, "end"+num);
        model.addSynonymTransition(first+"_E", "E", first+"_Ek-"+num, "Ek-"+num);
				model.addSynonymTransition(first+"_E", "end", first+"_Ek-"+num, "endk-"+num);
        
        model.setTransitionProbability(first+"_E", "E", 0.5);
       	model.setTransitionProbability(first+"_E", "end", 0.5);
				
				// SC K
				model.addSynonymTransition("AABase1", "AABase2", "ASC"+num, first+"_K"+num);
        model.addSynonymTransition("AABase1", "AABase2", "ASCk-"+num, first+"_Kk-"+num);

        model.addSynonymTransition(first+"_K", "K", first+"_K"+num, "K"+num);
				model.addSynonymTransition(first+"_K", "end", first+"_K"+num, "end"+num);
        model.addSynonymTransition(first+"_K", "K", first+"_Kk-"+num, "Kk-"+num);
				model.addSynonymTransition(first+"_K", "end", first+"_Kk-"+num, "endk-"+num);

        model.setTransitionProbability(first+"_K", "K", 0.5);
        model.setTransitionProbability(first+"_K", "end", 0.5);
			
				// SC H
				model.addSynonymTransition("AABase1", "AABase2", "ASC"+num, first+"_H"+num);
        model.addSynonymTransition("AABase1", "AABase2", "ASCk-"+num, first+"_Hk-"+num);

        model.addSynonymTransition(first+"_H", "H", first+"_H"+num, "H"+num);
				model.addSynonymTransition(first+"_H", "end", first+"_H"+num, "end"+num);
        model.addSynonymTransition(first+"_H", "H", first+"_Hk-"+num, "Hk-"+num);
				model.addSynonymTransition(first+"_H", "end", first+"_Hk-"+num, "endk-"+num);

        model.setTransitionProbability(first+"_H", "H", 0.5);
        model.setTransitionProbability(first+"_H", "end", 0.5);

				// SC R	
				model.addSynonymTransition("AABase1", "AABase2", "ASC"+num, first+"_R"+num);
        model.addSynonymTransition("AABase1", "AABase2", "ASCk-"+num, first+"_Rk-"+num);

        model.addSynonymTransition(first+"_R", "R", first+"_R"+num, "R"+num);
				model.addSynonymTransition(first+"_R", "end", first+"_R"+num, "end"+num);
        model.addSynonymTransition(first+"_R", "R", first+"_Rk-"+num, "Rk-"+num);
				model.addSynonymTransition(first+"_R", "end", first+"_Rk-"+num, "endk-"+num);

        model.setTransitionProbability(first+"_R", "R", 0.5);
        model.setTransitionProbability(first+"_R", "end", 0.5);
				

				
				for (set<const Residue*>::const_iterator jt = residues.begin(); jt != residues.end(); ++jt)
				{
					AASequence second_aa;
					second_aa += *jt;
					String second(second_aa.toString());

					model.addSynonymTransition("AABase1", "AABase2", "AA"+num, first+second+"_bxyz"+num);
					model.addSynonymTransition("AABase1", "AABase2", "AAk-"+num, first+second+"_bxyzk-"+num);

					model.addSynonymTransition("AABase1", "AABase2", "AA"+num, first+second+"_axyz"+num);
					model.addSynonymTransition("AABase1", "AABase2", "AAk-"+num, first+second+"_axyzk-"+num);
				
					if (i <= 2)
					{
						if (second == "P")
						{
            	model.setTransitionProbability(first+second+"_bxyz"+num, "bxyz"+num, 0.5);
							model.addSynonymTransition(first+second+"_bxyz"+num, "end", first+second+"_bxyz"+num, "end"+num);
							model.setTransitionProbability(first+second+"_bxyz"+num, "end", 0.5);

							model.setTransitionProbability(first+second+"_axyz"+num, "axyz"+num, 0.5);
							model.addSynonymTransition(first+second+"_axyz"+num, "end", first+second+"_axyz"+num, "end"+num);
							model.setTransitionProbability(first+second+"_axyz"+num, "end", 0.5);

						}
						else
						{
							model.setTransitionProbability(first+second+"_bxyz"+num, "bxyz"+num, 0.5);
              model.addSynonymTransition(first+second+"_bxyz"+num, "end", first+second+"_bxyz"+num, "end"+num);
              model.setTransitionProbability(first+second+"_bxyz"+num, "end", 0.5);

							model.setTransitionProbability(first+second+"_axyz"+num, "axyz"+num, 0.5);
							model.addSynonymTransition(first+second+"_axyz"+num, "end", first+second+"_axyz"+num, "end"+num);
							model.setTransitionProbability(first+second+"_axyz"+num, "end", 0.5);
						}
					}
					else
					{
					
						model.addSynonymTransition(first+second+"_bxyz", "bxyz", first+second+"_bxyz"+num, "bxyz"+num);
						model.addSynonymTransition(first+second+"_bxyz", "end", first+second+"_bxyz"+num, "end"+num);

						model.addSynonymTransition(first+second+"_axyz", "axyz", first+second+"_axyz"+num, "axyz"+num);
						model.addSynonymTransition(first+second+"_axyz", "end", first+second+"_axyz"+num, "end"+num);
					}
					
					model.addSynonymTransition(first+second+"_bxyz", "bxyz", first+second+"_bxyzk-"+num, "bxyzk-"+num);
					model.addSynonymTransition(first+second+"_bxyz", "end", first+second+"_bxyzk-"+num, "endk-"+num);
					model.setTransitionProbability(first+second+"_bxyz", "bxyz", 0.5);
					model.setTransitionProbability(first+second+"_bxyz", "end", 0.5);
					model.addSynonymTransition(first+second+"_bxyz", "end", first+second+"_bxyzk-"+num, "end"+num);

					model.addSynonymTransition(first+second+"_axyz", "axyz", first+second+"_axyzk-"+num, "axyzk-"+num);
					model.addSynonymTransition(first+second+"_axyz", "end", first+second+"_axyzk-"+num, "endk-"+num);
					model.setTransitionProbability(first+second+"_axyz", "axyz", 0.5);
					model.setTransitionProbability(first+second+"_axyz", "end", 0.5);
					model.addSynonymTransition(first+second+"_axyz", "end", first+second+"_axyzk-"+num, "end"+num);
				}
			}
		}

		model.disableTransitions();
	}

} // namespace OpenMS

