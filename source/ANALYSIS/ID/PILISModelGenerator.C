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

		for (UInt i = 1; i <= visible_model_depth; ++i)
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
			model.addNewState(new HMMState("b" + num + "+", false));
			model.addNewState(new HMMState("bk-" + num + "+", false));
			model.addNewState(new HMMState("y" + num + "+", false));
			model.addNewState(new HMMState("yk-" + num + "+", false));
			model.addNewState(new HMMState("a" + num + "+", false));
			model.addNewState(new HMMState("ak-" + num + "+", false));

			model.addNewState(new HMMState("b"+num+"++", false));
      model.addNewState(new HMMState("bk-" + num + "++", false));
      model.addNewState(new HMMState("y" + num + "++", false));
      model.addNewState(new HMMState("yk-" + num + "++", false));
      model.addNewState(new HMMState("a" + num + "++", false));
      model.addNewState(new HMMState("ak-" + num + "++", false));

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

    // CR(?) bk-1, bk-2
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
		for (UInt i = 1; i <= visible_model_depth; ++i)
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
		//model.buildSynonyms();

		//model.write(cerr);
	}

	/*
	void PILISModelGenerator::getPrecursorModel(HiddenMarkovModel& precursor_model)
	{
    set<const Residue*> residues(ResidueDB::getInstance()->getResidues(ResidueDB::NATURAL_20));
    StringList variable_modifications = param_.getValue("variable_modifications");

    for (StringList::const_iterator it = variable_modifications.begin(); it != variable_modifications.end(); ++it)
    {
      residues.insert(ResidueDB::getInstance()->getModifiedResidue(*it));
#ifdef PRECURSOR_MODEL_DEBUG
			AASequence aa;
			aa += ResidueDB::getInstance()->getModifiedResidue(*it);
			cerr << "AddingModifiedResidue: " << aa.toString() << endl;
#endif
    }

		set<String> losses;
		for (set<const Residue*>::const_iterator it = residues.begin(); it != residues.end(); ++it)
		{
			vector<EmpiricalFormula> res_losses = (*it)->getLossFormulas();
			for (vector<EmpiricalFormula>::const_iterator loss_it = res_losses.begin(); loss_it != res_losses.end(); ++loss_it)
			{
				String loss = loss_it->getString();
				losses.insert(loss);
#ifdef PRECURSOR_MODEL_DEBUG
				cerr << "Loss: " << loss << ", of residue: " << (*it)->getName() << endl;
#endif
			}
		}

#ifdef PRECURSOR_MODEL_DEBUG
		cerr << "Adding states..." << endl;
#endif
		
		// precursor is not fragmented
		precursor_model.addNewState(new HMMState("p", false));

		// emitting nodes for single losses
		for (set<String>::const_iterator it = losses.begin(); it != losses.end(); ++it)
		{
			precursor_model.addNewState(new HMMState("p-" + *it, false));
		}

		// emitting nodes for double losses
		set<String> double_losses;
		for (set<String>::const_iterator it1 = losses.begin(); it1 != losses.end(); ++it1)
		{
			set<String>::const_iterator it2 = it1;
			++it2;
			for (; it2 != losses.end(); ++it2)
			{
				double_losses.insert(*it1 + "-" + *it2);
			}
		}
		for (set<String>::const_iterator it1 = losses.begin(); it1 != losses.end(); ++it1)
		{
			for (set<String>::const_iterator it2 = it1; it2 != losses.end(); ++it2)
			{
				precursor_model.addNewState(new HMMState("p-" + *it1 + "-" + *it2, false));
			}
		}
		

		// H2O loss from the C-terminus
		String h2o(EmpiricalFormula("H2O").getString());
		precursor_model.addNewState(new HMMState("COOH-" + h2o, true));
		precursor_model.addNewState(new HMMState("start", true));

		// TODO put this into a parameter
		UInt num_explicit(4);
		
		// add double loss states
		// add edges from double loss states to single loss states and emitting states
#ifdef PRECURSOR_MODEL_DEBUG
		cerr << "Adding double loss states" << endl;
#endif

		for (UInt i = 0; i != num_explicit; ++i)
		{
			precursor_model.addNewState(new HMMState("COOH-" + h2o + "_" + String(i + 1)));
		}
		
		for (set<const Residue*>::const_iterator it1 = residues.begin(); it1 != residues.end(); ++it1)
		{
			AASequence aa1;
			aa1 += *it1;

			vector<EmpiricalFormula> res_losses1 = (*it1)->getLossFormulas();
      for (vector<EmpiricalFormula>::const_iterator loss_it1 = res_losses1.begin(); loss_it1 != res_losses1.end(); ++loss_it1)
      {
				String loss1 = loss_it1->getString();
				if (loss1 == "")
				{
					continue;
				}
			
				precursor_model.addNewState(new HMMState(aa1.toString() + "-" + loss1));
				precursor_model.addNewState(new HMMState(aa1.toString() + "-" + loss1 + "-next"));

				for (UInt i = 0; i != num_explicit; ++i)
				{
					precursor_model.addNewState(new HMMState(aa1.toString() + "-" + loss1 + "_" + String(i + 1)));
				}
			
				String losses;
				if (h2o < loss1)
				{
					losses = "-" + h2o  + "-" + loss1;
				}
				else
				{
					losses = "-" + loss1 + "-" + h2o;
				}
				String cooh_name = aa1.toString() + "COOH" + losses;
				precursor_model.addNewState(new HMMState(cooh_name));
				precursor_model.addNewState(new HMMState(cooh_name + "-next"));
				precursor_model.setTransitionProbability(cooh_name, "p" + losses, 0.25);
				precursor_model.setTransitionProbability(cooh_name, aa1.toString() + "-" + loss1, 0.25);
				precursor_model.setTransitionProbability(cooh_name, "COOH-" + h2o, 0.25);
				precursor_model.setTransitionProbability(cooh_name, cooh_name + "-next", 0.25);

				for (UInt i = 0; i != num_explicit; ++i)
				{
					String cooh_name_num = cooh_name + "_" + String(i + 1);
					precursor_model.addNewState(new HMMState(cooh_name_num));
					precursor_model.addSynonymTransition(cooh_name, "p" + losses, cooh_name_num, "p" + losses);
				}
			}
		}
		
			for (set<const Residue*>::const_iterator it1 = residues.begin(); it1 != residues.end(); ++it1)
			{
				AASequence aa1;
				aa1 += *it1;

				vector<EmpiricalFormula> res_losses1 = (*it1)->getLossFormulas();
	      for (vector<EmpiricalFormula>::const_iterator loss_it1 = res_losses1.begin(); loss_it1 != res_losses1.end(); ++loss_it1)
  	    {
    	    String loss1 = loss_it1->getString();				
				
					if (loss1 == "")
					{
						continue;
					}
				
					for (set<const Residue*>::const_iterator it2 = it1; it2 != residues.end(); ++it2)
					{
						AASequence aa2;
						aa2 += *it2;

						vector<EmpiricalFormula> res_losses2 = (*it2)->getLossFormulas();
			      for (vector<EmpiricalFormula>::const_iterator loss_it2 = res_losses2.begin(); loss_it2 != res_losses2.end(); ++loss_it2)
 	  	 		  {
							String loss2 = loss_it2->getString();
							if (loss2 != "")
							{
								String name;
								String losses;
								if (loss1 < loss2)
								{
									losses = "-" + loss1 + "-" + loss2;
								}
								else
								{
									losses = "-" + loss2 + "-" + loss1;
								}
								if (aa1 < aa2)
								{
									name = aa1.toString() + aa2.toString();
								}
								else
								{
									name = aa2.toString() + aa1.toString();
								}
								precursor_model.addNewState(new HMMState(name + losses, true));
								precursor_model.addNewState(new HMMState(name + losses + "-next", true));
			
								for (UInt i = 0; i != num_explicit; ++i)
								{
									precursor_model.addNewState(new HMMState(name + losses + "_" + String(i + 1), true));
									//precursor_model.addNewState(new HMMState(aa1.toString() + "COOH-" + loss2 + "-" + h2o + "_" +  String(i + 1), true));
								}
			
								precursor_model.setTransitionProbability(name + losses, "p" + losses, 0.25);
								precursor_model.setTransitionProbability(name + losses, aa1.toString() + "-" + loss1, 0.25);
								precursor_model.setTransitionProbability(name + losses, aa2.toString() + "-" + loss2, 0.25);
								precursor_model.setTransitionProbability(name + losses, name + losses + "-next", 0.25);
			
								for (UInt i = 0; i != num_explicit; ++i)
								{
									String state_name_num = name + losses + "_" + String(i + 1);
									String state_name = name + losses;
									precursor_model.addSynonymTransition(state_name, "p" + losses,  state_name_num, "p" + losses);
									precursor_model.addSynonymTransition(state_name, aa1.toString() + "-" + loss1, state_name_num, aa1.toString() + "-" + loss1 + "_" + String(i + 1));
									precursor_model.addSynonymTransition(state_name, aa2.toString() + "-" + loss2, state_name_num, aa2.toString() + "-" + loss2 + "_" + String(i + 1));
								}
							}
						}
					}
				}
		}

#ifdef PRECURSOR_MODEL_DEBUG
		cerr << "Adding single loss states" << endl;
#endif
		String cooh_name = "COOH-" + h2o;
		precursor_model.addNewState(new HMMState(cooh_name + "-next"));
		precursor_model.setTransitionProbability(cooh_name, "p-" + h2o, 0.25);
		precursor_model.setTransitionProbability(cooh_name, "p", 0.25);
		precursor_model.setTransitionProbability(cooh_name, "COOH-" + h2o + "-next", 0.5);
	
		precursor_model.addSynonymTransition(cooh_name, "p", cooh_name + "_1", "p");
		precursor_model.addSynonymTransition(cooh_name, "p-" + h2o, cooh_name + "_1", "p-" + h2o);
		
		for (set<const Residue*>::const_iterator it = residues.begin(); it != residues.end(); ++it)
		{
			AASequence aa;
      aa += *it;
			
			vector<EmpiricalFormula> res_losses = (*it)->getLossFormulas();
      for (vector<EmpiricalFormula>::const_iterator loss_it = res_losses.begin(); loss_it != res_losses.end(); ++loss_it)
			{
				String loss = loss_it->getString();
			
	      if (loss == "")
 	     	{
 	      	continue;
 	     	}

				precursor_model.setTransitionProbability(aa.toString() + "-" + loss, "p-" + loss, 0.25);
				precursor_model.setTransitionProbability(aa.toString() + "-" + loss, "p", 0.25);
				precursor_model.setTransitionProbability(aa.toString() + "-" + loss, aa.toString() + "-" + loss + "-next", 0.5);

				for (UInt i = 0; i != num_explicit; ++i)
				{
					String name_num = aa.toString() + "-" + loss + "_" + String(i + 1);
					String name = aa.toString() + "-" + loss;
					precursor_model.addSynonymTransition(name, "p-" + loss, name_num, "p-" + loss);
					precursor_model.addSynonymTransition(name, "p", name_num, "p");
				}
			}
		}

#ifdef PRECURSOR_MODEL_DEBUG
		cerr << "Finalizing HMM" << endl;
#endif
		precursor_model.disableTransitions();
		//precursor_model.buildSynonyms();


#ifdef PRECURSOR_MODEL_DEBUG
    cerr << "#States: " << precursor_model.getNumberOfStates() << endl;
    precursor_model.writeGraphMLFile("precursor_model.graphML");
#endif
		
		return;
	}
*/


} // namespace OpenMS

