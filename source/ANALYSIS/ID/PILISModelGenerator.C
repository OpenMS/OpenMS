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
		defaults_.setValue("modifications", StringList::create("MOD:00719,MOD:01214"), "Modifications which should be included in the model, represented by PSI-MOD accessions.");
		defaultsToParam_();
	}

	PILISModelGenerator::PILISModelGenerator(const PILISModelGenerator& rhs)
		:	DefaultParamHandler(rhs)
	{
		model_ = rhs.model_;
	}

	PILISModelGenerator& PILISModelGenerator::operator = (const PILISModelGenerator& rhs)
	{
		if (this != &rhs)
		{
			DefaultParamHandler::operator = (rhs);
			model_ = rhs.model_;
		}
		return *this;
	}

	PILISModelGenerator::~PILISModelGenerator()
	{
	}

	const PILISModel& PILISModelGenerator::getModel()
	{
		initModels_();
		model_.valid_ = true;
		return model_;
	}
	
	void PILISModelGenerator::initModels_()
	{
		initMainModel_();
		//initLossModels_();
		initPrecursorModel_();
		return;
	}

	void PILISModelGenerator::initMainModel_()
	{
		UInt visible_model_depth = (UInt)param_.getValue("visible_model_depth");
		UInt model_depth = (UInt)param_.getValue("model_depth");

		model_.hmm_.addNewState(new HMMState("endcenter", false));
		model_.hmm_.addNewState(new HMMState("end", false));
		
		model_.hmm_.addNewState("BBcenter");
		model_.hmm_.addNewState("AAcenter");
		model_.hmm_.addNewState("CRcenter");
		model_.hmm_.addNewState("Acenter");
		model_.hmm_.addNewState("SCcenter");
		model_.hmm_.addNewState("ASCcenter");
		
		model_.hmm_.addNewState("bxyz");
		model_.hmm_.addNewState("axyz");
		model_.hmm_.addNewState("D");
		model_.hmm_.addNewState("E");
		
		model_.hmm_.addNewState("AABase1");
		model_.hmm_.addNewState("AABase2");
		
		model_.hmm_.addNewState("K");
		model_.hmm_.addNewState("H");
		model_.hmm_.addNewState("R");

		// 
		set<const Residue*> residues(ResidueDB::getInstance()->getResidues(ResidueDB::NATURAL_20));
		StringList modifications = param_.getValue("modifications");
		for (StringList::const_iterator it = modifications.begin(); it != modifications.end(); ++it)
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
      
			model_.hmm_.addNewState(first + "_D");
			model_.hmm_.addNewState(first + "_E");

			model_.hmm_.addNewState(first + "_K");
			model_.hmm_.addNewState(first + "_H");
			model_.hmm_.addNewState(first + "_R");
			
			for (set<const Residue*>::const_iterator jt = residues.begin(); jt != residues.end(); ++jt)
			{
				AASequence second_aa;
				second_aa += *jt;
				String second(second_aa.toString());
      	
				model_.hmm_.addNewState(first+second +"_bxyz");
				model_.hmm_.addNewState(first+second +"_axyz");
			}
			model_.hmm_.addNewState(first + "_bk-1");
			model_.hmm_.addNewState(first + "_bk-2");
    }
		model_.hmm_.addNewState("bk-1");
		model_.hmm_.addNewState("bk-2");

		for (UInt i = 1; i <= visible_model_depth; ++i)
		{
			// these states are really created
			// charge states
			String num(i);
			model_.hmm_.addNewState("BB" + num);
			model_.hmm_.addNewState("BBk-" + num);
			model_.hmm_.addNewState("CR" + num); 
			model_.hmm_.addNewState("CRk-" + num);

			model_.hmm_.addNewState("SC" + num);
			model_.hmm_.addNewState("SCk-" + num);

			// states for trans mapping 
			model_.hmm_.addNewState("AA"+num);
			model_.hmm_.addNewState("AAk-"+num);

			model_.hmm_.addNewState("A"+num);
			model_.hmm_.addNewState("Ak-"+num);

			model_.hmm_.addNewState("ASC"+num);
			model_.hmm_.addNewState("ASCk-"+num);

			// emitting ion states
			model_.hmm_.addNewState(new HMMState("b" + num + "+", false));
			model_.hmm_.addNewState(new HMMState("bk-" + num + "+", false));
			model_.hmm_.addNewState(new HMMState("y" + num + "+", false));
			model_.hmm_.addNewState(new HMMState("yk-" + num + "+", false));
			model_.hmm_.addNewState(new HMMState("a" + num + "+", false));
			model_.hmm_.addNewState(new HMMState("ak-" + num + "+", false));

			model_.hmm_.addNewState(new HMMState("b"+num+"++", false));
      model_.hmm_.addNewState(new HMMState("bk-" + num + "++", false));
      model_.hmm_.addNewState(new HMMState("y" + num + "++", false));
      model_.hmm_.addNewState(new HMMState("yk-" + num + "++", false));
      model_.hmm_.addNewState(new HMMState("a" + num + "++", false));
      model_.hmm_.addNewState(new HMMState("ak-" + num + "++", false));

			model_.hmm_.addNewState(new HMMState("end"+num, false));
			model_.hmm_.addNewState(new HMMState("endk-"+num, false));

			// post AA collector states
			model_.hmm_.addNewState("bxyz"+num);
			model_.hmm_.addNewState("bxyzk-"+num);

			model_.hmm_.addNewState("axyz"+num);
			model_.hmm_.addNewState("axyzk-"+num);
			
			model_.hmm_.addNewState("D"+num);
			model_.hmm_.addNewState("Dk-"+num);
			model_.hmm_.addNewState("E"+num);
			model_.hmm_.addNewState("Ek-"+num);

			model_.hmm_.addNewState("K"+num);
			model_.hmm_.addNewState("Kk-"+num);
			model_.hmm_.addNewState("H"+num);
			model_.hmm_.addNewState("Hk-"+num);
			model_.hmm_.addNewState("R"+num);
			model_.hmm_.addNewState("Rk-"+num);

			// map the residue states
			for (set<const Residue*>::const_iterator it = residues.begin(); it != residues.end(); ++it)
			{
				AASequence first_aa;
				first_aa += *it;
				String first(first_aa.toString());
				
				model_.hmm_.addNewState(first + "_D" + num);
				model_.hmm_.addNewState(first + "_Dk-" + num);

        model_.hmm_.addNewState(first + "_E" + num);
        model_.hmm_.addNewState(first + "_Ek-" + num);

				model_.hmm_.addNewState(first + "_K" + num);
				model_.hmm_.addNewState(first + "_Kk-" + num);

				model_.hmm_.addNewState(first + "_H" + num);
        model_.hmm_.addNewState(first + "_Hk-" + num);

				model_.hmm_.addNewState(first + "_R" + num);
        model_.hmm_.addNewState(first + "_Rk-" + num);
				
				for (set<const Residue*>::const_iterator jt = residues.begin(); jt != residues.end(); ++jt)
				{
					AASequence second_aa;
					second_aa += *jt;
					String second(second_aa.toString());
					
					model_.hmm_.addNewState(first + second + "_bxyz" + num);
					model_.hmm_.addNewState(first + second + "_bxyzk-" + num);

					model_.hmm_.addNewState(first + second + "_axyz" + num);
					model_.hmm_.addNewState(first + second + "_axyzk-" + num);
				}
			}
		}

		model_.hmm_.setTransitionProbability("AABase1", "AABase2", 1);

    // CR(?) bk-1, bk-2
		for (set<const Residue*>::const_iterator it = residues.begin(); it != residues.end(); ++it)
		{
			AASequence first_aa;
			first_aa += *it;
			String first(first_aa.toString());

    	model_.hmm_.addSynonymTransition("AABase1", "AABase2", "Ak-1", first+"_bk-1");

    	model_.hmm_.setTransitionProbability(first+"_bk-1", "bk-1", 0.5);
    	model_.hmm_.setTransitionProbability(first+"_bk-1", "endk-1", 0.5);

			model_.hmm_.addSynonymTransition("AABase1", "AABase2", "Ak-2", first+"_bk-2");
			model_.hmm_.setTransitionProbability(first+"_bk-2", "bk-2", 0.5);
			model_.hmm_.setTransitionProbability(first+"_bk-2", "endk-2", 0.5);
			
		}

		// set the initial transitions
		for (UInt i = 1; i <= visible_model_depth; ++i)
		{
			String num(i);
			if (i <= model_depth)
			{
		
				model_.hmm_.setTransitionProbability("BB"+num, "end"+num, 0.5);
				model_.hmm_.setTransitionProbability("BBk-"+num, "endk-"+num, 0.5);
			
				model_.hmm_.setTransitionProbability("SC"+num, "end"+num, 0.5);
				model_.hmm_.setTransitionProbability("SCk-"+num, "endk-"+num, 0.5);

				model_.hmm_.setTransitionProbability("CR"+num, "end"+num, 0.5);
				model_.hmm_.setTransitionProbability("CRk-"+num, "endk-"+num, 0.5);
	
				model_.hmm_.setTransitionProbability("BB"+num, "AA"+num, 0.5);
				model_.hmm_.setTransitionProbability("BBk-"+num, "AAk-"+num, 0.5);

				model_.hmm_.setTransitionProbability("CR"+num, "A"+num, 0.5);
				model_.hmm_.setTransitionProbability("CRk-"+num, "Ak-"+num, 0.5);

        model_.hmm_.setTransitionProbability("SC"+num, "ASC"+num, 0.5);
        model_.hmm_.setTransitionProbability("SCk-"+num, "ASCk-"+num, 0.5);
			}
			else
			{
				model_.hmm_.addSynonymTransition("BBcenter", "endcenter", "BB"+num, "end"+num);
				model_.hmm_.addSynonymTransition("BBcenter", "endcenter", "BBk-"+num, "endk-"+num);
				model_.hmm_.setTransitionProbability("BBcenter", "endcenter", 0.5);
			
				model_.hmm_.addSynonymTransition("CRcenter", "endcenter", "CR"+num, "end"+num);
				model_.hmm_.addSynonymTransition("CRcenter", "endcenter", "CRk-"+num, "endk-"+num);
				model_.hmm_.setTransitionProbability("CRcenter", "endcenter", 0.5);

				model_.hmm_.addSynonymTransition("SCcenter", "endcenter", "SC"+num, "end"+num);
				model_.hmm_.addSynonymTransition("SCcenter", "endcenter", "SCk-"+num, "endk-"+num);
				model_.hmm_.setTransitionProbability("SCcenter", "endcenter", 0.5);

				model_.hmm_.addSynonymTransition("BBcenter", "AAcenter", "BB"+num, "AA"+num);
				model_.hmm_.addSynonymTransition("BBcenter", "AAcenter", "BBk-"+num, "AAk-"+num);
				model_.hmm_.setTransitionProbability("BBcenter", "AAcenter", 0.5);

				model_.hmm_.addSynonymTransition("CRcenter", "Acenter", "CR"+num, "A"+num);
				model_.hmm_.addSynonymTransition("CRcenter", "Acenter", "CRk-"+num, "Ak-"+num);
				model_.hmm_.setTransitionProbability("CRcenter", "Acenter", 0.5);

        model_.hmm_.addSynonymTransition("SCcenter", "ASCcenter", "SC"+num, "ASC"+num);
        model_.hmm_.addSynonymTransition("SCcenter", "ASCcenter", "SCk-"+num, "ASCk-"+num);
        model_.hmm_.setTransitionProbability("SCcenter", "ASCcenter", 0.5);
			}
			
			for (set<const Residue*>::const_iterator it = residues.begin(); it != residues.end(); ++it)
			{
				AASequence first_aa;
				first_aa += *it;
				String first(first_aa.toString());

				// CR D
				model_.hmm_.addSynonymTransition("AABase1", "AABase2", "A"+num, first+"_D"+num);
				model_.hmm_.addSynonymTransition("AABase1", "AABase2", "Ak-"+num, first+"_Dk-"+num);
			
				model_.hmm_.addSynonymTransition(first+"_D", "D", first+"_D"+num, "D"+num);
				model_.hmm_.addSynonymTransition(first+"_D", "end", first+"_D"+num, "end"+num);
				model_.hmm_.addSynonymTransition(first+"_D", "D", first+"_Dk-"+num, "Dk-"+num);
				model_.hmm_.addSynonymTransition(first+"_D", "end", first+"_Dk-"+num, "endk-"+num);

				model_.hmm_.setTransitionProbability(first+"_D", "D", 0.5);
				model_.hmm_.setTransitionProbability(first+"_D", "end", 0.5);
		
				// CR E
				model_.hmm_.addSynonymTransition("AABase1", "AABase2", "A"+num, first+"_E"+num);
        model_.hmm_.addSynonymTransition("AABase1", "AABase2", "Ak-"+num, first+"_Ek-"+num);
        
        model_.hmm_.addSynonymTransition(first+"_E", "E", first+"_E"+num, "E"+num);
				model_.hmm_.addSynonymTransition(first+"_E", "end", first+"_E"+num, "end"+num);
        model_.hmm_.addSynonymTransition(first+"_E", "E", first+"_Ek-"+num, "Ek-"+num);
				model_.hmm_.addSynonymTransition(first+"_E", "end", first+"_Ek-"+num, "endk-"+num);
        
        model_.hmm_.setTransitionProbability(first+"_E", "E", 0.5);
       	model_.hmm_.setTransitionProbability(first+"_E", "end", 0.5);
				
				// SC K
				model_.hmm_.addSynonymTransition("AABase1", "AABase2", "ASC"+num, first+"_K"+num);
        model_.hmm_.addSynonymTransition("AABase1", "AABase2", "ASCk-"+num, first+"_Kk-"+num);

        model_.hmm_.addSynonymTransition(first+"_K", "K", first+"_K"+num, "K"+num);
				model_.hmm_.addSynonymTransition(first+"_K", "end", first+"_K"+num, "end"+num);
        model_.hmm_.addSynonymTransition(first+"_K", "K", first+"_Kk-"+num, "Kk-"+num);
				model_.hmm_.addSynonymTransition(first+"_K", "end", first+"_Kk-"+num, "endk-"+num);

        model_.hmm_.setTransitionProbability(first+"_K", "K", 0.5);
        model_.hmm_.setTransitionProbability(first+"_K", "end", 0.5);
			
				// SC H
				model_.hmm_.addSynonymTransition("AABase1", "AABase2", "ASC"+num, first+"_H"+num);
        model_.hmm_.addSynonymTransition("AABase1", "AABase2", "ASCk-"+num, first+"_Hk-"+num);

        model_.hmm_.addSynonymTransition(first+"_H", "H", first+"_H"+num, "H"+num);
				model_.hmm_.addSynonymTransition(first+"_H", "end", first+"_H"+num, "end"+num);
        model_.hmm_.addSynonymTransition(first+"_H", "H", first+"_Hk-"+num, "Hk-"+num);
				model_.hmm_.addSynonymTransition(first+"_H", "end", first+"_Hk-"+num, "endk-"+num);

        model_.hmm_.setTransitionProbability(first+"_H", "H", 0.5);
        model_.hmm_.setTransitionProbability(first+"_H", "end", 0.5);

				// SC R	
				model_.hmm_.addSynonymTransition("AABase1", "AABase2", "ASC"+num, first+"_R"+num);
        model_.hmm_.addSynonymTransition("AABase1", "AABase2", "ASCk-"+num, first+"_Rk-"+num);

        model_.hmm_.addSynonymTransition(first+"_R", "R", first+"_R"+num, "R"+num);
				model_.hmm_.addSynonymTransition(first+"_R", "end", first+"_R"+num, "end"+num);
        model_.hmm_.addSynonymTransition(first+"_R", "R", first+"_Rk-"+num, "Rk-"+num);
				model_.hmm_.addSynonymTransition(first+"_R", "end", first+"_Rk-"+num, "endk-"+num);

        model_.hmm_.setTransitionProbability(first+"_R", "R", 0.5);
        model_.hmm_.setTransitionProbability(first+"_R", "end", 0.5);
				

				
				for (set<const Residue*>::const_iterator jt = residues.begin(); jt != residues.end(); ++jt)
				{
					AASequence second_aa;
					second_aa += *jt;
					String second(second_aa.toString());

					model_.hmm_.addSynonymTransition("AABase1", "AABase2", "AA"+num, first+second+"_bxyz"+num);
					model_.hmm_.addSynonymTransition("AABase1", "AABase2", "AAk-"+num, first+second+"_bxyzk-"+num);

					model_.hmm_.addSynonymTransition("AABase1", "AABase2", "AA"+num, first+second+"_axyz"+num);
					model_.hmm_.addSynonymTransition("AABase1", "AABase2", "AAk-"+num, first+second+"_axyzk-"+num);
				
					if (i <= 2)
					{
						if (second == "P")
						{
            	model_.hmm_.setTransitionProbability(first+second+"_bxyz"+num, "bxyz"+num, 0.5);
							model_.hmm_.addSynonymTransition(first+second+"_bxyz"+num, "end", first+second+"_bxyz"+num, "end"+num);
							model_.hmm_.setTransitionProbability(first+second+"_bxyz"+num, "end", 0.5);

							model_.hmm_.setTransitionProbability(first+second+"_axyz"+num, "axyz"+num, 0.5);
							model_.hmm_.addSynonymTransition(first+second+"_axyz"+num, "end", first+second+"_axyz"+num, "end"+num);
							model_.hmm_.setTransitionProbability(first+second+"_axyz"+num, "end", 0.5);

						}
						else
						{
							model_.hmm_.setTransitionProbability(first+second+"_bxyz"+num, "bxyz"+num, 0.5);
              model_.hmm_.addSynonymTransition(first+second+"_bxyz"+num, "end", first+second+"_bxyz"+num, "end"+num);
              model_.hmm_.setTransitionProbability(first+second+"_bxyz"+num, "end", 0.5);

							model_.hmm_.setTransitionProbability(first+second+"_axyz"+num, "axyz"+num, 0.5);
							model_.hmm_.addSynonymTransition(first+second+"_axyz"+num, "end", first+second+"_axyz"+num, "end"+num);
							model_.hmm_.setTransitionProbability(first+second+"_axyz"+num, "end", 0.5);
						}
					}
					else
					{
					
						model_.hmm_.addSynonymTransition(first+second+"_bxyz", "bxyz", first+second+"_bxyz"+num, "bxyz"+num);
						model_.hmm_.addSynonymTransition(first+second+"_bxyz", "end", first+second+"_bxyz"+num, "end"+num);

						model_.hmm_.addSynonymTransition(first+second+"_axyz", "axyz", first+second+"_axyz"+num, "axyz"+num);
						model_.hmm_.addSynonymTransition(first+second+"_axyz", "end", first+second+"_axyz"+num, "end"+num);
					}
					
					model_.hmm_.addSynonymTransition(first+second+"_bxyz", "bxyz", first+second+"_bxyzk-"+num, "bxyzk-"+num);
					model_.hmm_.addSynonymTransition(first+second+"_bxyz", "end", first+second+"_bxyzk-"+num, "endk-"+num);
					model_.hmm_.setTransitionProbability(first+second+"_bxyz", "bxyz", 0.5);
					model_.hmm_.setTransitionProbability(first+second+"_bxyz", "end", 0.5);
					model_.hmm_.addSynonymTransition(first+second+"_bxyz", "end", first+second+"_bxyzk-"+num, "end"+num);

					model_.hmm_.addSynonymTransition(first+second+"_axyz", "axyz", first+second+"_axyzk-"+num, "axyzk-"+num);
					model_.hmm_.addSynonymTransition(first+second+"_axyz", "end", first+second+"_axyzk-"+num, "endk-"+num);
					model_.hmm_.setTransitionProbability(first+second+"_axyz", "axyz", 0.5);
					model_.hmm_.setTransitionProbability(first+second+"_axyz", "end", 0.5);
					model_.hmm_.addSynonymTransition(first+second+"_axyz", "end", first+second+"_axyzk-"+num, "end"+num);
				}
			}
		}

		model_.hmm_.disableTransitions();
		model_.hmm_.buildSynonyms();

		//model_.hmm_.write(cerr);
	}

	/*
	void PILISModelGenerator::initLossModels_()
	{
		//# new states from this HMM
		//# format State <Name> [<Hidden?>]
		//# emitting states
		//State y-H2O false
		//State y-NH3 false
		//State y_loss_end false
		model_.hmms_losses_[PILISModel::YIon].addNewState(new HMMStateLight(PILISModel::Y_H2O, false));
		model_.hmms_losses_[PILISModel::YIon].addNewState(new HMMStateLight(PILISModel::Y_NH3, false));
		//model_.hmms_losses_[PILISModel::YIon].addNewState(new HMMStateLight(PILISModel::Y_LOSS_END, false));
		model_.hmms_losses_[PILISModel::YIon].addNewState(new HMMStateLight(PILISModel::Y_LOSS_END, false));

		//# base states
		//State yion
		//State y_Base1
		//State y_Base2
		model_.hmms_losses_[PILISModel::YIon].addNewState(new HMMStateLight(PILISModel::Y_ION));
		model_.hmms_losses_[PILISModel::YIon].addNewState(new HMMStateLight(PILISModel::Y_BASE1));
		model_.hmms_losses_[PILISModel::YIon].addNewState(new HMMStateLight(PILISModel::Y_BASE2));

		//# loss pathway states
		//State y_H2O_D
		//State y_H2O_E
		//State y_H2O_T
		//State y_H2O_S
		//State y_H2O_Q1
		//State y_H2O_Cterm
		//State y_NH3_K
		//State y_NH3_R
		model_.hmms_losses_[PILISModel::YIon].addNewState(new HMMStateLight(PILISModel::Y_H2O_D));
		model_.hmms_losses_[PILISModel::YIon].addNewState(new HMMStateLight(PILISModel::Y_H2O_E));
		model_.hmms_losses_[PILISModel::YIon].addNewState(new HMMStateLight(PILISModel::Y_H2O_S));
		model_.hmms_losses_[PILISModel::YIon].addNewState(new HMMStateLight(PILISModel::Y_H2O_T));
		model_.hmms_losses_[PILISModel::YIon].addNewState(new HMMStateLight(PILISModel::Y_H2O_Q1));
		model_.hmms_losses_[PILISModel::YIon].addNewState(new HMMStateLight(PILISModel::Y_H2O_CTERM));
		model_.hmms_losses_[PILISModel::YIon].addNewState(new HMMStateLight(PILISModel::Y_NH3_K));
		model_.hmms_losses_[PILISModel::YIon].addNewState(new HMMStateLight(PILISModel::Y_NH3_R));
		model_.hmms_losses_[PILISModel::YIon].addNewState(new HMMStateLight(PILISModel::Y_NH3_Q));
		model_.hmms_losses_[PILISModel::YIon].addNewState(new HMMStateLight(PILISModel::Y_NH3_N));
		
		//# synonym states are not needed in this state
		//# transitions
		//Transition y_Base1 y_Base2 1
		model_.hmms_losses_[PILISModel::YIon].setTransitionProbability(PILISModel::Y_BASE1, PILISModel::Y_BASE2, 1.0);

		//# H2O losses
		//Transition y_H2O_D y-H2O 0.1
		//Transition y_H2O_D y_loss_end 0.9
		//Transition y_H2O_E y-H2O 0.1
		//Transition y_H2O_E y_loss_end 0.9
		//Transition y_H2O_T y-H2O 0.1
		//Transition y_H2O_T y_loss_end 0.9
		//Transition y_H2O_S y-H2O 0.1
		//Transition y_H2O_S y_loss_end 0.9
		//Transition y_H2O_Q1 y-H2O 0.1
		//Transition y_H2O_Q1 y_loss_end 0.9
		//Transition y_H2O_Cterm y-H2O 0.1
		//Transition y_H2O_Cterm y_loss_end 0.9
		model_.hmms_losses_[PILISModel::YIon].setTransitionProbability(PILISModel::Y_H2O_D, PILISModel::Y_H2O, 0.1);
		model_.hmms_losses_[PILISModel::YIon].setTransitionProbability(PILISModel::Y_H2O_D, PILISModel::Y_LOSS_END, 0.9);
		model_.hmms_losses_[PILISModel::YIon].setTransitionProbability(PILISModel::Y_H2O_E, PILISModel::Y_H2O, 0.1);
		model_.hmms_losses_[PILISModel::YIon].setTransitionProbability(PILISModel::Y_H2O_E, PILISModel::Y_LOSS_END, 0.9);
		model_.hmms_losses_[PILISModel::YIon].setTransitionProbability(PILISModel::Y_H2O_S, PILISModel::Y_H2O, 0.1);
		model_.hmms_losses_[PILISModel::YIon].setTransitionProbability(PILISModel::Y_H2O_S, PILISModel::Y_LOSS_END, 0.9);
		model_.hmms_losses_[PILISModel::YIon].setTransitionProbability(PILISModel::Y_H2O_T, PILISModel::Y_H2O, 0.1);
		model_.hmms_losses_[PILISModel::YIon].setTransitionProbability(PILISModel::Y_H2O_T, PILISModel::Y_LOSS_END, 0.9);
		model_.hmms_losses_[PILISModel::YIon].setTransitionProbability(PILISModel::Y_H2O_Q1, PILISModel::Y_H2O, 0.1);
		model_.hmms_losses_[PILISModel::YIon].setTransitionProbability(PILISModel::Y_H2O_Q1, PILISModel::Y_LOSS_END, 0.9);
		model_.hmms_losses_[PILISModel::YIon].setTransitionProbability(PILISModel::Y_H2O_CTERM, PILISModel::Y_H2O, 0.1);
		model_.hmms_losses_[PILISModel::YIon].setTransitionProbability(PILISModel::Y_H2O_CTERM, PILISModel::Y_LOSS_END, 0.9);

		//# NH3 losses
		//Transition y_NH3_K y-NH3 0.1
		//Transition y_NH3_K y_loss_end 0.9
		//Transition y_NH3_R y-NH3 0.1
		//Transition y_NH3_R y_loss_end 0.9
		model_.hmms_losses_[PILISModel::YIon].setTransitionProbability(PILISModel::Y_NH3_K, PILISModel::Y_NH3, 0.1);
    model_.hmms_losses_[PILISModel::YIon].setTransitionProbability(PILISModel::Y_NH3_K, PILISModel::Y_LOSS_END, 0.9);
    model_.hmms_losses_[PILISModel::YIon].setTransitionProbability(PILISModel::Y_NH3_R, PILISModel::Y_NH3, 0.1);
    model_.hmms_losses_[PILISModel::YIon].setTransitionProbability(PILISModel::Y_NH3_R, PILISModel::Y_LOSS_END, 0.9);
		model_.hmms_losses_[PILISModel::YIon].setTransitionProbability(PILISModel::Y_NH3_Q, PILISModel::Y_NH3, 0.1);
    model_.hmms_losses_[PILISModel::YIon].setTransitionProbability(PILISModel::Y_NH3_Q, PILISModel::Y_LOSS_END, 0.9);
		model_.hmms_losses_[PILISModel::YIon].setTransitionProbability(PILISModel::Y_NH3_N, PILISModel::Y_NH3, 0.1);
		model_.hmms_losses_[PILISModel::YIon].setTransitionProbability(PILISModel::Y_NH3_N, PILISModel::Y_LOSS_END, 0.9);
				
		//# synonym transitions
		//Synonym yion y_H2O_D y_Base1 y_Base2
		//Synonym yion y_H2O_E y_Base1 y_Base2
		//Synonym yion y_H2O_T y_Base1 y_Base2
		//Synonym yion y_H2O_S y_Base1 y_Base2
		//Synonym yion y_H2O_Q1 y_Base1 y_Base2
		//Synonym yion y_H2O_Cterm y_Base1 y_Base2
		//Synonym yion y_NH3_K y_Base1 y_Base2
		//Synonym yion y_NH3_R y_Base1 y_Base2
		model_.hmms_losses_[PILISModel::YIon].addSynonymTransition(PILISModel::Y_BASE1, PILISModel::Y_BASE2, PILISModel::Y_ION, PILISModel::Y_H2O_D);
    model_.hmms_losses_[PILISModel::YIon].addSynonymTransition(PILISModel::Y_BASE1, PILISModel::Y_BASE2, PILISModel::Y_ION, PILISModel::Y_H2O_E);
    model_.hmms_losses_[PILISModel::YIon].addSynonymTransition(PILISModel::Y_BASE1, PILISModel::Y_BASE2, PILISModel::Y_ION, PILISModel::Y_H2O_S);
    model_.hmms_losses_[PILISModel::YIon].addSynonymTransition(PILISModel::Y_BASE1, PILISModel::Y_BASE2, PILISModel::Y_ION, PILISModel::Y_H2O_T);
    model_.hmms_losses_[PILISModel::YIon].addSynonymTransition(PILISModel::Y_BASE1, PILISModel::Y_BASE2, PILISModel::Y_ION, PILISModel::Y_H2O_Q1);
    model_.hmms_losses_[PILISModel::YIon].addSynonymTransition(PILISModel::Y_BASE1, PILISModel::Y_BASE2, PILISModel::Y_ION, PILISModel::Y_H2O_CTERM);
    model_.hmms_losses_[PILISModel::YIon].addSynonymTransition(PILISModel::Y_BASE1, PILISModel::Y_BASE2, PILISModel::Y_ION, PILISModel::Y_NH3_K);
    model_.hmms_losses_[PILISModel::YIon].addSynonymTransition(PILISModel::Y_BASE1, PILISModel::Y_BASE2, PILISModel::Y_ION, PILISModel::Y_NH3_R);
		model_.hmms_losses_[PILISModel::YIon].addSynonymTransition(PILISModel::Y_BASE1, PILISModel::Y_BASE2, PILISModel::Y_ION, PILISModel::Y_NH3_Q);
		model_.hmms_losses_[PILISModel::YIon].addSynonymTransition(PILISModel::Y_BASE1, PILISModel::Y_BASE2, PILISModel::Y_ION, PILISModel::Y_NH3_N);

		model_.hmms_losses_[PILISModel::YIon].disableTransitions();
		model_.hmms_losses_[PILISModel::YIon].buildSynonyms();


		// b-ions
		model_.hmms_losses_[PILISModel::BIon].addNewState(new HMMStateLight(PILISModel::B_H2O, false));
    model_.hmms_losses_[PILISModel::BIon].addNewState(new HMMStateLight(PILISModel::B_NH3, false));
		model_.hmms_losses_[PILISModel::BIon].addNewState(new HMMStateLight(PILISModel::A_ION, false));
    model_.hmms_losses_[PILISModel::BIon].addNewState(new HMMStateLight(PILISModel::B_LOSS_END, false));

    model_.hmms_losses_[PILISModel::BIon].addNewState(new HMMStateLight(PILISModel::B_ION));
    model_.hmms_losses_[PILISModel::BIon].addNewState(new HMMStateLight(PILISModel::B_BASE1));
    model_.hmms_losses_[PILISModel::BIon].addNewState(new HMMStateLight(PILISModel::B_BASE2));

    model_.hmms_losses_[PILISModel::BIon].addNewState(new HMMStateLight(PILISModel::B_H2O_D));
    model_.hmms_losses_[PILISModel::BIon].addNewState(new HMMStateLight(PILISModel::B_H2O_E));
    model_.hmms_losses_[PILISModel::BIon].addNewState(new HMMStateLight(PILISModel::B_H2O_S));
    model_.hmms_losses_[PILISModel::BIon].addNewState(new HMMStateLight(PILISModel::B_H2O_T));
		model_.hmms_losses_[PILISModel::BIon].addNewState(new HMMStateLight(PILISModel::B_H2O_Q1));
    model_.hmms_losses_[PILISModel::BIon].addNewState(new HMMStateLight(PILISModel::B_NH3_K));
    model_.hmms_losses_[PILISModel::BIon].addNewState(new HMMStateLight(PILISModel::B_NH3_R));
		model_.hmms_losses_[PILISModel::BIon].addNewState(new HMMStateLight(PILISModel::B_NH3_Q));
		model_.hmms_losses_[PILISModel::BIon].addNewState(new HMMStateLight(PILISModel::B_NH3_N));
		model_.hmms_losses_[PILISModel::BIon].addNewState(new HMMStateLight(PILISModel::B_CO));
		

    model_.hmms_losses_[PILISModel::BIon].setTransitionProbability(PILISModel::B_BASE1, PILISModel::B_BASE2, 1.0);

    model_.hmms_losses_[PILISModel::BIon].setTransitionProbability(PILISModel::B_H2O_D, PILISModel::B_H2O, 0.1);
    model_.hmms_losses_[PILISModel::BIon].setTransitionProbability(PILISModel::B_H2O_D, PILISModel::B_LOSS_END, 0.9);
    model_.hmms_losses_[PILISModel::BIon].setTransitionProbability(PILISModel::B_H2O_E, PILISModel::B_H2O, 0.1);
    model_.hmms_losses_[PILISModel::BIon].setTransitionProbability(PILISModel::B_H2O_E, PILISModel::B_LOSS_END, 0.9);
    model_.hmms_losses_[PILISModel::BIon].setTransitionProbability(PILISModel::B_H2O_S, PILISModel::B_H2O, 0.1);
    model_.hmms_losses_[PILISModel::BIon].setTransitionProbability(PILISModel::B_H2O_S, PILISModel::B_LOSS_END, 0.9);
    model_.hmms_losses_[PILISModel::BIon].setTransitionProbability(PILISModel::B_H2O_T, PILISModel::B_H2O, 0.1);
    model_.hmms_losses_[PILISModel::BIon].setTransitionProbability(PILISModel::B_H2O_T, PILISModel::B_LOSS_END, 0.9);
		model_.hmms_losses_[PILISModel::BIon].setTransitionProbability(PILISModel::B_H2O_Q1, PILISModel::B_H2O, 0.1);
		model_.hmms_losses_[PILISModel::BIon].setTransitionProbability(PILISModel::B_H2O_Q1, PILISModel::B_LOSS_END, 0.9);

		model_.hmms_losses_[PILISModel::BIon].setTransitionProbability(PILISModel::B_NH3_K, PILISModel::B_NH3, 0.1);
    model_.hmms_losses_[PILISModel::BIon].setTransitionProbability(PILISModel::B_NH3_K, PILISModel::B_LOSS_END, 0.9);
    model_.hmms_losses_[PILISModel::BIon].setTransitionProbability(PILISModel::B_NH3_R, PILISModel::B_NH3, 0.1);
    model_.hmms_losses_[PILISModel::BIon].setTransitionProbability(PILISModel::B_NH3_R, PILISModel::B_LOSS_END, 0.9);
		model_.hmms_losses_[PILISModel::BIon].setTransitionProbability(PILISModel::B_NH3_Q, PILISModel::B_NH3, 0.1);
		model_.hmms_losses_[PILISModel::BIon].setTransitionProbability(PILISModel::B_NH3_Q, PILISModel::B_LOSS_END, 0.9);
		model_.hmms_losses_[PILISModel::BIon].setTransitionProbability(PILISModel::B_NH3_N, PILISModel::B_NH3, 0.1);
		model_.hmms_losses_[PILISModel::BIon].setTransitionProbability(PILISModel::B_NH3_N, PILISModel::B_LOSS_END, 0.9);
						
		
		model_.hmms_losses_[PILISModel::BIon].setTransitionProbability(PILISModel::B_CO, PILISModel::A_ION, 0.1);
		model_.hmms_losses_[PILISModel::BIon].setTransitionProbability(PILISModel::B_CO, PILISModel::B_LOSS_END, 0.9);

		model_.hmms_losses_[PILISModel::BIon].addSynonymTransition(PILISModel::B_BASE1, PILISModel::B_BASE2, PILISModel::B_ION, PILISModel::B_H2O_D);
    model_.hmms_losses_[PILISModel::BIon].addSynonymTransition(PILISModel::B_BASE1, PILISModel::B_BASE2, PILISModel::B_ION, PILISModel::B_H2O_E);
    model_.hmms_losses_[PILISModel::BIon].addSynonymTransition(PILISModel::B_BASE1, PILISModel::B_BASE2, PILISModel::B_ION, PILISModel::B_H2O_S);
    model_.hmms_losses_[PILISModel::BIon].addSynonymTransition(PILISModel::B_BASE1, PILISModel::B_BASE2, PILISModel::B_ION, PILISModel::B_H2O_T);
		model_.hmms_losses_[PILISModel::BIon].addSynonymTransition(PILISModel::B_BASE1, PILISModel::B_BASE2, PILISModel::B_ION, PILISModel::B_H2O_Q1);
    model_.hmms_losses_[PILISModel::BIon].addSynonymTransition(PILISModel::B_BASE1, PILISModel::B_BASE2, PILISModel::B_ION, PILISModel::B_NH3_K);
    model_.hmms_losses_[PILISModel::BIon].addSynonymTransition(PILISModel::B_BASE1, PILISModel::B_BASE2, PILISModel::B_ION, PILISModel::B_NH3_R);
		model_.hmms_losses_[PILISModel::BIon].addSynonymTransition(PILISModel::B_BASE1, PILISModel::B_BASE2, PILISModel::B_ION, PILISModel::B_NH3_Q);
		model_.hmms_losses_[PILISModel::BIon].addSynonymTransition(PILISModel::B_BASE1, PILISModel::B_BASE2, PILISModel::B_ION, PILISModel::B_NH3_N);
		model_.hmms_losses_[PILISModel::BIon].addSynonymTransition(PILISModel::B_BASE1, PILISModel::B_BASE2, PILISModel::B_ION, PILISModel::B_CO);

		model_.hmms_losses_[PILISModel::BIon].disableTransitions();
    model_.hmms_losses_[PILISModel::BIon].buildSynonyms();

		// b2-ions
    model_.hmms_losses_[PILISModel::B2Ion].addNewState(new HMMStateLight(PILISModel::B_H2O, false));
    model_.hmms_losses_[PILISModel::B2Ion].addNewState(new HMMStateLight(PILISModel::B_NH3, false));
    model_.hmms_losses_[PILISModel::B2Ion].addNewState(new HMMStateLight(PILISModel::A_ION, false));
    model_.hmms_losses_[PILISModel::B2Ion].addNewState(new HMMStateLight(PILISModel::B_LOSS_END, false));

    model_.hmms_losses_[PILISModel::B2Ion].addNewState(new HMMStateLight(PILISModel::B_ION));
    model_.hmms_losses_[PILISModel::B2Ion].addNewState(new HMMStateLight(PILISModel::B_BASE1));
    model_.hmms_losses_[PILISModel::B2Ion].addNewState(new HMMStateLight(PILISModel::B_BASE2));

    model_.hmms_losses_[PILISModel::B2Ion].addNewState(new HMMStateLight(PILISModel::B_H2O_D));
    model_.hmms_losses_[PILISModel::B2Ion].addNewState(new HMMStateLight(PILISModel::B_H2O_E));
    model_.hmms_losses_[PILISModel::B2Ion].addNewState(new HMMStateLight(PILISModel::B_H2O_S));
    model_.hmms_losses_[PILISModel::B2Ion].addNewState(new HMMStateLight(PILISModel::B_H2O_T));
		model_.hmms_losses_[PILISModel::B2Ion].addNewState(new HMMStateLight(PILISModel::B_H2O_Q1));
    model_.hmms_losses_[PILISModel::B2Ion].addNewState(new HMMStateLight(PILISModel::B_NH3_K));
    model_.hmms_losses_[PILISModel::B2Ion].addNewState(new HMMStateLight(PILISModel::B_NH3_R));
    model_.hmms_losses_[PILISModel::B2Ion].addNewState(new HMMStateLight(PILISModel::B_NH3_Q));
    model_.hmms_losses_[PILISModel::B2Ion].addNewState(new HMMStateLight(PILISModel::B_NH3_N));
    model_.hmms_losses_[PILISModel::B2Ion].addNewState(new HMMStateLight(PILISModel::B_CO));


    model_.hmms_losses_[PILISModel::B2Ion].setTransitionProbability(PILISModel::B_BASE1, PILISModel::B_BASE2, 1.0);

    model_.hmms_losses_[PILISModel::B2Ion].setTransitionProbability(PILISModel::B_H2O_D, PILISModel::B_H2O, 0.1);
    model_.hmms_losses_[PILISModel::B2Ion].setTransitionProbability(PILISModel::B_H2O_D, PILISModel::B_LOSS_END, 0.9);
    model_.hmms_losses_[PILISModel::B2Ion].setTransitionProbability(PILISModel::B_H2O_E, PILISModel::B_H2O, 0.1);
    model_.hmms_losses_[PILISModel::B2Ion].setTransitionProbability(PILISModel::B_H2O_E, PILISModel::B_LOSS_END, 0.9);
    model_.hmms_losses_[PILISModel::B2Ion].setTransitionProbability(PILISModel::B_H2O_S, PILISModel::B_H2O, 0.1);
    model_.hmms_losses_[PILISModel::B2Ion].setTransitionProbability(PILISModel::B_H2O_S, PILISModel::B_LOSS_END, 0.9);
    model_.hmms_losses_[PILISModel::B2Ion].setTransitionProbability(PILISModel::B_H2O_T, PILISModel::B_H2O, 0.1);
    model_.hmms_losses_[PILISModel::B2Ion].setTransitionProbability(PILISModel::B_H2O_T, PILISModel::B_LOSS_END, 0.9);
		model_.hmms_losses_[PILISModel::B2Ion].setTransitionProbability(PILISModel::B_H2O_Q1, PILISModel::B_H2O, 0.1);
		model_.hmms_losses_[PILISModel::B2Ion].setTransitionProbability(PILISModel::B_H2O_Q1, PILISModel::B_H2O, 0.1);

    model_.hmms_losses_[PILISModel::B2Ion].setTransitionProbability(PILISModel::B_NH3_K, PILISModel::B_NH3, 0.1);
    model_.hmms_losses_[PILISModel::B2Ion].setTransitionProbability(PILISModel::B_NH3_K, PILISModel::B_LOSS_END, 0.9);
    model_.hmms_losses_[PILISModel::B2Ion].setTransitionProbability(PILISModel::B_NH3_R, PILISModel::B_NH3, 0.1);
    model_.hmms_losses_[PILISModel::B2Ion].setTransitionProbability(PILISModel::B_NH3_R, PILISModel::B_LOSS_END, 0.9);
    model_.hmms_losses_[PILISModel::B2Ion].setTransitionProbability(PILISModel::B_NH3_Q, PILISModel::B_NH3, 0.1);
    model_.hmms_losses_[PILISModel::B2Ion].setTransitionProbability(PILISModel::B_NH3_Q, PILISModel::B_LOSS_END, 0.9);
    model_.hmms_losses_[PILISModel::B2Ion].setTransitionProbability(PILISModel::B_NH3_N, PILISModel::B_NH3, 0.1);
    model_.hmms_losses_[PILISModel::B2Ion].setTransitionProbability(PILISModel::B_NH3_N, PILISModel::B_LOSS_END, 0.9);


    model_.hmms_losses_[PILISModel::B2Ion].setTransitionProbability(PILISModel::B_CO, PILISModel::A_ION, 0.1);
    model_.hmms_losses_[PILISModel::B2Ion].setTransitionProbability(PILISModel::B_CO, PILISModel::B_LOSS_END, 0.9);

    model_.hmms_losses_[PILISModel::B2Ion].addSynonymTransition(PILISModel::B_BASE1, PILISModel::B_BASE2, PILISModel::B_ION, PILISModel::B_H2O_D);
    model_.hmms_losses_[PILISModel::B2Ion].addSynonymTransition(PILISModel::B_BASE1, PILISModel::B_BASE2, PILISModel::B_ION, PILISModel::B_H2O_E);
    model_.hmms_losses_[PILISModel::B2Ion].addSynonymTransition(PILISModel::B_BASE1, PILISModel::B_BASE2, PILISModel::B_ION, PILISModel::B_H2O_S);
    model_.hmms_losses_[PILISModel::B2Ion].addSynonymTransition(PILISModel::B_BASE1, PILISModel::B_BASE2, PILISModel::B_ION, PILISModel::B_H2O_T);
		model_.hmms_losses_[PILISModel::B2Ion].addSynonymTransition(PILISModel::B_BASE1, PILISModel::B_BASE2, PILISModel::B_ION, PILISModel::B_H2O_Q1);
    model_.hmms_losses_[PILISModel::B2Ion].addSynonymTransition(PILISModel::B_BASE1, PILISModel::B_BASE2, PILISModel::B_ION, PILISModel::B_NH3_K);
    model_.hmms_losses_[PILISModel::B2Ion].addSynonymTransition(PILISModel::B_BASE1, PILISModel::B_BASE2, PILISModel::B_ION, PILISModel::B_NH3_R);
    model_.hmms_losses_[PILISModel::B2Ion].addSynonymTransition(PILISModel::B_BASE1, PILISModel::B_BASE2, PILISModel::B_ION, PILISModel::B_NH3_Q);
    model_.hmms_losses_[PILISModel::B2Ion].addSynonymTransition(PILISModel::B_BASE1, PILISModel::B_BASE2, PILISModel::B_ION, PILISModel::B_NH3_N);
    model_.hmms_losses_[PILISModel::B2Ion].addSynonymTransition(PILISModel::B_BASE1, PILISModel::B_BASE2, PILISModel::B_ION, PILISModel::B_CO);

    model_.hmms_losses_[PILISModel::B2Ion].disableTransitions();
    model_.hmms_losses_[PILISModel::B2Ion].buildSynonyms();
		return;
	}
*/

	void PILISModelGenerator::initPrecursorModel_()
	{
    set<const Residue*> residues(ResidueDB::getInstance()->getResidues(ResidueDB::NATURAL_20));
    StringList modifications = param_.getValue("modifications");

    for (StringList::const_iterator it = modifications.begin(); it != modifications.end(); ++it)
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
			String loss = (*it)->getLossFormula().getString();

			if (loss != "")
			{
				losses.insert(loss);
#ifdef PRECURSOR_MODEL_DEBUG
				cerr << "Loss: " << loss << " (name=" << (*it)->getLossName() << ") of residue: " << (*it)->getName() << endl;
#endif
			}
		}

		// precursor is not fragmented
		model_.hmm_precursor_.addNewState(new HMMState("p", false));

		// emitting nodes for single losses
		for (set<String>::const_iterator it = losses.begin(); it != losses.end(); ++it)
		{
			model_.hmm_precursor_.addNewState(new HMMState("p-" + *it, false));
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
				model_.hmm_precursor_.addNewState(new HMMState("p-" + *it1 + "-" + *it2, false));
			}
		}
		

		// H2O loss from the C-terminus
		String h2o(EmpiricalFormula("H2O").getString());
		model_.hmm_precursor_.addNewState(new HMMState("COOH-" + h2o, true));
		model_.hmm_precursor_.addNewState(new HMMState("start", true));

		// TODO put this into a parameter
		UInt num_explicit(4);
		
		// add double loss states
		// add edges from double loss states to single loss states and emitting states
		for (set<const Residue*>::const_iterator it1 = residues.begin(); it1 != residues.end(); ++it1)
		{
			AASequence aa1;
			aa1 += *it1;
			String loss1 = (*it1)->getLossFormula().getString();
			if (loss1 == "")
			{
				continue;
			}
			
			model_.hmm_precursor_.addNewState(new HMMState(aa1.toString() + "-" + loss1));
			model_.hmm_precursor_.addNewState(new HMMState(aa1.toString() + "-" + loss1 + "-next"));

			for (UInt i = 0; i != num_explicit; ++i)
			{
				model_.hmm_precursor_.addNewState(new HMMState(aa1.toString() + "-" + loss1 + "_" + String(i + 1)));
			}
			
			model_.hmm_precursor_.addNewState(new HMMState(aa1.toString()  + "COOH-" + loss1 + "-" + h2o));
			model_.hmm_precursor_.addNewState(new HMMState(aa1.toString()  + "COOH-" + loss1 + "-" + h2o + "-next"));

			String losses;
			if (h2o < loss1)
			{
				losses = "-" + h2o  + "-" + loss1;
			}
			else
			{
				losses = "-" + loss1 + "-" + h2o;
			}
			String cooh_name = aa1.toString() + "COOH-" + loss1 + "-" + h2o;
			model_.hmm_precursor_.setTransitionProbability(cooh_name, "p" + losses, 0.25);
			model_.hmm_precursor_.setTransitionProbability(cooh_name, aa1.toString() + "-" + loss1, 0.25);
			model_.hmm_precursor_.setTransitionProbability(cooh_name, "COOH-" + h2o, 0.25);
			model_.hmm_precursor_.setTransitionProbability(cooh_name, cooh_name + "-next", 0.25);
		}
		for (set<const Residue*>::const_iterator it1 = residues.begin(); it1 != residues.end(); ++it1)
		{
			AASequence aa1;
			aa1 += *it1;
			String loss1 = (*it1)->getLossFormula().getString();
			if (loss1 == "")
			{
				continue;
			}
			
			for (set<const Residue*>::const_iterator it2 = it1; it2 != residues.end(); ++it2)
			{
				AASequence aa2;
				aa2 += *it2;
				String loss2 = (*it2)->getLossFormula().getString();
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
					model_.hmm_precursor_.addNewState(new HMMState(name + losses, true));
					model_.hmm_precursor_.addNewState(new HMMState(name + losses + "-next", true));

					for (UInt i = 0; i != num_explicit; ++i)
					{
						model_.hmm_precursor_.addNewState(new HMMState(name + losses + "_" + String(i + 1), true));
					}

					model_.hmm_precursor_.setTransitionProbability(name + losses, "p" + losses, 0.25);
					model_.hmm_precursor_.setTransitionProbability(name + losses, aa1.toString() + "-" + loss1, 0.25);
					model_.hmm_precursor_.setTransitionProbability(name + losses, aa2.toString() + "-" + loss2, 0.25);
					model_.hmm_precursor_.setTransitionProbability(name + losses, name + losses + "-next", 0.25);

					for (UInt i = 0; i != num_explicit; ++i)
					{
						String state_name_num = name + losses + "_" + String(i + 1);
						String state_name = name + losses;
						model_.hmm_precursor_.addSynonymTransition(state_name, "p" + losses,  state_name_num, "p" + losses);
						model_.hmm_precursor_.addSynonymTransition(state_name, aa1.toString() + "-" + loss1, state_name_num, aa1.toString() + "-" + loss1 + "_" + String(i + 1));
						model_.hmm_precursor_.addSynonymTransition(state_name, aa2.toString() + "-" + loss2, state_name_num, aa2.toString() + "-" + loss2 + "_" + String(i + 1));
					}
				}
			}
		}

		for (set<const Residue*>::const_iterator it = residues.begin(); it != residues.end(); ++it)
		{
			AASequence aa;
      aa += *it;
      String loss = (*it)->getLossFormula().getString();
      if (loss == "")
      {
        continue;
      }

			model_.hmm_precursor_.setTransitionProbability(aa.toString() + "-" + loss, "p-" + loss, 0.25);
			model_.hmm_precursor_.setTransitionProbability(aa.toString() + "-" + loss, "p", 0.25);
			model_.hmm_precursor_.setTransitionProbability(aa.toString() + "-" + loss, aa.toString() + "-" + loss + "-next", 0.5);

			for (UInt i = 0; i != num_explicit; ++i)
			{
				String name_num = aa.toString() + "-" + loss + "_" + String(i + 1);
				String name = aa.toString() + "-" + loss;
				model_.hmm_precursor_.addSynonymTransition(name, "p-" + loss, name_num, "p-" + loss);
				model_.hmm_precursor_.addSynonymTransition(name, "p", name_num, "p");
			}
		}

		model_.hmm_precursor_.disableTransitions();
		model_.hmm_precursor_.buildSynonyms();


#ifdef PRECURSOR_MODEL_DEBUG
    cerr << "#States: " << model_.hmm_precursor_.getNumberOfStates() << endl;
    model_.hmm_precursor_.writeGraphMLFile("precursor_model.graphML");
#endif
		
		return;
	}

} // namespace OpenMS

