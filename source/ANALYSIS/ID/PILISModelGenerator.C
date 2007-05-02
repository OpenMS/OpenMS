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


#include <OpenMS/ANALYSIS/ID/PILISModelGenerator.h>

using namespace std;

namespace OpenMS 
{
	PILISModelGenerator::PILISModelGenerator()
		: DefaultParamHandler("PILISModelGenerator")
	{
		defaults_.setValue("model_depth", 4);
		defaults_.setValue("visible_model_depth", 30);
		defaultsToParam_();
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
		initLossModels_();
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
	
		const String residues("ACDEFGHIKLMNPQRSTVWY");
	
		// create the residue states states
		for (UInt ii = 0; ii != residues.size(); ++ii)
    {
      String first;
      first += residues[ii];
      model_.hmm_.addNewState(first+"D");
			model_.hmm_.addNewState(first+"E");

			model_.hmm_.addNewState(first+"K");
			model_.hmm_.addNewState(first+"H");
			model_.hmm_.addNewState(first+"RSC");
			
      for (UInt j = 0; j != residues.size(); ++j)
      {
        String second;
        second += residues[j];
        model_.hmm_.addNewState(first+second+"bxyz");
				model_.hmm_.addNewState(first+second+"axyz");
      }
    }

		for (UInt i = 1; i <= visible_model_depth; ++i)
		{
			// these states are really created
			// charge states
			String num(i);
			model_.hmm_.addNewState("BB"+num);
			model_.hmm_.addNewState("BBk-"+num);
			model_.hmm_.addNewState("CR"+num); 
			model_.hmm_.addNewState("CRk-"+num);

			model_.hmm_.addNewState("SC"+num);
			model_.hmm_.addNewState("SCk-"+num);

			// states for trans mapping 
			model_.hmm_.addNewState("AA"+num);
			model_.hmm_.addNewState("AAk-"+num);

			model_.hmm_.addNewState("A"+num);
			model_.hmm_.addNewState("Ak-"+num);

			model_.hmm_.addNewState("ASC"+num);
			model_.hmm_.addNewState("ASCk-"+num);

			// emitting ion states
			model_.hmm_.addNewState(new HMMState("b"+num+"+", false));
			model_.hmm_.addNewState(new HMMState("bk-"+num+"+", false));
			model_.hmm_.addNewState(new HMMState("y"+num+"+", false));
			model_.hmm_.addNewState(new HMMState("yk-"+num+"+", false));
			model_.hmm_.addNewState(new HMMState("a"+num+"+", false));
			model_.hmm_.addNewState(new HMMState("ak-"+num+"+", false));

			model_.hmm_.addNewState(new HMMState("b"+num+"++", false));
      model_.hmm_.addNewState(new HMMState("bk-"+num+"++", false));
      model_.hmm_.addNewState(new HMMState("y"+num+"++", false));
      model_.hmm_.addNewState(new HMMState("yk-"+num+"++", false));
      model_.hmm_.addNewState(new HMMState("a"+num+"++", false));
      model_.hmm_.addNewState(new HMMState("ak-"+num+"++", false));

			model_.hmm_.addNewState(new HMMState("end"+num, false));
			model_.hmm_.addNewState(new HMMState("endk-"+num, false));

			// emitting neutral loss states
		
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
			for (UInt ii = 0; ii != residues.size(); ++ii)
			{
				String first;
				first += residues[ii];
				model_.hmm_.addNewState(first+"D"+num);
				model_.hmm_.addNewState(first+"Dk-"+num);

        model_.hmm_.addNewState(first+"E"+num);
        model_.hmm_.addNewState(first+"Ek-"+num);

				model_.hmm_.addNewState(first+"K"+num);
				model_.hmm_.addNewState(first+"Kk-"+num);

				model_.hmm_.addNewState(first+"H"+num);
        model_.hmm_.addNewState(first+"Hk-"+num);

				model_.hmm_.addNewState(first+"RSC"+num);
        model_.hmm_.addNewState(first+"RSCk-"+num);
				
				for (UInt j = 0; j != residues.size(); ++j)
				{
					String second;
					second += residues[j];
					model_.hmm_.addNewState(first+second+"bxyz"+num);
					model_.hmm_.addNewState(first+second+"bxyzk-"+num);
					model_.hmm_.addNewState(first+second+"axyz"+num);
					model_.hmm_.addNewState(first+second+"axyzk-"+num);
				}
			}
		}

		model_.hmm_.setTransitionProbability("AABase1", "AABase2", 1);

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
			
			for (UInt ii = 0; ii != residues.size(); ++ii)
			{
				String first;
				first += residues[ii];

				// CR D
				model_.hmm_.addSynonymTransition("AABase1", "AABase2", "A"+num, first+"D"+num);
				model_.hmm_.addSynonymTransition("AABase1", "AABase2", "Ak-"+num, first+"Dk-"+num);
			
				model_.hmm_.addSynonymTransition(first+"D", "D", first+"D"+num, "D"+num);
				model_.hmm_.addSynonymTransition(first+"D", "end", first+"D"+num, "end"+num);
				model_.hmm_.addSynonymTransition(first+"D", "D", first+"Dk-"+num, "Dk-"+num);
				model_.hmm_.addSynonymTransition(first+"D", "end", first+"Dk-"+num, "endk-"+num);

				model_.hmm_.setTransitionProbability(first+"D", "D", 0.5);
				model_.hmm_.setTransitionProbability(first+"D", "end", 0.5);
		
				// CR E
				model_.hmm_.addSynonymTransition("AABase1", "AABase2", "A"+num, first+"E"+num);
        model_.hmm_.addSynonymTransition("AABase1", "AABase2", "Ak-"+num, first+"Ek-"+num);
        
        model_.hmm_.addSynonymTransition(first+"E", "E", first+"E"+num, "E"+num);
				model_.hmm_.addSynonymTransition(first+"E", "end", first+"E"+num, "end"+num);
        model_.hmm_.addSynonymTransition(first+"E", "E", first+"Ek-"+num, "Ek-"+num);
				model_.hmm_.addSynonymTransition(first+"E", "end", first+"Ek-"+num, "endk-"+num);
        
        model_.hmm_.setTransitionProbability(first+"E", "E", 0.5);
       	model_.hmm_.setTransitionProbability(first+"E", "end", 0.5);
		
				// SC K
				model_.hmm_.addSynonymTransition("AABase1", "AABase2", "ASC"+num, first+"K"+num);
        model_.hmm_.addSynonymTransition("AABase1", "AABase2", "ASCk-"+num, first+"Kk-"+num);

        model_.hmm_.addSynonymTransition(first+"K", "K", first+"K"+num, "K"+num);
				model_.hmm_.addSynonymTransition(first+"K", "end", first+"K"+num, "end"+num);
        model_.hmm_.addSynonymTransition(first+"K", "K", first+"Kk-"+num, "Kk-"+num);
				model_.hmm_.addSynonymTransition(first+"K", "end", first+"Kk-"+num, "endk-"+num);

        model_.hmm_.setTransitionProbability(first+"K", "K", 0.5);
        model_.hmm_.setTransitionProbability(first+"K", "end", 0.5);
			
				// SC H
				model_.hmm_.addSynonymTransition("AABase1", "AABase2", "ASC"+num, first+"H"+num);
        model_.hmm_.addSynonymTransition("AABase1", "AABase2", "ASCk-"+num, first+"Hk-"+num);

        model_.hmm_.addSynonymTransition(first+"H", "H", first+"H"+num, "H"+num);
				model_.hmm_.addSynonymTransition(first+"H", "end", first+"H"+num, "end"+num);
        model_.hmm_.addSynonymTransition(first+"H", "H", first+"Hk-"+num, "Hk-"+num);
				model_.hmm_.addSynonymTransition(first+"H", "end", first+"Hk-"+num, "endk-"+num);

        model_.hmm_.setTransitionProbability(first+"H", "H", 0.5);
        model_.hmm_.setTransitionProbability(first+"H", "end", 0.5);

				// SC R	
				model_.hmm_.addSynonymTransition("AABase1", "AABase2", "ASC"+num, first+"RSC"+num);
        model_.hmm_.addSynonymTransition("AABase1", "AABase2", "ASCk-"+num, first+"RSCk-"+num);

        model_.hmm_.addSynonymTransition(first+"RSC", "R", first+"RSC"+num, "R"+num);
				model_.hmm_.addSynonymTransition(first+"RSC", "end", first+"RSC"+num, "end"+num);
        model_.hmm_.addSynonymTransition(first+"RSC", "R", first+"RSCk-"+num, "Rk-"+num);
				model_.hmm_.addSynonymTransition(first+"RSC", "end", first+"RSCk-"+num, "endk-"+num);

        model_.hmm_.setTransitionProbability(first+"RSC", "R", 0.5);
        model_.hmm_.setTransitionProbability(first+"RSC", "end", 0.5);
				

				for (UInt j = 0; j != residues.size(); ++j)
				{
					String second;
					second += residues[j];

					model_.hmm_.addSynonymTransition("AABase1", "AABase2", "AA"+num, first+second+"bxyz"+num);
					model_.hmm_.addSynonymTransition("AABase1", "AABase2", "AA"+num, first+second+"axyz"+num);

					model_.hmm_.addSynonymTransition("AABase1", "AABase2", "AAk-"+num, first+second+"bxyzk-"+num);
					model_.hmm_.addSynonymTransition("AABase1", "AABase2", "AAk-"+num, first+second+"axyzk-"+num);
				
					if (/*(second == "P")&&*/ i <= 2)
					{
						if (second == "P")
						{
            	model_.hmm_.setTransitionProbability(first+second+"bxyz"+num, "bxyz"+num, 0.01);
							model_.hmm_.addSynonymTransition(first+second+"bxyz"+num, "end", first+second+"bxyz"+num, "end"+num);
							model_.hmm_.setTransitionProbability(first+second+"bxyz"+num, "end", 0.99);

            	model_.hmm_.setTransitionProbability(first+second+"axyz"+num, "axyz"+num, 0.01);
							model_.hmm_.addSynonymTransition(first+second+"axyz"+num, "end", first+second+"axyz"+num, "end"+num);
							model_.hmm_.setTransitionProbability(first+second+"axyz"+num, "end", 0.99);
						}
						else
						{
							model_.hmm_.setTransitionProbability(first+second+"bxyz"+num, "bxyz"+num, 0.5);
              model_.hmm_.addSynonymTransition(first+second+"bxyz"+num, "end", first+second+"bxyz"+num, "end"+num);
              model_.hmm_.setTransitionProbability(first+second+"bxyz"+num, "end", 0.5);

							model_.hmm_.setTransitionProbability(first+second+"axyz"+num, "axyz"+num, 0.01);
              model_.hmm_.addSynonymTransition(first+second+"axyz"+num, "end", first+second+"axyz"+num, "end"+num);
              model_.hmm_.setTransitionProbability(first+second+"axyz"+num, "end", 0.99);
						}
					}
					else
					{
					
						model_.hmm_.addSynonymTransition(first+second+"bxyz", "bxyz", first+second+"bxyz"+num, "bxyz"+num);
						model_.hmm_.addSynonymTransition(first+second+"bxyz", "end", first+second+"bxyz"+num, "end"+num);

						model_.hmm_.addSynonymTransition(first+second+"axyz", "axyz", first+second+"axyz"+num, "axyz"+num);
						model_.hmm_.addSynonymTransition(first+second+"axyz", "end", first+second+"axyz"+num, "end"+num); 

						model_.hmm_.addSynonymTransition(first+second+"bxyz", "end", first+second+"bxyz"+num, "end"+num);
						model_.hmm_.addSynonymTransition(first+second+"axyz", "end", first+second+"axyz"+num, "end"+num);
					}
					
					model_.hmm_.addSynonymTransition(first+second+"bxyz", "bxyz", first+second+"bxyzk-"+num, "bxyzk-"+num);
					model_.hmm_.addSynonymTransition(first+second+"bxyz", "end", first+second+"bxyzk-"+num, "endk-"+num);
					
					model_.hmm_.addSynonymTransition(first+second+"axyz", "axyz", first+second+"axyzk-"+num, "axyzk-"+num);
					model_.hmm_.addSynonymTransition(first+second+"axyz", "end", first+second+"axyzk-"+num, "endk-"+num);
					
					model_.hmm_.setTransitionProbability(first+second+"bxyz", "bxyz", 0.5);
					model_.hmm_.setTransitionProbability(first+second+"axyz", "axyz", 0.5);

					model_.hmm_.setTransitionProbability(first+second+"bxyz", "end", 0.5);
					model_.hmm_.setTransitionProbability(first+second+"axyz", "end", 0.5);
					model_.hmm_.addSynonymTransition(first+second+"bxyz", "end", first+second+"bxyzk-"+num, "end"+num);
          model_.hmm_.addSynonymTransition(first+second+"axyz", "end", first+second+"axyzk-"+num, "end"+num);
				}
			}
		}

		model_.hmm_.disableTransitions();
		model_.hmm_.buildSynonyms();
	}

	void PILISModelGenerator::initLossModels_()
	{
		//# new states from this HMM
		//# format State <Name> [<Hidden?>]
		//# emitting states
		//State y-H2O false
		//State y-NH3 false
		//State y_loss_end false
		model_.hmms_losses_[Residue::YIon].addNewState(new HMMStateLight(PILISModel::Y_H2O, false));
		model_.hmms_losses_[Residue::YIon].addNewState(new HMMStateLight(PILISModel::Y_NH3, false));
		model_.hmms_losses_[Residue::YIon].addNewState(new HMMStateLight(PILISModel::Y_LOSS_END, false));

		//# base states
		//State yion
		//State y_Base1
		//State y_Base2
		model_.hmms_losses_[Residue::YIon].addNewState(new HMMStateLight(PILISModel::Y_ION));
		model_.hmms_losses_[Residue::YIon].addNewState(new HMMStateLight(PILISModel::Y_BASE1));
		model_.hmms_losses_[Residue::YIon].addNewState(new HMMStateLight(PILISModel::Y_BASE2));

		//# loss pathway states
		//State y_H2O_D
		//State y_H2O_E
		//State y_H2O_T
		//State y_H2O_S
		//State y_H2O_Q1
		//State y_H2O_Cterm
		//State y_NH3_K
		//State y_NH3_R
		model_.hmms_losses_[Residue::YIon].addNewState(new HMMStateLight(PILISModel::Y_H2O_D));
		model_.hmms_losses_[Residue::YIon].addNewState(new HMMStateLight(PILISModel::Y_H2O_E));
		model_.hmms_losses_[Residue::YIon].addNewState(new HMMStateLight(PILISModel::Y_H2O_S));
		model_.hmms_losses_[Residue::YIon].addNewState(new HMMStateLight(PILISModel::Y_H2O_T));
		model_.hmms_losses_[Residue::YIon].addNewState(new HMMStateLight(PILISModel::Y_H2O_Q1));
		model_.hmms_losses_[Residue::YIon].addNewState(new HMMStateLight(PILISModel::Y_H2O_CTERM));
		model_.hmms_losses_[Residue::YIon].addNewState(new HMMStateLight(PILISModel::Y_NH3_K));
		model_.hmms_losses_[Residue::YIon].addNewState(new HMMStateLight(PILISModel::Y_NH3_R));
		
		//# synonym states are not needed in this state
		//# transitions
		//Transition y_Base1 y_Base2 1
		model_.hmms_losses_[Residue::YIon].setTransitionProbability(PILISModel::Y_BASE1, PILISModel::Y_BASE2, 1.0);

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
		model_.hmms_losses_[Residue::YIon].setTransitionProbability(PILISModel::Y_H2O_D, PILISModel::Y_H2O, 0.1);
		model_.hmms_losses_[Residue::YIon].setTransitionProbability(PILISModel::Y_H2O_D, PILISModel::Y_LOSS_END, 0.9);
		model_.hmms_losses_[Residue::YIon].setTransitionProbability(PILISModel::Y_H2O_E, PILISModel::Y_H2O, 0.1);
		model_.hmms_losses_[Residue::YIon].setTransitionProbability(PILISModel::Y_H2O_E, PILISModel::Y_LOSS_END, 0.9);
		model_.hmms_losses_[Residue::YIon].setTransitionProbability(PILISModel::Y_H2O_S, PILISModel::Y_H2O, 0.1);
		model_.hmms_losses_[Residue::YIon].setTransitionProbability(PILISModel::Y_H2O_S, PILISModel::Y_LOSS_END, 0.9);
		model_.hmms_losses_[Residue::YIon].setTransitionProbability(PILISModel::Y_H2O_T, PILISModel::Y_H2O, 0.1);
		model_.hmms_losses_[Residue::YIon].setTransitionProbability(PILISModel::Y_H2O_T, PILISModel::Y_LOSS_END, 0.9);
		model_.hmms_losses_[Residue::YIon].setTransitionProbability(PILISModel::Y_H2O_Q1, PILISModel::Y_H2O, 0.1);
		model_.hmms_losses_[Residue::YIon].setTransitionProbability(PILISModel::Y_H2O_Q1, PILISModel::Y_LOSS_END, 0.9);
		model_.hmms_losses_[Residue::YIon].setTransitionProbability(PILISModel::Y_H2O_CTERM, PILISModel::Y_H2O, 0.1);
		model_.hmms_losses_[Residue::YIon].setTransitionProbability(PILISModel::Y_H2O_CTERM, PILISModel::Y_LOSS_END, 0.9);

		//# NH3 losses
		//Transition y_NH3_K y-NH3 0.1
		//Transition y_NH3_K y_loss_end 0.9
		//Transition y_NH3_R y-NH3 0.1
		//Transition y_NH3_R y_loss_end 0.9
		model_.hmms_losses_[Residue::YIon].setTransitionProbability(PILISModel::Y_NH3_K, PILISModel::Y_NH3, 0.1);
    model_.hmms_losses_[Residue::YIon].setTransitionProbability(PILISModel::Y_NH3_K, PILISModel::Y_LOSS_END, 0.9);
    model_.hmms_losses_[Residue::YIon].setTransitionProbability(PILISModel::Y_NH3_R, PILISModel::Y_NH3, 0.1);
    model_.hmms_losses_[Residue::YIon].setTransitionProbability(PILISModel::Y_NH3_R, PILISModel::Y_LOSS_END, 0.9);
	
		//# synonym transitions
		//Synonym yion y_H2O_D y_Base1 y_Base2
		//Synonym yion y_H2O_E y_Base1 y_Base2
		//Synonym yion y_H2O_T y_Base1 y_Base2
		//Synonym yion y_H2O_S y_Base1 y_Base2
		//Synonym yion y_H2O_Q1 y_Base1 y_Base2
		//Synonym yion y_H2O_Cterm y_Base1 y_Base2
		//Synonym yion y_NH3_K y_Base1 y_Base2
		//Synonym yion y_NH3_R y_Base1 y_Base2
		model_.hmms_losses_[Residue::YIon].addSynonymTransition(PILISModel::Y_BASE1, PILISModel::Y_BASE2, PILISModel::Y_ION, PILISModel::Y_H2O_D);
    model_.hmms_losses_[Residue::YIon].addSynonymTransition(PILISModel::Y_BASE1, PILISModel::Y_BASE2, PILISModel::Y_ION, PILISModel::Y_H2O_E);
    model_.hmms_losses_[Residue::YIon].addSynonymTransition(PILISModel::Y_BASE1, PILISModel::Y_BASE2, PILISModel::Y_ION, PILISModel::Y_H2O_S);
    model_.hmms_losses_[Residue::YIon].addSynonymTransition(PILISModel::Y_BASE1, PILISModel::Y_BASE2, PILISModel::Y_ION, PILISModel::Y_H2O_T);
    model_.hmms_losses_[Residue::YIon].addSynonymTransition(PILISModel::Y_BASE1, PILISModel::Y_BASE2, PILISModel::Y_ION, PILISModel::Y_H2O_Q1);
    model_.hmms_losses_[Residue::YIon].addSynonymTransition(PILISModel::Y_BASE1, PILISModel::Y_BASE2, PILISModel::Y_ION, PILISModel::Y_H2O_CTERM);
    model_.hmms_losses_[Residue::YIon].addSynonymTransition(PILISModel::Y_BASE1, PILISModel::Y_BASE2, PILISModel::Y_ION, PILISModel::Y_NH3_K);
    model_.hmms_losses_[Residue::YIon].addSynonymTransition(PILISModel::Y_BASE1, PILISModel::Y_BASE2, PILISModel::Y_ION, PILISModel::Y_NH3_R);

		model_.hmms_losses_[Residue::YIon].disableTransitions();
		model_.hmms_losses_[Residue::YIon].buildSynonyms();


		// b-ions
		model_.hmms_losses_[Residue::BIon].addNewState(new HMMStateLight(PILISModel::B_H2O, false));
    model_.hmms_losses_[Residue::BIon].addNewState(new HMMStateLight(PILISModel::B_NH3, false));
    model_.hmms_losses_[Residue::BIon].addNewState(new HMMStateLight(PILISModel::B_LOSS_END, false));

    model_.hmms_losses_[Residue::BIon].addNewState(new HMMStateLight(PILISModel::B_ION));
    model_.hmms_losses_[Residue::BIon].addNewState(new HMMStateLight(PILISModel::B_BASE1));
    model_.hmms_losses_[Residue::BIon].addNewState(new HMMStateLight(PILISModel::B_BASE2));

    model_.hmms_losses_[Residue::BIon].addNewState(new HMMStateLight(PILISModel::B_H2O_D));
    model_.hmms_losses_[Residue::BIon].addNewState(new HMMStateLight(PILISModel::B_H2O_E));
    model_.hmms_losses_[Residue::BIon].addNewState(new HMMStateLight(PILISModel::B_H2O_S));
    model_.hmms_losses_[Residue::BIon].addNewState(new HMMStateLight(PILISModel::B_H2O_T));
    model_.hmms_losses_[Residue::BIon].addNewState(new HMMStateLight(PILISModel::B_NH3_K));
    model_.hmms_losses_[Residue::BIon].addNewState(new HMMStateLight(PILISModel::B_NH3_R));

    model_.hmms_losses_[Residue::BIon].setTransitionProbability(PILISModel::B_BASE1, PILISModel::B_BASE2, 1.0);

    model_.hmms_losses_[Residue::BIon].setTransitionProbability(PILISModel::B_H2O_D, PILISModel::B_H2O, 0.1);
    model_.hmms_losses_[Residue::BIon].setTransitionProbability(PILISModel::B_H2O_D, PILISModel::B_LOSS_END, 0.9);
    model_.hmms_losses_[Residue::BIon].setTransitionProbability(PILISModel::B_H2O_E, PILISModel::B_H2O, 0.1);
    model_.hmms_losses_[Residue::BIon].setTransitionProbability(PILISModel::B_H2O_E, PILISModel::B_LOSS_END, 0.9);
    model_.hmms_losses_[Residue::BIon].setTransitionProbability(PILISModel::B_H2O_S, PILISModel::B_H2O, 0.1);
    model_.hmms_losses_[Residue::BIon].setTransitionProbability(PILISModel::B_H2O_S, PILISModel::B_LOSS_END, 0.9);
    model_.hmms_losses_[Residue::BIon].setTransitionProbability(PILISModel::B_H2O_T, PILISModel::B_H2O, 0.1);
    model_.hmms_losses_[Residue::BIon].setTransitionProbability(PILISModel::B_H2O_T, PILISModel::B_LOSS_END, 0.9);

		model_.hmms_losses_[Residue::BIon].setTransitionProbability(PILISModel::B_NH3_K, PILISModel::B_NH3, 0.1);
    model_.hmms_losses_[Residue::BIon].setTransitionProbability(PILISModel::B_NH3_K, PILISModel::B_LOSS_END, 0.9);
    model_.hmms_losses_[Residue::BIon].setTransitionProbability(PILISModel::B_NH3_R, PILISModel::B_NH3, 0.1);
    model_.hmms_losses_[Residue::BIon].setTransitionProbability(PILISModel::B_NH3_R, PILISModel::B_LOSS_END, 0.9);

		model_.hmms_losses_[Residue::BIon].addSynonymTransition(PILISModel::B_BASE1, PILISModel::B_BASE2, PILISModel::B_ION, PILISModel::B_H2O_D);
    model_.hmms_losses_[Residue::BIon].addSynonymTransition(PILISModel::B_BASE1, PILISModel::B_BASE2, PILISModel::B_ION, PILISModel::B_H2O_E);
    model_.hmms_losses_[Residue::BIon].addSynonymTransition(PILISModel::B_BASE1, PILISModel::B_BASE2, PILISModel::B_ION, PILISModel::B_H2O_S);
    model_.hmms_losses_[Residue::BIon].addSynonymTransition(PILISModel::B_BASE1, PILISModel::B_BASE2, PILISModel::B_ION, PILISModel::B_H2O_T);
    model_.hmms_losses_[Residue::BIon].addSynonymTransition(PILISModel::B_BASE1, PILISModel::B_BASE2, PILISModel::B_ION, PILISModel::B_NH3_K);
    model_.hmms_losses_[Residue::BIon].addSynonymTransition(PILISModel::B_BASE1, PILISModel::B_BASE2, PILISModel::B_ION, PILISModel::B_NH3_R);

		model_.hmms_losses_[Residue::BIon].disableTransitions();
    model_.hmms_losses_[Residue::BIon].buildSynonyms();

		return;
	}

	void PILISModelGenerator::initPrecursorModel_()
	{
		//# new states from this HMM
		//# format State <Name> [<Hidden?>]
		//# emitting states
		//State [M+H] false
		//State [M+H]-NH3 false
		//State [M+H]-H2O false
		//State Preend false
		model_.hmm_precursor_.addNewState(new HMMStateLight(PILISModel::PRE_MH, false));
		model_.hmm_precursor_.addNewState(new HMMStateLight(PILISModel::PRE_MH_NH3, false));
		model_.hmm_precursor_.addNewState(new HMMStateLight(PILISModel::PRE_MH_H2O, false));
		model_.hmm_precursor_.addNewState(new HMMStateLight(PILISModel::PRE_END, false));

		//# base states
		//State Pre
		//State Pre_Base1
		//State Pre_Base2
		model_.hmm_precursor_.addNewState(new HMMStateLight(PILISModel::PRE_ION));
		model_.hmm_precursor_.addNewState(new HMMStateLight(PILISModel::PRE_BASE1));
		model_.hmm_precursor_.addNewState(new HMMStateLight(PILISModel::PRE_BASE2));

		//# loss pathway states
		//State Pre_H2O_D
		//State Pre_H2O_E
		//State Pre_H2O_T
		//State Pre_H2O_S
		//State Pre_H2O_Q1
		//State Pre_H2O_Cterm
		//State Pre_NH3_K
		//State Pre_NH3_R
		model_.hmm_precursor_.addNewState(new HMMStateLight(PILISModel::PRE_H2O_S));
		model_.hmm_precursor_.addNewState(new HMMStateLight(PILISModel::PRE_H2O_T));
		model_.hmm_precursor_.addNewState(new HMMStateLight(PILISModel::PRE_H2O_E));
		model_.hmm_precursor_.addNewState(new HMMStateLight(PILISModel::PRE_H2O_D));
		model_.hmm_precursor_.addNewState(new HMMStateLight(PILISModel::PRE_H2O_Q1));
		model_.hmm_precursor_.addNewState(new HMMStateLight(PILISModel::PRE_H2O_CTERM));
		model_.hmm_precursor_.addNewState(new HMMStateLight(PILISModel::PRE_NH3_K));
		model_.hmm_precursor_.addNewState(new HMMStateLight(PILISModel::PRE_NH3_R));

		//# synonym states are not needed in this state
		//# transitions
		//Transition Pre_Base1 Pre_Base2 1
		model_.hmm_precursor_.setTransitionProbability(PILISModel::PRE_BASE1, PILISModel::PRE_BASE2, 1);

		//# H2O losses
		//Transition Pre_H2O_D [M+H]-H2O 0.1
		//Transition Pre_H2O_D [M+H] 0.1
		//Transition Pre_H2O_D Preend 0.8
		//Transition Pre_H2O_E [M+H]-H2O 0.1
		//Transition Pre_H2O_E [M+H] 0.1
		//Transition Pre_H2O_E Preend 0.8
		//Transition Pre_H2O_T [M+H]-H2O 0.1
		//Transition Pre_H2O_T [M+H] 0.1
		//Transition Pre_H2O_T Preend 0.8
		//Transition Pre_H2O_S [M+H]-H2O 0.1
		//Transition Pre_H2O_S [M+H] 0.1
		//Transition Pre_H2O_S Preend 0.8
		//Transition Pre_H2O_Q1 [M+H]-H2O 0.1
		//Transition Pre_H2O_Q1 [M+H] 0.1
		//Transition Pre_H2O_Q1 Preend 0.8
		//Transition Pre_H2O_Cterm [M+H]-H2O 0.1
		//Transition Pre_H2O_Cterm [M+H] 0.1
		//Transition Pre_H2O_Cterm Preend 0.8
		model_.hmm_precursor_.setTransitionProbability(PILISModel::PRE_H2O_D, PILISModel::PRE_MH_H2O, 0.1);
		model_.hmm_precursor_.setTransitionProbability(PILISModel::PRE_H2O_D, PILISModel::PRE_MH, 0.1);
		model_.hmm_precursor_.setTransitionProbability(PILISModel::PRE_H2O_D, PILISModel::PRE_END, 0.8);

		model_.hmm_precursor_.setTransitionProbability(PILISModel::PRE_H2O_E, PILISModel::PRE_MH_H2O, 0.1);
		model_.hmm_precursor_.setTransitionProbability(PILISModel::PRE_H2O_E, PILISModel::PRE_MH, 0.1);
		model_.hmm_precursor_.setTransitionProbability(PILISModel::PRE_H2O_E, PILISModel::PRE_END, 0.8);

		model_.hmm_precursor_.setTransitionProbability(PILISModel::PRE_H2O_S, PILISModel::PRE_MH_H2O, 0.1);
		model_.hmm_precursor_.setTransitionProbability(PILISModel::PRE_H2O_S, PILISModel::PRE_MH, 0.1);
		model_.hmm_precursor_.setTransitionProbability(PILISModel::PRE_H2O_S, PILISModel::PRE_END, 0.8);

		model_.hmm_precursor_.setTransitionProbability(PILISModel::PRE_H2O_T, PILISModel::PRE_MH_H2O, 0.1);
		model_.hmm_precursor_.setTransitionProbability(PILISModel::PRE_H2O_T, PILISModel::PRE_MH, 0.1);
		model_.hmm_precursor_.setTransitionProbability(PILISModel::PRE_H2O_T, PILISModel::PRE_END, 0.8);

		model_.hmm_precursor_.setTransitionProbability(PILISModel::PRE_H2O_Q1, PILISModel::PRE_MH_H2O, 0.1);
		model_.hmm_precursor_.setTransitionProbability(PILISModel::PRE_H2O_Q1, PILISModel::PRE_MH, 0.1);
		model_.hmm_precursor_.setTransitionProbability(PILISModel::PRE_H2O_Q1, PILISModel::PRE_END, 0.8);

		model_.hmm_precursor_.setTransitionProbability(PILISModel::PRE_H2O_CTERM, PILISModel::PRE_MH_H2O, 0.1);
		model_.hmm_precursor_.setTransitionProbability(PILISModel::PRE_H2O_CTERM, PILISModel::PRE_MH, 0.1);
		model_.hmm_precursor_.setTransitionProbability(PILISModel::PRE_H2O_CTERM, PILISModel::PRE_END, 0.8);

		//# NH3 losses
		//Transition Pre_NH3_K [M+H]-NH3 0.1
		//Transition Pre_NH3_K [M+H] 0.1
		//Transition Pre_NH3_K Preend 0.8
		//Transition Pre_NH3_R [M+H]-NH3 0.1
		//Transition Pre_NH3_R [M+H] 0.1
		//Transition Pre_NH3_R Preend 0.8
		model_.hmm_precursor_.setTransitionProbability(PILISModel::PRE_NH3_K, PILISModel::PRE_MH_NH3, 0.1);
		model_.hmm_precursor_.setTransitionProbability(PILISModel::PRE_NH3_K, PILISModel::PRE_MH, 0.1);
		model_.hmm_precursor_.setTransitionProbability(PILISModel::PRE_NH3_K, PILISModel::PRE_END, 0.8);

		model_.hmm_precursor_.setTransitionProbability(PILISModel::PRE_NH3_R, PILISModel::PRE_MH_NH3, 0.1);
		model_.hmm_precursor_.setTransitionProbability(PILISModel::PRE_NH3_R, PILISModel::PRE_MH, 0.1);
		model_.hmm_precursor_.setTransitionProbability(PILISModel::PRE_NH3_R, PILISModel::PRE_END, 0.8);

		//# synonym transitions
		//Synonym Pre Pre_H2O_D Pre_Base1 Pre_Base2
		//Synonym Pre Pre_H2O_E Pre_Base1 Pre_Base2
		//Synonym Pre Pre_H2O_T Pre_Base1 Pre_Base2
		//Synonym Pre Pre_H2O_S Pre_Base1 Pre_Base2
		//Synonym Pre Pre_H2O_Q1 Pre_Base1 Pre_Base2
		//Synonym Pre Pre_H2O_Cterm Pre_Base1 Pre_Base2
		//Synonym Pre Pre_NH3_K Pre_Base1 Pre_Base2
		//Synonym Pre Pre_NH3_R Pre_Base1 Pre_Base2
		model_.hmm_precursor_.addSynonymTransition(PILISModel::PRE_BASE1, PILISModel::PRE_BASE2, PILISModel::PRE_ION, PILISModel::PRE_H2O_D);
		model_.hmm_precursor_.addSynonymTransition(PILISModel::PRE_BASE1, PILISModel::PRE_BASE2, PILISModel::PRE_ION, PILISModel::PRE_H2O_E);
		model_.hmm_precursor_.addSynonymTransition(PILISModel::PRE_BASE1, PILISModel::PRE_BASE2, PILISModel::PRE_ION, PILISModel::PRE_H2O_S);
		model_.hmm_precursor_.addSynonymTransition(PILISModel::PRE_BASE1, PILISModel::PRE_BASE2, PILISModel::PRE_ION, PILISModel::PRE_H2O_T);
		model_.hmm_precursor_.addSynonymTransition(PILISModel::PRE_BASE1, PILISModel::PRE_BASE2, PILISModel::PRE_ION, PILISModel::PRE_H2O_Q1);
		model_.hmm_precursor_.addSynonymTransition(PILISModel::PRE_BASE1, PILISModel::PRE_BASE2, PILISModel::PRE_ION, PILISModel::PRE_H2O_CTERM);
		model_.hmm_precursor_.addSynonymTransition(PILISModel::PRE_BASE1, PILISModel::PRE_BASE2, PILISModel::PRE_ION, PILISModel::PRE_NH3_R);
		model_.hmm_precursor_.addSynonymTransition(PILISModel::PRE_BASE1, PILISModel::PRE_BASE2, PILISModel::PRE_ION, PILISModel::PRE_NH3_K);

		model_.hmm_precursor_.disableTransitions();
		model_.hmm_precursor_.buildSynonyms();


		return;
	}

} // namespace OpenMS

