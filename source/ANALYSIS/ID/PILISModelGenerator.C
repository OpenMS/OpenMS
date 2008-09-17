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

using namespace std;

namespace OpenMS 
{
	PILISModelGenerator::PILISModelGenerator()
		: DefaultParamHandler("PILISModelGenerator")
	{
		defaults_.setValue("model_depth", 4, "The number of explicitly modeled backbone cleavages from N-terminus and C-terminus, would be 9 for the default value");
		defaults_.setValue("visible_model_depth", 30, "The maximal possible size of a peptide to be modeled");
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
			model_.hmm_.addNewState(first+"bk-1");
			model_.hmm_.addNewState(first+"bk-2");
    }
		model_.hmm_.addNewState("bk-1");
		model_.hmm_.addNewState("bk-2");

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

    // CR(?) bk-1, bk-2
		for (UInt i = 0; i != residues.size(); ++i)
		{
			String first;
			first += residues[i];
    	model_.hmm_.addSynonymTransition("AABase1", "AABase2", "Ak-1", first+"bk-1");

    	model_.hmm_.setTransitionProbability(first+"bk-1", "bk-1", 0.5);
    	model_.hmm_.setTransitionProbability(first+"bk-1", "endk-1", 0.5);

			model_.hmm_.addSynonymTransition("AABase1", "AABase2", "Ak-2", first+"bk-2");
			model_.hmm_.setTransitionProbability(first+"bk-2", "bk-2", 0.5);
			model_.hmm_.setTransitionProbability(first+"bk-2", "endk-2", 0.5);
			
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
					model_.hmm_.addSynonymTransition("AABase1", "AABase2", "AAk-"+num, first+second+"bxyzk-"+num);

					model_.hmm_.addSynonymTransition("AABase1", "AABase2", "AA"+num, first+second+"axyz"+num);
					model_.hmm_.addSynonymTransition("AABase1", "AABase2", "AAk-"+num, first+second+"axyzk-"+num);
				
					if (i <= 2)
					{
						if (second == "P")
						{
            	model_.hmm_.setTransitionProbability(first+second+"bxyz"+num, "bxyz"+num, 0.5);
							model_.hmm_.addSynonymTransition(first+second+"bxyz"+num, "end", first+second+"bxyz"+num, "end"+num);
							model_.hmm_.setTransitionProbability(first+second+"bxyz"+num, "end", 0.5);

							model_.hmm_.setTransitionProbability(first+second+"axyz"+num, "axyz"+num, 0.5);
							model_.hmm_.addSynonymTransition(first+second+"axyz"+num, "end", first+second+"axyz"+num, "end"+num);
							model_.hmm_.setTransitionProbability(first+second+"axyz"+num, "end", 0.5);

						}
						else
						{
							model_.hmm_.setTransitionProbability(first+second+"bxyz"+num, "bxyz"+num, 0.5);
              model_.hmm_.addSynonymTransition(first+second+"bxyz"+num, "end", first+second+"bxyz"+num, "end"+num);
              model_.hmm_.setTransitionProbability(first+second+"bxyz"+num, "end", 0.5);

							model_.hmm_.setTransitionProbability(first+second+"axyz"+num, "axyz"+num, 0.5);
							model_.hmm_.addSynonymTransition(first+second+"axyz"+num, "end", first+second+"axyz"+num, "end"+num);
							model_.hmm_.setTransitionProbability(first+second+"axyz"+num, "end", 0.5);
						}
					}
					else
					{
					
						model_.hmm_.addSynonymTransition(first+second+"bxyz", "bxyz", first+second+"bxyz"+num, "bxyz"+num);
						model_.hmm_.addSynonymTransition(first+second+"bxyz", "end", first+second+"bxyz"+num, "end"+num);

						model_.hmm_.addSynonymTransition(first+second+"axyz", "axyz", first+second+"axyz"+num, "axyz"+num);
						model_.hmm_.addSynonymTransition(first+second+"axyz", "end", first+second+"axyz"+num, "end"+num);
					}
					
					model_.hmm_.addSynonymTransition(first+second+"bxyz", "bxyz", first+second+"bxyzk-"+num, "bxyzk-"+num);
					model_.hmm_.addSynonymTransition(first+second+"bxyz", "end", first+second+"bxyzk-"+num, "endk-"+num);
					model_.hmm_.setTransitionProbability(first+second+"bxyz", "bxyz", 0.5);
					model_.hmm_.setTransitionProbability(first+second+"bxyz", "end", 0.5);
					model_.hmm_.addSynonymTransition(first+second+"bxyz", "end", first+second+"bxyzk-"+num, "end"+num);

					model_.hmm_.addSynonymTransition(first+second+"axyz", "axyz", first+second+"axyzk-"+num, "axyzk-"+num);
					model_.hmm_.addSynonymTransition(first+second+"axyz", "end", first+second+"axyzk-"+num, "endk-"+num);
					model_.hmm_.setTransitionProbability(first+second+"axyz", "axyz", 0.5);
					model_.hmm_.setTransitionProbability(first+second+"axyz", "end", 0.5);
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

	void PILISModelGenerator::initPrecursorModel_()
	{
		model_.hmm_pre_loss_.addNewState(new HMMState("PRE", true));
		model_.hmm_pre_loss_.addNewState(new HMMState("S1", true));
		model_.hmm_pre_loss_.addNewState(new HMMState("T1", true));
		model_.hmm_pre_loss_.addNewState(new HMMState("E1", true));
		model_.hmm_pre_loss_.addNewState(new HMMState("D1", true));
		model_.hmm_pre_loss_.addNewState(new HMMState("Q11", true));
		model_.hmm_pre_loss_.addNewState(new HMMState("COOH1", true));
		model_.hmm_pre_loss_.addNewState(new HMMState("S2", true));
		model_.hmm_pre_loss_.addNewState(new HMMState("T2", true));
		model_.hmm_pre_loss_.addNewState(new HMMState("E2", true));
		model_.hmm_pre_loss_.addNewState(new HMMState("D2", true));
		model_.hmm_pre_loss_.addNewState(new HMMState("Q12", true));
		model_.hmm_pre_loss_.addNewState(new HMMState("COOH2", true));
		model_.hmm_pre_loss_.addNewState(new HMMState("NH2CHNH", true));
		model_.hmm_pre_loss_.addNewState(new HMMState("M-NH2CHNH", false));
		model_.hmm_pre_loss_.addNewState(new HMMState("M-H2O", false));
		model_.hmm_pre_loss_.addNewState(new HMMState("M-NH3", false));
		model_.hmm_pre_loss_.addNewState(new HMMState("end", false));
		model_.hmm_pre_loss_.addNewState(new HMMState("M", false));
		model_.hmm_pre_loss_.addNewState(new HMMState("M-H2O-H2O", false));
		model_.hmm_pre_loss_.addNewState(new HMMState("M-H2O-NH3", false));
		model_.hmm_pre_loss_.addNewState(new HMMState("R1", true));
		model_.hmm_pre_loss_.addNewState(new HMMState("K1", true));
		model_.hmm_pre_loss_.addNewState(new HMMState("Q1", true));
		model_.hmm_pre_loss_.addNewState(new HMMState("N1", true));
		model_.hmm_pre_loss_.addNewState(new HMMState("R2", true));
		model_.hmm_pre_loss_.addNewState(new HMMState("K2", true));
		model_.hmm_pre_loss_.addNewState(new HMMState("N2", true));
		model_.hmm_pre_loss_.addNewState(new HMMState("Q2", true));
		model_.hmm_pre_loss_.addNewState(new HMMState("BASE1", true));
		model_.hmm_pre_loss_.addNewState(new HMMState("BASE2", true));
		model_.hmm_pre_loss_.setTransitionProbability("BASE1", "BASE2", 1);

		// first layer y to residue specific losses (delta 0-1 function)
		model_.hmm_pre_loss_.addSynonymTransition("BASE1", "BASE2", "PRE", "S1");
		model_.hmm_pre_loss_.addSynonymTransition("BASE1", "BASE2", "PRE", "T1");
		model_.hmm_pre_loss_.addSynonymTransition("BASE1", "BASE2", "PRE", "E1");
		model_.hmm_pre_loss_.addSynonymTransition("BASE1", "BASE2", "PRE", "D1");
		model_.hmm_pre_loss_.addSynonymTransition("BASE1", "BASE2", "PRE", "Q11");
		model_.hmm_pre_loss_.addSynonymTransition("BASE1", "BASE2", "PRE", "COOH1");

		model_.hmm_pre_loss_.addSynonymTransition("BASE1", "BASE2", "PRE", "NH2CHNH");
		
		// residue specific loss probabilities
		model_.hmm_pre_loss_.setTransitionProbability("S1", "M", 0.1);
		model_.hmm_pre_loss_.setTransitionProbability("S1", "end", 0.5);	
		model_.hmm_pre_loss_.setTransitionProbability("T1", "M", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("T1", "end", 0.5);
		model_.hmm_pre_loss_.setTransitionProbability("E1", "M", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("E1", "end", 0.5);
		model_.hmm_pre_loss_.setTransitionProbability("D1", "M", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("D1", "end", 0.5);
		model_.hmm_pre_loss_.setTransitionProbability("Q11", "M", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("Q11", "end", 0.5);
		model_.hmm_pre_loss_.setTransitionProbability("COOH1", "M", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("COOH1", "end", 0.5);

		model_.hmm_pre_loss_.setTransitionProbability("S1", "S2", 0.1);
		model_.hmm_pre_loss_.setTransitionProbability("S1", "T2", 0.1);
		model_.hmm_pre_loss_.setTransitionProbability("S1", "E2", 0.1);
		model_.hmm_pre_loss_.setTransitionProbability("S1", "D2", 0.1);
		model_.hmm_pre_loss_.setTransitionProbability("S1", "Q12", 0.1);
		model_.hmm_pre_loss_.setTransitionProbability("S1", "COOH2", 0.1);

    model_.hmm_pre_loss_.setTransitionProbability("T1", "S2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("T1", "T2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("T1", "E2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("T1", "D2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("T1", "Q12", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("T1", "COOH2", 0.1);

    model_.hmm_pre_loss_.setTransitionProbability("E1", "S2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("E1", "T2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("E1", "E2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("E1", "D2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("E1", "Q12", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("E1", "COOH2", 0.1);

    model_.hmm_pre_loss_.setTransitionProbability("D1", "S2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("D1", "T2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("D1", "E2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("D1", "D2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("D1", "Q12", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("D1", "COOH2", 0.1);

		model_.hmm_pre_loss_.setTransitionProbability("Q11", "S2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("Q11", "T2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("Q11", "E2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("Q11", "D2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("Q11", "Q12", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("Q11", "COOH2", 0.1);

    model_.hmm_pre_loss_.setTransitionProbability("COOH1", "S2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("COOH1", "T2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("COOH1", "E2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("COOH1", "D2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("COOH1", "Q12", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("COOH1", "COOH2", 0.1);

		// .. to NH3 loss states
    model_.hmm_pre_loss_.setTransitionProbability("S1", "K2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("S1", "N2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("S1", "Q2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("S1", "R2", 0.1);

    model_.hmm_pre_loss_.setTransitionProbability("T1", "K2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("T1", "N2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("T1", "Q2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("T1", "R2", 0.1);

    model_.hmm_pre_loss_.setTransitionProbability("E1", "K2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("E1", "N2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("E1", "Q2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("E1", "R2", 0.1);

    model_.hmm_pre_loss_.setTransitionProbability("D1", "K2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("D1", "N2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("D1", "Q2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("D1", "R2", 0.1);

    model_.hmm_pre_loss_.setTransitionProbability("Q11", "K2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("Q11", "N2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("Q11", "Q2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("Q11", "R2", 0.1);

		model_.hmm_pre_loss_.setTransitionProbability("COOH1", "K2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("COOH1", "N2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("COOH1", "Q2", 0.1);
    model_.hmm_pre_loss_.setTransitionProbability("COOH1", "R2", 0.1);
	
		// from second loss states to y-H2O-H2O
    model_.hmm_pre_loss_.setTransitionProbability("S2", "M-H2O-H2O", 0.4);
		model_.hmm_pre_loss_.setTransitionProbability("S2", "M-H2O", 0.4);
    model_.hmm_pre_loss_.setTransitionProbability("S2", "end", 0.2);
    model_.hmm_pre_loss_.setTransitionProbability("T2", "M-H2O-H2O", 0.4);
		model_.hmm_pre_loss_.setTransitionProbability("T2", "M-H2O", 0.4);
    model_.hmm_pre_loss_.setTransitionProbability("T2", "end", 0.2);
    model_.hmm_pre_loss_.setTransitionProbability("E2", "M-H2O-H2O", 0.4);
		model_.hmm_pre_loss_.setTransitionProbability("E2", "M-H2O", 0.4);
    model_.hmm_pre_loss_.setTransitionProbability("E2", "end", 0.2);
    model_.hmm_pre_loss_.setTransitionProbability("D2", "M-H2O-H2O", 0.4);
		model_.hmm_pre_loss_.setTransitionProbability("D2", "M-H2O", 0.4);
    model_.hmm_pre_loss_.setTransitionProbability("D2", "end", 0.2);
    model_.hmm_pre_loss_.setTransitionProbability("Q12", "M-H2O-H2O", 0.4);
		model_.hmm_pre_loss_.setTransitionProbability("Q12", "M-H2O", 0.4);
    model_.hmm_pre_loss_.setTransitionProbability("Q12", "end", 0.2);
    model_.hmm_pre_loss_.setTransitionProbability("COOH2", "M-H2O-H2O", 0.4);
		model_.hmm_pre_loss_.setTransitionProbability("COOH2", "M-H2O", 0.4);
    model_.hmm_pre_loss_.setTransitionProbability("COOH2", "end", 0.2);

		// single NH3 loss delta functions
    model_.hmm_pre_loss_.addSynonymTransition("BASE1", "BASE2", "PRE", "K1");
    model_.hmm_pre_loss_.addSynonymTransition("BASE1", "BASE2", "PRE", "N1");
    model_.hmm_pre_loss_.addSynonymTransition("BASE1", "BASE2", "PRE", "Q1");
    model_.hmm_pre_loss_.addSynonymTransition("BASE1", "BASE2", "PRE", "R1");

    model_.hmm_pre_loss_.setTransitionProbability("K1", "M-NH3", 0.4);
    model_.hmm_pre_loss_.setTransitionProbability("K1", "M", 0.4);
		model_.hmm_pre_loss_.setTransitionProbability("K1", "end", 0.2);
    model_.hmm_pre_loss_.setTransitionProbability("N1", "M-NH3", 0.4);
    model_.hmm_pre_loss_.setTransitionProbability("N1", "M", 0.4);
		model_.hmm_pre_loss_.setTransitionProbability("N1", "end", 0.2);
    model_.hmm_pre_loss_.setTransitionProbability("Q1", "M-NH3", 0.4);
		model_.hmm_pre_loss_.setTransitionProbability("Q1", "M", 0.4);
    model_.hmm_pre_loss_.setTransitionProbability("Q1", "end", 0.2);		
    model_.hmm_pre_loss_.setTransitionProbability("R1", "M-NH3", 0.4);
		model_.hmm_pre_loss_.setTransitionProbability("R1", "M", 0.4);
    model_.hmm_pre_loss_.setTransitionProbability("R1", "end", 0.2);
	
		model_.hmm_pre_loss_.setTransitionProbability("NH2CHNH", "M-NH2CHNH", 0.4);
		model_.hmm_pre_loss_.setTransitionProbability("NH2CHNH", "M", 0.4);
		model_.hmm_pre_loss_.setTransitionProbability("NH2CHNH", "end", 0.2);
		
    model_.hmm_pre_loss_.setTransitionProbability("K2", "M-H2O-NH3", 0.4);
    model_.hmm_pre_loss_.setTransitionProbability("K2", "M-H2O", 0.4);
		model_.hmm_pre_loss_.setTransitionProbability("K2", "end", 0.2);
    model_.hmm_pre_loss_.setTransitionProbability("N2", "M-H2O-NH3", 0.4);
    model_.hmm_pre_loss_.setTransitionProbability("N2", "M-H2O", 0.4);
		model_.hmm_pre_loss_.setTransitionProbability("N2", "end", 0.2);
    model_.hmm_pre_loss_.setTransitionProbability("Q2", "M-H2O-NH3", 0.4);
		model_.hmm_pre_loss_.setTransitionProbability("Q2", "M-H2O", 0.4);
    model_.hmm_pre_loss_.setTransitionProbability("Q2", "end", 0.2);	
    model_.hmm_pre_loss_.setTransitionProbability("R2", "M-H2O-NH3", 0.4);
		model_.hmm_pre_loss_.setTransitionProbability("R2", "M-H2O", 0.4);
    model_.hmm_pre_loss_.setTransitionProbability("R2", "end", 0.2);

		
		model_.hmm_pre_loss_.disableTransitions();
		model_.hmm_pre_loss_.buildSynonyms();

		return;
	}

} // namespace OpenMS

