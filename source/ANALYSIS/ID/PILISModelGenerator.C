// -*- Mode: C++; tab-width: 2; -*-
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


#include <OpenMS/ANALYSIS/ID/PILISModelGenerator.h>
#include <OpenMS/ANALYSIS/ID/HiddenMarkovModel.h>
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
	
		model.addNewState(new HMMState("axyz_center-prefix-ions", false));
    model.addNewState(new HMMState("axyz_center-suffix-ions", false));
    model.addNewState(new HMMState("axyz_center-prefix"));
    model.addNewState(new HMMState("axyz_center-suffix"));
    model.addNewState(new HMMState("bxyz_center-prefix-ions", false));
    model.addNewState(new HMMState("bxyz_center-suffix-ions", false));
    model.addNewState(new HMMState("bxyz_center-prefix"));
    model.addNewState(new HMMState("bxyz_center-suffix"));
    model.addNewState(new HMMState("D_center-prefix-ions", false));
    model.addNewState(new HMMState("D_center-suffix-ions", false));
    model.addNewState(new HMMState("D_center-prefix"));
    model.addNewState(new HMMState("D_center-suffix"));
    model.addNewState(new HMMState("E_center-prefix-ions", false));
    model.addNewState(new HMMState("E_center-suffix-ions", false));
    model.addNewState(new HMMState("E_center-prefix"));
    model.addNewState(new HMMState("E_center-suffix"));
    model.addNewState(new HMMState("K_center-prefix-ions", false));
    model.addNewState(new HMMState("K_center-suffix-ions", false));
    model.addNewState(new HMMState("K_center-prefix"));
    model.addNewState(new HMMState("K_center-suffix"));
    model.addNewState(new HMMState("R_center-prefix-ions", false));
    model.addNewState(new HMMState("R_center-suffix-ions", false));
    model.addNewState(new HMMState("R_center-prefix"));
    model.addNewState(new HMMState("R_center-suffix"));
    model.addNewState(new HMMState("H_center-prefix-ions", false));
    model.addNewState(new HMMState("H_center-suffix-ions", false));
    model.addNewState(new HMMState("H_center-prefix"));
    model.addNewState(new HMMState("H_center-suffix"));

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
		set<const Residue*> residues(ResidueDB::getInstance()->getResidues("Natural20"));
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

    model.addNewState(new HMMState("bk-1_-prefix-ions", false));
    model.addNewState(new HMMState("bk-1_-suffix-ions", false));
    model.addNewState(new HMMState("bk-1_-prefix"));
    model.addNewState(new HMMState("bk-1_-suffix"));
    model.addNewState(new HMMState("bk-2_-prefix-ions", false));
    model.addNewState(new HMMState("bk-2_-suffix-ions", false));
    model.addNewState(new HMMState("bk-2_-prefix"));
    model.addNewState(new HMMState("bk-2_-suffix"));

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
			model.addNewState(new HMMState("axyz_" + num + "-prefix-ions", false));
			model.addNewState(new HMMState("axyz_" + num + "-suffix-ions", false));
      model.addNewState(new HMMState("axyz_" + num + "-prefix"));
      model.addNewState(new HMMState("axyz_" + num + "-suffix"));
			model.addNewState(new HMMState("axyz_k-" + num + "-prefix-ions", false));
			model.addNewState(new HMMState("axyz_k-" + num + "-suffix-ions", false));
      model.addNewState(new HMMState("axyz_k-" + num + "-prefix"));
      model.addNewState(new HMMState("axyz_k-" + num + "-suffix"));
			
			model.addNewState(new HMMState("bxyz_" + num + "-prefix-ions", false));
			model.addNewState(new HMMState("bxyz_" + num + "-suffix-ions", false));
      model.addNewState(new HMMState("bxyz_" + num + "-prefix"));
      model.addNewState(new HMMState("bxyz_" + num + "-suffix"));
			model.addNewState(new HMMState("bxyz_k-" + num + "-prefix-ions", false));
			model.addNewState(new HMMState("bxyz_k-" + num + "-suffix-ions", false));
      model.addNewState(new HMMState("bxyz_k-" + num + "-prefix"));
      model.addNewState(new HMMState("bxyz_k-" + num + "-suffix"));
		
			String sc_and_cr("DEHKR");
			for (String::ConstIterator it = sc_and_cr.begin(); it != sc_and_cr.end(); ++it)
			{
				String aa(*it);
				model.addNewState(new HMMState(aa + "_" + num + "-prefix-ions", false));
				model.addNewState(new HMMState(aa + "_" + num + "-suffix-ions", false));
 	     	model.addNewState(new HMMState(aa + "_" + num + "-prefix"));
 	     	model.addNewState(new HMMState(aa + "_" + num + "-suffix"));
				model.addNewState(new HMMState(aa + "_k-" + num + "-prefix-ions", false));
				model.addNewState(new HMMState(aa + "_k-" + num + "-suffix-ions", false));
      	model.addNewState(new HMMState(aa + "_k-" + num + "-prefix"));
      	model.addNewState(new HMMState(aa + "_k-" + num + "-suffix"));
			}

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

    model.setTransitionProbability("bk-1_-prefix", "bk-1_-prefix-ions", 0.9);
    model.setTransitionProbability("bk-1_-prefix", "endk-1", 0.1);
    model.setTransitionProbability("bk-1_-suffix", "bk-1_-suffix-ions", 0.9);
		model.setTransitionProbability("bk-1_-suffix", "endk-1", 0.1);

    model.setTransitionProbability("bk-2_-prefix", "bk-2_-prefix-ions", 0.9);
    model.setTransitionProbability("bk-2_-prefix", "endk-1", 0.1);
    model.setTransitionProbability("bk-2_-suffix", "bk-2_-suffix-ions", 0.9);
    model.setTransitionProbability("bk-2_-suffix", "endk-1", 0.1);



		// set the initial transitions
		for (Size i = 1; i <= visible_model_depth; ++i)
		{
			String num(i);
			if (i <= model_depth)
			{
				model.setTransitionProbability("axyz_" + num + "-prefix", "axyz_" + num + "-prefix-ions", 0.9);
				model.setTransitionProbability("axyz_" + num + "-prefix", "end" + num, 0.1);
				model.setTransitionProbability("axyz_" + num + "-suffix", "axyz_" + num + "-suffix-ions", 0.9);
				model.setTransitionProbability("axyz_" + num + "-suffix", "end" + num, 0.1);
        model.setTransitionProbability("axyz_k-" + num + "-prefix", "axyz_k-" + num + "-prefix-ions", 0.9);
        model.setTransitionProbability("axyz_k-" + num + "-prefix", "endk-" + num, 0.1);
        model.setTransitionProbability("axyz_k-" + num + "-suffix", "axyz_k-" + num + "-suffix-ions", 0.9);
        model.setTransitionProbability("axyz_k-" + num + "-suffix", "endk-" + num, 0.1);

        model.setTransitionProbability("bxyz_" + num + "-prefix", "bxyz_" + num + "-prefix-ions", 0.9);
        model.setTransitionProbability("bxyz_" + num + "-prefix", "end" + num, 0.1);
        model.setTransitionProbability("bxyz_" + num + "-suffix", "bxyz_" + num + "-suffix-ions", 0.9);
        model.setTransitionProbability("bxyz_" + num + "-suffix", "end" + num, 0.1);
        model.setTransitionProbability("bxyz_k-" + num + "-prefix", "bxyz_k-" + num + "-prefix-ions", 0.9);
        model.setTransitionProbability("bxyz_k-" + num + "-prefix", "endk-" + num, 0.1);
        model.setTransitionProbability("bxyz_k-" + num + "-suffix", "bxyz_k-" + num + "-suffix-ions", 0.9);
        model.setTransitionProbability("bxyz_k-" + num + "-suffix", "endk-" + num, 0.1);

				String sc_and_cr("DEHRK");

				for (String::ConstIterator it = sc_and_cr.begin(); it != sc_and_cr.end(); ++it)
				{
					String aa(*it);
        	model.setTransitionProbability(aa + "_" + num + "-prefix", aa + "_" + num + "-prefix-ions", 0.9);
        	model.setTransitionProbability(aa + "_" + num + "-prefix", "end" + num, 0.1);
        	model.setTransitionProbability(aa + "_" + num + "-suffix", aa + "_" + num + "-suffix-ions", 0.9);
        	model.setTransitionProbability(aa + "_" + num + "-suffix", "end" + num, 0.1);
        	model.setTransitionProbability(aa + "_k-" + num + "-prefix", aa + "_k-" + num + "-prefix-ions", 0.9);
        	model.setTransitionProbability(aa + "_k-" + num + "-prefix", "endk-" + num, 0.1);
        	model.setTransitionProbability(aa + "_k-" + num + "-suffix", aa + "_k-" + num + "-suffix-ions", 0.9);
        	model.setTransitionProbability(aa + "_k-" + num + "-suffix", "endk-" + num, 0.1);
				}


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
				model.addSynonymTransition("axyz_center-prefix", "axyz_center-prefix-ions", "axyz_" + num + "-prefix", "axyz_" + num + "-prefix-ions");
        model.addSynonymTransition("axyz_center-prefix", "endcenter", "axyz_" + num + "-prefix", "end" + num);
				model.addSynonymTransition("axyz_center-suffix", "axyz_center-suffix-ions", "axyz_" + num + "-suffix", "axyz_" + num + "-suffix-ions");
        model.addSynonymTransition("axyz_center-suffix", "endcenter", "axyz_" + num + "-suffix", "end" + num);
        model.addSynonymTransition("axyz_center-prefix", "axyz_center-prefix-ions", "axyz_k-" + num + "-prefix", "axyz_k-" + num + "-prefix-ions");
        model.addSynonymTransition("axyz_center-prefix", "endcenter", "axyz_k-" + num + "-prefix", "endk-" + num);
        model.addSynonymTransition("axyz_center-suffix", "axyz_center-suffix-ions", "axyz_k-" + num + "-suffix", "axyz_k-" + num + "-suffix-ions");
        model.addSynonymTransition("axyz_center-suffix", "endcenter", "axyz_k-" + num + "-suffix", "endk-" + num);

        model.addSynonymTransition("bxyz_center-prefix", "bxyz_center-prefix-ions", "bxyz_" + num + "-prefix", "bxyz_" + num + "-prefix-ions");
        model.addSynonymTransition("bxyz_center-prefix", "endcenter", "bxyz_" + num + "-prefix", "end" + num);
        model.addSynonymTransition("bxyz_center-suffix", "bxyz_center-suffix-ions", "bxyz_" + num + "-suffix", "bxyz_" + num + "-suffix-ions");
        model.addSynonymTransition("bxyz_center-suffix", "endcenter", "bxyz_" + num + "-suffix", "end" + num);
        model.addSynonymTransition("bxyz_center-prefix", "bxyz_center-prefix-ions", "bxyz_k-" + num + "-prefix", "bxyz_k-" + num + "-prefix-ions");
        model.addSynonymTransition("bxyz_center-prefix", "endcenter", "bxyz_k-" + num + "-prefix", "endk-" + num);
        model.addSynonymTransition("bxyz_center-suffix", "bxyz_center-suffix-ions", "bxyz_k-" + num + "-suffix", "bxyz_k-" + num + "-suffix-ions");
        model.addSynonymTransition("bxyz_center-suffix", "endcenter", "bxyz_k-" + num + "-suffix", "endk-" + num);

        String sc_and_cr("DEHRK");
        for (String::ConstIterator it = sc_and_cr.begin(); it != sc_and_cr.end(); ++it)
        {
          String aa(*it);
          model.addSynonymTransition(aa + "_center-prefix", aa + "_center-prefix-ions", aa + "_" + num + "-prefix", aa + "_" + num + "-prefix-ions");
          model.addSynonymTransition(aa + "_center-prefix", "endcenter", aa + "_" + num + "-prefix", "end" + num);
					model.addSynonymTransition(aa + "_center-prefix", aa + "_center-prefix-ions", aa + "_k-" + num + "-prefix", aa + "_k-" + num + "-prefix-ions");
          model.addSynonymTransition(aa + "_center-prefix", "endcenter", aa + "_k-" + num + "-prefix", "endk-" + num);

          model.addSynonymTransition(aa + "_center-suffix", aa + "_center-suffix-ions", aa + "_" + num + "-suffix", aa + "_" + num + "-suffix-ions");
          model.addSynonymTransition(aa + "_center-suffix", "endcenter", aa + "_" + num + "-suffix", "end" + num);
          model.addSynonymTransition(aa + "_center-suffix", aa + "_center-suffix-ions", aa + "_k-" + num + "-suffix", aa + "_k-" + num + "-suffix-ions");
          model.addSynonymTransition(aa + "_center-suffix", "endcenter", aa + "_k-" + num + "-suffix", "endk-" + num);

        }

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

