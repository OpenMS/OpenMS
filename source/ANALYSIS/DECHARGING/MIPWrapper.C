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
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------
#include <OpenMS/ANALYSIS/DECHARGING/MIPWrapper.h>

#include <iostream>
#include <ctime>
#include <fstream>
#include <vector>

#include <OpenMS/DATASTRUCTURES/MassExplainer.h>
#include <OpenMS/CONCEPT/Exception.h>

#ifdef _MSC_VER // disable some COIN-OR warnings that distract from ours
#	pragma warning( push ) // save warning state
#	pragma warning( disable : 4267 )
#else
# pragma GCC diagnostic ignored "-Wunused-parameter"
#endif
// useful docu: https://projects.coin-or.org/Cbc
// useful example: https://projects.coin-or.org/Cbc/browser/trunk/Cbc/examples/sample5.cpp
// Cuts

#include "coin/CglGomory.hpp"
#include "coin/CglProbing.hpp"
#include "coin/CglKnapsackCover.hpp"
#include "coin/CglOddHole.hpp"
#include "coin/CglClique.hpp"
#include "coin/CglFlowCover.hpp"
#include "coin/CglMixedIntegerRounding.hpp"

// Heuristics
#include "coin/CbcHeuristic.hpp"
#include "coin/CbcHeuristicLocal.hpp"
#include "coin/CbcConfig.h"
#include "coin/CbcModel.hpp"
#include "coin/CoinModel.hpp"
#include "coin/OsiClpSolverInterface.hpp"
#include "coin/CoinTime.hpp"
#ifdef _MSC_VER
#	pragma warning( pop )  // restore old warning state
#else
# pragma GCC diagnostic warning "-Wunused-parameter"
#endif


namespace OpenMS {

	MIPWrapper::MIPWrapper(){}

	MIPWrapper::~MIPWrapper(){}

	DoubleReal MIPWrapper::compute(const MassExplainer& me, PairsType& pairs)
	{
    DoubleReal score = 0;
    #if 0
    //@TODO this is OpenMP capable
		// warning: this relies on RT-sorted Pairs!!!	
		UInt allowedPairs = 1000;	
		// split problem into slices and have each one solved by the ILPS
		DoubleReal scorel;
		for	(PairsIndex i=0;i+allowedPairs-1<pairs.size();i+=allowedPairs)
		{
			scorel = compute_slice(me, pairs, i, i+allowedPairs,0,0);
			score += scorel;
			//TODO take care of overlapping borders!!
			std::cout << "slice (" << i << "-" << (i+allowedPairs) << ")score: " << scorel << std::endl; 
		}
    #endif
    score = compute_slice_(me, pairs, 0, pairs.size());

		return score;
	}
		
	DoubleReal MIPWrapper::compute_slice_(const MassExplainer& me, 
																				PairsType& pairs, 
																				const PairsIndex margin_left, 
																				const PairsIndex margin_right)
	{

		std::cout << "compute .." << std::endl;
#ifdef COIN_HAS_CLP
	  OsiClpSolverInterface solver;
#elif COIN_HAS_OSL
	  OsiOslSolverInterface solver;
#endif
	  /* From now on we can build model in a solver independent way.
	     You can add rows one at a time but for large problems this is slow so
	     this example uses CoinBuild or CoinModel
	  */
		//OsiSolverInterface * solver = &solver1;

		std::cout << "compute 2.." << std::endl;
		CoinModel build;
		std::cout << "model created..." << std::endl;

		//------------------------------------objective function-----------------------------------------------

		// find maximal objective value
		float score = 0, score_offs=0;
		for (PairsIndex i=margin_left; i<margin_right; ++i)
		{
				score = getScore(i, pairs, me);
				if (score < score_offs) score_offs = score;
		}
		score_offs = fabs(score_offs) + 1;

		// fill in objective values 
		std::ostringstream namebuf;

		std::cout << "adding variables...\n";
		for (PairsIndex i=margin_left; i<margin_right; ++i)
		{

     if(i==2797)
     {
        ChargePair t =  pairs[i];
        std::cout << "found";
      }
      score = getScore(i, pairs, me);
			score += score_offs;
			namebuf.str("");
			namebuf<<"x#"<<i;
			// create the new variable object
			build.setColumnBounds(int(i-margin_left),0,1);
			build.setInteger(int(i-margin_left));
			build.setObjective(int(i-margin_left),score);
			//std::cout << "added #"<< i << "\n";
		}
		std::cout << "DONE adding variables..."<< std::endl;

		//------------------------------------adding constraints--------------------------------------------------

		bool is_conflicting;
		std::vector<int> conflict_idx(4);

		Adduct implicit_one(1, 1, 1, "H1", 0.9f);

		for (PairsIndex i=margin_left; i<margin_right; ++i)
		{
			for (PairsIndex j=i+1; j<margin_right; ++j)
			{
				is_conflicting = false;
				// add pairwise constraints (one for each two conflicting ChargePairs)
				// if features are identical they must have identical charges (because any single
				// feature can only have one unique charge)

				const Compomer& ci = me.getCompomerById(pairs[i].getCompomerId());
				const Compomer& cj = me.getCompomerById(pairs[j].getCompomerId());
				
				// getCharge should always be positive. getNegativeCharges() as well
				int hc_i_0 = pairs[i].getCharge(0) - ci.getNegativeCharges(); // this should always be positive! check!!
				int hc_j_0 = pairs[j].getCharge(0) - cj.getNegativeCharges(); // this should always be positive! check!!
				int hc_i_1 = pairs[i].getCharge(1) - ci.getPositiveCharges(); // this should always be positive! check!!
				int hc_j_1 = pairs[j].getCharge(1) - cj.getPositiveCharges(); // this should always be positive! check!!

				if (hc_i_0 < 0 || hc_i_1 < 0 || hc_j_0 < 0 || hc_j_1 < 0 )
				{
					std::cout << "WARNING!!! implicit number of H+ is negative!!!\nCP_i: " << pairs[i] << "\n" << ci << "\nCP_j: " << pairs[j] << "\n" << cj  << std::endl;
					exit(0);
				}


        if(i==2796 && j==2797) // every conflict involving the missing edge...
        {
          std::cout << "debug point";
        }

				//outgoing edges (from one feature)
				if (pairs[i].getElementIndex(0) == pairs[j].getElementIndex(0))
				{
					if ( (pairs[i].getCharge(0) != pairs[j].getCharge(0))  ||
							 ci.isConflicting(cj, true, true, implicit_one * hc_i_0, implicit_one * hc_j_0)  )
					{
						is_conflicting = true;
						++conflict_idx[0];
					}
				}

				//incoming edges (from one feature)
				if (pairs[i].getElementIndex(1) == pairs[j].getElementIndex(1))
				{
					if ( (pairs[i].getCharge(1) != pairs[j].getCharge(1))  ||
							 ci.isConflicting(cj, false, false, implicit_one * hc_i_1, implicit_one * hc_j_1) )
					{
						is_conflicting = true;
						++conflict_idx[1];
					}
				}

				//incoming/outgoing edge (from one feature) 
				if (pairs[i].getElementIndex(1) == pairs[j].getElementIndex(0))
				{
					if ( (pairs[i].getCharge(1) != pairs[j].getCharge(0))  ||
							 ci.isConflicting(cj, false, true , implicit_one * hc_i_1, implicit_one * hc_j_0) )
					{
						is_conflicting = true;
						++conflict_idx[2];
					}
				}
				
				//incoming/outgoing edge (from one feature) -- this should not happen
				if (pairs[i].getElementIndex(0) == pairs[j].getElementIndex(1))
				{
					if ( (pairs[i].getCharge(0) != pairs[j].getCharge(1))  ||
							 ci.isConflicting(cj, true, false, implicit_one * hc_i_0, implicit_one * hc_j_1) )
					{
						is_conflicting = true;
						++conflict_idx[3];
					}
				}
				
								
				if (is_conflicting)
				{
          if(i==2797 || j==2797) // every conflict involving the missing edge...
          {
            ChargePair ti =  pairs[i];
            ChargePair tj =  pairs[j];
            
            std::cout << "conflicting egde!";
          }

					namebuf.str("");
					//namebuf <<"cp"<<i<<"("<<pairs[i].getElementIndex(0)<<"[q"<<pairs[i].getCharge(0)<<"]-"<<pairs[i].getElementIndex(1)<<"[q"<<pairs[i].getCharge(1)<<"]"<<") vs. "
					//				<<"cp"<<j<<"("<<pairs[j].getElementIndex(0)<<"[q"<<pairs[j].getCharge(0)<<"]-"<<pairs[j].getElementIndex(1)<<"[q"<<pairs[j].getCharge(1)<<"]"<<")";
					namebuf << "C" << i << "." << j;
					String s = namebuf.str();
					//std::cout << "<-- conflict between " << s << "\n\n";


					double element[] = {1.0,1.0};
					int columns[] = {int(i-margin_left),int(j-margin_left)};
					// Now build rows: two variables, with indices 'columns', factors '1', and 0-1 bounds.

					build.addRow(2, columns, element, 0, 1, s.c_str());

					//std::cout << "Added row#" << i << " " << j << std::endl;
				}
			}
		}
		// add rows into solver
		solver.loadFromCoinModel(build);

		std::cout << " CONSTRAINTS count: " << conflict_idx[0] << " + " << conflict_idx[1] << " + " << conflict_idx[2] << " + " << conflict_idx[3] << "(0!) \n";

		// write the model (for debug)
		//build.writeMps ("Y:/datasets/simulated/coinor.mps");

		//---------------------------------------------------------------------------------------------------------
		//----------------------------------------Solving and querying result--------------------------------------
		//---------------------------------------------------------------------------------------------------------

		/* Now let MIP calculate a solution */
		// Pass to solver
		CbcModel model(solver);
		model.setObjSense(-1); // -1 = maximize, 1=minimize
		model.solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);

		// Output details
		model.messageHandler()->setLogLevel(2);
		model.solver()->messageHandler()->setLogLevel(1);
		
		
		CglProbing generator1;
		generator1.setUsingObjective(true);
		CglGomory generator2;
		generator2.setLimit(300);
		CglKnapsackCover generator3;
		CglOddHole generator4;
		generator4.setMinimumViolation(0.005);
		generator4.setMinimumViolationPer(0.00002);
		generator4.setMaximumEntries(200);
		CglClique generator5;
		generator5.setStarCliqueReport(false);
		generator5.setRowCliqueReport(false);
		CglMixedIntegerRounding mixedGen;
		CglFlowCover flowGen;
		
		// Add in generators (you should prefer the ones used often and disable the others as they increase solution time)
		//model.addCutGenerator(&generator1,-1,"Probing");
		//model.addCutGenerator(&generator2,-1,"Gomory");
		//model.addCutGenerator(&generator3,-1,"Knapsack");
		//model.addCutGenerator(&generator4,-1,"OddHole");
		model.addCutGenerator(&generator5,-10,"Clique");
		//model.addCutGenerator(&flowGen,-1,"FlowCover");
		//model.addCutGenerator(&mixedGen,-1,"MixedIntegerRounding");

		// Heuristics
		CbcRounding heuristic1(model);
		model.addHeuristic(&heuristic1);
		CbcHeuristicLocal heuristic2(model);
		model.addHeuristic(&heuristic2);

		// set maximum allowed CPU time before forced stop (dangerous!)
		//model.setDblParam(CbcModel::CbcMaximumSeconds,60.0*1);

		// Do initial solve to continuous
		model.initialSolve();
		
		
		// solve
		double time1 = CoinCpuTime();
		std::cout << "starting to solve..." << std::endl;
		model.branchAndBound();
		std::cout<<" Branch and cut took "<<CoinCpuTime()-time1<<" seconds, "
						 <<model.getNodeCount()<<" nodes with objective "
						 <<model.getObjValue()
						 <<(!model.status() ? " Finished" : " Not finished")
						 <<std::endl;



		/* variable values */
		UInt active_edges = 0;
		const double * solution = model.solver()->getColSolution();
		for (int iColumn=0; iColumn<model.solver()->getNumCols(); ++iColumn)
		{
			double value=solution[iColumn];
			if (fabs(value)>0.5 && model.solver()->isInteger(iColumn))
			{
				++active_edges;
				pairs[margin_left+iColumn].setActive(true);
			}
		}
		std::cout << "active edges: " << active_edges << " of overall " << pairs.size() << std::endl;

		DoubleReal opt_value = model.getObjValue();

		//objective function value of optimal(?) solution
		return opt_value;
	} // !compute

	float MIPWrapper::getScore(const PairsIndex& i, const PairsType& pairs, const MassExplainer& me)
	{
		// TODO think of something better here!
		return me.getCompomerById(pairs[i].getCompomerId()).getLogP();
	}


	MIPWrapper default_scipw;
}

