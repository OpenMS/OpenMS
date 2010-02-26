// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/ANALYSIS/DECHARGING/ILPDCWrapper.h>

#include <iostream>
#include <ctime>
#include <cmath>
#include <climits>
#include <fstream>
#include <vector>
#include <algorithm>
#include <utility>

#include <OpenMS/DATASTRUCTURES/MassExplainer.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FORMAT/TextFile.h>

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

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

	ILPDCWrapper::ILPDCWrapper(){}

	ILPDCWrapper::~ILPDCWrapper(){}

	DoubleReal ILPDCWrapper::compute(const MassExplainer& me, const FeatureMap<> fm, PairsType& pairs)
	{
    DoubleReal score = 0;

    // check number of components:
    using namespace boost;
    {
      typedef adjacency_list <vecS, vecS, undirectedS> Graph;

      Graph g;
      //add_edge(0, 1, G);
      for	(Size i=0;i<pairs.size();++i)
		  {
        add_edge(pairs[i].getElementIndex(0), pairs[i].getElementIndex(1), g);
      }
      std::vector<int> component(num_vertices(g));
      int num = connected_components(g, &component[0]);
      
      std::cout << "Total number of components: " << num << std::endl;
      // ** check size of components **
      OpenMS::Map <int,int> hist_components, hist_component_sum;
      // size of each component
      for (Size i = 0; i != component.size(); ++i) ++hist_components[component[i]];
      // summary of component size
      for (OpenMS::Map <int,int>::const_iterator it=hist_components.begin();it!=hist_components.end();++it)
      {
        ++hist_component_sum[it->second]; // e.g. component 2 has size 4; thus increase count for size 4 
      }
      std::cout << "components of size>2:\n";
      for (OpenMS::Map <int,int>::const_iterator it=hist_component_sum.begin();it!=hist_component_sum.end();++it)
      {
        if ((it->first)>2) std::cout << "  size " << it->first << " occurs " << it->second << "x\n";
      }
      //for (Size i = 0; i != component.size(); ++i)
      //  cout << "Vertex " << i <<" is in component " << component[i] << endl;
    }

    #if 0
    //@TODO this is OpenMP capable
		// warning: this relies on RT-sorted Pairs!!!	
		UInt allowedPairs = 1000;	
		// split problem into slices and have each one solved by the ILPS
		DoubleReal scorel;
		for	(PairsIndex i=0;i+allowedPairs-1<pairs.size();i+=allowedPairs)
		{
			scorel = computeSlice_(me, pairs, i, i+allowedPairs,0,0);
			score += scorel;
			//TODO take care of overlapping borders!!
			std::cout << "slice (" << i << "-" << (i+allowedPairs) << ")score: " << scorel << std::endl; 
		}
    #endif
    score = computeSlice_(me, fm, pairs, 0, pairs.size());

		return score;
	}
		
	DoubleReal ILPDCWrapper::computeSlice_(const MassExplainer& /*me*/,
																			 const FeatureMap<> fm,
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
		CoinModel build;
		std::cout << "model created..." << std::endl;

		//------------------------------------objective function-----------------------------------------------

		// find maximal objective value
		DoubleReal score = 0;
		DoubleReal score_min=10e10f, score_max = -10e10f;

		// fill in objective values 
		std::ostringstream namebuf;

		std::cout << "adding variables...\n";
		for (PairsIndex i=margin_left; i<margin_right; ++i)
		{

     //if(i==2797)
     //{
     //   ChargePair t =  pairs[i];
     //   std::cout << "found";
     // }
			
			// log scores are good for addition in ILP - but they are < 0, thus not suitable for maximizing
			// ... so we just add normal probabilities...
      score = exp(getLogScore_(i, pairs, fm));
			pairs[i].setEdgeScore(score * pairs[i].getEdgeScore()); // multiply with preset score
			namebuf.str("");
			namebuf<<"x#"<<i;
			//std::cout << "score: " << score << "\n";
			// create the new variable object
			build.setColumnBounds(int(i-margin_left),0,1);
			build.setInteger(int(i-margin_left));
			build.setObjective(int(i-margin_left), pairs[i].getEdgeScore());
			if (score_min > score ) score_min = score;
			if (score_max < score ) score_max = score;
			
			// DEBUG:
			//std::cerr << "MIP: egde#"<< i << " score: " << pairs[i].getEdgeScore() << " adduct:" << pairs[i].getCompomer().getAdductsAsString() << "\n";
		}
		std::cout << "score_min: " << score_min << " score_max: " << score_max << "\n";

		//------------------------------------adding constraints--------------------------------------------------

		bool is_conflicting;
		std::vector<int> conflict_idx(4);

		Map< Size, std::vector<Size> > conflict_map;

		for (PairsIndex i=margin_left; i<margin_right; ++i)
		{
			const Compomer& ci = pairs[i].getCompomer();
	
			for (PairsIndex j=i+1; j<margin_right; ++j)
			{
				const Compomer& cj = pairs[j].getCompomer();
				
				is_conflicting = false;
				// add pairwise constraints (one for each two conflicting ChargePairs)
				// if features are identical they must have identical charges (because any single
				// feature can only have one unique charge)

        if(i==2796 && j==2797) // every conflict involving the missing edge...
        {
          std::cout << "debug point";
        }

				//outgoing edges (from one feature)
				if (pairs[i].getElementIndex(0) == pairs[j].getElementIndex(0))
				{
					if ( (pairs[i].getCharge(0) != pairs[j].getCharge(0))  ||
							 ci.isConflicting(cj, Compomer::LEFT, Compomer::LEFT) )
					{
						is_conflicting = true;
						++conflict_idx[0];
					}
				}

				//incoming edges (into one feature)
				if (pairs[i].getElementIndex(1) == pairs[j].getElementIndex(1))
				{
					if ( (pairs[i].getCharge(1) != pairs[j].getCharge(1))  ||
							 ci.isConflicting(cj, Compomer::RIGHT, Compomer::RIGHT) )
					{
						is_conflicting = true;
						++conflict_idx[1];
					}
				}

				//incoming/outgoing edge (from one feature) 
				if (pairs[i].getElementIndex(1) == pairs[j].getElementIndex(0))
				{
					if ( (pairs[i].getCharge(1) != pairs[j].getCharge(0))  ||
							 ci.isConflicting(cj, Compomer::RIGHT, Compomer::LEFT) )
					{
						is_conflicting = true;
						++conflict_idx[2];
					}
				}
				
				//incoming/outgoing edge (from one feature) -- this should only happen to addiditonally inferred egdes
				if (pairs[i].getElementIndex(0) == pairs[j].getElementIndex(1))
				{
					if ( (pairs[i].getCharge(0) != pairs[j].getCharge(1))  ||
							 ci.isConflicting(cj, Compomer::LEFT, Compomer::RIGHT) )
					{
						is_conflicting = true;
						++conflict_idx[3];
					}
				}
				
								
				if (is_conflicting)
				{
				
					conflict_map[i].push_back(j);
					conflict_map[j].push_back(i);
					
          if(i==150 || j==2797) // every conflict involving the missing edge...
          {
            ChargePair ti =  pairs[i];
            ChargePair tj =  pairs[j];
            
            std::cout << "conflicting edge!";
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

		// DEBUG:
		{
		std::cout << "node count: " << fm.size() << "\n";
		std::cout << "edge count: " << pairs.size() << "\n";
		std::cout << "constraint count: " << (conflict_idx[0]+conflict_idx[1]+conflict_idx[2]+conflict_idx[3]) << " = " << conflict_idx[0] << " + " << conflict_idx[1] << " + " << conflict_idx[2] << " + " << conflict_idx[3] << "(0 or inferred)" << std::endl;
		TextFile conflict_map_out;
		for (Map<Size, std::vector <Size> >::const_iterator it = conflict_map.begin(); it!= conflict_map.end(); ++it)
		{
			String s;
			s = String(it->first) + ":";
			for (Size i = 0; i<it->second.size(); ++i)
			{
				s+= " " + String(it->second[i]);
			}
			conflict_map_out.push_back(s);
		}
		//conflict_map_out.store("c:/conflict_map.txt");
		//get rid of memory blockers:
		Map< Size, std::vector<Size> > tmp_map;
		conflict_map.swap(tmp_map);
		}
		

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
		model.addCutGenerator(&generator1,-1,"Probing");
		model.addCutGenerator(&generator2,-1,"Gomory");
		model.addCutGenerator(&generator3,-1,"Knapsack");
		//model.addCutGenerator(&generator4,-1,"OddHole"); // seg faults...
		model.addCutGenerator(&generator5,-10,"Clique");
		model.addCutGenerator(&flowGen,-1,"FlowCover");
		model.addCutGenerator(&mixedGen,-1,"MixedIntegerRounding");

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
		Map < String, Size > count_cmp;
		const double * solution = model.solver()->getColSolution();
		for (int iColumn=0; iColumn<model.solver()->getNumCols(); ++iColumn)
		{
			double value=solution[iColumn];
			if (fabs(value)>0.5 && model.solver()->isInteger(iColumn))
			{
				++active_edges;
				pairs[margin_left+iColumn].setActive(true);
				// for statistical purposes: collect compomer distribution
				String cmp = pairs[margin_left+iColumn].getCompomer().getAdductsAsString();
				count_cmp[cmp] = count_cmp[cmp] + 1;
				//std::cerr << " edge " << iColumn << " with " << value << "(active)\n";
			}
			else
			{
				// DEBUG
				//std::cerr << " edge " << iColumn << " with " << value << "\n";
			}
		}
		std::cout << "active edges: " << active_edges << " of overall " << pairs.size() << std::endl;

		for (Map < String, Size >::const_iterator it=count_cmp.begin(); it != count_cmp.end(); ++it)
		{
			//std::cout << "Cmp " << it->first << " x " << it->second << "\n";
		}

		DoubleReal opt_value = model.getObjValue();

		//objective function value of optimal(?) solution
		return opt_value;
	} // !compute

	DoubleReal ILPDCWrapper::getLogScore_(const PairsIndex& i, const PairsType& pairs, const FeatureMap<>& fm)
	{
		// TODO think of something better here!
		DoubleReal score;
		String e;
		if (getenv ("M") != 0) e = String(getenv ("M"));
		if (e == "")
		{
			//std::cout << "1";
			score = pairs[i].getCompomer().getLogP();
			/*DoubleReal charge_enhance = 0;
			
			if (pairs[i].getCharge(0) == fm[pairs[i].getElementIndex(0)].getCharge()) 
				charge_enhance += log(0.9); else charge_enhance += log(0.1);
				
			if (pairs[i].getCharge(1) == fm[pairs[i].getElementIndex(1)].getCharge()) 
				charge_enhance += log(0.9); else charge_enhance += log(0.1);
		
			score += charge_enhance;
			*/

		}
		else 
		{
			//std::cout << "2";
			DoubleReal rt_diff =  fabs(fm[pairs[i].getElementIndex(0)].getRT() - fm[pairs[i].getElementIndex(1)].getRT());
			// enhance correct charge
			DoubleReal charge_enhance = ( (pairs[i].getCharge(0) == fm[pairs[i].getElementIndex(0)].getCharge())
																	&&
																		(pairs[i].getCharge(1) == fm[pairs[i].getElementIndex(1)].getCharge()) )
																	? 100 : 1;
			score = charge_enhance * (1 / (pairs[i].getMassDiff()+1) + 1 / (rt_diff+1));
		}
		
		//std::cout << "logscore: " << score << "\n";
		
		return score;
	}


}

