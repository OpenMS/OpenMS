// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Chris Bielow $
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

#ifdef _MSC_VER // disable some COIN-OR warnings that distract from ours
#	pragma warning( push ) // save warning state
#	pragma warning( disable : 4267 )
#else
# pragma GCC diagnostic ignored "-Wunused-parameter"
#endif
// useful doc: https://projects.coin-or.org/Cbc
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

	DoubleReal ILPDCWrapper::compute(const MassExplainer& me, const FeatureMap<> fm, PairsType& pairs, Size verbose_level)
	{
    DoubleReal score = 0;

    if (fm.size()==0)
    {
      LOG_INFO << "ILPDC wrapper received empty feature list. Nothing to compute! Exiting..." << std::endl;
      return -1;
    }
    
    PairsType pairs_clique_ordered;
    pairs_clique_ordered.reserve(pairs.size());
    typedef std::vector < std::pair <Size, Size> > BinType;
    BinType bins;
    // check number of components for complete putative edge graph (usually not all will be set to 'active' during ILP):
    {
      //
      // find groups of edges
      //
      Size group_count(0);
      Map<Size, Size> f2g; // feature id to connected group
      Map<Size, std::set<Size> > g2pairs; // group id to all pairs involved
      Map<Size, std::set<Size> > g2f;     // group id to all features involved

      for	(Size i=0;i<pairs.size();++i)
      {
        Size f1 = pairs[i].getElementIndex(0);
        Size f2 = pairs[i].getElementIndex(1);
        if (f2g.has(f1) && f2g.has(f2)) // edge connects two distinct groups
        {
          Size group1 = f2g[f1];
          Size group2 = f2g[f2];
          if (group2 != group1)
          {
            // point group2 to group1
            g2pairs[group1].insert(g2pairs[group2].begin(), g2pairs[group2].end());
            g2pairs.erase(group2);
            for (std::set<Size>::const_iterator its=g2f[group2].begin(); its != g2f[group2].end(); ++its)
            {
              g2f[group1].insert(*its);
              f2g[*its] = group1; // reassign features of group2 to group1
            }
            g2f.erase(group2);
          }
        }
        else if (f2g.has(f1) && !f2g.has(f2)) // only f1 is part of a group
        {
          Size group1 = f2g[f1];
          f2g[f2] = group1;
          g2f[group1].insert(f2);
        }
        else if (!f2g.has(f1) && f2g.has(f2)) // only f2 is part of a group
        {
          Size group2 = f2g[f2];
          f2g[f1] = group2;
          g2f[group2].insert(f1);
        }
        else // neither feature has a group: make a new one
        {
          Size group = ++group_count;
          f2g[f1] = group;
          f2g[f2] = group;
          g2f[group].insert(f1);
          g2f[group].insert(f2);
        }
        // append current edge to common group
        g2pairs[f2g[f1]].insert(i);
      }

      if (g2pairs.size() != g2f.size())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Clique construction failed! Unequal number of groups produced!", String(g2pairs.size()) + "!=" + String(g2f.size()));
      }

      Map <Size,Size> hist_component_sum;
      // now walk though groups and see the size:
      for ( Map<Size, std::set<Size> >::const_iterator it=g2f.begin();it!=g2f.end();++it)
      {
        ++hist_component_sum[it->second.size()]; // e.g. component 2 has size 4; thus increase count for size 4 
      }
      if (verbose_level > 1)
      {
        LOG_INFO << "components:\n";
        LOG_INFO << "  size 1 occurs ?x\n";
        for (OpenMS::Map <Size,Size>::const_iterator it=hist_component_sum.begin();it!=hist_component_sum.end();++it)
        {
          LOG_INFO << "  size " << it->first << " occurs " << it->second << "x\n";
        }
      }

      /* partition the cliques into bins, one given to the ILP at a time */
      UInt pairs_per_bin = 1000;
      UInt big_clique_bin_threshold = 200;

      Size start(0);
      Size count(0);
      for (Map<Size, std::set<Size> >::ConstIterator it = g2pairs.begin(); it!=g2pairs.end(); ++it)
      {
        Size clique_size = it->second.size();
        if (count > pairs_per_bin || clique_size > big_clique_bin_threshold)
        {
          if (count > 0) // either bin is full or we have to close it due to big clique
          {
            if (verbose_level > 2) LOG_INFO << "Overstepping border of " << pairs_per_bin << " by " << SignedSize(count - pairs_per_bin) << " elements!\n";
            bins.push_back(std::make_pair(start, pairs_clique_ordered.size()));
            start = pairs_clique_ordered.size();
            count=0;
          }
          if (clique_size > big_clique_bin_threshold) // extra bin for this big clique
          {
            for (std::set<Size>::const_iterator i_p = it->second.begin(); i_p != it->second.end(); ++i_p)
            {
              pairs_clique_ordered.push_back(pairs[*i_p]);
            }
            if (verbose_level > 2) LOG_INFO << "Extra bin for big clique (" << clique_size << ") prepended to schedule\n";
            bins.insert(bins.begin(), std::make_pair(start, pairs_clique_ordered.size()));
            start = pairs_clique_ordered.size();
            continue; // next clique (this one is already processed)
          }
        }
        count += clique_size;
        for (std::set<Size>::const_iterator i_p = it->second.begin(); i_p != it->second.end(); ++i_p)
        {
          pairs_clique_ordered.push_back(pairs[*i_p]);
        }
      }
      if (count>0) bins.push_back(std::make_pair(start, pairs_clique_ordered.size()));
    }

    if (pairs_clique_ordered.size() != pairs.size())
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, __PRETTY_FUNCTION__, pairs_clique_ordered.size() - pairs.size());
    }
    /* swap pairs, such that edges are order by cliques (so we can make clean cuts) */
    pairs.swap(pairs_clique_ordered);

		// split problem into slices and have each one solved by the ILPS
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for	(SignedSize i=0; i<(SignedSize)bins.size(); ++i)
		{
			score += computeSlice_(me, fm, pairs, bins[i].first, bins[i].second, verbose_level);
		}
    //score = computeSlice_(me, fm, pairs, 0, pairs.size()); // all at once - no bins

		return score;
	}
		
	DoubleReal ILPDCWrapper::computeSlice_(const MassExplainer& /*me*/,
																			   const FeatureMap<> fm,
																			   PairsType& pairs, 
																			   const PairsIndex margin_left, 
																			   const PairsIndex margin_right,
                                         Size verbose_level)
	{

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

		//------------------------------------objective function-----------------------------------------------

		// find maximal objective value
		DoubleReal score = 0;
		DoubleReal score_min=10e10f, score_max = -10e10f;

		// fill in objective values 
		std::ostringstream namebuf;

		for (PairsIndex i=margin_left; i<margin_right; ++i)
		{

     //if(i==2797)
     //{
     //   ChargePair t =  pairs[i];
     //   std::cout << "found";
     // }
			
			// log scores are good for addition in ILP - but they are < 0, thus not suitable for maximizing
			// ... so we just add normal probabilities...
      score = exp(getLogScore_(pairs[i], fm));
			pairs[i].setEdgeScore(score * pairs[i].getEdgeScore()); // multiply with preset score
			namebuf.str("");
			namebuf<<"x#"<<i;
			// create the new variable object
			build.setColumnBounds(int(i-margin_left),0,1);
			build.setInteger(int(i-margin_left));
			build.setObjective(int(i-margin_left), pairs[i].getEdgeScore());
			if (score_min > score ) score_min = score;
			if (score_max < score ) score_max = score;
			
			// DEBUG:
			//std::cerr << "MIP: egde#"<< i << " score: " << pairs[i].getEdgeScore() << " adduct:" << pairs[i].getCompomer().getAdductsAsString() << "\n";
		}
		if (verbose_level > 2) LOG_INFO << "score_min: " << score_min << " score_max: " << score_max << "\n";

		//------------------------------------adding constraints--------------------------------------------------

		bool is_conflicting;
		std::vector<int> conflict_idx(4);

		for (PairsIndex i=margin_left; i<margin_right; ++i)
		{
			const Compomer& ci = pairs[i].getCompomer();
	
      // TODO: only go until next clique...
			for (PairsIndex j=i+1; j<margin_right; ++j)
			{
				const Compomer& cj = pairs[j].getCompomer();
				
				is_conflicting = false;
				// add pairwise constraints (one for each two conflicting ChargePairs)
				// if features are identical they must have identical charges (because any single
				// feature can only have one unique charge)

        /*if(i==2796 && j==2797) // every conflict involving the missing edge...
        {
          std::cout << "debug point";
        }*/

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
				
				//incoming/outgoing edge (from one feature) -- this should only happen to additionally inferred edges
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
          /*if(i==150 || j==2797) // every conflict involving the missing edge...
          {
            ChargePair ti =  pairs[i];
            ChargePair tj =  pairs[j];
            
            std::cout << "conflicting edge!";
          }*/

					String s = String("C") + i + "." + j;

          // Now build rows: two variables, with indices 'columns', factors '1', and 0-1 bounds.
					double element[] = {1.0, 1.0};
					int columns[] = {int(i-margin_left),int(j-margin_left)};
					build.addRow(2, columns, element, 0, 1, s.c_str());
				}
			}
		}
		// add rows into solver
		solver.loadFromCoinModel(build);

		if (verbose_level > 2) LOG_INFO << "node count: " << fm.size() << "\n";
		if (verbose_level > 2) LOG_INFO << "edge count: " << pairs.size() << "\n";
		if (verbose_level > 2) LOG_INFO << "constraint count: " << (conflict_idx[0]+conflict_idx[1]+conflict_idx[2]+conflict_idx[3]) << " = " << conflict_idx[0] << " + " << conflict_idx[1] << " + " << conflict_idx[2] << " + " << conflict_idx[3] << "(0 or inferred)" << std::endl;

    // DEBUG:
		/*
    {
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
		conflict_map_out.store("c:/conflict_map.txt");
		//get rid of memory blockers:
		Map< Size, std::vector<Size> > tmp_map;
		conflict_map.swap(tmp_map);

    // write the model (for debug)
    //build.writeMps ("Y:/datasets/simulated/coinor.mps");
    }
    */
		
		//---------------------------------------------------------------------------------------------------------
		//----------------------------------------Solving and querying result--------------------------------------
		//---------------------------------------------------------------------------------------------------------

		/* Now let MIP calculate a solution */
		// Pass to solver
		CbcModel model(solver);
		model.setObjSense(-1); // -1 = maximize, 1=minimize
		model.solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);

		// Output details
		model.messageHandler()->setLogLevel( verbose_level > 1 ? 2 : 0);
		model.solver()->messageHandler()->setLogLevel( verbose_level > 1 ? 1 : 0);
		
		//CglProbing generator1;
		//generator1.setUsingObjective(true);
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
		//CglFlowCover flowGen;
    CglMixedIntegerRounding mixedGen;
		
		// Add in generators (you should prefer the ones used often and disable the others as they increase solution time)
		//model.addCutGenerator(&generator1,-1,"Probing");
		model.addCutGenerator(&generator2,-1,"Gomory");
		model.addCutGenerator(&generator3,-1,"Knapsack");
		//model.addCutGenerator(&generator4,-1,"OddHole"); // seg faults...
		model.addCutGenerator(&generator5,-10,"Clique");
		//model.addCutGenerator(&flowGen,-1,"FlowCover");
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
		if (verbose_level > 0) LOG_INFO << "Starting to solve..." << std::endl;
		model.branchAndBound();
		if (verbose_level > 0) LOG_INFO << " Branch and cut took " << CoinCpuTime()-time1 << " seconds, "
						                        << model.getNodeCount()<<" nodes with objective "
						                        << model.getObjValue()
						                        << (!model.status() ? " Finished" : " Not finished")
						                        << std::endl;



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
				++count_cmp[cmp];
			}
			else
			{
				// DEBUG
				//std::cerr << " edge " << iColumn << " with " << value << "\n";
			}
		}
		if (verbose_level > 2) LOG_INFO << "active edges: " << active_edges << " of overall " << pairs.size() << std::endl;

		for (Map < String, Size >::const_iterator it=count_cmp.begin(); it != count_cmp.end(); ++it)
		{
			//std::cout << "Cmp " << it->first << " x " << it->second << "\n";
		}

		DoubleReal opt_value = model.getObjValue();

		//objective function value of optimal(?) solution
		return opt_value;
	} // !compute

	DoubleReal ILPDCWrapper::getLogScore_(const PairsType::value_type& pair, const FeatureMap<>& fm)
	{
		// TODO think of something better here!
		DoubleReal score;
		String e;
		if (getenv ("M") != 0) e = String(getenv ("M"));
		if (e == "")
		{
			//std::cout << "1";
			score = pair.getCompomer().getLogP();
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
			DoubleReal rt_diff =  fabs(fm[pair.getElementIndex(0)].getRT() - fm[pair.getElementIndex(1)].getRT());
			// enhance correct charge
			DoubleReal charge_enhance = ( (pair.getCharge(0) == fm[pair.getElementIndex(0)].getCharge())
																	&&
																		(pair.getCharge(1) == fm[pair.getElementIndex(1)].getCharge()) )
																	? 100 : 1;
			score = charge_enhance * (1 / (pair.getMassDiff()+1) + 1 / (rt_diff+1));
		}
		
		//std::cout << "logscore: " << score << "\n";
		
		return score;
	}


}

