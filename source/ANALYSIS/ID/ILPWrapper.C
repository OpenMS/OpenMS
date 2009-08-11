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
// $Maintainer: Alexandra Zerck $
// $Authors: $
// --------------------------------------------------------------------------
#include <OpenMS/ANALYSIS/ID/ILPWrapper.h>

#ifdef _MSC_VER // disable some COIN-OR warnings that distract from ours
#	pragma warning( push ) // save warning state
#	pragma warning( disable : 4267 )
#else
# pragma GCC diagnostic ignored "-Wunused-parameter"
#endif
// useful docu: https://projects.coin-or.org/Cbc
// useful example: https://projects.coin-or.org/Cbc/browser/trunk/Cbc/examples/sample5.cpp
// Cuts

// #include "coin/CglGomory.hpp"
// #include "coin/CglProbing.hpp"
// #include "coin/CglKnapsackCover.hpp"
// #include "coin/CglOddHole.hpp"
// #include "coin/CglClique.hpp"
// #include "coin/CglFlowCover.hpp"
// #include "coin/CglMixedIntegerRounding.hpp"

// Heuristics
#include "coin/CbcHeuristic.hpp"
#include "coin/CbcHeuristicLocal.hpp"
#include "coin/CbcConfig.h"
#include "coin/CbcModel.hpp"
//#include "coin/CoinModel.hpp"
#include "coin/OsiClpSolverInterface.hpp"
#include "coin/CoinTime.hpp"
#ifdef _MSC_VER
#	pragma warning( pop )  // restore old warning state
#else
# pragma GCC diagnostic warning "-Wunused-parameter"
#endif

using namespace OpenMS;


ILPWrapper::ILPWrapper():DefaultParamHandler("ILPWrapper")
{
}

ILPWrapper::~ILPWrapper()
{
  
}




void ILPWrapper::solve(std::vector<int>& solution_indices)
{

  std::cout << "compute .." << std::endl;
#ifdef COIN_HAS_CLP
  OsiClpSolverInterface solver;
#elif COIN_HAS_OSL
  OsiOslSolverInterface solver;
#endif

	// test if model exists
	if(model_.numberElements()==0)
		{
			std::cout << "Model is empty." <<std::endl;
			return;
		}
	
  // add rows into solver
  solver.loadFromCoinModel(model_);

  /* Now let MIP calculate a solution */
  // Pass to solver
  CbcModel model(solver);

  model.setObjSense(model_.optimizationDirection()); // -1 = maximize, 1=minimize
  model.solver()->setHintParam(OsiDoReducePrint,true,OsiHintTry);

  // Output details
  model.messageHandler()->setLogLevel(2);
  model.solver()->messageHandler()->setLogLevel(1);
		
		
//   CglProbing generator1;
//   generator1.setUsingObjective(true);
//   CglGomory generator2;
//   // try larger limit
//   generator2.setLimit(300);
//   CglKnapsackCover generator3;
//   CglOddHole generator4;
//   generator4.setMinimumViolation(0.005);
//   generator4.setMinimumViolationPer(0.00002);
//   // try larger limit
//   generator4.setMaximumEntries(200);
//   CglClique generator5;
//   generator5.setStarCliqueReport(false);
//   generator5.setRowCliqueReport(false);
//   CglMixedIntegerRounding mixedGen;
//   CglFlowCover flowGen;
//   // Add in generators
//   model.addCutGenerator(&generator1,-1,"Probing");
//   model.addCutGenerator(&generator2,-1,"Gomory");
//   model.addCutGenerator(&generator3,-1,"Knapsack");
//   model.addCutGenerator(&generator4,-1,"OddHole");
//   model.addCutGenerator(&generator5,-1,"Clique");
//   model.addCutGenerator(&flowGen,-1,"FlowCover");
//   model.addCutGenerator(&mixedGen,-1,"MixedIntegerRounding");

//   CbcRounding heuristic1(model);
//   model.addHeuristic(&heuristic1);

//   // And local search when new solution found
//   CbcHeuristicLocal heuristic2(model);
//   model.addHeuristic(&heuristic2);
  // Redundant definition of default branching (as Default == User)
  //		CbcBranchUserDecision branch;
  //		model.setBranchingMethod(&branch);
  // Definition of node choice
  //		CbcCompareUser compare;
  //		model.setNodeComparison(compare);

  model.setDblParam(CbcModel::CbcMaximumSeconds,60.0*1);

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

	// best_solution
	std::cout << model.solver()->getNumCols()<<" columns has solution"<<std::endl;
	const double * solution = model.solver()->getColSolution();

	for (int column=0; column<model.solver()->getNumCols(); ++column)
		{
			double value=solution[column];
			if (fabs(value)>0.5 && model.solver()->isInteger(column))
				{
					std::cout << model_.getColumnName(column)<<" is in optimal solution"<<std::endl;
					solution_indices.push_back(column);
				}
		}
		
		
  
}
