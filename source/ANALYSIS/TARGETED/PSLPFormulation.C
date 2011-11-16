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
// $Maintainer: Alexandra Zerck $
// $Authors: Alexandra Zerck $
// --------------------------------------------------------------------------
#include <OpenMS/ANALYSIS/TARGETED/PSLPFormulation.h>

namespace OpenMS
{


PSLPFormulation::PSLPFormulation()
{
  //model_ = new LPWrapper();
}

PSLPFormulation::~PSLPFormulation()
{
  //delete model_;
}

void PSLPFormulation::createAndSolveILP_(const FeatureMap<>& features,std::vector<std::vector<DoubleReal> >& intensity_weights,
																		std::set<Int>& charges_set,std::vector<std::vector<std::pair<Size,Size> > > & mass_ranges,
																		std::vector<IndexTriple>& variable_indices,std::vector<int>& solution_indices,
																		UInt ms2_spectra_per_rt_bin,
																		Size number_of_scans)
{
	Int counter = 0;
  model_ = new LPWrapper();
  //#define DEBUG_OPS	
	std::cout << "Feature Based: Build model: first objective"<<std::endl;
	///////////////////////////////////////////////////////////////////////
	// add objective function
	///////////////////////////////////////////////////////////////////////
	model_->setObjectiveSense(LPWrapper::MAX); // maximize
	// max \sum_j x_jk * signal_jk
	//                    column_index, feature_index,scan
	
	//

	for(Size i = 0; i < features.size(); ++i)
		{
			// first check if charge state is allowed
			// charge not in "ChargeFilter" list
#ifdef DEBUG_OPS
			std::cout << "feat: "<<i <<" charge "<<features[i].getCharge() << std::endl;
#endif
			if (charges_set.count(features[i].getCharge())<1) continue;
			if(mass_ranges[i].size()==0) continue;
#ifdef DEBUG_OPS
			if(mass_ranges[i].size() > 0)
				{
					std::cout << "start_scan "<< mass_ranges[i][0].first << " ?= "<<features[i].getQuality(0)
										<< " stop scan "<< (mass_ranges[i].end()-1)->first<< " ?= "	<< features[i].getQuality(1)-1<<std::endl;
				}
#endif

		
			Size c = 0;
			// go through all rts of the current feature
			for(Size s_idx = 0; s_idx < mass_ranges[i].size();s_idx+=2) // walk in size two steps
				{
					Size s = mass_ranges[i][s_idx].first;


					////////////////////////////////////////////////////////////////////////

						
					
#ifdef DEBUG_OPS
					std::cout << "add column "<<counter << std::endl;
#endif
					IndexTriple triple;
					triple.feature = i;
					triple.scan = s;
          Int index = model_->addColumn();
          triple.variable = index;
					variable_indices.push_back(triple);
          
          std::cout << index << " variable index"<<std::endl;
          model_->setColumnBounds(index,0,1,LPWrapper::DOUBLE_BOUNDED);
          model_->setColumnType(index,LPWrapper::BINARY); // binary variable
					model_->setColumnName(index,(String("x_")+i+","+s));
          //#ifdef DEBUG_OPS	
					std::cout << "feat "<<i << " scan "<< s << " intensity_weight "
										<< intensity_weights[i][c] <<std::endl;
          //#endif
					model_->setObjective(index,intensity_weights[i][c]);
					++counter;
					++c;
				}
		}
	
	///////////////////////////////////////////////////////////////////////
	// add constraints
	///////////////////////////////////////////////////////////////////////
#ifdef DEBUG_OPS	
	std::cout << "and now the constraints:"<<std::endl;
#endif
	///////////////////////////////////////////////////////////////////////
	// 1: ensure that each precursor is acquired maximally once
	///////////////////////////////////////////////////////////////////////

#ifdef DEBUG_OPS	
	std::cout << "first the number of times a precursors is acquired"<<std::endl;
#endif
	Size j = 0;
	for(Size i = 0; i < features.size();++i)
		{
			Size start = j;
			while(j < variable_indices.size() && variable_indices[j].feature == i)
				{
#ifdef DEBUG_OPS
					std::cout << j << " "<<variable_indices[j].variable << " "
										<< variable_indices[j].feature << " "
										<< variable_indices[j].scan<<std::endl;
#endif
					++j;
				}

			Size stop = j;
      std::vector<double> entries(stop-start);
      std::vector<int> indices(stop-start);
#ifdef DEBUG_OPS
			std::cout << "feature "<<i <<" "<<features[i].getMZ() <<" "<<features[i].getRT()<<" ";
			std::cout << stop-start<<"variables in equation\n";
#endif
			Size c = 0;
			for(Size k = start; k < stop; ++k)
				{
					entries[c] = 1.;
					indices[c] = (int) variable_indices[k].variable;
          std::cout << "indices["<<c<<"]= "<<indices[c]<<std::endl;
					++c;
				}
#ifdef DEBUG_OPS
			std::cout << "\nadd row "<<std::endl;
#endif
      //			String name = "PREC_ACQU_LIMIT_" + String(i);
			
			model_->addRow(indices,entries,(String("PREC_ACQU_LIMIT_")+i),0,1,LPWrapper::UPPER_BOUND_ONLY); // only upper bounded problem -> lower bound is ignored

#ifdef DEBUG_OPS
			std::cout << stop-start << " "<<name<<std::endl;
			std::cout << "added row"<<std::endl;
#endif
      //#undef DEBUG_OPS
			
		}



	///////////////////////////////////////////////////////////////////////
	// 2: do not exceed rt bin capacity
	///////////////////////////////////////////////////////////////////////
#ifdef DEBUG_OPS	
	std::cout << "and now the rt bin capacity"<<std::endl;
	std::cout << ms2_spectra_per_rt_bin << " rt bin capacity"<<std::endl;
#endif
	// sort variable_indices according to their scan number
	sort(variable_indices.begin(),variable_indices.end(),ScanLess());
	j = 0;
	for(Size i = 0; i < number_of_scans;++i)
		{
			// first determine number of indices:
			Size start = j;
			while(j < variable_indices.size() && variable_indices[j].scan == i)
				{
					++j;
				}
			// no feature occurring in this scan
			if(start == j) continue;

			Size stop = j;
			Size c = 0;			
      std::vector<double> entries(stop-start);
      std::vector<int> indices(stop-start);
			for(Size s = start; s < stop; ++s)
				{
					entries[c] = 1.;
					indices[c] = (int)  variable_indices[s].variable;
          std::cout << "indices["<<c<<"]= "<<indices[c]<<std::endl;
					++c;
				}
#ifdef DEBUG_OPS
			std::cout << "\nadd row "<<std::endl;
#endif
			model_->addRow(indices,entries,(String("RT_CAP")+i),0,ms2_spectra_per_rt_bin,LPWrapper::UPPER_BOUND_ONLY);// only upper bounded problem -> lower bound is ignored
#ifdef DEBUG_OPS
			std::cout << "added row"<<std::endl;
#endif
			
		}



#ifdef DEBUG_OPS	
	model_->writeProblem("/home/zerck/data/tmp/test_pis_problem.mps","MPS");
#endif
	
	solveILP_(solution_indices);
}



void PSLPFormulation::solveILP_(std::vector<int>& solution_indices)
{
#ifdef DEBUG_OPS	
  std::cout << "compute .." << std::endl;
#endif
	
	// test if model exists
	if(model_->getNumberOfColumns()==0)
		{
			std::cout << "Model is empty." <<std::endl;
			return;
		}
	
  // // add rows into solver
  // solver.loadFromCoinModel(*cmodel_);

  // /* Now let MIP calculate a solution */
  // // Pass to solver
  // CbcModel model(solver);

  // model.setObjSense(cmodel_->optimizationDirection()); // -1 = maximize, 1=minimize
  // model.solver()->setHintParam(OsiDoReducePrint,true,OsiHintTry);

  // // Output details
  // model.messageHandler()->setLogLevel(2);
  // model.solver()->messageHandler()->setLogLevel(1);
		
		
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

  // model.setDblParam(CbcModel::CbcMaximumSeconds,60.0*1);

  // // Do initial solve to continuous
  // model.initialSolve();
		
		
  // // solve
// #ifdef DEBUG_OPS	
//   double time1 = CoinCpuTime();
// 	std::cout << "starting to solve..." << std::endl;
// #endif
//   model.branchAndBound();
// #ifdef DEBUG_OPS	
//   std::cout<<" Branch and cut took "<<CoinCpuTime()-time1<<" seconds, "
// 	   <<model.getNodeCount()<<" nodes with objective "
// 	   <<model.getObjValue()
// 	   <<(!model.status() ? " Finished" : " Not finished")
// 	   <<std::endl;

// 	// best_solution
// 	std::cout << model.solver()->getNumCols()<<" columns has solution"<<std::endl;
// #endif
  LPWrapper::SolverParam param;
  model_->solve(param);

  for (Int column = 0; column < model_->getNumberOfColumns(); ++column)
  {
    double value = model_->getColumnValue(column);
    std::cout << value << " "<< model_->getColumnType(column) << std::endl;
    if ((fabs(value) > 0.5 && model_->getColumnType(column) == LPWrapper::BINARY) ||
      (fabs(value) > 0.5 && model_->getColumnType(column) == LPWrapper::INTEGER))
    {
#ifdef DEBUG_OPS	
      std::cout << model_->getColumnName(column) << " is in optimal solution" << std::endl;
#endif
      solution_indices.push_back((int) column);
    }
  }
		
		
  
}

} // namespace OpenMS
