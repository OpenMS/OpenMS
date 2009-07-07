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
	defaults_.setValue("ms2_spectra_per_rt_bin",5,"Number of allowed MS/MS spectra in a retention time bin.");
	defaults_.setMinInt("ms2_spectra_per_rt_bin",1);
	defaults_.setValue("peptides_per_protein",2,"Minimal number of peptides selected for each protein.");
	defaults_.setMinInt("peptides_per_protein",1);
}

ILPWrapper::~ILPWrapper()
{
  
}


void ILPWrapper::getXIC_(std::vector<std::pair<Size,Size> >& end_points,
												 std::vector<DoubleReal>& weights,MSExperiment<>& experiment,bool normalize)
{
	DoubleReal max_weight = 0.;
	weights.clear();
	for(Size i = 0; i < end_points.size();i+=2)
		{
			DoubleReal weight = 0.;
			for(Size j = end_points[i].second;j <= end_points[i+1].second;++j)
				{
					weight += experiment[end_points[i].first][j].getIntensity();
					//					std::cout << " add "<<experiment[end_points[i].first][j].getIntensity()<<std::endl;
				}
			if(weight > max_weight)  max_weight = weight;
			weights.push_back(weight);
		}

	if(normalize)
		{
			// normalize weights
			for(Size i = 0; i < weights.size();++i)
				{
#ifdef DEBUG_OPS
					if(end_points.size()>=i)
						{
							std::cout << "scan "<< end_points[i].first << " "<<weights[i] << " "<<max_weight
												<< " " << weights[i] / max_weight << std::endl;
						}
#endif
					weights[i] /= max_weight;
				}
		}
}



void ILPWrapper::encodeModelForKnownLCMSMapFeatureBased(FeatureMap<>& features, MSExperiment<>& experiment,
																												std::vector<IndexTriple>& variable_indices,
																												std::vector<std::vector<std::pair<Size,Size> > > & mass_ranges,
																												std::set<Int>& charges_set,UInt ms2_spectra_per_rt_bin)
{

	std::cout << "Feature Based: Build model: first objective"<<std::endl;
	///////////////////////////////////////////////////////////////////////
	// add objective function
	///////////////////////////////////////////////////////////////////////
	model_.setOptimizationDirection(-1); // maximize
	// max \sum_j x_jk * signal_jk
	//                    column_index, feature_index,scan
	
	//
	Int counter = 0;
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

			std::vector<DoubleReal> intensity_weights;
			getXIC_(mass_ranges[i],intensity_weights,experiment,true);
#ifdef DEBUG_OPS
			std::cout << "got xic"<<std::endl;
#endif

// 			if(intensity_weights.size() != (mass_ranges[i].end()-1)->first - mass_ranges[i][0].first)
// 				{
// 					std::cout << "attention: "<<intensity_weights.size() << " != "
// 										<<(mass_ranges[i].end()-1)->first - mass_ranges[i][0].first<<std::endl;
// 				}
			Size c = 0;
			// go through all rts of the current feature
			for(Size s = mass_ranges[i][0].first; s <= (mass_ranges[i].end()-1)->first;++s)
				{ 
					model_.setColumnName(counter,(String("x_")+i+","+s).c_str());
#ifdef DEBUG_OPS
					std::cout << "add column "<<counter << std::endl;
#endif
					IndexTriple triple;
					triple.feature = i;
					triple.scan = s;
					triple.variable = counter;
					variable_indices.push_back(triple);
					model_.setColumnUpper(counter,1.);
					model_.setColumnLower(counter,0.);
					model_.setColumnIsInteger(counter,true);
					
#ifdef DEBUG_OPS	
					std::cout << "feat "<<i << " scan "<< s << " intensity_weight "
										<< intensity_weights[c] <<std::endl;
#endif
					model_.setObjective(counter,intensity_weights[c]);
					++counter;
					++c;
				}
		}
	
	///////////////////////////////////////////////////////////////////////
	// add constraints
	///////////////////////////////////////////////////////////////////////
	std::cout << "and now the constraints:"<<std::endl;

	///////////////////////////////////////////////////////////////////////
	// 1: ensure that each precursor is acquired maximally once
	///////////////////////////////////////////////////////////////////////
	std::cout << "first the number of times a precursors is acquired"<<std::endl;
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
			double* entries = new double[stop-start];
			int* indices = new int[stop-start];
#ifdef DEBUG_OPS
			std::cout << "feature "<<i <<" "<<features[i].getMZ() <<" "<<features[i].getRT()<<" ";
			std::cout << stop-start<<"variables in equation\n";
#endif
			Size c = 0;
			for(Size k = start; k < stop; ++k)
				{
					entries[c] = 1.;
					indices[c] = variable_indices[k].variable;
					//					std::cout << j<<" "<<indices[j]<<std::endl;
					++c;
				}
#ifdef DEBUG_OPS
			std::cout << "\nadd row "<<std::endl;
#endif
			String name = "PREC_ACQU_LIMIT_" + String(i);
			
			model_.addRow((int)(stop-start),indices,entries,-COIN_DBL_MAX,1,name.c_str());
#ifdef DEBUG_OPS
			std::cout << stop-start << " "<<name<<std::endl;
			std::cout << "added row"<<std::endl;
#endif
			delete entries;
			delete indices;
			
		}



	///////////////////////////////////////////////////////////////////////
	// 2: do not exceed rt bin capacity
	///////////////////////////////////////////////////////////////////////
	std::cout << "and now the rt bin capacity"<<std::endl;
	std::cout << ms2_spectra_per_rt_bin << " rt bin capacity"<<std::endl;
	// sort variable_indices according to their scan number
	sort(variable_indices.begin(),variable_indices.end(),OfflinePrecursorIonSelection::ScanLess());
	j = 0;
	for(Size i = 0; i < experiment.size();++i)
		{
			// first determine number of indices:
			Size start = j;
			while(j < variable_indices.size() && variable_indices[j].scan == i)
				{
					++j;
				}
			// no feature occuring in this scan
			if(start == j) continue;

			Size stop = j;
			Size c = 0;			
			double* entries = new double[stop-start];
			int* indices = new int[stop-start];
			for(Size s = start; s < stop; ++s)
				{
					entries[c] = 1.;
					indices[c] = variable_indices[s].variable;
					++c;
				}
#ifdef DEBUG_OPS
			std::cout << "\nadd row "<<std::endl;
#endif
			model_.addRow((int)(stop-start),indices,entries,-COIN_DBL_MAX,ms2_spectra_per_rt_bin,(String("RT_CAP")+i).c_str());
#ifdef DEBUG_OPS
			std::cout << "added row"<<std::endl;
#endif
			delete entries;
			delete indices;
			
		}



#ifdef DEBUG_OPS	
	model_.writeMps("/home/zerck/data/tmp/test_pis_problem.mps",0,0,2,true);
#endif

}

void ILPWrapper::encodeModelForOptimalSolution(FeatureMap<>& features,
																							 MSExperiment<>& experiment,
																							 std::vector<std::vector<std::pair<Size,Size> > > & mass_ranges,
																							 std::map<String,std::vector<Size> >& protein_precursor_map,
																							 std::vector<IndexTriple>& variable_indices,
																							 UInt ms2_spectra_per_rt_bin)
{
 	std::cout << "Find optimal solution: Build model: first objective"<<std::endl;
	///////////////////////////////////////////////////////////////////////
	// add objective function
	///////////////////////////////////////////////////////////////////////
	model_.setOptimizationDirection(-1); // maximize
	// max \sum_j x_jk * signal_jk
	//                    column_index, feature_index,scan
	
	//
	
	Int counter = 0;
	for(Size i = 0; i < features.size(); ++i)
		{
			if(features[i].getPeptideIdentifications().size()==0 || mass_ranges[i].size()==0) continue;
#ifdef DEBUG_OPS
			if(mass_ranges[i].size() > 0)
				{
					std::cout << "start_scan "<< mass_ranges[i][0].first << " ?= "<<features[i].getQuality(0)
										<< " stop scan "<< (mass_ranges[i].end()-1)->first<< " ?= "	<< features[i].getQuality(1)-1<<std::endl;
				}
#endif

// 			std::vector<DoubleReal> intensity_weights;
// 			getXIC_(mass_ranges[i],intensity_weights,experiment,true);
#ifdef DEBUG_OPS
			std::cout << "got xic"<<std::endl;
#endif
			Size c = 0;
			// go through all rts of the current feature
			for(Size s = mass_ranges[i][0].first; s <= (mass_ranges[i].end()-1)->first;++s)
				{ 
					model_.setColumnName(counter,(String("x_")+i+","+s).c_str());
#ifdef DEBUG_OPS
					std::cout << "add column "<<counter << std::endl;
#endif
					IndexTriple triple;
					triple.feature = i;
					triple.scan = s;
					triple.variable = counter;
					variable_indices.push_back(triple);
					model_.setColumnUpper(counter,1.);
					model_.setColumnLower(counter,0.);
					model_.setColumnIsInteger(counter,true);
					
#ifdef DEBUG_OPS	
					std::cout << "feat "<<i << " scan "<< s << " intensity_weight "
										<< intensity_weights[c] <<std::endl;
#endif
					//					model_.setObjective(counter,intensity_weights[c]);
					model_.setObjective(counter,1);
					++counter;
					++c;
				}
		}
	
	///////////////////////////////////////////////////////////////////////
	// add constraints
	///////////////////////////////////////////////////////////////////////
	std::cout << "and now the constraints:"<<std::endl;

	///////////////////////////////////////////////////////////////////////
	// 1: ensure that each precursor is acquired maximally once
	///////////////////////////////////////////////////////////////////////
	std::cout << "first the number of times a precursors is acquired"<<std::endl;
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
			double* entries = new double[stop-start];
			int* indices = new int[stop-start];
#ifdef DEBUG_OPS
			std::cout << "feature "<<i <<" "<<features[i].getMZ() <<" "<<features[i].getRT()<<" ";
			std::cout << stop-start<<"variables in equation\n";
#endif
			Size c = 0;
			for(Size k = start; k < stop; ++k)
				{
					entries[c] = 1.;
					indices[c] = variable_indices[k].variable;
					//					std::cout << j<<" "<<indices[j]<<std::endl;
					++c;
				}
#ifdef DEBUG_OPS
			std::cout << "\nadd row "<<std::endl;
#endif
			String name = "PREC_ACQU_LIMIT_" + String(i);
			
			model_.addRow((int)(stop-start),indices,entries,-COIN_DBL_MAX,1,name.c_str());
#ifdef DEBUG_OPS
			std::cout << stop-start << " "<<name<<std::endl;
			std::cout << "added row"<<std::endl;
#endif
			delete entries;
			delete indices;
			
		}



	///////////////////////////////////////////////////////////////////////
	// 2: do not exceed rt bin capacity
	///////////////////////////////////////////////////////////////////////
	std::cout << "and now the rt bin capacity"<<std::endl;
	std::cout << ms2_spectra_per_rt_bin << " rt bin capacity"<<std::endl;
	// sort variable_indices according to their scan number
	sort(variable_indices.begin(),variable_indices.end(),OfflinePrecursorIonSelection::ScanLess());
	j = 0;
	for(Size i = 0; i < experiment.size();++i)
		{
			// first determine number of indices:
			Size start = j;
			while(j < variable_indices.size() && variable_indices[j].scan == i)
				{
					++j;
				}
			// no feature occuring in this scan
			if(start == j) continue;

			Size stop = j;
			Size c = 0;			
			double* entries = new double[stop-start];
			int* indices = new int[stop-start];
			for(Size s = start; s < stop; ++s)
				{
					entries[c] = 1.;
					indices[c] = variable_indices[s].variable;
					++c;
				}
#ifdef DEBUG_OPS
			std::cout << "\nadd row "<<std::endl;
#endif
			model_.addRow((int)(stop-start),indices,entries,-COIN_DBL_MAX,ms2_spectra_per_rt_bin,(String("RT_CAP")+i).c_str());
#ifdef DEBUG_OPS
			std::cout << "added row"<<std::endl;
#endif
			delete entries;
			delete indices;
			
		}


	///////////////////////////////////////////////////////////////////////
	// 3: protein coverage
	///////////////////////////////////////////////////////////////////////
	std::cout << "and now the protein coverage"<<std::endl;
	sort(variable_indices.begin(),variable_indices.end(),OfflinePrecursorIonSelection::VariableIndexLess());	
	std::map<String,std::vector<Size> >::iterator map_iter = protein_precursor_map.begin();	
	for(; map_iter != protein_precursor_map.end();++map_iter)
		{
			std::vector<Int> indices_vec;
			std::vector<Size>::iterator f_index_iter = map_iter->second.begin();
			// go through all feature that have ids belonging to this protein
			for(; f_index_iter != map_iter->second.end(); ++f_index_iter)
				{
					// now go through all x_variable for this feature
					for(Size v = 0; v < variable_indices.size();++v)
						{
							if(variable_indices[v].feature == *f_index_iter)
								{
									// if there are duplicates in the vector CoinModel will abort
									if(find(indices_vec.begin(),indices_vec.end(),variable_indices[v].variable) == indices_vec.end())
										{
											indices_vec.push_back(variable_indices[v].variable);
										}
								}
							else if(variable_indices[v].feature > *f_index_iter)  break; // the indices are sorted, hence if the current index is larger, we are finished
						}
				}
			if(indices_vec.size() == 0) continue;
			if(indices_vec.size() < 2)
				{
					std::cout << "too few features with ids for this protein, skipping protein"<<std::endl;
					continue;
				}
			std::vector<double> entries(indices_vec.size(),1.);
#ifdef DEBUG_OPS
			std::cout << "\nadd row "<<std::endl;
#endif
			Int i = distance(protein_precursor_map.begin(),map_iter);
			std::cout << (String("PROT_COV_")+i) <<"\t"<<(String("PROT_COV_")+i).c_str() <<std::endl;
			std::cout << indices_vec.size()<< " "<<&(indices_vec[0])<< " "<<&(entries[0])<<std::endl;
			std::cout << (indices_vec[0])<< " "<<(entries[0])<<std::endl;
			// at the moment we want maximally 2 precursors for each protein
			model_.addRow((int)indices_vec.size(),&(indices_vec[0]),&(entries[0]),-COIN_DBL_MAX,2,(String("PROT_COV_")+i).c_str());
#ifdef DEBUG_OPS
			std::cout << "added row"<<std::endl;
#endif
			
		}
	std::cout << "model was built" <<std::endl;
//#ifdef DEBUG_OPS	
	model_.writeMps("/home/zerck/data/tmp/test_pis_problem.mps",0,0,2,true);
//#endif
 
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
