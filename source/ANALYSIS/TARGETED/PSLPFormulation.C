// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/ANALYSIS/TARGETED/PrecursorIonSelectionPreprocessing.h>

namespace OpenMS
{


PSLPFormulation::PSLPFormulation():DefaultParamHandler("PSLPFormulation"), solver_(LPWrapper::SOLVER_GLPK)
{
  //model_ = new LPWrapper();
  
  defaults_.setValue("rt:min_rt",960.,"Minimal rt in seconds.");
	defaults_.setMinFloat("rt:min_rt",0.);

  defaults_.setValue("rt:max_rt",3840.,"Maximal rt in seconds.");
	defaults_.setMinFloat("rt:max_rt",0.);

  defaults_.setValue("rt:rt_step_size",30.,"rt step size in seconds.");
	defaults_.setMinFloat("rt:rt_step_size",1.);

  defaults_.setValue("rt:rt_window_size",100,"rt window size in seconds.");
  defaults_.setMinInt("rt:rt_window_size",1);
  
  // defaults_.setValue("thresholds:min_protein_probability",0.2,"Minimal protein probability for a protein to be considered in the ILP");
	// defaults_.setMinFloat("thresholds:min_protein_probability",0.);
	// defaults_.setMaxFloat("thresholds:min_protein_probability",1.);

  // defaults_.setValue("thresholds:min_protein_id_probability",0.95,"Minimal protein probability for a protein to be considered identified.");
  // defaults_.setMinFloat("thresholds:min_protein_id_probability",0.);
	// defaults_.setMaxFloat("thresholds:min_protein_id_probability",1.);

  defaults_.setValue("thresholds:min_pt_weight",0.5,"Minimal pt weight of a precursor");
	defaults_.setMinFloat("thresholds:min_pt_weight",0.);
	defaults_.setMaxFloat("thresholds:min_pt_weight",1.);

  defaults_.setValue("thresholds:min_mz",500.,"Minimal mz to be considered in protein based LP formulation.");
  defaults_.setMinFloat("thresholds:min_mz",0.);

  defaults_.setValue("thresholds:max_mz",5000.,"Minimal mz to be considered in protein based LP formulation.");
  defaults_.setMinFloat("thresholds:max_mz",0.);
  
  // defaults_.setValue("thresholds:min_pred_pep_prob",0.5,"Minimal predicted peptide probability of a precursor");
	// defaults_.setMinFloat("thresholds:min_pred_pep_prob",0.);
	// defaults_.setMaxFloat("thresholds:min_pred_pep_prob",1.);

	// defaults_.setValue("thresholds:min_rt_weight",0.5,"Minimal rt weight of a precursor");
	// defaults_.setMinFloat("thresholds:min_rt_weight",0.);
	// defaults_.setMaxFloat("thresholds:min_rt_weight",1.);


	// defaults_.setValue("thresholds:use_peptide_rule","false","Use peptide rule instead of minimal protein id probability");
	// defaults_.setValidStrings("thresholds:use_peptide_rule",StringList::create("true,false"));

  // defaults_.setValue("thresholds:min_peptide_ids",2,"If use_peptide_rule is true, this parameter sets the minimal number of peptide ids for a protein id");
  // defaults_.setMinInt("thresholds:min_peptide_ids",1);
  
  // defaults_.setValue("thresholds:min_peptide_probability",0.95,"If use_peptide_rule is true, this parameter sets the minimal probability for a peptide to be safely identified");
  // defaults_.setMinFloat("thresholds:min_peptide_probability",0.);
  // defaults_.setMaxFloat("thresholds:min_peptide_probability",1.);  

	// defaults_.setValue("mz_tolerance",25.,"Allowed precursor mass error tolerance in ppm.");
	// defaults_.setMinFloat("mz_tolerance",0.);
  
  // defaults_.setValue("combined_ilp:k1",0.2,"combined ilp: weight for z_i");
  // defaults_.setMinFloat("combined_ilp:k1",0.);
  // //	defaults_.setMaxFloat("combined_ilp:k1",1.);
  // defaults_.setValue("combined_ilp:k2",0.2,"combined ilp: weight for x_j,s*int_j,s");
  // defaults_.setMinFloat("combined_ilp:k2",0.);
  // //	defaults_.setMaxFloat("combined_ilp:k1",1.);
  // defaults_.setValue("combined_ilp:k3",0.4,"combined ilp: weight for -x_j,s*w_j,s");
  // defaults_.setMinFloat("combined_ilp:k3",0.);
  // //	defaults_.setMaxFloat("combined_ilp:k1",1.);
  // defaults_.setValue("combined_ilp:scale_matching_probs","true","flag if detectability * rt_weight shall be scaled to cover all [0,1]");
  // defaults_.setValidStrings("combined_ilp:scale_matching_probs",StringList::create("true,false"));
	defaultsToParam_();

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
#ifdef DEBUG_OPS
	std::cout << "Feature Based: Build model: first objective"<<std::endl;
#endif
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
			if(mass_ranges[i].empty()) continue;
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
#ifdef DEBUG_OPS          
          std::cout << index << " variable index"<<std::endl;
#endif
          model_->setColumnBounds(index,0,1,LPWrapper::DOUBLE_BOUNDED);
          model_->setColumnType(index,LPWrapper::BINARY); // binary variable
					model_->setColumnName(index,(String("x_")+i+","+s));
#ifdef DEBUG_OPS	
					std::cout << "feat "<<i << " scan "<< s << " intensity_weight "
										<< intensity_weights[i][c] <<std::endl;
#endif
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
#ifdef DEBUG_OPS
          std::cout << "indices["<<c<<"]= "<<indices[c]<<std::endl;
#endif
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
			while(j < variable_indices.size() && (Size)variable_indices[j].scan == i)
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
#ifdef DEBUG_OPS
          std::cout << "indices["<<c<<"]= "<<indices[c]<<std::endl;
#endif
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
#ifdef DEBUG_OPS
    std::cout << value << " "<< model_->getColumnType(column) << std::endl;
#endif
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


  void PSLPFormulation::createAndSolveILPForInclusionListCreation(PrecursorIonSelectionPreprocessing& preprocessing,
                                                                UInt ms2_spectra_per_rt_bin,
                                                                UInt max_list_size,
                                                                FeatureMap<>& precursors,
                                                                bool solve_ILP)
{
  const std::map<String, std::vector<DoubleReal> >& pt_prot_map = preprocessing.getProteinPTMap();
  std::map<String, std::vector<DoubleReal> >::const_iterator map_iter = pt_prot_map.begin();
  
  model_ = new LPWrapper();
  model_->setSolver(solver_);
  model_->setObjectiveSense(LPWrapper::MAX); // maximize
  
  DoubleReal min_rt = param_.getValue("rt:min_rt");
  DoubleReal max_rt = param_.getValue("rt:max_rt");
  DoubleReal rt_step_size = param_.getValue("rt:rt_step_size");
    
  Size max_index = (Size)ceil((max_rt-min_rt) / rt_step_size);
  Size counter = 0;
  Size feature_counter = 0;
  std::vector<IndexTriple> variable_indices;
  std::map<String,Size> protein_penalty_index_map;
  Size pep_counter = 0;
#ifdef DEBUG_OPS
  FeatureMap<> test_map;
#endif
  //  std::cout << "add proteins to LP"<<std::endl;
  for(;map_iter!=pt_prot_map.end();++map_iter)
    {
      addProteinToILP_(preprocessing,map_iter,counter,pep_counter,feature_counter,
                       variable_indices,protein_penalty_index_map,precursors);
    } //  for(;map_iter!=pt_prot_map.end();++map_iter)

  ///////////////////////////////////////////////////////////////////////
  // constraint 1: ensure that the maximal total number of precursors is not exceeded
  ///////////////////////////////////////////////////////////////////////
  addMaxInclusionListSizeConstraints_(variable_indices,/*feature_counter, */max_list_size);
  //////////////////////////////////////////////////////////////////////////
  // constraint 2: rt_bin_capacity
  //////////////////////////////////////////////////////////////////////////
  addRTBinCapacityConstraint_(variable_indices,max_index,ms2_spectra_per_rt_bin);
  ///////////////////////////////////////////////////////////////////////
  // constraint 3: protein coverage
  ///////////////////////////////////////////////////////////////////////
  addProteinCoverageConstraint_(variable_indices,preprocessing,protein_penalty_index_map);
  //  model_->writeMps("/project/maldi/data/projects/data_openms/results/inclusion_lists/20110502/test_pis_problem.mps",0,0,2,true);
#ifdef DEBUG_OPS
  //cmodel_->writeMps("/home/zerck/data/tmp/test_pis_problem.mps",0,0,2,true);
  //FeatureXMLFile().store("/home/zerck/data_project_maldi/results/inclusion_lists/all_feats.featureXML",test_map);
#endif
  if(solve_ILP)
    {
      precursors.clear(true);
      std::vector<int> solution_indices;
      solveILP_(solution_indices);
      assembleInclusionListForProteinBasedLP_(variable_indices,precursors,solution_indices,preprocessing);
   }
  else
    {

    }
  //  std::cout <<"ILPwrapper is finished"<<std::endl;

}


void PSLPFormulation::addProteinToILP_(PrecursorIonSelectionPreprocessing& preprocessing,
                                       std::map<String,std::vector<DoubleReal> >::const_iterator map_iter,
                                       Size& counter,Size& pep_counter,Size& feature_counter,
                                       std::vector<IndexTriple>& variable_indices,
                                       std::map<String,Size>& protein_variable_index_map,
                                       FeatureMap<>& precursors)
{
  DoubleReal min_rt = param_.getValue("rt:min_rt");
  DoubleReal max_rt = param_.getValue("rt:max_rt");
  DoubleReal rt_step_size = param_.getValue("rt:rt_step_size");
  DoubleReal min_pt = param_.getValue("thresholds:min_pt_weight");
  DoubleReal min_mz = param_.getValue("thresholds:min_mz");
  DoubleReal max_mz = param_.getValue("thresholds:max_mz");
  Size rt_window_size = param_.getValue("rt:rt_window_size");
  Size max_index = (Size)ceil((max_rt-min_rt) / rt_step_size);

  const std::vector<DoubleReal>& masses =  preprocessing.getMasses(map_iter->first);
  // insert protein variable
  Int index = model_->addColumn();
  model_->setColumnName(index,(String("y_")+map_iter->first).c_str());
  model_->setColumnBounds(index,0.,0.,LPWrapper::LOWER_BOUND_ONLY);
  //model_->setColumnBounds(index,0.,1.,LPWrapper::DOUBLE_BOUNDED);
  //  cmodel_->setColumnUpper(counter,1.); // test for inlcusion list protein based

 //  cmodel_->setColumnIsInteger(counter,true); // testweise abgeschaltet, da er sonst, wenn nicht ausreichend
  // Peptide da waren, um das Protein auf 1 zu setzen, gar keine Variablen f?r das Protein ausw?hlt
  model_->setObjective(index,1.);
  protein_variable_index_map.insert(make_pair(map_iter->first,counter));
  ++counter;

  std::cout << "protein "<<map_iter->first <<" with "<<map_iter->second.size()<<" peptides."<<std::endl;
  std::cout << min_pt  <<std::endl;
  for(Size p = 0;p < map_iter->second.size();++p)
    {
      // check if peptide is in desired mass range
      if(masses[p] < min_mz || masses[p] > max_mz) continue;
      if(map_iter->second[p] > min_pt)
        {
          //          std::cout <<"feature_counter "<<feature_counter<<std::endl;
          std::cout << "peptide "<<p<<"\tmz "<<masses[p]<<"\t";
          ++pep_counter;
          Feature test_feature;
          DoubleReal rt = preprocessing.getRT(map_iter->first,p);
          std::cout<<"predicted rt "<<rt;
          // DoubleReal mz_weight = preprocessing.getWeight(masses[p]);
          // if(mz_weight == 0.) mz_weight = 1.;
          // first find closest rt
          Size rt_index = std::min((Size)std::max(0.,ceil( (rt-min_rt)/rt_step_size)) ,max_index);
          DoubleReal curr_rt = min_rt + rt_index * rt_step_size;
          
          std::cout << "\tnearest rt "<<curr_rt<<"\t";
          Size curr_rt_index = rt_index;
          std::cout << "dt "<<map_iter->second[p]<<std::endl;
#ifdef DEBUG_OPS
          test_feature.setMZ(masses[p]);
          test_feature.setRT(curr_rt);
          std::map<String, std::vector<String> >::const_iterator prot_pep_iter = prot_pep_seq_map.find(map_iter->first);
          if(prot_pep_iter != prot_pep_seq_map.end() )
            {
              const std::vector<String>& pep_seqs =  prot_pep_iter->second;
              test_feature.setMetaValue("sequence",pep_seqs[p]);
            }
          
          
#endif
          String name = map_iter->first+"_"+String(p);
          Int last_index=0,current_index=0;
          Int rt_center_variable_index=0;
          // then go to the left as long as possible
          // with a fixed window this means until we are not more than rt_window_size from rt_index
          //          std::cout << curr_rt << " >=? "<< min_rt << "\t" << rt_index << " - "<< curr_rt_index <<" <=?" <<  rt_window_size << std::endl;
          while(curr_rt >= min_rt && fabs(curr_rt - rt) <= rt_window_size)
            {
              // create variable
              index = model_->addColumn();
              model_->setColumnName(index,(String("x_")+map_iter->first+"_"+String(p)+","+String(curr_rt)).c_str());
              //#ifdef DEBUG_OPS
              std::cout << "add column "<<index<< " "
                        << (String("x_")+map_iter->first+"_"+String(p)+","+String(curr_rt)) << " "
                        << std::endl;
              //#endif
              //  std::cout <<"adding "<<counter<<std::endl;
              IndexTriple triple;
              triple.prot_acc = map_iter->first;
              triple.feature = p;
              triple.scan = curr_rt_index;
              triple.variable = index;
              variable_indices.push_back(triple);
              model_->setColumnBounds(index,0.,1.,LPWrapper::DOUBLE_BOUNDED);
              model_->setColumnType(index,LPWrapper::BINARY);
              ConvexHull2D hull;
              hull.addPoint(DPosition<2>(masses[p],curr_rt));
              hull.addPoint(DPosition<2>(masses[p]+ 3.,curr_rt)); // TODO: different charge states!!!
              test_feature.getConvexHulls().push_back(hull);

#ifdef DEBUG_OPS
              std::cout << "protein "<<map_iter->first <<" peptide "<< p

                        << " scan "<< curr_rt_index << " weight "
                        << curr_rt_weight*map_iter->second[p]*mz_weight
                        << " = "<< curr_rt_weight<< " * "<< map_iter->second[p] << std::endl;
              //                            << " * "<<mz_weight<<std::endl;
#endif
              model_->setObjective(index,0.0);
              //cmodel_->setObjective(counter,-0.05);
              
              ++counter;

              if(rt_index == curr_rt_index)
                {
                  //                  last_index = index;
                  current_index = index;
                  rt_center_variable_index = index;
                }
              else
                {
                  last_index = current_index;
                  current_index = index;
                }
              // enter constraint that ensures only consecutive scans are chosen for a feature
              if(rt_index != curr_rt_index)
                {
                  std::vector<DoubleReal> entries(2);
                  std::vector<Int> indices(2);
                  indices[0] = last_index;
                  indices[1] = current_index;
                  entries[0] = -1.;
                  entries[1] = 1.;                  
                  //#ifdef DEBUG_OPS
                  std::cout << "\nadd row "<<std::endl;
                  std::cout << "indices: "<< indices[0] << " "<< indices[1]<<std::endl;
                  //#endif
                  model_->addRow(indices,entries,String("rt_cons_") +name+String("_")+String(curr_rt_index)+String("_") +String(curr_rt_index+1),
                                 0,0,LPWrapper::UPPER_BOUND_ONLY);// only upper bounded problem -> lower bound is ignored
                  
                }
              
              // update curr_rt
              if(curr_rt_index == 0) break;
              --curr_rt_index;
              curr_rt -= rt_step_size;
            }
#ifdef DEBUG_OPS
          std::cout << "links fertig----nun nach rechts"<<std::endl;
#endif
          curr_rt_index = rt_index+1;
          curr_rt = min_rt + curr_rt_index * rt_step_size;
          std::cout << curr_rt<<" <= "<<max_rt << " "<< curr_rt_index - rt_index <<" <= "<< rt_window_size<<std::endl;
          // then to the right
          // here the same applies as before
          // if we use a fixed rt window size we continue until we are not more than rt_window_size away from our origin
          while(curr_rt <= max_rt && fabs(curr_rt - rt) <= rt_window_size)
            {
              // create variable
              index = model_->addColumn();
              model_->setColumnName(index,(String("x_")+map_iter->first+"_"+String(p)+","+String(curr_rt)).c_str());
              //#ifdef DEBUG_OPS
              std::cout << "add column "<<index << " "
                        << (String("x_")+map_iter->first+"_"+String(p)+","+String(curr_rt)) << " "
                        << std::endl;
              //#endif
              //std::cout <<"adding "<<counter<<std::endl;
              IndexTriple triple;
              triple.prot_acc = map_iter->first;
              triple.feature = p;
              triple.scan = curr_rt_index;
              triple.variable = index;
              variable_indices.push_back(triple);
              model_->setColumnBounds(index,0.,1.,LPWrapper::DOUBLE_BOUNDED);
              model_->setColumnType(index,LPWrapper::BINARY);

              ConvexHull2D hull;
              hull.addPoint(DPosition<2>(masses[p],curr_rt));
              hull.addPoint(DPosition<2>(masses[p]+ 3.,curr_rt)); // TODO: different charge states!!!
              test_feature.getConvexHulls().push_back(hull);

#ifdef DEBUG_OPS
              std::cout << "protein "<<map_iter->first <<" peptide "<< p
                        << " scan "<< curr_rt_index << " weight "
                        << curr_rt_weight*map_iter->second[p]*mz_weight
                        << " = "<< curr_rt_weight<< " * "<< map_iter->second[p]<<std::endl;
              //  << " * "<<mz_weight <<std::endl;
#endif
              //cmodel_->setObjective(counter,1.-curr_rt_weight*map_iter->second[p]);
              //cmodel_->setObjective(counter,-0.05);
              model_->setObjective(index,0.0);
              ++counter;

              if(curr_rt_index == rt_index +1)  last_index = rt_center_variable_index;
              else last_index = current_index;
              current_index = index;
              
              // enter constraint that ensures only consecutive scans are chosen for a feature
              std::vector<DoubleReal> entries(2);
              std::vector<Int> indices(2);
              indices[0] = current_index;
              indices[1] = last_index;
              entries[0] = 1.;
              entries[1] = -1.;                  
#ifdef DEBUG_OPS
              std::cout << "\nadd row "<<std::endl;
#endif
              model_->addRow(indices,entries,String("rt_cons_") +name+String("_")+String(curr_rt_index)+String("_") +String(curr_rt_index+1),
                             0,0,LPWrapper::UPPER_BOUND_ONLY);// only upper bounded problem -> lower bound is ignored

              // update curr_rt
              if(curr_rt_index >= max_index) break;
              ++curr_rt_index;
              curr_rt += rt_step_size;
            }
          precursors.push_back(test_feature);
        }// if(map_iter->second[p] > min_pt)

      ++feature_counter;
    } // for(Size p = 0;p < map_iter->second.size();++p)

}

void PSLPFormulation::addMaxInclusionListSizeConstraints_(std::vector<IndexTriple>& variable_indices,
                                                          /*Size number_of_features,*/UInt max_list_size)
{
	///////////////////////////////////////////////////////////////////////
	// ensure that each precursor is counted maximally once, so create new variable
	///////////////////////////////////////////////////////////////////////
  //	std::cout << "now the number of times a precursors is acquired"<<std::endl;
  Size old_num = model_->getNumberOfColumns();
  //	Size j = 0;
  Size old_size = variable_indices.size();
	for(Size i = 0; i < old_size;++i)
		{

      // get protein accession and feature counter
      String acc = variable_indices[i].prot_acc;
      Size feature_count = variable_indices[i].feature;
      
      std::cout << "i:"<<std::endl;
      // get start and end point of this feature in variable_indices-vector
      Size start = i;
      Size end = i;
      //      std::cout << acc<< " "<< variable_indices[end+1].prot_acc << " "<< feature_count << " "<<variable_indices[end+1].feature << std::endl;
      while(end+1 < variable_indices.size() && variable_indices[end+1].prot_acc == acc && variable_indices[end+1].feature == feature_count)
        {
          ++end;
          //          std::cout << acc<< " "<< variable_indices[end+1].prot_acc << " "<< feature_count << " "<<variable_indices[end+1].feature << std::endl;
        }
      //std::cout << "start "<<start << " end "<<end <<std::endl;
      // enter new x_j variable, which indicates if a (hypothetical) peptide is part of the solution, independent from its RT region (which is given by the x_js-variables)
      Int index = model_->addColumn();
      IndexTriple triple;
      triple.prot_acc = acc;
      triple.feature = feature_count;
      triple.scan = -1; // to show this is not a x_j,s-variable
      triple.variable = index;
      variable_indices.push_back(triple);
      model_->setColumnBounds(index,0.,1.,LPWrapper::DOUBLE_BOUNDED);
      model_->setColumnType(index,LPWrapper::BINARY);
      model_->setColumnName(index,String("x_")+acc+String("_")+String(feature_count));
      // ensure that x_j is 1 iff at least one of the x_js-variables is 1
      for(Size f = start; f <= end;++f)
        {
          std::vector<DoubleReal> entries(2);
          std::vector<Int> indices(2);
          entries[0] = 1.;
          entries[1] = -1.;
          indices[0] = variable_indices[f].variable;
          indices[1] = triple.variable;
          String name = "x_" + String(i) + String(",") + String(variable_indices[f].scan) + "_x_" + String(i);
          model_->addRow(indices,entries,name,0,0,LPWrapper::UPPER_BOUND_ONLY);
        }
      std::cout << "added variable and constraint, now set i to "<<end<<std::endl;
      i = end;
		}
  std::cout << "now add the actual max list size constraint\nmax_list_size:"<<max_list_size<< "\n"<<model_->getNumberOfColumns() << " - "<< old_num<<std::endl;
  // now add constraints that controls the size of the inclusion list
  std::vector<DoubleReal> entries(model_->getNumberOfColumns() - old_num);
  std::vector<Int> indices(model_->getNumberOfColumns() - old_num);
  for(Int i = old_num; i < model_->getNumberOfColumns(); ++i)
    {
      entries[i - old_num] = 1.;
      indices[i - old_num] = i;
      //      std::cout << "indices["<<i - old_num<<"]="<<i<<std::endl;
    }
  model_->addRow(indices,entries,"max_list_size",0,max_list_size,LPWrapper::UPPER_BOUND_ONLY);
}


void PSLPFormulation::addRTBinCapacityConstraint_(std::vector<IndexTriple>& variable_indices,
                                                  Size max_rt_index,UInt ms2_spectra_per_rt_bin,bool sequential_order)
{
	//////////////////////////////////////////////////////////////////////////
	// constraint : rt_bin_capacity
	//////////////////////////////////////////////////////////////////////////
	// sort variable_indices according to their scan number
	sort(variable_indices.begin(),variable_indices.end(),ScanLess());
	Size j = 0;
	for(Size i = 0; i < max_rt_index;++i)
		{
      //      std::cout << "RT cap, scan "<< i << std::endl;
			// first determine number of indices:
			Size start = j;
			while(j < variable_indices.size())
				{
          if(variable_indices[j].scan <0)
            {
              ++j;
              continue;
            }
          if((Size)variable_indices[j].scan != i) break;
					++j;
				}
			// no feature occuring in this scan
			if(start == j) continue;

			Size stop = j;
			Size c = 0;			
      std::vector<double> entries(stop-start);
      std::vector<int> indices(stop-start);
			for(Size s = start; s < stop; ++s)
				{
					entries[c] = 1.;
					indices[c] = (int)  variable_indices[s].variable;
#ifdef DEBUG_OPS
          std::cout << "indices["<<c<<"]= "<<indices[c]<<std::endl;
#endif
					++c;
				}
      //#ifdef DEBUG_OPS
      std::cout << "add row with "<< entries.size()<< " indices"<<std::endl;
      //#endif
      if(sequential_order && i!=0)
        {
          model_->addRow(indices,entries,(String("RT_CAP")+i),0,0,LPWrapper::FIXED);// fix capacities of all other spectra than the 1st to 0
        }
      else
        {
          model_->addRow(indices,entries,(String("RT_CAP")+i),0,ms2_spectra_per_rt_bin,LPWrapper::UPPER_BOUND_ONLY);// only upper bounded problem -> lower bound is ignored
        }
#ifdef DEBUG_OPS
			std::cout << "added row"<<std::endl;
#endif
			
		}

}


void PSLPFormulation::addProteinCoverageConstraint_(std::vector<IndexTriple>& variable_indices,
                                                    PrecursorIonSelectionPreprocessing& preprocessing,
                                                    std::map<String,Size> protein_variable_index_map)
{
  //  DoubleReal min_prot_prob = param_.getValue("thresholds:min_protein_id_probability");
    
  std::cout << "and now the protein coverage"<<std::endl;
  const std::map<String, std::vector<DoubleReal> >& pt_prot_map = preprocessing.getProteinPTMap();
  std::map<String, std::vector<DoubleReal> >::const_iterator map_iter = pt_prot_map.begin();
  sort(variable_indices.begin(),variable_indices.end(),VariableIndexLess());
  Size feature_counter = 0;
  for(; map_iter != pt_prot_map.end();++map_iter)
    {
      std::cout << "protein: "<<map_iter->first<<std::endl;
      std::vector<Int> indices_vec;
      std::vector<double> entries;
      std::vector<DoubleReal>::const_iterator f_index_iter = map_iter->second.begin();
      //const std::vector<DoubleReal>& masses =  preprocessing.getMasses(map_iter->first);
      // go through all feature that have ids belonging to this protein
      for(; f_index_iter != map_iter->second.end(); ++f_index_iter)
        {
          // now go through all x_variable for this feature
          for(Size v = 0; v < variable_indices.size();++v)
            {
              if(variable_indices[v].prot_acc == map_iter->first && (Int)variable_indices[v].feature == distance(map_iter->second.begin(),f_index_iter) )
                {
                  // if there are duplicates in the vector CoinModel will abort
                  if(find(indices_vec.begin(),indices_vec.end(),variable_indices[v].variable) == indices_vec.end())
                    {
                      
                      indices_vec.push_back((Int)variable_indices[v].variable);
                      DoubleReal dt = map_iter->second[distance(map_iter->second.begin(),f_index_iter)];
                      if(fabs(1.-dt) <  0.000001)
                        {
                          entries.push_back(log(0.000001));
                          //entries.push_back(-log(0.000001) / log(1.-min_prot_prob));
                        }
                      else
                        {
                          //entries.push_back(-log(1.-dt) / log(1.-min_prot_prob));
                          entries.push_back(log(1.-dt));
                        }
                    
                    }
                }
              else if(variable_indices[v].feature > feature_counter)  break; // the indices are sorted, hence if the current index is larger, we are finished
            }
          ++feature_counter;
        }
 //  if(indices_vec.size() == 0) continue;
      // if(indices_vec.size() < 2)
      //            {
      //              std::cout << "too few features with ids for this protein, skipping protein"<<std::endl;
      //              continue;
      //            }

      // enter protein variable
      indices_vec.push_back(protein_variable_index_map[map_iter->first]);
      entries.push_back(1.);
      
      Int i = distance(pt_prot_map.begin(),map_iter);
#ifdef DEBUG_OPS
      std::cout << "\nadd row "<<std::endl;
      std::cout << (String("PROT_COV_")+map_iter->first) <<"\t"<<(String("PROT_COV_")+i).c_str() <<std::endl;
      std::cout << indices_vec.size()<< " "<<(indices_vec[0])<< " "<<(entries[0])
                << " min "<<-COIN_DBL_MAX<<", max "<<log(1.-min_prot_coverage)<<std::endl;
#endif
      // at the moment we want a certain coverage for each protein
      model_->addRow(indices_vec,entries,(String("PROT_COV_")+i),0.,0.,LPWrapper::UPPER_BOUND_ONLY);
      // test for inclusion list protein based
//       cmodel_->addRow((int)indices_vec.size(),&(indices_vec[0]),&(entries[0]),
//                      -COIN_DBL_MAX,COIN_DBL_MAX,(String("PROT_COV_")+i).c_str());
      std::cout << "\nadded row "<<std::endl;
#ifdef DEBUG_OPS
      std::cout << "added row"<<std::endl;
#endif

    }
}

void PSLPFormulation::assembleInclusionListForProteinBasedLP_(std::vector<IndexTriple>& variable_indices,FeatureMap<>& precursors,
                                                              std::vector<int>& solution_indices,
                                                              PrecursorIonSelectionPreprocessing& preprocessing)
{
  DoubleReal min_rt = param_.getValue("rt:min_rt");
  DoubleReal max_rt = param_.getValue("rt:max_rt");
  DoubleReal rt_step_size = param_.getValue("rt:rt_step_size");
  Size max_index = (Size)ceil((max_rt-min_rt) / rt_step_size);
  //  const std::map<String, std::vector<DoubleReal> >& pt_prot_map = preprocessing.getProteinPTMap();
  const std::map<String, std::vector<String> >& prot_pep_seq_map = preprocessing.getProteinPeptideSequenceMap();
  std::cout << solution_indices.size()<< " out of "<<variable_indices.size()<<" possible precursors were chosen."<<std::endl;
  sort(variable_indices.begin(),variable_indices.end(),VariableIndexLess());
  sort(solution_indices.begin(),solution_indices.end());
  for(Size idx = 0; idx < solution_indices.size();++idx)
    {
      std::cout << idx<<" "<< variable_indices.size()<<" "<<solution_indices[idx]<<std::endl;
      Size var_counter = 0;
      while(var_counter < variable_indices.size() && (Int)variable_indices[var_counter].variable != solution_indices[idx])
        {
          ++var_counter;
        }
      std::cout <<var_counter << " "<<variable_indices[var_counter].feature<<std::endl;
      if(var_counter ==  variable_indices.size() ) {std::cout << "varaible index not found..skipping"<<std::endl;continue;}
      std::cout << variable_indices[var_counter]<<std::endl;

      if(variable_indices[var_counter].scan == -1) continue; // we don't enter x_j variables

      // now get all variables belonging to the same peptide
      String acc = variable_indices[var_counter].prot_acc;
      Size feature = variable_indices[var_counter].feature;
      Size start = var_counter;
      Size end = var_counter+1;
      // we start at var_counter as the solution indices are sorted in increasing order
      for(; end < variable_indices.size();++end)
        {
          if(!(variable_indices[end].prot_acc == acc && variable_indices[end].feature == feature))
            {
              break;
            }
        }
      Feature tmp_feat;
      tmp_feat.setMZ(preprocessing.getMasses(acc)[feature]);
      // set RT to nearest spectrum to the predicted rt
      DoubleReal rt = preprocessing.getRT(acc,feature);
      // first find closest rt
      Size rt_index = std::min((Size)std::max(0.,ceil( (rt-min_rt)/rt_step_size)) ,max_index);
      DoubleReal curr_rt = min_rt + rt_index * rt_step_size;
      tmp_feat.setRT(curr_rt);
      std::cout << "now add the convex hulls"<<std::endl;
      std::vector< DPosition<2> > points;
      Size max_sol_index = idx;
      for(Size i = start; i < end; ++i)
        {
          if(variable_indices[i].scan == -1) continue; // we don't enter x_j variables
          // now check if this index is in solution
          std::vector<Int>::iterator iter = find(solution_indices.begin(),solution_indices.end(),variable_indices[i].variable);
          if( iter == solution_indices.end()) continue;
          // we need to remember the index in the solution_indices
          else if(distance(solution_indices.begin(),iter) > (Int)max_sol_index) max_sol_index = distance(solution_indices.begin(),iter);
          points.push_back(DPosition<2>(min_rt + variable_indices[i].scan * rt_step_size,tmp_feat.getMZ()-0.1));
          points.push_back(DPosition<2>(min_rt + variable_indices[i].scan * rt_step_size,tmp_feat.getMZ()+3.));
          
        }
      ConvexHull2D hull;
      hull.addPoints(points);
      tmp_feat.getConvexHulls().push_back(hull);
      std::cout << "now add sequence and protein acc"<<std::endl;
      tmp_feat.setMetaValue("protein",acc);
      tmp_feat.ensureUniqueId();
      std::map<String, std::vector<String> >::const_iterator prot_pep_seq_map_iter = prot_pep_seq_map.find(acc);
      if( prot_pep_seq_map_iter != prot_pep_seq_map.end())
        {
          tmp_feat.setMetaValue("sequence",prot_pep_seq_map_iter->second[feature]);
        }
      std::cout <<"added"<<std::endl;
      precursors.push_back(tmp_feat);
          
      idx = max_sol_index;
      std::cout << tmp_feat.getRT() <<" "<<tmp_feat.getMZ() <<" "<<precursors.size()<<std::endl;
    }
  std::cout << "all solution indices added "<<precursors.size()<<std::endl;
  precursors.ensureUniqueId();
}


} // namespace OpenMS
