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

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/DATASTRUCTURES/LPWrapper.h>
#include <OpenMS/DATASTRUCTURES/MassExplainer.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/SYSTEM/StopWatch.h>



namespace OpenMS {

	ILPDCWrapper::ILPDCWrapper(){}

	ILPDCWrapper::~ILPDCWrapper(){}

	DoubleReal ILPDCWrapper::compute(const FeatureMap<> fm, PairsType& pairs, Size verbose_level)
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
    for	(SignedSize i=0; i < (SignedSize)bins.size(); ++i)
		{
      /*PairsType p2 = pairs;
      writeProblem_(fm, p2, bins[i].first, bins[i].second, verbose_level, String("Z:\\myDocuments\\analysis\\DC_Lagrange\\output\\p") + String(i) + ".txt");
      StopWatch t1,t2;
      t1.start();*/
  		DoubleReal s = computeSlice_(fm, pairs, bins[i].first, bins[i].second, verbose_level);
      score += s;
      //t1.stop();
      //std::cerr << "CoinOr Slice: " << i << "  time: "<< 	t1.getClockTime () << "  score: " << s <<std::endl;
		}
    //score = computeSlice_(fm, pairs, 0, pairs.size()); // all at once - no bins

		return score;
	}

  

  void ILPDCWrapper::writeProblem_(const FeatureMap<> fm,
                                    PairsType& pairs, 
                                    const PairsIndex margin_left, 
                                    const PairsIndex margin_right,
                                    Size verbose_level,
                                    String filename)
  {
    // feature idx --> rotamers set 
    typedef std::map<Size, std::set<String> > r_type;
    r_type residues;

    for (PairsIndex i = margin_left; i < margin_right; ++i)
    {
      DoubleReal score = exp(getLogScore_(pairs[i], fm));
      pairs[i].setEdgeScore(score * pairs[i].getEdgeScore()); // multiply with preset score

      String rota_l =  String(pairs[i].getElementIndex(0)) + pairs[i].getCompomer().getAdductsAsString(0) + "_" + pairs[i].getCharge(0);
      residues[pairs[i].getElementIndex(0)].insert(rota_l);
      String rota_r =  String(pairs[i].getElementIndex(1)) + pairs[i].getCompomer().getAdductsAsString(1) + "_" + pairs[i].getCharge(1);
      residues[pairs[i].getElementIndex(1)].insert(rota_r);
    }

    // map featurenumber to Residue index:
    std::map<Size, Size> f_2_r;
    Size count(0);
    for (r_type::iterator it=residues.begin(); it!= residues.end(); ++it)
    {
      f_2_r[it->first] = count++;
    }

    TextFile f;
    f.push_back("[NUMBER_OF_RESIDUES]");
    f.push_back(String(residues.size()) + "\n");
    f.push_back("\n[NUMBER_OF_ROTAMERS]");
    for (r_type::iterator it=residues.begin(); it!= residues.end(); ++it)
    {
      f.push_back(String(f_2_r[it->first]) + " " + String(it->second.size()));
    }
    f.push_back("\n[SELF_ENERGIES]");
    for (r_type::iterator it=residues.begin(); it!= residues.end(); ++it)
    {
      for (Size i=0;i<it->second.size();++i)
      {
        f.push_back(String(f_2_r[it->first]) + " " + i + " 0");
      }
    }
    f.push_back("\n[PAIRWISE_ENERGIES]");

    // assign index to "Rotamer"(=adduct composition) to find it later for pairwise energies
    typedef std::map<String, Size> Cmp2Idx;
    Cmp2Idx cmpidx;
    for (r_type::iterator it=residues.begin(); it!= residues.end(); ++it)
    {
      Size count(0);
      for (std::set<String>::iterator it2=it->second.begin(); it2!=it->second.end(); ++it2)
      {
        cmpidx[*it2] = count++;
      }
    }
    // fill data
    for (PairsIndex i=margin_left; i<margin_right; ++i)
    {
      String rota_l =  String(pairs[i].getElementIndex(0)) + pairs[i].getCompomer().getAdductsAsString(0) + "_" + pairs[i].getCharge(0);
      String rota_r =  String(pairs[i].getElementIndex(1)) + pairs[i].getCompomer().getAdductsAsString(1) + "_" + pairs[i].getCharge(1);
      f.push_back(String(f_2_r[pairs[i].getElementIndex(0)]) + " " + cmpidx[rota_l] + " " + String(f_2_r[pairs[i].getElementIndex(1)]) + " " + cmpidx[rota_r] + " " + -pairs[i].getEdgeScore());
    }

    f.store(filename);

  }


	DoubleReal ILPDCWrapper::computeSlice_(const FeatureMap<> fm,
																			   PairsType& pairs, 
																			   const PairsIndex margin_left, 
																			   const PairsIndex margin_right,
                                         Size verbose_level)
	{
		LPWrapper build;
    //build.setSolver(LPWrapper::SOLVER_GLPK);
    build.setObjectiveSense(LPWrapper::MAX); // maximize

		//------------------------------------objective function-----------------------------------------------
		// find maximal objective value
		DoubleReal score = 0;
		DoubleReal score_min=10e10f, score_max = -10e10f;

		// fill in objective values 
		std::ostringstream namebuf;

		for (PairsIndex i=margin_left; i<margin_right; ++i)
		{
		
			// log scores are good for addition in ILP - but they are < 0, thus not suitable for maximizing
			// ... so we just add normal probabilities...
      score = exp(getLogScore_(pairs[i], fm));
			pairs[i].setEdgeScore(score * pairs[i].getEdgeScore()); // multiply with preset score
			namebuf.str("");
			namebuf<<"x#"<<i;
			// create the new variable object
      Size index = build.addColumn();
      build.setColumnBounds(index,0,1,LPWrapper::DOUBLE_BOUNDED);
			build.setColumnType(index,LPWrapper::INTEGER); // integer variable
			build.setObjective(index, pairs[i].getEdgeScore());
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
					String s = String("C") + i + "." + j;

          // Now build rows: two variables, with indices 'columns', factors '1', and 0-1 bounds.
          std::vector<double> element(2,1.0);
					std::vector<int> columns;
          columns.push_back(int(i-margin_left));
          columns.push_back(int(j-margin_left));
					build.addRow(columns, element,s, 0., 1.,LPWrapper::DOUBLE_BOUNDED);
				}
			}
		}
    
		if (verbose_level > 2) LOG_INFO << "node count: " << fm.size() << "\n";
		if (verbose_level > 2) LOG_INFO << "edge count: " << pairs.size() << "\n";
		if (verbose_level > 2) LOG_INFO << "constraint count: " << (conflict_idx[0]+conflict_idx[1]+conflict_idx[2]+conflict_idx[3]) << " = " << conflict_idx[0] << " + " << conflict_idx[1] << " + " << conflict_idx[2] << " + " << conflict_idx[3] << "(0 or inferred)" << std::endl;
	
		//---------------------------------------------------------------------------------------------------------
		//----------------------------------------Solving and querying result--------------------------------------
		//---------------------------------------------------------------------------------------------------------

    if (verbose_level > 0) LOG_INFO << "Starting to solve..." << std::endl;
    LPWrapper::SolverParam param;
    param.enable_mir_cuts= true;
    param.enable_cov_cuts= true;
    param.enable_feas_pump_heuristic=true;
    param.enable_binarization=false;
    param.enable_clq_cuts =true;
    param.enable_gmi_cuts = true;
    param.enable_presolve = true;
    StopWatch time1;
    time1.start();
    build.solve(param);
    time1.stop();
		if (verbose_level > 0) LOG_INFO << " Branch and cut took " << time1.getClockTime() << " seconds, "
						                        << " with objective value: " << build.getObjectiveValue() << "."
						                        << " Status: " << (!build.getStatus() ? " Finished" : " Not finished")
						                        << std::endl;

		/* variable values */
		UInt active_edges = 0;
		Map <String, Size> count_cmp;

      for (UInt iColumn = 0; iColumn < build.getNumberOfColumns(); ++iColumn)
        {
          double value=build.getColumnValue(iColumn);
          if (fabs(value)>0.5)
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
		if (verbose_level > 2) LOG_INFO << "Active edges: " << active_edges << " of overall " << pairs.size() << std::endl;

		for (Map < String, Size >::const_iterator it=count_cmp.begin(); it != count_cmp.end(); ++it)
		{
			//std::cout << "Cmp " << it->first << " x " << it->second << "\n";
		}

		DoubleReal opt_value = build.getObjectiveValue();

		//objective function value of optimal(?) solution
		return opt_value;
	} // !compute

	DoubleReal ILPDCWrapper::getLogScore_(const PairsType::value_type& pair, const FeatureMap<>& fm)
	{
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

