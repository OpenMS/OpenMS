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

#ifndef OPENMS_ANALYSIS_DECHARGING_FEATUREDECONVOLUTION_H
#define OPENMS_ANALYSIS_DECHARGING_FEATUREDECONVOLUTION_H

// OpenMS
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/DATASTRUCTURES/ChargePair.h>
#include <OpenMS/DATASTRUCTURES/Compomer.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/MassExplainer.h>
#include <OpenMS/ANALYSIS/DECHARGING/MIPWrapper.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

#include <OpenMS/FORMAT/ConsensusXMLFile.h> // tmp

// STL
#include <vector>

#include <iostream>
#include <fstream>

namespace OpenMS
{
  /** 
    @brief An algorithm to decharge features (i.e. as found by FeatureFinder).
    
		@ref FeatureDeconvolution_Parameters are explained on a separate page.
    
    @ingroup Analysis
  */
  
  class FeatureDeconvolution : public DefaultParamHandler
  {
    public:
    
      typedef FeatureMap<> FeatureMapType;
      typedef Feature FeatureType;
      typedef DPosition<2> ClusterPointType;
			typedef FeatureMapType::FeatureType::CoordinateType CoordinateType;
			typedef MIPWrapper::PairsType PairsType;
			//typedef std::vector<ChargePair> PairsType;
			
      /** @name Constructors and Destructors
      */
      //@{
      /// default constructor
      FeatureDeconvolution()
      : DefaultParamHandler("FeatureDeconvolution"),
        potential_adducts_()
      {
        defaults_.setValue("charge_min", 1, "minimal possible charge");
        defaults_.setValue("charge_max", 10, "maximal possible charge");
				// TODO how to enforce this later after pairs merged? build pairs only for span=2?
        defaults_.setValue("charge_span_max", 4, "maximal range of charges for a single analyte, i.e. observing q1=[5,6,7] implies span=3");
        defaults_.setValue("mass_max", 25000.0, "maximal mass (not m/z!) expected", StringList::create("advanced"));
        defaults_.setValue("retention_max_diff", 0.1, "maximum allowed RT difference between between two co-features");
				/// TODO should be m/z difference?!
        defaults_.setValue("mass_max_diff", 0.5, "maximum allowed mass difference between between two co-features");
        defaults_.setValue("potential_adducts", StringList::create("H+:0.7,K+:0.1,NH4+:0.1,Na+:0.1"), "Adducts used to explain mass differences in format: 'Element(+)*:Probablity', i.e. the number of '+' indicate the charge, e.g. 'Ca++:0.5' indicates +2. Probabilites are normalized to one.");
        
        defaults_.setValue("max_minority_bound", 2, "maximum count of the least probable adduct (according to 'potential_adducts' param) within a charge variant. E.g. setting this to 2 will not allow an adduct composition of '1(H+),3(Na+)' if Na+ is the least probable adduct");
        defaults_.setMinInt("max_minority_bound", 0);
        
        defaults_.setValue("min_rt_overlap", 0.66, "Minimum overlap of the convex hull' RT intersection measured against the union from two features (if CHs are given)");
        defaults_.setMinFloat("min_rt_overlap", 0);
        defaults_.setMaxFloat("min_rt_overlap", 1);
        
        defaultsToParam_();
      }


      void updateMembers_()
      {
        StringList potential_adducts_s = param_.getValue("potential_adducts");
        potential_adducts_.clear();
				
				float prob_sum=0;
				
        for (StringList::const_iterator it=potential_adducts_s.begin(); it != potential_adducts_s.end(); ++it)
        {
          StringList adduct;
          it->split(':', adduct);
          if (adduct.size()!=2)
          {
            String error = "FeatureDeconvolution::potential_adducts (" + (*it) + ") does not have two entries ('Element:Probability')!";
            throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, error);
          }
          String adduct_formula = adduct[0];
					// determine charge of adduct (by # of '+')
					Size l_charge = adduct[0].size();
					l_charge -= adduct[0].remove('+').size();
					// determine probability
          float prob = adduct[1].toFloat();
          if (prob > 1.0 || prob<0.0)
					{
            String error = "FeatureDeconvolution::potential_adducts (" + (*it) + ") does not have a proper probablity (" + String(prob) + ") in [0,1]!";
            throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, error);
          }
          EmpiricalFormula ef(adduct_formula);
          Adduct a((Int)l_charge, 1, ef.getMonoWeight(), ef.getString(), log(prob));
          prob_sum += prob;
          //std::cout << "FeatureDeconvolution: inserting potential adduct " << ef.getString() << "[q:" << l_charge << ", pr:" << prob << "(" << a.getLogProb() << ")" << "]\n";
          potential_adducts_.push_back(a);
        }
        
        // normalize probablity to 1
        for (Size i=0; i<potential_adducts_.size(); ++i)
        {
					potential_adducts_[i].setLogProb(log(exp(potential_adducts_[i].getLogProb()) / prob_sum));
        }
        
      }
      /// Copy constructor
      inline FeatureDeconvolution(const FeatureDeconvolution& source)
          : DefaultParamHandler(source),
            potential_adducts_(source.potential_adducts_)
      {}

      /// Assignment operator
      inline FeatureDeconvolution& operator=(const FeatureDeconvolution& source)
      {
        if (&source==this)
        {
          return *this;
        }

        DefaultParamHandler::operator=(source);
        potential_adducts_ = source.potential_adducts_;

        return *this;          
      };

        

      /// destructor
      virtual ~FeatureDeconvolution() 
      {};
      //@}    

      /// compute a zero-charge feature map from a set of charged features (@p map)
			/// i.e. find putative ChargePairs, then score them and hand over to ILP
			/// @param map Input feature-map
			/// @param map Output feature-map (sorted by position and augmented with user params)
			/// @param cons_map [out] Output of grouped features belonging to a charge group
      void compute(const FeatureMapType &fm_in, FeatureMapType &fm_out, ConsensusMap &cons_map, ConsensusMap &cons_map_p) 
      {
      
				ConsensusMap cons_map_p_neg; // tmp
				
        Int q_min = param_.getValue("charge_min");
				Int q_max = param_.getValue("charge_max");
				Int q_span = param_.getValue("charge_span_max");

				DoubleReal rt_diff_max = param_.getValue("retention_max_diff");
				DoubleReal mz_diff_max = param_.getValue("mass_max_diff");
				
				DoubleReal rt_min_overlap = param_.getValue("min_rt_overlap");
				
				// sort by RT and then m/z
				fm_out = fm_in;
				fm_out.sortByPosition();

        
        // search for most & least probable adduct to fix p threshold
        DoubleReal adduct_lowest_log_p = log(1.0);
        DoubleReal adduct_highest_log_p = log(0.0000000001);
        for (Size i=0; i<potential_adducts_.size(); ++i)
        {
					adduct_lowest_log_p  = std::min(adduct_lowest_log_p , potential_adducts_[i].getLogProb() );
					adduct_highest_log_p = std::max(adduct_highest_log_p, potential_adducts_[i].getLogProb() );
        }
        Int max_minority_bound = param_.getValue("max_minority_bound");
				DoubleReal thresh_logp = adduct_lowest_log_p * max_minority_bound +
																 adduct_highest_log_p * std::max(q_max - max_minority_bound, 0);


				// create mass difference list

				std::cout << "Generating Masses with threshold: " << thresh_logp << " ...\n";
        MassExplainer me(potential_adducts_, q_min, q_max, q_span, thresh_logp);
				me.compute();
				std::cout << "done\n";
				
				// holds query results for a mass difference
				MassExplainer::CompomerIterator md_s, md_e;
				Compomer null_compomer(0,0,-std::numeric_limits<DoubleReal>::max());
				SignedSize hits = 0;
				
				CoordinateType mz1,mz2, m1;
				
				Size possibleEdges = 0, overallHits = 0;
				
				PairsType feature_relation;
				
        for (Size i_RT = 0; i_RT < fm_out.size(); ++i_RT)
        { // ** RT-sweepline
					
					mz1 = fm_out[i_RT].getMZ();
					
					for (Size i_RT_window = i_RT + 1
							 ; (i_RT_window < fm_out.size())
							   &&((fm_out[i_RT_window].getRT() - fm_out[i_RT].getRT()) < rt_diff_max) 
							 ; ++i_RT_window)
					{ // ** RT-window
						
						// knock-out criterion first: RT overlap
						// use sorted structure and use 2nd start--1stend / 1st start--2ndend
						const Feature& f1 = fm_out[i_RT];
						const Feature& f2 = fm_out[i_RT_window];
						
						if (! ( f1.getConvexHull().getBoundingBox().isEmpty() || f2.getConvexHull().getBoundingBox().isEmpty() ) )
						{
							DoubleReal f_start1 = std::min(f1.getConvexHull().getBoundingBox().minX(),f2.getConvexHull().getBoundingBox().minX());
							DoubleReal f_start2 = std::max(f1.getConvexHull().getBoundingBox().minX(),f2.getConvexHull().getBoundingBox().minX());
							DoubleReal f_end1 = std::min(f1.getConvexHull().getBoundingBox().maxX(),f2.getConvexHull().getBoundingBox().maxX());
							DoubleReal f_end2 = std::max(f1.getConvexHull().getBoundingBox().maxX(),f2.getConvexHull().getBoundingBox().maxX());
							
							DoubleReal union_length = f_end2 - f_start1;
							DoubleReal intersect_length = std::max(0., f_end1 - f_start2);
							
							if (intersect_length / union_length < rt_min_overlap) continue;
						}
	
						// start guessing charges ...
	
						mz2 = fm_out[i_RT_window].getMZ();
						
						//@improvement TODO only check for charges which approx have the same quotient as mz1/mz2
						for (Int q1=q_min; q1<=q_max; ++q1)
						{ // ** q1

              //DEBUG:
              if (fm_out[i_RT_window].getRT()>1930.08 && fm_out[i_RT_window].getRT()<1931.2 && mz1>1443 && mz2>1443 && mz1<2848 && mz2<2848)
              {
                std::cout << "we are at debug location\n" << fm_out[i_RT_window].getRT() <<"   : " << mz1 << "; " << mz2 << "\n";
              }
              if (i_RT == 930 && i_RT_window == 931)
              {
                std::cout << "we are at debug location\n" << fm_out[i_RT_window].getRT() <<"   : " << mz1 << "; " << mz2 << "\n";
              }
              // \DEBUG

								m1 = mz1*q1;
								// additionally: forbid q1 and q2 with distance greater than q_span
								for (Int q2= std::max(q_min, q1-q_span+1)
								     ; (q2<=q_max) && (q2<=q1+q_span-1)
										 ; ++q2)
								{ // ** q2

									++possibleEdges; // internal count, not vital
									
									// find possible adduct combinations
									hits = me.query(q2-q1, mz2*q2-m1, mz_diff_max, thresh_logp, md_s, md_e);
									OPENMS_PRECONDITION(hits>=0,"FeatureDeconvolution querying #hits got negative result!");

									overallHits+=hits;
									// choose most probable hit (TODO think of something clever here)
									// for now, we take the one that has highest p in terms of the compomer structure
									if (hits>0)
									{
										std::cout << "#hits: " << hits << "\n";
										//std::cout << "q2: " << q2 << " (mass to explain: " << (mz2*q2-m1) << ")\n";
										CoordinateType naive_mass_diff = mz2*q2-m1;
										Compomer best_hit = null_compomer;
										for (; md_s!=md_e; ++md_s)
										{
											//std::cout << "neg: " << md_s->getNegativeCharges() << " pos: " << md_s->getPositiveCharges() << " p: " << md_s->getLogP() << " \n";
											if (
													// compomer has better probability
													(best_hit.getLogP() < md_s->getLogP())
												 	&&
												 	// compomer fits charge assignment of left & right feature
												 	(q1 >= abs(md_s->getNegativeCharges()) && q2 >= abs(md_s->getPositiveCharges()))
												 )
											{
												best_hit = *md_s;
												/** testing: we just add every explaining edge 
													- a first estimate shows that 90% of hits are of |1|
													- the remaining 10% have |2|, so the additional overhead is minimal
												**/
												#if 1
												ChargePair cp(i_RT, i_RT_window, q1, q2, best_hit.getID(), naive_mass_diff - best_hit.getMass(), false);
												//std::cout << "CP # "<< feature_relation.size() << " :" << i_RT << " " << i_RT_window<< " " << q1<< " " << q2 << "\n";
												feature_relation.push_back(cp);
												#endif
											}
										}
										if (best_hit == null_compomer)
										{	
											std::cout << "FeatureDeconvolution.h:: could not find a compomer which complies with assumed q1 and q2 values!\n with q1: " << q1 << " q2: " << q2 << "\n";
										}
										else
										{
											// disabled while we add every hit (and not only the best - see above)
											#if 0
											ChargePair cp(i_RT, i_RT_window, q1, q2, best_hit.getID(), naive_mass_diff - best_hit.getMass(), false);
											//std::cout << "CP # "<< feature_relation.size() << " :" << i_RT << " " << i_RT_window<< " " << q1<< " " << q2 << "\n";
											feature_relation.push_back(cp);
											#endif
										}
									}
									
								} // q2
						} // q1
					} // RT-window
        } // RT sweepline


				std::cout << "found " << feature_relation.size() << " putative edges (of " << possibleEdges << ")"
									<< " and avg hit-size of " << (overallHits/feature_relation.size())
									<< std::endl;
				
				std::cout << "creating wrapper..." << std::endl;
				// forward set of putative edges to ILP
				MIPWrapper lp_wrapper;
				std::cout << "Done" << std::endl;
				// compute best solution
				DoubleReal ilp_score = lp_wrapper.compute(me, fm_out, feature_relation);
				std::cout << "score is: " << ilp_score << std::endl;
				
				// prepare output consensusMaps
				cons_map.setProteinIdentifications( fm_out.getProteinIdentifications() );
				cons_map_p.setProteinIdentifications( fm_out.getProteinIdentifications() );
				
				// some statistics about compomers used
				std::map < ChargePair, Int > compomer_stats;
				
				UInt agreeing_fcharge = 0;
				std::vector<Size> f_idx_v(2);
				Size aedges=0;

				String scores_clean_egde, scores_dirty_edge;
				EmpiricalFormula ef_clean_edge, ef_dirty_edge;

				// write related features
				for (Size i=0; i<feature_relation.size(); ++i)
				{
          if (!feature_relation[i].isActive()) continue;
					
					//std::cout << "active charge_pair #" << i << "\n";

          f_idx_v[0] = feature_relation[i].getElementIndex(0);
          f_idx_v[1] = feature_relation[i].getElementIndex(1);
          ++aedges;
					
					bool dirty = false;

					for (Size f_idx=0;f_idx<2;++f_idx)
					{
						// check if the local feature charges agree
						if (fm_out[f_idx_v[f_idx]].getCharge() == feature_relation[i].getCharge((UInt)f_idx))
						{
							++agreeing_fcharge;
						}
						else
						{
							DoubleReal rt_diff =  fabs(fm_out[feature_relation[i].getElementIndex(0)].getRT() - fm_out[feature_relation[i].getElementIndex(1)].getRT());
							std::cout << "conflict in Q:: RT:" << fm_out[f_idx_v[f_idx]].getRT() << " MZ:" << fm_out[f_idx_v[f_idx]].getMZ() << " Q:" << fm_out[f_idx_v[f_idx]].getCharge() << " PredictedQ:" << feature_relation[i].getCharge((UInt)f_idx) << "[[ dRT: " << rt_diff << " dMZ: " << feature_relation[i].getMassDiff() << " score:" << feature_relation[i].getEdgeScore() << " f#:" << f_idx_v[f_idx] << " ]]\n";
							dirty = true;
						}
					}
					
					Compomer c = me.getCompomerById( feature_relation[i].getCompomerId());
					EmpiricalFormula ef = (c.getAdductsAsString(-1)) + (c.getAdductsAsString(+1));
					
					
					// store score distribution:
					if (!dirty)
					{
						scores_clean_egde += String(" ") + feature_relation[i].getEdgeScore();
						ef_clean_edge += ef;
					}
					else
					{
						scores_dirty_edge += String(" ") + feature_relation[i].getEdgeScore();
						ef_dirty_edge += ef;
					}
					
				}
				std::cout << "agreeing charges: " << agreeing_fcharge << "/" << (aedges*2) << std::endl;

				std::cout << "Edge score distribution (clean):\n" + scores_clean_egde + "\n(dirty)\n" + scores_dirty_edge + "\n\n";

				std::cout << "Edge emprirical formula (clean):\n" + ef_clean_edge.getString() + "\n(dirty)\n" + ef_dirty_edge.getString() + "\n\n";

				// fresh start for meta annotation
				for (Size i=0;i<fm_out.size(); ++i)
				{
					if (fm_out[i].metaValueExists("dc_charge_adducts")) fm_out[i].removeMetaValue("dc_charge_adducts");
				}
				
				// write groups to consensusXML (with type="charge_groups")

				// **find cliques from pairs
				// find which featureIdx maps to which consensusFeatureIdx
				// if no mapping is found, make a new CF.
				// if new pair spans two existing CFs -> merge CFs
				typedef std::map<Size, Size> Map;
				typedef Map::const_iterator MapCI;
				Map clique_register; 
				MapCI clique_it;
				
				for (Size i=0; i<feature_relation.size(); ++i)
				{
          Size f0_idx = feature_relation[i].getElementIndex(0);
          Size f1_idx = feature_relation[i].getElementIndex(1);
          
          Int old_q0 = fm_out[f0_idx].getCharge();
          Int old_q1 = fm_out[f1_idx].getCharge();

          if (feature_relation[i].isActive())
          {
            //std::cout << "feature #" << f0_idx << " #" << f1_idx << " ACTIVE with RT: " << fm_out[f1_idx].getRT() << "\n";

            //
						// annotate the affected features
						// ... and check consistency
						//
						Compomer c = me.getCompomerById( feature_relation[i].getCompomerId());
						// - left
						EmpiricalFormula ef_l(c.getAdductsAsString(-1));
						ef_l += "H" + String(feature_relation[i].getCharge(0) - c.getNegativeCharges());
						if (fm_out[f0_idx].metaValueExists("dc_charge_adducts"))
    				{
    					OPENMS_POSTCONDITION (ef_l.getString() == fm_out[f0_idx].getMetaValue("dc_charge_adducts"), "Decharging produced inconsistent adduct annotation!");
    				}
    				else
    				{
    					fm_out[f0_idx].setMetaValue("dc_charge_adducts", ef_l.getString());
    				}
						fm_out[f0_idx].setMetaValue("dc_charge_adduct_mass", ef_l.getMonoWeight());
						fm_out[f0_idx].setCharge(feature_relation[i].getCharge(0));
						
						// - right
						EmpiricalFormula ef_r(c.getAdductsAsString(+1));
						ef_r += "H" + String(feature_relation[i].getCharge(1) - c.getPositiveCharges());
						if (fm_out[f1_idx].metaValueExists("dc_charge_adducts"))
    				{
    					OPENMS_POSTCONDITION (ef_r.getString() == fm_out[f1_idx].getMetaValue("dc_charge_adducts"), "Decharging produced inconsistent adduct annotation!");
    				}
    				else
    				{
    					fm_out[f1_idx].setMetaValue("dc_charge_adducts", ef_r.getString());
    				}
						fm_out[f1_idx].setMetaValue("dc_charge_adduct_mass", ef_r.getMonoWeight());
						fm_out[f1_idx].setCharge(feature_relation[i].getCharge(1));
						
						//
						// create cliques
						//
            SignedSize target_cf0 = -1, target_cf1 = -1;

            // find the index of the ConsensusFeatures for the current pair
            if (clique_register.count(f0_idx) > 0)
            {
              target_cf0 = clique_register[f0_idx];
            }
            if (clique_register.count(f1_idx) > 0)
            {
              target_cf1 = clique_register[f1_idx];
            }

            ConsensusFeature cf(fm_out[f0_idx]);
            cf.insert(0,f0_idx, fm_out[f0_idx]);
            cf.insert(0,f1_idx, fm_out[f1_idx]);
            cf.setMetaValue("Local", String(old_q0)+":"+String(old_q1));
            cf.setMetaValue("CP", String(fm_out[f0_idx].getCharge())+"("+ String(fm_out[f0_idx].getMetaValue("dc_charge_adducts")) +"):"
																 +String(fm_out[f1_idx].getCharge())+"("+ String(fm_out[f1_idx].getMetaValue("dc_charge_adducts")) +") "
																 +String("Score: ") + feature_relation[i].getEdgeScore());
            //cf.computeDechargeConsensus(fm_out);
            #if 1
            // print pairs only
                cons_map_p.push_back(cf);
                cons_map_p.getFileDescriptions()[0].size = fm_out.size();
                cons_map_p.getFileDescriptions()[0].label = "charged features pairs";
            #endif

            // seen both features for the first time
            if ((target_cf0 == -1) &&
                (target_cf1 == -1))
            {//** new ConsensusFeature required
                cons_map.push_back(cf);
                clique_register[f0_idx] = cons_map.size()-1;
                clique_register[f1_idx] = cons_map.size()-1;
                //std::cout << "new: F" << f0_idx << " + F" << f1_idx << " are " << (cons_map.size()-1) << "\n";
            }
            else if (target_cf0 != target_cf1)
            {
              if (target_cf0 == -1)
              {//** add f0 to the already existing cf of f1
                cons_map[target_cf1].insert(0,f0_idx, fm_out[f0_idx]);
                clique_register[f0_idx] = target_cf1;
                //std::cout << "add: F" << f0_idx << " to " <<target_cf1 << " dueto F" << f1_idx << "\n";
              }
              else if (target_cf1 == -1)
              {//** add f1 to the already existing cf of f0
                cons_map[target_cf0].insert(0,f1_idx, fm_out[f1_idx]);
                clique_register[f1_idx] = target_cf0;
                //std::cout << "add: F" << f1_idx << " to " <<target_cf0 << " dueto F" << f0_idx << "\n";
              } else
              { //** conflict: the two elements of the pair already have separate CF´s --> merge
                // take every feature from second CF and: #1 put into first CF, #2 change registration with map
                ConsensusFeature::HandleSetType hst = cons_map[target_cf1].getFeatures();
                for (ConsensusFeature::HandleSetType::const_iterator it=hst.begin(); it!=hst.end();++it)
                { //** update cf_index
                  clique_register[it->getElementIndex()] = target_cf0;
                }
                // insert features from cf1 to cf0
                cons_map[target_cf0].insert(hst);
                // clear cf1; do NOT delete cf1 (will invalidate higher indices) - do that afterwards
                cons_map[target_cf1].clear();
                //std::cout << "conflict: F" << f0_idx << " + F" << f1_idx << " --> "<< target_cf0 << "(" << target_cf1 << " killed)" << "\n";
              }
            }

          }// active
          else
          {
            ConsensusFeature cf(fm_out[f0_idx]);
            cf.insert(0,f0_idx, fm_out[f0_idx]);
            cf.insert(0,f1_idx, fm_out[f1_idx]);
            cf.setMetaValue("Local", String(old_q0)+":"+String(old_q1));
            cf.setMetaValue("CP", String(fm_out[f0_idx].getCharge())+":"
																 +String(fm_out[f1_idx].getCharge()));
            //cf.computeDechargeConsensus(fm_out);
            #if 1
            // print pairs only
                cons_map_p_neg.push_back(cf);
                cons_map_p_neg.getFileDescriptions()[0].size = fm_out.size();
                cons_map_p_neg.getFileDescriptions()[0].label = "charged features pairs";
            #endif          
            //std::cout << "Not active: cp#" << i <<"\n";
          }

        } // !for

				// tmp
				ConsensusXMLFile cf_neg;
				//cf_neg.store("dc_pairs_neg.consensusXML", cons_map_p_neg);

				// remove empty ConsensusFeatures from map
				ConsensusMap cons_map_tmp(cons_map);
				cons_map_tmp.clear(); // keep other meta information (like ProteinIDs & Map)
				for (ConsensusMap::ConstIterator it = cons_map.begin(); it!=cons_map.end(); ++it)
				{
					// skip if empty
					if (it->getFeatures().size()==0) continue;
					cons_map_tmp.push_back(*it);					
					// set a centroid
					cons_map_tmp.back().computeDechargeConsensus(fm_out);

				}
				cons_map_tmp.swap(cons_map);

				Size singletons_count = 0;
        // include single features without a buddy!
        for (Size i=0; i<fm_out.size(); ++i)
				{
          // find the index of the ConsensusFeature for the current feature
          if (clique_register.count(i) > 0) continue;

          ConsensusFeature cf(fm_out[i]);
          cf.insert(0, i, fm_out[i]);

          cons_map.push_back(cf);
          cons_map.back().computeDechargeConsensus(fm_out);
          ++singletons_count;
        }
				
				std::cout << "Single features without charge ladder: " << singletons_count << " of " << fm_out.size() << "\n";
				
				
				// fill the header
				//cons_map.setMetaValue("meta",String("value"));
				//cons_map.setIdentifier("some lsid");
				cons_map.getFileDescriptions()[0].filename = "TODO - take from FeatureMAP.getLoadedFilePath () ";
				cons_map.getFileDescriptions()[0].size = fm_out.size();
				cons_map.getFileDescriptions()[0].label = "charged features";
				//cons_map.getFileDescriptions()[0].setMetaValue("meta",String("meta"));
				
				// TODO count Mass^->Charge violations, ie if charge is higher, mass^ should be higher as well
				// (mass^ is the peptide mass + adduct mass)
				
				
				
        return;
      }

    protected:
      /// List of adducts used to explain mass differences
      MassExplainer::AdductsType potential_adducts_;

  };
} // namespace OpenMS

#endif // OPENMS_ANALYSIS_DECHARGING_FEATUREDECONVOLUTION_H

