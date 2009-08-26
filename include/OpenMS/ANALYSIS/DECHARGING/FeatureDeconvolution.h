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
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/CONCEPT/Constants.h>

#include <OpenMS/FORMAT/ConsensusXMLFile.h> // tmp

// STL
#include <vector>

#include <iostream>
#include <fstream>

namespace OpenMS
{
  /** 
    @brief An algorithm to decharge features (i.e. as found by FeatureFinder).
    
    @htmlinclude OpenMS_FeatureDeconvolution.parameters	
    
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
          EmpiricalFormula ef(adduct[0].remove('+'));
          ef -= "H"+String(l_charge); ef.setCharge(l_charge); // effectively substract electron masses

          Adduct a((Int)l_charge, 1, ef.getMonoWeight(), adduct[0].remove('+'), log(prob));
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

				// edges				
				PairsType feature_relation;
				// for each feature, hold the explicit adduct type induced by egdes
				Map<Size, std::set<CmpInfo_> > feature_adducts;
						
				// # compomer results that either passed or failed the feature charge constraints
				Size no_cmp_hit=0, cmp_hit=0;
				
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
	
						// implicit adducts dont't cost anything
						Adduct proton(1, 1, Constants::PROTON_MASS_U, "H1", log(1.0));

						// start guessing charges ...
						mz2 = fm_out[i_RT_window].getMZ();
						
						//@improvement TODO only check for charges which approx have the same quotient as mz1/mz2
						for (Int q1=q_min; q1<=q_max; ++q1)
						{ // ** q1

              //DEBUG:
              /**if (fm_out[i_RT_window].getRT()>1930.08 && fm_out[i_RT_window].getRT()<1931.2 && mz1>1443 && mz2>1443 && mz1<2848 && mz2<2848)
              {
                std::cout << "we are at debug location\n" << fm_out[i_RT_window].getRT() <<"   : " << mz1 << "; " << mz2 << "\n";
              }
              if (i_RT == 930 && i_RT_window == 931)
              {
                std::cout << "we are at debug location\n" << fm_out[i_RT_window].getRT() <<"   : " << mz1 << "; " << mz2 << "\n";
              }*/
              // \DEBUG

								m1 = mz1*q1;
								// additionally: forbid q1 and q2 with distance greater than q_span
								for (Int q2= std::max(q_min, q1-q_span+1)
								     ; (q2<=q_max) && (q2<=q1+q_span-1)
										 ; ++q2)
								{ // ** q2

									++possibleEdges; // internal count, not vital
									
									//if (i_RT==439 && i_RT_window==442 && q1==4 && q2==3)
									//{
									//	std::cout << "DEBUG reached\n";
									//}
									
									// find possible adduct combinations
									CoordinateType naive_mass_diff = mz2*q2-m1;
									hits = me.query(q2-q1, naive_mass_diff, mz_diff_max, thresh_logp, md_s, md_e);
									OPENMS_PRECONDITION(hits>=0,"FeatureDeconvolution querying #hits got negative result!");

									overallHits+=hits;
									// choose most probable hit (TODO think of something clever here)
									// for now, we take the one that has highest p in terms of the compomer structure
									if (hits>0)
									{
										//std::cout << "#hits: " << hits << "\n";
										//std::cout << "q2: " << q2 << " (mass to explain: " << (mz2*q2-m1) << ")\n";
										
										Compomer best_hit = null_compomer;
										for (; md_s!=md_e; ++md_s)
										{
											//std::cout << "neg: " << md_s->getNegativeCharges() << " pos: " << md_s->getPositiveCharges() << " p: " << md_s->getLogP() << " \n";
											if (// compomer fits charge assignment of left & right feature
												 	(q1 >= md_s->getNegativeCharges()) && (q2 >= md_s->getPositiveCharges())
												 )
											{
												// compomer has better probability
												if (best_hit.getLogP() < md_s->getLogP())	best_hit = *md_s;

											
												/** testing: we just add every explaining edge 
													- a first estimate shows that 90% of hits are of |1|
													- the remaining 10% have |2|, so the additional overhead is minimal
												**/
												#if 1
												Compomer cmp = me.getCompomerById(md_s->getID());
												if ( ((q1 - cmp.getNegativeCharges()) % proton.getCharge() != 0) ||
														 ((q2 - cmp.getPositiveCharges()) % proton.getCharge() != 0))
												{
													std::cerr << "Cannot add enough default adduct (" << proton.getFormula() << ") to exactly fit feature charge! Next...)\n";
													continue;		
												}

												int hc_left  = (q1 - cmp.getNegativeCharges()) / proton.getCharge(); // this should always be positive! check!!
												int hc_right = (q2 - cmp.getPositiveCharges()) / proton.getCharge(); // this should always be positive! check!!
												
												
												if (hc_left < 0 || hc_right < 0)
												{
													throw Exception::Postcondition(__FILE__,__LINE__,__PRETTY_FUNCTION__,"WARNING!!! implicit number of H+ is negative!!! left:" + String(hc_left) + " right: " + String(hc_right) + "\n");
												}
												
												//if (i_RT==272 && i_RT_window==271)
												//{ std::cout << "DEBUG\n";}
												
												// get non-default adducts of this edge
												Compomer cmp_stripped(cmp.removeAdduct(proton));
												
												// save new adduct candidate
												if (cmp_stripped.getComponent()[Compomer::LEFT].size()>0)
												{
													String tmp = cmp_stripped.getAdductsAsString(Compomer::LEFT);
													CmpInfo_ cmp_left(tmp, feature_relation.size(), Compomer::LEFT);
													feature_adducts[i_RT].insert( cmp_left );
												}
												if (cmp_stripped.getComponent()[Compomer::RIGHT].size()>0)
												{
													String tmp = cmp_stripped.getAdductsAsString(Compomer::RIGHT);
													CmpInfo_ cmp_right(tmp, feature_relation.size(), Compomer::RIGHT);
													feature_adducts[i_RT_window].insert( cmp_right );
												}
												
												// add implicit H+ (if != 0)
												if (hc_left>0) cmp.add(proton*hc_left,Compomer::LEFT);
												if (hc_right>0)cmp.add(proton*hc_right,Compomer::RIGHT);
												
												ChargePair cp(i_RT, i_RT_window, q1, q2, cmp, naive_mass_diff - md_s->getMass(), false);
												//std::cout << "CP # "<< feature_relation.size() << " :" << i_RT << " " << i_RT_window<< " " << q1<< " " << q2 << " score: " << cp.getCompomer().getLogP() << "\n";
												feature_relation.push_back(cp);
												#endif
											}
										} // ! hits loop
										
										if (best_hit == null_compomer)
										{	
											//std::cout << "FeatureDeconvolution.h:: could not find a compomer which complies with assumed q1 and q2 values!\n with q1: " << q1 << " q2: " << q2 << "\n";
											++no_cmp_hit;
										}
										else
										{
											++cmp_hit;
											// disabled while we add every hit (and not only the best - see above)
											#if 0
											TODO if reactivated: add implicits (see above)
											ChargePair cp(i_RT, i_RT_window, q1, q2, me.getCompomerById(best_hit.getID()), naive_mass_diff - best_hit.getMass(), false);
											//std::cout << "CP # "<< feature_relation.size() << " :" << i_RT << " " << i_RT_window<< " " << q1<< " " << q2 << "\n";
											feature_relation.push_back(cp);
											#endif
										}
									}
									
								} // q2
						} // q1
					} // RT-window
        } // RT sweepline
        
        std::cout << no_cmp_hit << " of " << (no_cmp_hit+cmp_hit) << " valid net charge compomer results did not pass the feature charge constraints\n";
        
        inferMoreEdges_(feature_relation, feature_adducts);

				// DEBUG:
				//printEdgesOfConnectedFeatures_(271, 272, feature_relation);

				std::cout << "found " << feature_relation.size() << " putative edges (of " << possibleEdges << ")"
									<< " and avg hit-size of " << (overallHits/feature_relation.size())
									<< std::endl;

				// -------------------------- //
				// ** compute ILP solution ** //
				// -------------------------- //
				
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

				// -------------------------- //
				// **       DEBUG          ** //
				// -------------------------- //

				Map <Size, Size> features_aes, features_des; // count of adjacent active and dead edges
				UInt agreeing_fcharge = 0;
				std::vector<Size> f_idx_v(2);
				Size aedges=0;
				StringList scores_clean_edge, scores_dirty_edge;
				StringList scores_clean_edge_idx, scores_dirty_edge_idx;
				EmpiricalFormula ef_clean_edge, ef_dirty_edge;
				// find # edges (active and dead) for each feature
				for (Size i=0; i<feature_relation.size(); ++i)
				{
          Size f0_idx = feature_relation[i].getElementIndex(0);
          Size f1_idx = feature_relation[i].getElementIndex(1);
          if (feature_relation[i].isActive())
          {
						++features_aes[f0_idx];
						++features_aes[f1_idx];
          }
          else
          {
						++features_des[f0_idx];
						++features_des[f1_idx];
          }
				}
				TextFile out_dead;
				for (Size i=0; i<feature_relation.size(); ++i)
				{
          f_idx_v[0] = feature_relation[i].getElementIndex(0);
          f_idx_v[1] = feature_relation[i].getElementIndex(1);
					
					Compomer c = feature_relation[i].getCompomer();
          
          if (!feature_relation[i].isActive())
          {
						out_dead.push_back(String("dead e") + i + " (" + (c.getAdductsAsString(Compomer::LEFT)) + " -> " + (c.getAdductsAsString(Compomer::RIGHT)) + "): " 
							+ f_idx_v[0] + " (q_ff:" + fm_out[f_idx_v[0]].getCharge() + " q_de:" + feature_relation[i].getCharge(0) + ")"
							+ f_idx_v[1] + " (q_ff:" + fm_out[f_idx_v[1]].getCharge() + " q_de:" + feature_relation[i].getCharge(1) + ")"
							);
						continue;
          }
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
							//ConsensusFeature cf_in_conflict = cons_map[clique_register[f_idx_v[f_idx]]];
							//cf_in_conflict.computeDechargeConsensus(fm_out);
							std::cout << "conflict in f_Q! f_RT:" << fm_out[f_idx_v[f_idx]].getRT() << " f_MZ:" << fm_out[f_idx_v[f_idx]].getMZ() << " f_int:" << fm_out[f_idx_v[f_idx]].getIntensity() << " Q:" << fm_out[f_idx_v[f_idx]].getCharge() << " PredictedQ:" << feature_relation[i].getCharge((UInt)f_idx) << "[[ dRT: " << rt_diff << " dMZ: " << feature_relation[i].getMassDiff() << " score[" << i << "]:" << feature_relation[i].getEdgeScore() << " f#:" << f_idx_v[f_idx] << "(a" << features_aes[f_idx_v[f_idx]] << ":d" << features_des[f_idx_v[f_idx]] << ") ]]\n";
							dirty = true;
						}
					}
					
					EmpiricalFormula ef(c.getAdductsAsString(Compomer::LEFT) + (c.getAdductsAsString(Compomer::RIGHT)) );
					
					// store score distribution:
					if (!dirty)
					{
						scores_clean_edge.push_back(String(feature_relation[i].getEdgeScore()));
						scores_clean_edge_idx.push_back(String(i));
						ef_clean_edge += ef;
					}
					else
					{
						scores_dirty_edge.push_back(String(feature_relation[i].getEdgeScore()));
						scores_dirty_edge_idx.push_back(String(i));
						ef_dirty_edge += ef;
					}
					
				}
				
				//TODO out_dead.store("dead_edges.txt");
				std::cout << "agreeing charges: " << agreeing_fcharge << "/" << (aedges*2) << std::endl;
				std::cout << "Edge score distribution (clean):\n" + scores_clean_edge.concatenate(" ") + "\n(dirty)\n" + scores_dirty_edge.concatenate(" ") + "\n\n";
				std::cout << "Edge emprirical formula (clean):\n" + ef_clean_edge.getString() + "\n(dirty)\n" + ef_dirty_edge.getString() + "\n\n";


				// END DEBUG
				
				// ---------------------------- //
				// ** write related features ** //
				// ---------------------------- //

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
				typedef std::map<Size, Size> CliqueMap;
				typedef CliqueMap::const_iterator MapCI;
				CliqueMap clique_register; 
				MapCI clique_it;
				
				StringList scores;
				StringList scores_e_inactive_idx, scores_e_active_idx;
				
				for (Size i=0; i<feature_relation.size(); ++i)
				{
          Size f0_idx = feature_relation[i].getElementIndex(0);
          Size f1_idx = feature_relation[i].getElementIndex(1);
          
          Int old_q0 = fm_out[f0_idx].getCharge();
          Int old_q1 = fm_out[f1_idx].getCharge();

					Int new_q0 = feature_relation[i].getCharge(0);
					Int new_q1 = feature_relation[i].getCharge(1);

					scores.push_back(String(feature_relation[i].getEdgeScore()));
					
          if (feature_relation[i].isActive())
          {
            //std::cout << "feature #" << f0_idx << " #" << f1_idx << " ACTIVE with RT: " << fm_out[f1_idx].getRT() << "\n";

            //
						// annotate the affected features
						// ... and check consistency
						//
						Compomer c = feature_relation[i].getCompomer();
						// - left
						EmpiricalFormula ef_l(c.getAdductsAsString(Compomer::LEFT));
						if (fm_out[f0_idx].metaValueExists("dc_charge_adducts"))
    				{
    					if (ef_l.getString() != fm_out[f0_idx].getMetaValue("dc_charge_adducts")) throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("Decharging produced inconsistent adduct annotation! [expected: ") + String(fm_out[f0_idx].getMetaValue("dc_charge_adducts")) + "]", ef_l.getString());
    				}
    				else
    				{
    					fm_out[f0_idx].setMetaValue("dc_charge_adducts", ef_l.getString());
    				}
						fm_out[f0_idx].setMetaValue("dc_charge_adduct_mass", ef_l.getMonoWeight());
						fm_out[f0_idx].setCharge(new_q0);
						
						// - right
						EmpiricalFormula ef_r(c.getAdductsAsString(Compomer::RIGHT));
						if (fm_out[f1_idx].metaValueExists("dc_charge_adducts"))
    				{
    					if (ef_r.getString() != fm_out[f1_idx].getMetaValue("dc_charge_adducts")) throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("Decharging produced inconsistent adduct annotation! [expected: ") + String(fm_out[f1_idx].getMetaValue("dc_charge_adducts")) + "]", ef_r.getString());
    				}
    				else
    				{
    					fm_out[f1_idx].setMetaValue("dc_charge_adducts", ef_r.getString());
    				}
						fm_out[f1_idx].setMetaValue("dc_charge_adduct_mass", ef_r.getMonoWeight());
						fm_out[f1_idx].setCharge(new_q1);
						
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

						scores_e_active_idx.push_back(String(i));
          }
          else // inactive edges
          {
						scores_e_inactive_idx.push_back(String(i));
          
						// DEBUG:
						ConsensusFeature cf(fm_out[f0_idx]);
						cf.insert(0,f0_idx, fm_out[f0_idx]);
						cf.insert(0,f1_idx, fm_out[f1_idx]);
            cf.setMetaValue("Local", String(old_q0)+":"+String(old_q1));
            cf.setMetaValue("CP", String(fm_out[f0_idx].getCharge())+"("+ String(fm_out[f0_idx].getMetaValue("dc_charge_adducts")) +"):"
																 +String(fm_out[f1_idx].getCharge())+"("+ String(fm_out[f1_idx].getMetaValue("dc_charge_adducts")) +") "
																 +String("Score: ") + feature_relation[i].getEdgeScore());
																 
						//cf.computeDechargeConsensus(fm_out);
						// print pairs only
            cons_map_p_neg.push_back(cf);
            cons_map_p_neg.getFileDescriptions()[0].size = fm_out.size();
            cons_map_p_neg.getFileDescriptions()[0].label = "charged features pairs";
          }

        } // !for feature_releation (i.e. edges)


				//  DEBUG 			
				ConsensusXMLFile cf_neg;
				//cf_neg.store("dc_pairs_neg.consensusXML", cons_map_p_neg);
				// DEBUG print scores
				TextFile tf;
				tf.push_back("scr = c(" + scores.concatenate(", ") + ")");
				tf.push_back("s_ia_idx = c(" + scores_e_inactive_idx.concatenate(", ") + ")+1");
				tf.push_back("s_a_idx =   c(" +   scores_e_active_idx.concatenate(", ") + ")+1");
				tf.push_back("s_a_idx_clean = c(" + scores_clean_edge_idx.concatenate(", ") + ")+1");
				tf.push_back("s_a_idx_dirty = c(" +   scores_dirty_edge_idx.concatenate(", ") + ")+1");
				
				tf.push_back("plot( density(scr[s_ia_idx]), xlim=range( scr ), main="", xlab="" )");
				tf.push_back("lines(density(scr[s_a_idx_dirty]), col=2)");
				tf.push_back("lines(density(scr[s_a_idx_clean]), col=3)");
				tf.push_back("legend(x=\"topright\",c(\"dead\", \"active_dirty\", \"active_clean\"), text.col=c(1,2,3))");
				//tf.store("plot_scores.r");



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
				// Warning: from here on cons_map indices have changes --> clique_register[]'s values are not reliable any longer (keys are still goood)

        // include single features without a buddy!
				Size singletons_count = 0;
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
				
        return;
      }

    protected:

			struct CmpInfo_
			{
				String s_comp;
				Size idx_cp;
				UInt side_cp;
				
				// C'tor
				CmpInfo_(): s_comp(), idx_cp(), side_cp() {}
				
				// C'tor
				CmpInfo_(String & s, Size idx, UInt side): s_comp(s), idx_cp(idx), side_cp(side) {}
				
				// Copy C'tor
				CmpInfo_(const CmpInfo_& rhs): s_comp(rhs.s_comp), idx_cp(rhs.idx_cp), side_cp(rhs.side_cp) {}
				
				// Assignment
				CmpInfo_& operator = (const CmpInfo_& rhs)
				{
					if (&rhs == this) return *this;
					s_comp = rhs.s_comp;
					idx_cp = rhs.idx_cp;
					side_cp= rhs.side_cp;
					return *this;					
				}
				
				// Comparator
				bool operator < (const CmpInfo_& other) const
				{
					if (s_comp < other.s_comp) return true; else return false;
				}
				bool operator == (const CmpInfo_& other) const
				{
					if (s_comp == other.s_comp) return true; else return false;
				}				
				
			};

			/// test if "simple" edges have alternative
			/// (more difficult explanation) supported by neighbouring edges
			/// e.g. (.)   -> (H+) might be augmented to
			///      (Na+) -> (H+Na+)
			void inferMoreEdges_(PairsType& edges, Map<Size, std::set<CmpInfo_> >& feature_adducts)
			{
				Adduct default_adduct(1,1,Constants::PROTON_MASS_U,"H1", log(1.0));
			
				Size edges_size = edges.size();
				
				for (Size i=0; i<edges_size; ++i)
				{
					Size idx0 = edges[i].getElementIndex(0);
					Size idx1 = edges[i].getElementIndex(1);
					//if (idx0 == 272 && idx1==271)
					//{ std::cout << "DEBUG reached";}
					
					std::set<CmpInfo_> result;
					
					// look at the intersection of the two adjacent features
					// (we thus require the non-H adduct to be present on both sides of the edge,
					//  if one is deemed enough just change to union)
					std::set_intersection(feature_adducts[idx0].begin(), feature_adducts[idx0].end(),
																feature_adducts[idx1].begin(), feature_adducts[idx1].end(),
																std::inserter(result, result.begin()) );
					std::set<CmpInfo_>::iterator result_it = result.begin();
					
					// add new edge with each adduct of the intersection
					while (result_it!=result.end())
					{
						Compomer::CompomerSide to_add = edges[result_it->idx_cp].getCompomer().removeAdduct(default_adduct).getComponent()[result_it->side_cp];
						// we currently do not punish additional two-side adducts
						for (Compomer::CompomerSide::iterator it=to_add.begin(); it!=to_add.end(); ++it)
						{
							it->second.setLogProb(0);
						}
						ChargePair cp(edges[i]); // make a copy
						Compomer new_cmp = cp.getCompomer().removeAdduct(default_adduct);
						new_cmp.add(to_add, Compomer::LEFT);
						new_cmp.add(to_add, Compomer::RIGHT);
						// refill with default adducts (usually H+):
						if ( ((cp.getCharge(0) - new_cmp.getNegativeCharges()) % default_adduct.getCharge() == 0) &&
								 ((cp.getCharge(1) - new_cmp.getPositiveCharges()) % default_adduct.getCharge() == 0)) // for singly charged default_adducts this should always be true
						{
							int hc_left  = (cp.getCharge(0) - new_cmp.getNegativeCharges()) / default_adduct.getCharge(); // this should always be positive! check!!
							int hc_right = (cp.getCharge(1) - new_cmp.getPositiveCharges()) / default_adduct.getCharge(); // this should always be positive! check!!
							
							// we have not stepped over the charge capacity of the features
							if (hc_left >= 0 && hc_right >= 0)
							{
								// fill up with defaults:
								if (hc_left>0)  new_cmp.add( default_adduct*hc_left, Compomer::LEFT);
								if (hc_right>0) new_cmp.add( default_adduct*hc_right, Compomer::RIGHT);
												
								// charge constraints of feature still fulfilled?
								if ( (new_cmp.getNegativeCharges() == cp.getCharge(0)) &&
										 (new_cmp.getPositiveCharges() == cp.getCharge(1)) )
								{
									cp.setCompomer(new_cmp);
									cp.setEdgeScore(0.99); //TODO how to score this new edge?
									edges.push_back(cp);	  // add edge
									//std::cout << "adding infer CMP with log-score: " << new_cmp.getLogP() << "\n";
								}
								else
								{
									throw Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__, "FeatureDeconvolution::inferMoreEdges_(): discovered internal error!", String(new_cmp.getNegativeCharges()) );
								}
							}
						
						}
						
						++result_it;
					}
				} // edge for
				
				std::cout << "inferMoreEdges_ raised edge count from " << edges_size << " to " << edges.size() << "\n";
			}

			void printEdgesOfConnectedFeatures_(Size idx_1, Size idx_2, const PairsType& feature_relation)
			{
				std::cout << " +++++ printEdgesOfConnectedFeatures_ +++++\n";
				for (Size i=0; i<feature_relation.size(); ++i)
				{
					if (
							((feature_relation[i].getElementIndex(0) == idx_1) && (feature_relation[i].getElementIndex(1)==idx_2))
							 ||
							((feature_relation[i].getElementIndex(0) == idx_2) && (feature_relation[i].getElementIndex(1)==idx_1))
						 )
					{
						std::cout << "\n" << feature_relation[i].getCompomer() << "\n";
					}
				}
				std::cout << " ----- printEdgesOfConnectedFeatures_ -----\n";
				return;
			}

      /// List of adducts used to explain mass differences
      MassExplainer::AdductsType potential_adducts_;

  };
} // namespace OpenMS

#endif // OPENMS_ANALYSIS_DECHARGING_FEATUREDECONVOLUTION_H

