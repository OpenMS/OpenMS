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
#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/MassExplainer.h>
#include <OpenMS/ANALYSIS/DECHARGING/MIPWrapper.h>

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
        defaults_.setValue("potential_adducts", StringList::create("H;1;0.7,K;1;0.1,NH4;1;0.1,Na;1;0.1"), "default adducts used to explain mass differences");
        
        defaultsToParam_();
      }


      void updateMembers_()
      {
        StringList potential_adducts_s = param_.getValue("potential_adducts");
				MassExplainer::AdductsType potential_adducts_;
        for (StringList::const_iterator it=potential_adducts_s.begin(); it != potential_adducts_s.end(); ++it)
        {
          StringList adduct;
          it->split(';', adduct);
          if (adduct.size()!=3)
          {
            String error = "FeatureDeconvolution::potential_adducts (" + (*it) + ") does not have three entries!";
            throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, error);
          }
          String adduct_formula = adduct[0];
          Int charge = adduct[1].toInt();
          float prob = adduct[2].toFloat();
          if (prob > 1.0 || prob<0.0)
					{
            String error = "FeatureDeconvolution::potential_adducts (" + (*it) + ") does not have a proper probablity (" + String(prob) + ") in [0,1]!";
            throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, error);
          }
          EmpiricalFormula ef(adduct_formula);
          Adduct a(charge, 1, ef.getMonoWeight(), ef.getString(), log(prob));
          //std::cout << "FeatureDeconvolution: inserting potential adduct " << ef.getString() << "[q:" << charge << ", pr:" << prob << "]\n";
          potential_adducts_.push_back(a);
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
			/// @param map [in] Input feature-map
			/// @param cons_map [out] Output of grouped features belonging to a charge group
      void compute(FeatureMapType &map, ConsensusMap &cons_map, ConsensusMap &cons_map_p) 
      {
        Int q_min = param_.getValue("charge_min");
				Int q_max = param_.getValue("charge_max");
				Int q_span = param_.getValue("charge_span_max");

				DoubleReal rt_diff_max = param_.getValue("retention_max_diff");
				DoubleReal mz_diff_max = param_.getValue("mass_max_diff");
				
				// sort by RT and then m/z
				map.sortByPosition();

        //TODO: make this accessibile from INI file
				float thresh_logp = log(0.3)*2 + log(0.7)*std::max(q_max-2, 0);
				//create mass difference list



				std::cout << "Generating Masses with threshold: " << thresh_logp << " ...\n";
				//MassExplainer me(adduct_base, q_min, q_max, q_span, thresh_logp);
        MassExplainer me(potential_adducts_, q_min, q_max, q_span, thresh_logp);
				me.compute();
				std::cout << "done\n";
				
				//holds query results for a mass difference
				MassExplainer::CompomerIterator md_s, md_e, md_tmp;
				SignedSize hits = 0;
				
				CoordinateType mz1,mz2, m1;
				
				Size possibleEdges = 0, overallHits = 0;
				
				PairsType feature_relation;
				
        for (UInt i_RT = 0; i_RT < map.size(); ++i_RT)
        { // ** RT-sweepline
					
					mz1 = map[i_RT].getMZ();
					
					for (UInt i_RT_window = i_RT + 1
							 ; (i_RT_window < map.size())
							   &&((map[i_RT_window].getRT() - map[i_RT].getRT()) < rt_diff_max) 
							 ; ++i_RT_window)
					{ // ** RT-window
							
						mz2 = map[i_RT_window].getMZ();
						
						//@improvement TODO only check for charges which approx have the same quotient as mz1/mz2
						for (Int q1=q_min; q1<=q_max; ++q1)
						{ // ** q1

                //DEBUG:
              if (map[i_RT_window].getRT()>1777.08 && map[i_RT_window].getRT()<1777.2 && mz1>216 && mz2>432)
              {
                std::cout << "we are at debug location\n" << map[i_RT_window].getRT() <<"   : " << mz1 << "; " << mz2 << "\n";
              }

								m1 = mz1*q1;
								// additionally: forbid q1 and q2 with distance greater than q_span
								for (Int q2= std::max(q_min, q1-q_span+1)
								     ; (q2<=q_max) && (q2<=q1+q_span-1)
										 ; ++q2)
								{ // ** q2

									++possibleEdges; //internal count, not vital
									
									// find possible adduct combinations
									hits = me.query(q2-q1, mz2*q2-m1, mz_diff_max, thresh_logp, md_s, md_e);
									OPENMS_PRECONDITION(hits>=0,"FeatureDeconvolution querying #hits got negative result!");

									overallHits+=hits;
									// choose most probable hit (TODO think of something clever here)
									// for now, we take the one that has highest p in terms of the compomer structure
									if (hits>0)
									{
											CoordinateType naive_mass_diff = mz2*q2-m1;
											md_tmp = md_s;
											for (; md_s!=md_e; ++md_s)
											{
													if (md_tmp->getLogP() < md_s->getLogP()) md_tmp = md_s;
											}
											ChargePair cp(i_RT, i_RT_window, q1, q2, md_tmp->getID(), naive_mass_diff - md_tmp->getMass(), false);
											if (q1 < abs(md_tmp->getNegativeCharges()) || q2 < abs(md_tmp->getPositiveCharges()))
											{	
												// this should not happen!!! because we assume the compomer to fit the charge assumption
												std::cout << "FeatureDeconvolution.h! compomer does not fit chargePair\n" << (*md_tmp) << "\n" << cp << "\n";
												exit(0);
											}
											//std::cout << "CP # "<< feature_relation.size() << " :" << i_RT << " " << i_RT_window<< " " << q1<< " " << q2 << "\n";
											feature_relation.push_back(cp);
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
				DoubleReal ilp_score = lp_wrapper.compute(me, feature_relation);
				std::cout << "score is: " << ilp_score << std::endl;
				
				
				// print all compomers that were assigned to a feature:
//				std::ofstream myfile;
//				myfile.open ("compomer_association.info");
//				UInt i0,i1;
//        //          feature ID   --> (CompomerID, left?)
//				std::map < UInt, std::set < std::pair < UInt, bool > > > feat_to_comp;
//				for (UInt i=0; i<feature_relation.size(); ++i)
//				{
//          if (!feature_relation[i].isActive()) continue;
//
//          if(i==962)
//          {
//            ChargePair t =  feature_relation[i];
//          }
//
//          i0 = feature_relation[i].getElementIndex(0);
//          i1 = feature_relation[i].getElementIndex(1);
//
//          feat_to_comp[i0].insert( std::make_pair<UInt, bool>(feature_relation[i].getCompomerId() , true) );
//          feat_to_comp[i1].insert( std::make_pair<UInt, bool>(feature_relation[i].getCompomerId() , false) );
//
//          //THINK!!! hack: store assigned feature charge in feature
//          //map[i0].setCharge(feature_relation[i].getCharge(0));
//          //map[i1].setCharge(feature_relation[i].getCharge(1));
//				}
//
//				// print each feature's compomers
//				for (std::map < UInt, std::set < std::pair < UInt, bool > > >::iterator it= feat_to_comp.begin();
//						 it!= feat_to_comp.end(); ++it)
//				{
//
//					myfile << "--------------\nFeature #" << it->first << " q=" << map[it->first].getCharge() << std::endl;
//					for (std::set < std::pair < UInt, bool > >::iterator it_set = it->second.begin();
//							 it_set != it->second.end(); ++it_set)
//					{
//						myfile << "    left? " << it_set->second << "\n";
//						myfile << "    " << me.getCompomerById(it_set->first) << "\n";
//					}
//
//				}
//
//				myfile.close();
				
				//some statistics about compomers used
				std::map < ChargePair, Int > compomer_stats;
				
				UInt agreeing_fcharge = 0;
				UInt f0,f1, aedges=0;
//				myfile.open ("pairs.info");
				//write related features
				for (UInt i=0; i<feature_relation.size(); ++i)
				{
          if (!feature_relation[i].isActive()) continue;

          f0 = feature_relation[i].getElementIndex(0);
          f1 = feature_relation[i].getElementIndex(1);
          ++aedges;

//          myfile << map[f0].getRT() << " " << map[f1].getRT() << " "
//                 << map[f0].getMZ() << " " << map[f1].getMZ() << " "
//                 << map[f0].getCharge() << " " << map[f1].getCharge()  << " "
//                 << me.getCompomerById(feature_relation[i].getCompomerId()).getMass() << "\n";


          // check if the local feature charges agree
          if (map[f0].getCharge() == feature_relation[i].getCharge(0))
          {
            ++agreeing_fcharge;
          }
          else
          {
            std::cout << "conflict in Q:: RT:" << map[f0].getRT() << " MZ:" << map[f0].getMZ() << " Q:" << map[f0].getCharge() << " PredictedQ:" << feature_relation[i].getCharge(0) << "\n";
          }
          if (map[f1].getCharge() == feature_relation[i].getCharge(1))
          {
            ++agreeing_fcharge;
          }
          else
          {
            std::cout << "conflict in Q:: RT:" << map[f1].getRT() << " MZ:" << map[f1].getMZ() << " Q:" << map[f1].getCharge() << " PredictedQ:" << feature_relation[i].getCharge(1) << "\n";
          }
				}
				std::cout << "agreeing charges: " << agreeing_fcharge << "/" << (aedges*2) << std::endl;

//			  myfile.close();
//				myfile.open ("feature.info");
//				for (UInt i=0; i< map.size(); ++i)
//				{
//						myfile << map[i].getRT() << " " << map[i].getMZ() << " " << map[i].getCharge() << "\n";
//				}
//			  myfile.close();
				
				
				
				
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
          if (feature_relation[i].isActive())
          {
            SignedSize target_cf0 = -1, target_cf1 = -1;
            Size f0_idx = feature_relation[i].getElementIndex(0);
            Size f1_idx = feature_relation[i].getElementIndex(1);

            //std::cout << "feature #" << f0_idx << " #" << f1_idx << " ACTIVE with RT: " << map[f1_idx].getRT() << "\n";

            // find the index of the ConsensusFeatures for the current pair
            if (clique_register.count(f0_idx) > 0)
            {
              target_cf0 = clique_register[f0_idx];
            }
            if (clique_register.count(f1_idx) > 0)
            {
              target_cf1 = clique_register[f1_idx];
            }

            Peak2D peak;
            peak.setRT( (map[f0_idx].getRT() +  map[f1_idx].getRT() )/2);
            peak.setMZ( (map[f0_idx].getMZ() +  map[f1_idx].getMZ() )/2);
            peak.setIntensity( (map[f0_idx].getIntensity() +  map[f1_idx].getIntensity() ));
            ConsensusFeature cf(peak);
            cf.insert(0,f0_idx, map[f0_idx]);
            cf.insert(0,f1_idx, map[f1_idx]);
            #if 1
            // print pairs only
                cons_map_p.push_back(cf);
                cons_map_p.getFileDescriptions()[0].size = map.size();
                cons_map_p.getFileDescriptions()[0].label = "charged features pairs";
            #endif

            //seen both features for the first time
            if ((target_cf0 == target_cf1) &&
                (target_cf0 == -1))
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
                cons_map[target_cf1].insert(0,f0_idx, map[f0_idx]);
                clique_register[f0_idx] = target_cf1;
                //std::cout << "add: F" << f0_idx << " to " <<target_cf1 << " dueto F" << f1_idx << "\n";
              }
              else if (target_cf1 == -1)
              {//** add f1 to the already existing cf of f0
                cons_map[target_cf0].insert(0,f1_idx, map[f1_idx]);
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
            //std::cout << "Not active: cp#" << i <<"\n";
          }

        } // !for

				// remove empty ConsensusFeatures from map
				ConsensusMap cons_map_tmp;
				cons_map_tmp.reserve(cons_map.size());
				for(ConsensusMap::ConstIterator it = cons_map.begin(); it!=cons_map.end(); ++it)
				{
					// skip if empty
					if (it->getFeatures().size()==0) continue;
					cons_map_tmp.push_back(*it);					
					// set a centroid
					cons_map_tmp.back().setRT(map[it->getFeatures().begin()->getElementIndex()].getRT());
					cons_map_tmp.back().setMZ(map[it->getFeatures().begin()->getElementIndex()].getMZ());
					cons_map_tmp.back().setIntensity(map[it->getFeatures().begin()->getElementIndex()].getIntensity());

				}
				cons_map_tmp.swap(cons_map);

        // include single features without a buddy!
        for (UInt i=0; i<map.size(); ++i)
				{
          // find the index of the ConsensusFeature for the current feature
          if (clique_register.count(i) > 0) continue;

          Peak2D peak;
          peak.setRT( map[i].getRT() );
          peak.setMZ( map[i].getMZ() );
          peak.setIntensity( map[i].getIntensity() );
          ConsensusFeature cf(peak);
          cf.insert(0, i, map[i]);

          cons_map.push_back(cf);
        }
				
				// fill the header
				//cons_map.setMetaValue("meta",String("value"));
				//cons_map.setIdentifier("some lsid");
				cons_map.getFileDescriptions()[0].filename = "TODO - take from FeatureMAP.getLoadedFilePath () ";
				cons_map.getFileDescriptions()[0].size = map.size();
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

