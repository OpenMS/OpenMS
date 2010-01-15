// -*- Mode: C++; tab-width: 2; -*-
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
//

#include <OpenMS/ANALYSIS/ID/PrecursorIonSelection.h>
#include <OpenMS/ANALYSIS/ID/PrecursorIonSelectionPreprocessing.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>

using namespace std;
//#define PIS_DEBUG
#undef PIS_DEBUG
namespace OpenMS
{
	PrecursorIonSelection::PrecursorIonSelection()
		: DefaultParamHandler("PrecursorIonSelection"),
			max_score_(0.)
	{
	  defaults_.setValue("type","IPS","Strategy for precursor ion selection.");
		defaults_.setValidStrings("type",StringList::create("IPS,SPS,Upshift,Downshift,DEX"));
		defaults_.setValue("precursor_mass_tolerance", 10., "Precursor mass tolerance which is used to query the peptide database for peptides");
		defaults_.setMinFloat("precursor_mass_tolerance",0.);
		defaults_.setValue("precursor_mass_tolerance_unit", "ppm", "Precursor mass tolerance unit.");
		defaults_.setValidStrings("precursor_mass_tolerance_unit",StringList::create("ppm,Da"));
		defaults_.setValue("missed_cleavages",1,"Number of allowed missed cleavages.");
		defaults_.setMinInt("missed_cleavages",0);
		defaults_.setValue("min_pep_ids",2,"Minimal number of identified peptides required for a protein identification.");
		defaults_.setMinInt("min_pep_ids",1);
		defaults_.setValue("max_iteration",100,"Maximal number of iterations.");
		defaults_.setMinInt("max_iteration",1);

		
		defaultsToParam_();
		updateMembers_();
	}

	PrecursorIonSelection::PrecursorIonSelection(const PrecursorIonSelection& source)
		: DefaultParamHandler(source),
			min_pep_ids_(source.min_pep_ids_),
			max_score_(source.max_score_),
			type_(source.type_)
	{
		updateMembers_();
	}

	PrecursorIonSelection::~PrecursorIonSelection()
	{

	}


	const DoubleReal& PrecursorIonSelection::getMaxScore() const 
	{
	  return max_score_; 
	}
	
	void PrecursorIonSelection::setMaxScore(const DoubleReal& max_score)
	{
	  max_score_ = max_score; 
	}


	void PrecursorIonSelection::getNextPrecursors(FeatureMap<>& features,FeatureMap<>& next_features,UInt number)
	{
		sortByTotalScore(features);
		UInt count = 0;
		FeatureMap<>::Iterator iter = features.begin();
		while(iter != features.end() && count < number)
			{
				if((iter->metaValueExists("fragmented") && iter->getMetaValue("fragmented")!= "true")
					 || !iter->metaValueExists("fragmented"))
					{
						if(type_==DEX)
							{
								// if it matching a mass of an identified protein, continue
								if(iter->metaValueExists("shifted") && iter->getMetaValue("shifted")== "down")
									{
										++iter;
										continue;
									}
							}
						iter->setMetaValue("fragmented",(String)"true");
						// store them 
						next_features.push_back(*iter);
						++count;
					}
				++iter;
			}
	}

	void PrecursorIonSelection::rescore_(FeatureMap<>& features,std::vector<PeptideIdentification>& new_pep_ids,
 																			 PrecursorIonSelectionPreprocessing& preprocessed_db)
	{
		// get maximal score in the list
		FeatureMap<>::Iterator iter = features.begin();
		for(;iter != features.end();++iter)
			{
				if((DoubleReal)iter->getMetaValue("msms_score")>max_score_)
					{
						max_score_ = (DoubleReal)iter->getMetaValue("msms_score");
					}
			}
		std::set<String> status_changed;
		// for all new peptide ids
		for(UInt i=0;i<new_pep_ids.size();++i)
			{
#ifdef PIS_DEBUG
				std::cout << new_pep_ids[i].getHits().size() << " hits"<<std::endl;
#endif
				const std::vector<PeptideHit>& hits = new_pep_ids[i].getHits();
				// for all peptide hits
				for(UInt h=0;h<hits.size();++h)
					{
#ifdef PIS_DEBUG
						std::cout << hits[h].getScore() << " >= "<< new_pep_ids[i].getSignificanceThreshold()<<" "
											<< hits[h].getMetaValue("Rank") << std::endl;
#endif
								std::vector< String >::const_iterator acc_it = hits[h].getProteinAccessions().begin();
								for(;acc_it != hits[h].getProteinAccessions().end();++acc_it)
									{
#ifdef PIS_DEBUG
										std::cout << *acc_it << std::endl;
#endif
										if(prot_id_counter_.find(*acc_it) != prot_id_counter_.end())
											{
#ifdef PIS_DEBUG
												std::cout << *acc_it << " hat "<<prot_id_counter_[*acc_it].size() << " hits"<<std::endl;
#endif
												if(prot_id_counter_[*acc_it].size()>=min_pep_ids_) // protein is already identified
													{
#ifdef PIS_DEBUG
														std::cout << prot_id_counter_[*acc_it].size() << " >= "<<min_pep_ids_ << std::endl;
#endif
														// enter peptide id
														prot_id_counter_[*acc_it].insert(hits[h].getSequence().toString());
													}
												else
													{
#ifdef PIS_DEBUG
														std::cout << prot_id_counter_[*acc_it].size() << " < "<<min_pep_ids_ << std::endl;
#endif
														String seq = hits[h].getSequence().toString();
														// if peptide isn't already present
														if(prot_id_counter_[*acc_it].count(seq)==0)
															{
																// enter sequence
																prot_id_counter_[*acc_it].insert(hits[h].getSequence().toString());
																status_changed.insert(*acc_it);
#ifdef PIS_DEBUG
																std::cout <<prot_id_counter_[*acc_it].size()
																					<< " >= "<<min_pep_ids_ << std::endl;
#endif
															}
													}
											}
										else
											{
#ifdef PIS_DEBUG
												std::cout << *acc_it << " ist nich im prot_counter_ drin"<< std::endl;
#endif
												std::set<String> seq_set;
												seq_set.insert(hits[h].getSequence().toString());
												prot_id_counter_.insert(make_pair(*acc_it,seq_set));
												status_changed.insert(*acc_it);
										}
									}
							}
			}

		std::set<String>::iterator it=status_changed.begin();
		// go through all proteins whose status changed
		for( ;it!=status_changed.end();++it)
			{
				if(prot_id_counter_[*it].size() >= min_pep_ids_)
					{
#ifdef PIS_DEBUG
						std::cout << *it << "\t"<<prot_id_counter_[*it].size()<<" hits \n";
#endif
						//						change ordering of corresponding features ->down
						if(type_ == UPSHIFT || type_ == SPS) continue;
						shiftDown_(features,preprocessed_db,*it);
					}
				else
					{
						if(type_ != DOWNSHIFT && type_ != DEX && type_ != SPS)
							{
								// change ordering of corresponding features ->up
#ifdef PIS_DEBUG
						std::cout << *it << "\t"<<prot_id_counter_[*it].size()<<" hits \n";
#endif

								shiftUp_(features,preprocessed_db,*it);
							}
					}
			}

	}


	void PrecursorIonSelection::rescore(FeatureMap<>& features,std::vector<PeptideIdentification>& new_pep_ids,
																			std::vector<ProteinIdentification>& prot_ids,
																			PrecursorIonSelectionPreprocessing& preprocessed_db,bool check_meta_values)
	{
		// check for required MetaValues in FeatureMap
		if(check_meta_values) 		checkForRequiredUserParams_(features);
#ifdef PIS_DEBUG
		std::cout << "checked for user params"<<std::endl;
#endif
		// filter significant peptide ids
		std::vector<PeptideIdentification> filtered_pep_ids = filterPeptideIds_(new_pep_ids);
#ifdef PIS_DEBUG
		std::cout << "filtered peptides ids"<<std::endl;
#endif
		// map ids on features 
		IDMapper mapper;
		Param p = mapper.getParameters();
		p.setValue("rt_delta", 0.2);
		p.setValue("mz_delta", 0.05);
		p.setValue("mz_measure","Da");
		mapper.setParameters(p);

		mapper.annotate(features,filtered_pep_ids,prot_ids);
#ifdef PIS_DEBUG
		std::cout << "mapped ids"<<std::endl;
#endif
		// make the rescoring
		rescore_(features,filtered_pep_ids,preprocessed_db);
		
	}
	
	

	void PrecursorIonSelection::shiftDown_(FeatureMap<>& features,PrecursorIonSelectionPreprocessing& preprocessed_db,
																				 String protein_acc)
	{
		const std::vector<DoubleReal>& masses =preprocessed_db.getMasses(protein_acc);
#ifdef PIS_DEBUG
		std::cout << protein_acc << "  shift down  "<<masses.size()<< " peptides"<<std::endl;
#endif
		// go through all tryptic peptide masses of this protein
		std::vector<DoubleReal>::const_iterator aa_vec_iter = masses.begin();
		for(;aa_vec_iter != masses.end(); ++aa_vec_iter)
			{
				FeatureMap<>::Iterator f_iter = features.begin();
				for(;f_iter != features.end();++f_iter)
					{
						if((DoubleReal) f_iter->getMetaValue("msms_score")>0
							 && f_iter->getMetaValue("fragmented") == "false"
							 && f_iter->getMetaValue("shifted") != "down"
							 && f_iter->getMetaValue("shifted") != "both")
							{
								DoubleReal weight = preprocessed_db.getWeight(*aa_vec_iter);
								if(mz_tolerance_unit_ == "ppm")
									{
																										
										if(fabs(f_iter->getMZ() - *aa_vec_iter) < ((f_iter->getMZ()*mz_tolerance_) / 1e06))
											{
												// matching mass
												DoubleReal score = (DoubleReal)f_iter->getMetaValue("msms_score");
												f_iter->setMetaValue("msms_score",(max( score-(score/2.)*(1-weight),0.)));
#ifdef PIS_DEBUG
										std::cout << f_iter->getRT() << "\t" <<f_iter->getMZ() << "\t"
															<< score << " --> "<<(DoubleReal)f_iter->getMetaValue("msms_score")
															<< std::endl;
#endif
												if(f_iter->getMetaValue("shifted") == "up")
													{
														f_iter->setMetaValue("shifted",(String)"both");
													}
												else
													{
														f_iter->setMetaValue("shifted",(String)"down");
													}
											}
									}
								else if(fabs(f_iter->getMZ() - *aa_vec_iter) < mz_tolerance_)
									{
										// matching mass
										DoubleReal score = (DoubleReal)f_iter->getMetaValue("msms_score");
										f_iter->setMetaValue("msms_score",(max( score-(score/2.)*(1-weight),0.)));
#ifdef PIS_DEBUG
										std::cout << f_iter->getRT() << "\t" <<f_iter->getMZ() << "\t"
															<< score << " --> "<<(DoubleReal)f_iter->getMetaValue("msms_score")
															<< std::endl;
#endif
										if(f_iter->getMetaValue("shifted") == "up")
											{
												f_iter->setMetaValue("shifted",(String)"both");
											}
										else
											{
												f_iter->setMetaValue("shifted",(String)"down");
											}
																										
									}
							}
					}
			}
	}

	void PrecursorIonSelection::shiftUp_(FeatureMap<>& features,PrecursorIonSelectionPreprocessing& preprocessed_db,
																			 String protein_acc)
	{
		const std::vector<DoubleReal>& masses =preprocessed_db.getMasses(protein_acc);
#ifdef PIS_DEBUG
		std::cout << protein_acc << "  shift up  "<<masses.size()<< " peptides"<<std::endl;
#endif
		// go through all tryptic peptide masses of this protein
		std::vector<DoubleReal>::const_iterator aa_vec_iter = masses.begin();
		for(;aa_vec_iter != masses.end(); ++aa_vec_iter)
			{
				FeatureMap<>::Iterator f_iter = features.begin();
				for(;f_iter != features.end();++f_iter)
					{
						if((DoubleReal) f_iter->getMetaValue("msms_score")>0
							 && f_iter->getMetaValue("fragmented") == "false"
							 && f_iter->getMetaValue("shifted") != "up"
							 && f_iter->getMetaValue("shifted") != "both")
							{
								DoubleReal weight = preprocessed_db.getWeight(*aa_vec_iter);
								if(mz_tolerance_unit_ == "ppm")
									{
																										
										if(fabs(f_iter->getMZ() - *aa_vec_iter) < ((f_iter->getMZ()*mz_tolerance_) / 1e06)
											 && f_iter->getMetaValue("shifted") != "up"
											 && f_iter->getMetaValue("shifted") != "both")
											{
												// matching mass
												DoubleReal score = (DoubleReal)f_iter->getMetaValue("msms_score");
												f_iter->setMetaValue("msms_score",score+(max_score_-score)*(1-weight));
#ifdef PIS_DEBUG
										std::cout << f_iter->getRT() << "\t" <<f_iter->getMZ() << "\t"
															<< score << " --> "<<(DoubleReal)f_iter->getMetaValue("msms_score")
															<< std::endl;
#endif
												if(f_iter->getMetaValue("shifted") == "down")
													{
														f_iter->setMetaValue("shifted",(String)"both");
													}
												else
													{
														f_iter->setMetaValue("shifted",(String)"up");
													}
											}
									}
								else if(fabs(f_iter->getMZ() - *aa_vec_iter) < mz_tolerance_
												&& f_iter->getMetaValue("shifted") != "up"
												&& f_iter->getMetaValue("shifted") != "both")
									{
										// matching mass
										DoubleReal score = (DoubleReal)f_iter->getMetaValue("msms_score");
										f_iter->setMetaValue("msms_score",score+(max_score_-score)*(1-weight));
#ifdef PIS_DEBUG
										std::cout << f_iter->getRT() << "\t" <<f_iter->getMZ() << "\t"
															<< score << " --> "<<(DoubleReal)f_iter->getMetaValue("msms_score")
															<< std::endl;
#endif
										if(f_iter->getMetaValue("shifted") == "down")
											{
												f_iter->setMetaValue("shifted",(String)"both");
											}
										else
											{
												f_iter->setMetaValue("shifted",(String)"up");
											}
									}
							}
					}
			}
	}

	void PrecursorIonSelection::checkForRequiredUserParams_(FeatureMap<>& features)
	{
#ifdef PIS_DEBUG
		std::cout << "check for required metadata"<<std::endl;
#endif
		for(UInt i=0;i<features.size();++i)
			{
				if(!features[i].metaValueExists("shifted")) features[i].setMetaValue("shifted",(String)"false");
				if(!features[i].metaValueExists("fragmented")) features[i].setMetaValue("fragmented",(String)"false");
				// default value for msms_score???
				if(!features[i].metaValueExists("msms_score")) features[i].setMetaValue("msms_score",features[i].getIntensity());
				if(!features[i].metaValueExists("init_msms_score")) features[i].setMetaValue("init_msms_score",features[i].getIntensity());
			}
#ifdef PIS_DEBUG
		std::cout << "...finished"<<std::endl;
#endif
		
	}

	std::vector<PeptideIdentification> PrecursorIonSelection::filterPeptideIds_(std::vector<PeptideIdentification>& pep_ids)
	{
		std::vector<PeptideIdentification> filtered_pep_ids;
		
		for(UInt id_c = 0; id_c < pep_ids.size();++id_c)
			{
				std::vector<PeptideHit> tmp_hits;
				if(pep_ids[id_c].getHits()[0].metaValueExists("Rank"))
					{
						for(UInt hit_c=0; hit_c < pep_ids[id_c].getHits().size();++hit_c)
							{
								if(pep_ids[id_c].getHits()[hit_c].getScore() >= pep_ids[id_c].getSignificanceThreshold() && 
									 (Int)pep_ids[id_c].getHits()[hit_c].getMetaValue("Rank") == 1 )
									{
										tmp_hits.push_back(pep_ids[id_c].getHits()[hit_c]);
									}
							}
					}
				else // if meta value rank doesn't exist, take highest scoring peptide hit
					{
						if(pep_ids[id_c].getHits().size() == 1 &&
							 pep_ids[id_c].getHits()[0].getScore() >= pep_ids[id_c].getSignificanceThreshold())
							{
								tmp_hits.push_back(pep_ids[id_c].getHits()[0]);
							}
						else if(pep_ids[id_c].getHits().size()>1)
							{
								UInt max_score_idx=0;
								for(UInt hit_c=1; hit_c < pep_ids[id_c].getHits().size();++hit_c)
									{
										if(pep_ids[id_c].getHits()[hit_c].getScore() >
											 pep_ids[id_c].getHits()[max_score_idx].getScore())
											{
												max_score_idx = hit_c;
											}
									}
								// check if highest scoring peptide hit is significant
								if(pep_ids[id_c].getHits()[max_score_idx].getScore() >= pep_ids[id_c].getSignificanceThreshold())
									{
										tmp_hits.push_back(pep_ids[id_c].getHits()[max_score_idx]);
									}
							}
					}
				
				if(!tmp_hits.empty()) // if there were significant hits save them
					{
						PeptideIdentification tmp_id = pep_ids[id_c];
						tmp_id.setHits(tmp_hits);
						filtered_pep_ids.push_back(tmp_id);
					}  
			}
	
	return filtered_pep_ids;
}


	void PrecursorIonSelection::reset()
	{
		prot_id_counter_.clear();
	}
	
	void PrecursorIonSelection::simulateRun(FeatureMap<>& features,std::vector<PeptideIdentification>& pep_ids,
																					std::vector<ProteinIdentification>& prot_ids,
																					PrecursorIonSelectionPreprocessing& preprocessed_db,
																					UInt step_size, String path)
	{

		sortByTotalScore(features);
		std::ofstream outf(path.c_str());
#ifdef PIS_DEBUG
		std::cout << "type "<<type_<<std::endl;
		std::cout << "upshift "<<UPSHIFT << "\n"
							<< "downshift "<<DOWNSHIFT << "\n"
							<< "ips "<<IPS << "\n"
							<< "sps "<<SPS << "\n"
							<< "dex "<<DEX << "\n";
#endif
		if(features.size() == 0)  return;
		// check if feature map has required user_params-> else add them
		checkForRequiredUserParams_(features);

		std::vector<PeptideIdentification> filtered_pep_ids = filterPeptideIds_(pep_ids);
		
		// annotate map with ids
		// TODO: wirklich mit deltas? oder lieber ueber convex hulls? Anm v. Chris: IDMapper benutzt CH's + Deltas wenn CH vorhanden sind
		IDMapper mapper;
		Param p = mapper.getParameters();
		p.setValue("rt_delta", 0.2);
		p.setValue("mz_delta", 0.05);
		p.setValue("mz_measure","Da");
		mapper.setParameters(p);
		mapper.annotate(features,filtered_pep_ids,prot_ids,true);
		
		// get first precursors
		FeatureMap<> new_features;
		getNextPrecursors(features,new_features,step_size);	

		Size precursors = 0;
		UInt iteration = 0;
		UInt pep_id_number = 0;
		std::vector<PeptideIdentification> curr_pep_ids,all_pep_ids;
		std::vector<ProteinIdentification> curr_prot_ids,all_prot_ids;
#ifdef PIS_DEBUG
		std::cout << max_iteration_ << std::endl;
#endif
		// while there are precursors left and the maximal number of iterations isn't arrived
		while((new_features.size()  > 0 && iteration < max_iteration_) )
			{
					
				++iteration;
#ifdef PIS_DEBUG
				std::cout << "================================ iteration "<<iteration<<std::endl;
#endif
				curr_pep_ids.clear();
				curr_prot_ids.clear();
					
				// go through the new compounds
				for(UInt c=0;c<new_features.size();++c)
					{
#ifdef PIS_DEBUG
						// print info
						std::cout << "rt "<< new_features[c].getRT()
											<< " mz: "<<  new_features[c].getMZ()
											<< " shifted: "<<new_features[c].getMetaValue("shifted") << " "
											<< " init_msms_score " <<  new_features[c].getMetaValue("init_msms_score") << " -> "
											<<  new_features[c].getMetaValue("msms_score") ;
#endif
		    

		    
						// get their peptide ids
						std::vector<PeptideIdentification> & pep_ids = new_features[c].getPeptideIdentifications();
#ifdef PIS_DEBUG
						if(pep_ids.size() > 0)
							{
								std::cout << " ids "	<< std::endl;
								std::cout << pep_ids[0].getHits()[0].getSequence().toString()<<std::endl;
								std::cout << pep_ids[0].getHits()[0].getScore() << " "
													<< pep_ids[0].getSignificanceThreshold() << " "
													<< pep_ids[0].getHits()[0].getMetaValue("Rank")<<std::endl;
							}
						else std::cout << std::endl;
						if((DoubleReal)new_features[c].getMetaValue("init_msms_score")==0 && pep_ids.size()>0)		    
							{
								std::cout << "Attention: score was 0, but spectrum led to an id!!! "
													<< std::endl;
								std::cout << "Score : "
													<<pep_ids[0].getHits()[0].getScore()<<std::endl;
							}
#endif
						for(UInt pep_id=0;pep_id<pep_ids.size();++pep_id)
							{
								// save peptide id
								all_pep_ids.push_back(pep_ids[pep_id]);
								curr_pep_ids.push_back(pep_ids[pep_id]);
								// go through peptide hits
								const std::vector<PeptideHit> & pep_hits = pep_ids[pep_id].getHits();
								for(UInt pep_hit=0;pep_hit<pep_hits.size();++pep_hit)
									{
										// get their accessions
										const std::vector<String>& accs = pep_hits[pep_hit].getProteinAccessions();
										//std::cout << accs.size() << std::endl;
										const std::vector<ProteinIdentification> & prot_ids = features.getProteinIdentifications();
										// get ProteinIds for accession and save them
										for(UInt prot_id=0;prot_id < prot_ids.size();++prot_id)
											{
												const std::vector<ProteinHit> & prot_hits = prot_ids[prot_id].getHits();
												for(UInt prot_hit=0;prot_hit < prot_hits.size();++prot_hit)
													{
														if(find(accs.begin(),accs.end(),prot_hits[prot_hit].getAccession()) != accs.end())
															{
																//std::cout << "found "<<prot_hits[prot_hit].getAccession() << std::endl;
																// check if protein is already in all_prot_ids
																bool exists = false;
																for(UInt s_prot_id = 0; s_prot_id < all_prot_ids.size();++s_prot_id) // should be at most one
																	{
																		for(UInt s_prot_hit = 0; s_prot_hit < all_prot_ids[s_prot_id].getHits().size();++s_prot_hit)
																			{
																				if(all_prot_ids[s_prot_id].getHits()[s_prot_hit].getAccession()
																					 == prot_hits[prot_hit].getAccession())
																					{
																						exists = true;
																						break;
																					}
																			}
																	} // for(UInt s_prot_id = 0; s_prot_id <...
																// add only if this protein doesn't exist as a hit yet 
																if(!exists)
																	{ 
																		if(all_prot_ids.size()>0)
																			{
																				all_prot_ids[0].insertHit(prot_hits[prot_hit]);
																				//std::cout << "enter prot id "<< prot_hits[prot_hit].getAccession()<<std::endl;
																			}
																		else
																			{
																				ProteinIdentification prot_identification;
																				all_prot_ids.push_back(prot_identification);
																				all_prot_ids[0].insertHit(prot_hits[prot_hit]);
																								//	std::cout << "enter prot id "<< prot_hits[prot_hit].getAccession()<<std::endl;
																			}
																		if(curr_prot_ids.size()>0)
																			{
																				curr_prot_ids[0].insertHit(prot_hits[prot_hit]);
																			}
																		else
																			{
																				ProteinIdentification prot_identification;
																				curr_prot_ids.push_back(prot_identification);
																				curr_prot_ids[0].insertHit(prot_hits[prot_hit]);
																			}
																		
																	}//if(!exists)
																
															} // if(find(accs.begin()...
													}//for(UInt prot_hit=0;prot_hit <...
											} // for(UInt prot_id=0;prot_id < ...
									}//for(UInt pep_hit=0;pep_hit<...
							}//for(UInt pep_id=0;pep_id<...
					}//for(UInt c=0;c<new_cl.size();++c)
				
				precursors += new_features.size();
#ifdef PIS_DEBUG
				std::cout << "new features "<<new_features.size() << std::endl;
				std::cout << "curr_pep_ids "<<curr_pep_ids.size() << std::endl;
#endif
				if(curr_pep_ids.size() == 0)
					{
						// necessary for making the figures later
						UInt num_prots = filterProtIds_(all_prot_ids);
						outf << iteration << "\t\t"  
								 << num_prots << "\t\t" 
								 << precursors << "\t\t"
								 << pep_id_number
								 << std::endl;

						new_features.clear(true);
						getNextPrecursors(features,new_features,step_size);	
						continue;
					}


				rescore_(features,curr_pep_ids,preprocessed_db);
		
				
				UInt num_prots = filterProtIds_(all_prot_ids);

#ifdef PIS_DEBUG
				for(UInt p=0;p<curr_pep_ids.size();++p)
					{
						for(UInt h=0;h<curr_pep_ids[p].getHits().size();++h)
							{
								for(UInt a=0; a<curr_pep_ids[p].getHits()[h].getProteinAccessions().size();++a)
									{
										std::cout << curr_pep_ids[p].getHits()[h].getProteinAccessions()[a] << "\t"
															<< prot_id_counter_[curr_pep_ids[p].getHits()[h].getProteinAccessions()[a]].size()
															<< std::endl;
										std::set<String>::iterator str_it= prot_id_counter_[curr_pep_ids[p].getHits()[h].getProteinAccessions()[a]].begin();
										for(;str_it!=prot_id_counter_[curr_pep_ids[p].getHits()[h].getProteinAccessions()[a]].end();++str_it)
											{
												std::cout << *str_it<<std::endl;
												
												
											}
									}
							}
						
					}
#endif
				
				// write results of current iteration				
				outf << iteration << "\t\t"
						 << num_prots << "\t\t" 
						 << precursors << "\t\t"
						 << pep_id_number
						 << std::endl;
				new_features.clear(true);
				getNextPrecursors(features,new_features,step_size);	
#ifdef PIS_DEBUG
				std::cout << new_features.size() << " compounds for msms"<< std::endl;
#endif
			}//while(new_features.size() > 0 && iteration < max_iteration)
		
#ifdef PIS_DEBUG		
		for(UInt p=0;p<all_pep_ids.size();++p)
			{
				for(UInt h=0;h<all_pep_ids[p].getHits().size();++h)
					{
						if(all_pep_ids[p].getHits()[h].getProteinAccessions().size() > 1)
							{
								std::cout << all_pep_ids[p].getHits()[h].getSequence() << " in "
													<<all_pep_ids[p].getHits()[h].getProteinAccessions().size() 
													<< " proteins:\n";
								for(UInt a=0; a<all_pep_ids[p].getHits()[h].getProteinAccessions().size();++a)
									{
										std::cout << all_pep_ids[p].getHits()[h].getProteinAccessions()[a] << "\t"
															<< prot_id_counter_[all_pep_ids[p].getHits()[h].getProteinAccessions()[a]].size()
															<< std::endl;
										std::set<String>::iterator str_it= prot_id_counter_[all_pep_ids[p].getHits()[h].getProteinAccessions()[a]].begin();
										for(;str_it!=prot_id_counter_[all_pep_ids[p].getHits()[h].getProteinAccessions()[a]].end();++str_it)
											{
												std::cout << *str_it<<std::endl;
												
												
											}
									}
							}
					}
				
			}
#endif
		
		
	}



	UInt PrecursorIonSelection::filterProtIds_(std::vector<ProteinIdentification>& prot_ids)
	{
		std::vector<UInt> not_count;
		UInt prot_count = 0;
		UInt prot_id_count = 0;
		for(UInt k = 0; k <  prot_ids.size();++k)
			{

				std::vector<ProteinHit> hits = prot_ids[k].getHits();
				for(UInt i=0;i<hits.size();++i)
					{
						if(prot_id_counter_[hits[i].getAccession()].size()<min_pep_ids_) continue; // if this protein has less than 3 peptides it isn't identified
						// for debugging purposes
						++prot_id_count;
						// is this already combined?
						if(find(not_count.begin(),not_count.end(),i) == not_count.end())
							{
								++prot_count;
								// go through all prot_hits and search for same score
								for(UInt j=i+1;j< hits.size();++j)
									{
										if(prot_id_counter_[hits[j].getAccession()].size()<min_pep_ids_) continue;
										// get referencing hits
										std::set<String> pep_hits_i = prot_id_counter_[hits[i].getAccession()];
										std::set<String> pep_hits_j = prot_id_counter_[hits[j].getAccession()];
												
										// if one includes the other -> do not count as separate hits
										if(includes(pep_hits_i.begin(),pep_hits_i.end(),
																pep_hits_j.begin(),pep_hits_j.end()))
											{
												not_count.push_back(j);
											}
										else if(includes(pep_hits_j.begin(),pep_hits_j.end(),pep_hits_i.begin(),pep_hits_i.end()))
											{
												not_count.push_back(i);
											}

									}
							 
							}

					}
			}

		sort(not_count.begin(),not_count.end());
		std::vector<UInt>::iterator new_end = unique(not_count.begin(),not_count.end());
		prot_count =  prot_id_count - distance(not_count.begin(),new_end);
#ifdef PIS_DEBUG
		std::cout <<"\tnew num id. prots: "<< prot_count<<std::endl;
#endif
		
		return prot_count;
	}

	


	
	
	void PrecursorIonSelection::updateMembers_()
	{
#ifdef PIS_DEBUG
		std::cout << "update members"<<std::endl;
#endif
		
		if(param_.getValue("type") == "IPS") type_ = IPS;
		else if(param_.getValue("type") == "Upshift") type_ = UPSHIFT;
		else if(param_.getValue("type") == "Downshift") type_ = DOWNSHIFT;
		else if(param_.getValue("type") == "SPS") type_ = SPS;
		else type_ = DEX;
		min_pep_ids_ = (UInt)param_.getValue("min_pep_ids");
		mz_tolerance_unit_ = (String)param_.getValue("precursor_mass_tolerance_unit");
		mz_tolerance_ = (DoubleReal)param_.getValue("precursor_mass_tolerance");
		max_iteration_ = (UInt) param_.getValue("max_iteration");
	}

	
} //namespace

