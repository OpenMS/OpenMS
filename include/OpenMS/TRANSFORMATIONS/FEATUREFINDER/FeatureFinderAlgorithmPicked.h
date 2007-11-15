// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm$
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMPICKED_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMPICKED_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithm.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>

namespace OpenMS
{
	/** 
		@brief FeatureFinderAlgorithm for picked peaks.
		
		@ingroup FeatureFinder
	*/
	template<class PeakType, class FeatureType> class FeatureFinderAlgorithmPicked 
		: public FeatureFinderAlgorithm<PeakType, FeatureType>,
			public FeatureFinderDefs
	{
		public:
			///@name Type definitions
			//@{
			typedef typename FeatureFinderAlgorithm<PeakType, FeatureType>::MapType MapType;
			typedef typename MapType::SpectrumType SpectrumType;
			typedef typename PeakType::CoordinateType CoordinateType;
			typedef typename PeakType::IntensityType IntensityType;
			typedef MSExperiment<Peak1D> FilteredMapType;
			//@}
			
			using FeatureFinderAlgorithm<PeakType, FeatureType>::map_;
			using FeatureFinderAlgorithm<PeakType, FeatureType>::param_;
				
		protected:
			///Helper structure for seeds used in FeatureFinderAlgorithmPicked
			struct Seed
			{
				///Spectrum index
				UInt spectrum;
				///Peak index
				UInt peak;
				///Intensity
				Real intensity;
				
				/// Comparison operator
				bool operator<(const Seed& rhs) const
				{
					return intensity<rhs.intensity;
				}
			};
			
			///Helper struct for mass traces used in FeatureFinderAlgorithmPicked
			struct MassTrace
			{
				///Maximum peak pointer
				const FilteredMapType::PeakType* max_peak;
				///RT of maximum peak
				DoubleReal max_rt;
				
				///Contained peaks (pair of RT and pointer to peak)
				std::vector<std::pair<DoubleReal, const FilteredMapType::PeakType*> > peaks;
				
				///determines the charge estimate (voting of all contained peaks)
				UInt getChargeEstimate() const
				{
					UInt charge = 0;
					
					std::map<UInt,UInt> charges;
					for (UInt i=0; i<peaks.size(); ++i)
					{
						charges[peaks[i].second->getMetaValue("charge")]++;
					}
	
					if (charges.size()!=0)
					{
						UInt max_charge = charges.begin()->first;
						UInt max_votes = charges.begin()->second;
						for (std::map<UInt,UInt>::const_iterator it=charges.begin()++; it!=charges.end(); ++it)
						{
							if (it->second>max_votes)
							{
								max_votes = it->second;
								max_charge = it->first;
							}
						}
						charge = max_charge;
					}
					
					return charge;
				}
				
				///determindes the convex hull of the trace
				ConvexHull2D getConvexhull() const
				{
					ConvexHull2D::PointArrayType hull_points(peaks.size());
					for (UInt i=0; i<peaks.size(); ++i)
					{
						hull_points[i][0] = peaks[i].first;
						hull_points[i][1] = peaks[i].second->getMZ();
					}
					return hull_points;
				}
			};
			
		public:			
			/// default constructor 
			FeatureFinderAlgorithmPicked() 
				: FeatureFinderAlgorithm<PeakType,FeatureType>()
			{
				//debugging
				this->defaults_.setValue("debug",0,"If not 0 debug mode is activated. Then several files with intermediate results are written.");
				//intensity
				this->defaults_.setValue("intensity:percentage",15,"Percentage of most intense peaks per spectrum that are probably part of a feature.", false);
				this->defaults_.setValue("intensity:soft_margin",15,"Percentage next lower intense peaks that might be part of a spectrum.", false);
				this->defaults_.setSectionDescription("intensity","Settings for the calculation of a score indicating if a peak is part of a mass trace (between 0 and 1)");
				//mass trace search parameters
				this->defaults_.setValue("mass_trace:mz_tolerance",0.1,"m/z difference tolerance of peaks belonging to the same mass trace.", false);
				this->defaults_.setValue("mass_trace:min_spectra",4,"Number of spectra upstream and downstream of the current spectrum that have to show the same peak mass.", false);
				this->defaults_.setValue("mass_trace:max_missing",1,"Number of spectra where a high mass deviation or missing peak is acceptable", false);
				this->defaults_.setValue("mass_trace:average",1,"Calulation of the mass trace score:\n1: Average the m/z deviations between spectra\n0: multipy the m/z deviations between spectra");
				this->defaults_.setValue("mass_trace:slope_bound",0.1,"The maximum slope of mass trace intensities when extending from the highest peak");
				this->defaults_.setSectionDescription("mass_trace","Settings for the calculation of a score indicating if a peak is part of a mass trace (between 0 and 1).");
				//Isotopic pattern search paramters
				this->defaults_.setValue("isotopic_pattern:charge_low",1,"Lowest charge to search for.", false);
				this->defaults_.setValue("isotopic_pattern:charge_high",4,"Highest charge to search for.", false);
				this->defaults_.setValue("isotopic_pattern:mz_tolerance",0.1,"Tolerated mass deviation from the theoretical isotopic pattern.", false);		
				this->defaults_.setSectionDescription("isotopic_pattern","Settings for the calculation of a score indicating if a peak is part of a isotoipic pattern.");
				//Quality assessment
				this->defaults_.setValue("quality:min_isotope_fit",0.5,"Minimum isotope fit quality.");
				this->defaults_.setValue("quality:mass_trace_intensity_decrease",0.5, "Factor how much intensity the border peaks of a mass trace are allowed to have in comarison to the maximum.");
				this->defaults_.setValue("quality:overall_threshold",0.5, "Overall quality threshold for a feature to be reported.", false);
				this->defaults_.setSectionDescription("quality","Parametes of feature Quality assesment.");
				
				this->defaultsToParam_();
			}
			
			/// Main method for actual FeatureFinder
			virtual void run()
			{
				//-------------------------------------------------------------------------
				//Initialization
				//-------------------------------------------------------------------------
				//TODO Magic constant alert?
				const DoubleReal pattern_bound = 4.5;

				//Initialization for step 2
				UInt min_spectra = param_.getValue("mass_trace:min_spectra");
				bool average = (UInt)(param_.getValue("mass_trace:average"));

				//Initialization for step 3
				UInt charge_low = param_.getValue("isotopic_pattern:charge_low");
				UInt charge_high = param_.getValue("isotopic_pattern:charge_high");
				DoubleReal min_isotope_fit = param_.getValue("quality:min_isotope_fit");
				//Initialization for step 4
				std::vector< Seed > seeds;
				
				//Initialization for step 5
				DoubleReal mass_trace_intensity_decrease = param_.getValue("quality:mass_trace_intensity_decrease");
				
				//Initialization of high score map
				high_score_map_.resize(map_->size());
				
				this->features_->reserve(1000);
				//---------------------------------------------------------------------------
				//Step 1:
				//For each spectrum find the intensity threshold of the x% most intense peaks
				//---------------------------------------------------------------------------
				std::vector< std::pair<IntensityType,IntensityType> > intensity_thresholds;
				intensity_thresholds.reserve(map_->size());
				//Open a new scope to make the local variables used in this step disappear afterwards
				{
					UInt percentage = param_.getValue("intensity:percentage");
					UInt margin = param_.getValue("intensity:soft_margin");
					
					std::vector<CoordinateType> tmp;
					for (typename MapType::const_iterator it=map_->begin(); it!=map_->end(); ++it)
					{
						if (it->size()!=0)
						{
							tmp.resize(it->size());
							for (UInt j=0; j< it->size(); ++j)
							{
								tmp[j] = (*it)[j].getIntensity();
							}
							std::sort(tmp.begin(), tmp.end());
							UInt index = (UInt)std::ceil(tmp.size()*(100.0-percentage)/100.0);
							UInt margin_index = (UInt)std::ceil(tmp.size()*(100.0-percentage-margin)/100.0);
							intensity_thresholds.push_back( std::make_pair(tmp[index],tmp[margin_index]) );
						}
						else
						{
							intensity_thresholds.push_back( std::make_pair(0.0,0.0) );
						}
					}
				}

				//-------------------------------------------------------------------------
				//debugging
				bool debug = ((UInt)(param_.getValue("debug"))!=0);				
				MapType int_map;
				MapType trace_map;
				MapType pattern_map;
				MapType charge_map;
				FeatureMap<> seed_map;

				if (debug)
				{
					int_map.resize(map_->size());
					trace_map.resize(map_->size());
					pattern_map.resize(map_->size());
					charge_map.resize(map_->size());
					
					CoordinateType rt = (*map_)[0].getRT();
					for (UInt s=0; s<min_spectra; ++s)
					{
						int_map[s].setRT(rt);
						trace_map[s].setRT(rt);
						pattern_map[s].setRT(rt);
						charge_map[s].setRT(rt);
						int_map[map_->size()-1-s].setRT(rt);
						trace_map[map_->size()-1-s].setRT(rt);
						pattern_map[map_->size()-1-s].setRT(rt);
						charge_map[map_->size()-1-s].setRT(rt);
					}
				}
				
				//MAIN LOOP FOR SPECTRA
				for (UInt s=0; s<map_->size(); ++s)
				{
					const SpectrumType& spectrum = (*map_)[s];
					high_score_map_[s].setRT(spectrum.getRT());
					
					//do nothing for the first few and last few spectra as
					//the scans required to search for traces are missing
					if (s<min_spectra || s>=map_->size()-min_spectra)
					{
						continue;
					}
					
					//-----------------------------------------------------------
					//Step 3: Calculate IsotopePattern score for all peaks
					// - Test all charges using m/z and intensity fit
					// - Take the charge with the highest score
					//-----------------------------------------------------------
					std::vector<DoubleReal> pattern_scores(spectrum.size(), 0.0);
					std::vector<UInt> pattern_charges(spectrum.size(), 0);
					//Intermediate storage of all charge scores
					std::vector<std::map<UInt,DoubleReal> > all_pattern_scores(spectrum.size());
					for (UInt p=0; p<spectrum.size(); ++p)
					{
						for (UInt c=charge_low; c<=charge_high; ++c)
						{
							DoubleReal mz = spectrum[p].getMZ();
							//get isotope distribution for this mass
							const std::vector<DoubleReal>& isotopes = getIsotopeDistribution_(mz*c);
							//determine highest peak in isopope distribution
							UInt max_isotope = std::max_element(isotopes.begin(), isotopes.end()) - isotopes.begin();
							//Look up expected isotopic peaks
							UInt peak_index = 0; //TODO Optimierung: Startwert schneller finden
							std::vector<DoubleReal> peaks(isotopes.size());
							std::vector<Int> peak_indices(isotopes.size());
							DoubleReal mz_score = 0.0;
							UInt missing_peaks = 0;
							for (UInt i=0; i<isotopes.size(); ++i)
							{
								DoubleReal isotope_pos = mz + ((DoubleReal)i-max_isotope)/c;
								peak_index = nearest_(isotope_pos, spectrum, peak_index);
								DoubleReal single_isotope_score = positionScore_(isotope_pos, spectrum[peak_index].getMZ(), pattern_tolerance_);
								if (single_isotope_score>0.0)
								{
									peaks[i] = spectrum[peak_index].getIntensity();
									mz_score += single_isotope_score;
									peak_indices[i] = peak_index;
								}
								else
								{
									peaks[i] = 0.0;
									peak_indices[i] = -1;
									++missing_peaks;
								}
							}
							//normalize mz_score with the number of possible isotope hits
							mz_score /= isotopes.size();
							//scale mz_score to zero if all but one isotope is missing
							if (missing_peaks==isotopes.size()-1)
							{
								mz_score=0.0;
							}

							//TODO Optimierung: Wenn mz_score gleich 0, braucht man die Korrelation nicht mehr auszurechnen
							//calculate correlation of peak and isotope intensities
							DoubleReal int_score = isotopeScore_(peaks, isotopes);


							//calculate overall score from m/z and intensity score
							DoubleReal pattern_score = mz_score * int_score;
							
							for (std::vector<Int>::const_iterator it=peak_indices.begin(); it!=peak_indices.end();++it)
							{
								if (*it!=-1)
								{
									if (pattern_score>pattern_scores[*it])
									{
										pattern_scores[*it]=pattern_score;
										pattern_charges[*it]=c;
									}
									if (pattern_score>all_pattern_scores[*it][c])
									{
										all_pattern_scores[*it][c] = pattern_score;
									}
								}
							}
						}
					}
					//-----------------------------------------------------------
					//Step 2:
					//Find the mass traces
					//-----------------------------------------------------------					
					if (debug)
					{
						int_map[s].setRT(spectrum.getRT());
						trace_map[s].setRT(spectrum.getRT());
						pattern_map[s].setRT(spectrum.getRT());
						charge_map[s].setRT(spectrum.getRT());
					}
					std::vector<UInt> indices_after(min_spectra+1, 0);
					std::vector<UInt> indices_before(min_spectra+1, 0);
					while( indices_after[0] < spectrum.size() )
					{
						DoubleReal trace_score = 0.0;
						
						//--------------------------------------------------------------
						//Calculate the distances to the nearest peaks in adjacent spectra
						std::vector<DoubleReal> scores;
						
						CoordinateType pos = spectrum[indices_after[0]].getMZ();
						IntensityType inte = spectrum[indices_after[0]].getIntensity();
						
						//std::cout << std::endl << "Peak: " << pos << std::endl;
						bool is_max_peak = true; //checking the maximum intensity peaks -> use them later as feature seeds.
						for (UInt i=1; i<=min_spectra; ++i)
						{
							const SpectrumType& spec = (*map_)[s+i];
							indices_after[i] = nearest_(pos, spec, indices_after[i]);
							DoubleReal position_score = positionScore_( pos, spec[indices_after[i]].getMZ(), trace_tolerance_);
							if (position_score >0 && spec[indices_after[i]].getIntensity()>inte) is_max_peak = false;
							scores.push_back(position_score);
						}
						indices_before[0] = indices_after[0];
						for (UInt i=1; i<=min_spectra; ++i)
						{
							const SpectrumType& spec = (*map_)[s-i];
							indices_before[i] = nearest_(pos, spec, indices_before[i]);
							DoubleReal position_score = positionScore_( pos, spec[indices_before[i]].getMZ(), trace_tolerance_);
							if (position_score>0 && spec[indices_before[i]].getIntensity()>inte) is_max_peak = false;
							scores.push_back(position_score);
						}
						
						//--------------------------------------------------------------
						//Calculate a consensus score out of the scores calculated before
						std::sort(scores.begin(), scores.end());
						
						if (average)
						{
							for(UInt i=max_missing_trace_peaks_; i<scores.size(); ++i)
							{
								trace_score += scores[i];
							}
							trace_score /= 2*min_spectra;
						}
						else
						{
							trace_score = 1.0;
							for(UInt i=max_missing_trace_peaks_; i<scores.size(); ++i)
							{
								trace_score *= scores[i];
							}
						}
						

						//------------------------------------------------------------------
						//Look up precalculated isotope pattern score
						DoubleReal pattern_score = pattern_scores[indices_after[0]];
						
						//------------------------------------------------------------------
						//Calculate intensity score and final score. Determine seeds
						DoubleReal intensity_score = 0.0;
						if (inte>=intensity_thresholds[s].first)
						{
							intensity_score = 1.0;
						}
						else if (inte>=intensity_thresholds[s].second)
						{
							intensity_score = (inte-intensity_thresholds[s].second)/(intensity_thresholds[s].first-intensity_thresholds[s].second);
						}

						//------------------------------------------------------------------
						//Calculate final score. Determine seeds
						DoubleReal final_score = intensity_score*trace_score*pattern_score;
						if (final_score>=0.01)
						{
							//peak index
							UInt p = indices_after[0];
							
							FilteredMapType::PeakType feature_peak;
							feature_peak.setIntensity(spectrum[p].getIntensity());
							feature_peak.setMZ(spectrum[p].getMZ());
							high_score_map_[s].push_back(feature_peak);
							//TODO Optimierung: Charge nicht in MetaInfo speichern
							high_score_map_[s].back().setMetaValue("charge",(Int)pattern_charges[p]);
							//local maximum peaks are considered seeds
							if (is_max_peak)
							{
								Seed seed;
								seed.spectrum = s;
								seed.peak = high_score_map_[s].size()-1;
								seed.intensity = inte;
								seeds.push_back(seed);
								
								Feature tmp;
								tmp.setIntensity(inte);
								tmp.setRT(spectrum.getRT());
								tmp.setMZ(high_score_map_[s].back().getMZ());
								tmp.setCharge(pattern_charges[p]);
								tmp.setOverallQuality(final_score);
								for (std::map<UInt,DoubleReal>::const_iterator it=all_pattern_scores[p].begin(); it!=all_pattern_scores[p].end(); ++it)
								{
									tmp.setMetaValue(String("charge_")+it->first,it->second);
								}
								seed_map.push_back(tmp);
							}
						}
						
						//-------------------------------------------------------------------------
						//debug output
						if (debug)
						{
							PeakType tmp;
							tmp.setPos(pos);
							
							tmp.setIntensity(trace_score);
							trace_map[s].push_back(tmp);

							tmp.setIntensity(pattern_score);
							pattern_map[s].push_back(tmp);
							
							tmp.setIntensity(intensity_score);
							int_map[s].push_back(tmp);
							
							tmp.setIntensity((Int)pattern_charges[indices_after[0]]);
							charge_map[s].push_back(tmp);
						}
						++indices_after[0];
					}//for peaks
				}//for spectra

				//------------------------------------------------------------------
				//Step 4:
				//Extension from a seed
				//------------------------------------------------------------------
				std::sort(seeds.rbegin(),seeds.rend());
				for (UInt i=0; i<seeds.size(); ++i)
				{
					std::cout << std::endl << "Seed " << i << " ~ " << 100.0*i/seeds.size()<< "%" << std::endl;
					//If the intensity is zero this seed is already uses in another feature
					if (seeds[i].intensity == 0.0)
					{
						std::cout << "Abort: Seed intensity is zero." << std::endl;
						continue;
					}
					const FilteredMapType::SpectrumType& spectrum = high_score_map_[seeds[i].spectrum];
					const FilteredMapType::PeakType peak = spectrum[seeds[i].peak];
					std::cout << " - Int: " << peak.getIntensity() << std::endl;
					std::cout << " - RT: " << spectrum.getRT() << std::endl;
					std::cout << " - MZ: " << peak.getMZ() << std::endl;
											
					//----------------------------------------------------------------
					//determine charge of the seed (voting of the whole trace)
					std::vector<MassTrace> traces(1);
					traces.reserve(8);
					extendSeedToMassTrace_(seeds[i], traces.back());
					if (traces[0].peaks.size()<2)
					{
						std::cout << "Abort: Could not extend seed trace." << std::endl;
						continue;
					}
					UInt charge = traces[0].getChargeEstimate();
					std::cout << "Charge (from trace vote): " << charge << std::endl;
					//Abort if no charge could be determined
					if (charge==0)
					{
						std::cout << "Abort: Could not determine a charge." << std::endl;
						continue;
					}

					//----------------------------------------------------------------
					//Find m/z boundaries
					//search end
					UInt end = seeds[i].peak;
					while(end<spectrum.size() && spectrum[end].getMZ()<peak.getMZ()+pattern_bound/(DoubleReal)charge)
					{
						++end;
					}
					--end;
					//search begin
					Int tmp = seeds[i].peak;
					while(tmp>=0 && spectrum[tmp].getMZ()>peak.getMZ()-pattern_bound/(DoubleReal)charge)
					{
						--tmp;
					}
					++tmp;
					UInt begin = tmp;
					
					//----------------------------------------------------------------
					//Find best fitting isotope pattern for this charge (using averagene)
					std::vector<std::pair<Int,DoubleReal> > isotopic_peaks; //pair or index and position (if index is -1)
					DoubleReal isotope_fit_quality = findBestIsotopeFit_(seeds[i], begin, end, isotopic_peaks, charge);
					if (isotope_fit_quality<=min_isotope_fit)
					{
						std::cout << "Abort: Isotope pattern correlation too low. (threshold: " << min_isotope_fit << ")" << std::endl;
						continue;
					}
					
					//extend the convex hull in m/z dimension (starting from the trace peaks)
					//missing traces (index is -1) are simply skipped
					std::cout << "Collecting mass traces" << std::endl;
					for (UInt p=0; p<isotopic_peaks.size(); ++p)
					{
						std::cout << " - Trace " << p << std::endl;
						Seed starting_peak;
						starting_peak.spectrum = seeds[i].spectrum;
						starting_peak.peak = isotopic_peaks[p].first;
						if (isotopic_peaks[p].first==-1)
						{
							DoubleReal missing_isotope_mass = isotopic_peaks[p].second;
							std::cout << "   - missing (" << missing_isotope_mass << ")" << std::endl;
							//try to extend from previous spectrum
							//TODO Optimierung: Besseren Startwert finden
							const FilteredMapType::SpectrumType& spectrum_before = high_score_map_[seeds[i].spectrum-1];
							UInt index_before = nearest_(missing_isotope_mass, spectrum_before, 0);
							if (positionScore_(missing_isotope_mass, spectrum_before[index_before].getMZ(), pattern_tolerance_)!=0.0)
							{
								std::cout << "   - found peak in previous spectrum" << std::endl;
								starting_peak.spectrum = seeds[i].spectrum-1;
								starting_peak.peak = index_before;
							}
							else
							{
								//try to extend from next spectrum
								//TODO Optimierung: Besseren Startwert finden
								const FilteredMapType::SpectrumType& spectrum_after = high_score_map_[seeds[i].spectrum+1];
								UInt index_after = nearest_(missing_isotope_mass, spectrum_after, 0);
								if (positionScore_(missing_isotope_mass, spectrum_after[index_after].getMZ(), pattern_tolerance_)!=0.0)
								{
									std::cout << "   - found peak in next spectrum" << std::endl;
									starting_peak.spectrum = seeds[i].spectrum+1;
									starting_peak.peak = index_after;
								}
								else
								{
									std::cout << "   - could not find peaks in neigboring spectra" << std::endl;
									continue;
								}							
							}							
						}
						else if (starting_peak.peak==seeds[i].peak)
						{
							std::cout << "   - skipped (seed trace)" << std::endl;
							continue;
						}
						starting_peak.intensity = spectrum[starting_peak.peak].getIntensity();
						std::cout << "   - extending from " << high_score_map_[starting_peak.spectrum].getRT() << " / " << high_score_map_[starting_peak.spectrum][starting_peak.peak].getMZ() << std::endl;
						
						//search for nearby maximum of the mass trace
						//as the extension assumes that it starts at the maximum
						for(UInt s=i+1; s<seeds.size(); ++s)
						{
							//the scan is nearby
							if (std::abs((Int)seeds[s].spectrum-(Int)starting_peak.spectrum)<(Int)min_spectra)
							{
								//the peak mass fits
								if (std::fabs(high_score_map_[starting_peak.spectrum][starting_peak.peak].getMZ()-high_score_map_[seeds[s].spectrum][seeds[s].peak].getMZ())<pattern_tolerance_)
								{
									starting_peak.spectrum = seeds[s].spectrum;
									starting_peak.peak = seeds[s].peak;
									std::cout << "   - found nearby seed to extend from at " << high_score_map_[starting_peak.spectrum].getRT() << " / " << high_score_map_[starting_peak.spectrum][starting_peak.peak].getMZ() << std::endl;
								}
							} 
						}
						//extend
						MassTrace tmp;
						extendSeedToMassTrace_(starting_peak, tmp);
						if (tmp.peaks.size()<2)
						{
							std::cout << "   - could not extend trace " << std::endl;
							continue;
						}
						traces.push_back(tmp);
					}
					if (traces.size()<2)
					{
						std::cout << "Abort: Found less than two traces!" << std::endl;
						continue;
					}
					
					//------------------------------------------------------------------
					//Step 5:
					//Quality estimation (isotope fit, convexity, mass trace shape)
					//------------------------------------------------------------------
					std::cout << "Quality estimation" << std::endl;
					Feature f;
					DoubleReal overall_quality = 0.0;
					
					//------------------------------------------------------------------
					//isotope fit of the collected mass trace maxima
					
					overall_quality += isotope_fit_quality;
					if (debug) f.setMetaValue("quality_iso_topefit",isotope_fit_quality);
					std::cout << " - Isotope fit: " << isotope_fit_quality << std::endl;
					
					//------------------------------------------------------------------
					//convexity 
					//calculate convex hull of trace maxima
					ConvexHull2D::PointArrayType border_points;
					for (UInt j=0; j<traces.size(); ++j)
					{
						if (traces[j].peaks.size()!=0)
						{
							border_points.push_back(ConvexHull2D::PointType(traces[j].peaks[0].first, traces[j].peaks[0].second->getMZ()));
							border_points.push_back(ConvexHull2D::PointType(traces[j].peaks.back().first, traces[j].peaks.back().second->getMZ()));
						}
					}
					ConvexHull2D hull = border_points;
					//Score: fraction of removed points is a measure for convexity
					DoubleReal convexity_quality = (hull.getPoints().size()-2.0)/(traces.size()*2.0-2.0);
					overall_quality += convexity_quality;
					std::cout << " - Convexity: " << convexity_quality << " ("<<border_points.size() << " -> " << hull.getPoints().size() << ")" << std::endl;
					if (debug) f.setMetaValue("quality_convexity",convexity_quality);

					//------------------------------------------------------------------
					//trace intensity goes down towards the border
					UInt error_count = 0;
					for (UInt j=0; j<traces.size(); ++j)
					{
						UInt size = traces[j].peaks.size();
						if (size>=2)
						{
							DoubleReal low_int = (traces[j].peaks[0].second->getIntensity() + traces[j].peaks[1].second->getIntensity()) / 2.0;
							DoubleReal high_int = (traces[j].peaks[size-2].second->getIntensity() + traces[j].peaks[size-1].second->getIntensity()) / 2.0;
							if (low_int > mass_trace_intensity_decrease*traces[j].max_peak->getIntensity() || high_int > mass_trace_intensity_decrease*traces[j].max_peak->getIntensity())
							{
								++error_count;
							}
						}
					}
					//Score: fraction of mal-formed rt profiles
					DoubleReal rt_shape_quality = 1.0 - (DoubleReal)error_count / (DoubleReal)traces.size();					
					overall_quality += rt_shape_quality;
					std::cout << " - RT shape: " << rt_shape_quality << std::endl;
					if (debug) f.setMetaValue("quality_rt_shape",rt_shape_quality);
					
					//----------------------------------------------------------------
					//abort if quality too low
					overall_quality /= 3.0;
					if (overall_quality<overall_quality_)
					{
						std::cout << "Abort: feature quality too low: " << overall_quality << " (threshold: " << overall_quality_ << ")" << std::endl;
						continue;
					}

					//------------------------------------------------------------------
					//Step 6:
					//Feature creation
					//------------------------------------------------------------------
					f.setCharge(charge);
					f.setOverallQuality(overall_quality);
					//set feature position and intensity
					for (UInt j=0; j<traces.size(); ++j)
					{
						if (traces[j].max_peak->getIntensity()>f.getIntensity())
						{
							f.setIntensity(traces[j].max_peak->getIntensity());
							f.setRT(traces[j].max_rt);
							f.setMZ(traces[j].max_peak->getMZ());
						}						
					}
					//add convex hulls of mass traces
					for (UInt j=0; j<traces.size(); ++j)
					{
						f.getConvexHulls().push_back(traces[j].getConvexhull());
					}
					//add feature to feature list
					this->features_->push_back(f);
					std::cout << "Feature number: " << this->features_->size() << std::endl;


					//----------------------------------------------------------------
					//Remove all seeds that lie inside the convex hull of the new feature
					DBoundingBox<2> bb = f.getConvexHull().getBoundingBox();
					for (UInt j=i+1; j<seeds.size(); ++j)
					{
						DoubleReal rt = high_score_map_[seeds[j].spectrum].getRT();
						DoubleReal mz = high_score_map_[seeds[j].spectrum][seeds[j].peak].getMZ();
						if (bb.encloses(rt,mz) && f.encloses(rt,mz))
						{
							//set intensity to zero => the peak will be skipped!
							seeds[j].intensity = 0.0;
						}
					}
				}

				//------------------------------------------------------------------
				//store debug info
				if (debug)
				{
					MzDataFile file;
					file.store("trace_scores.mzData",trace_map);
					file.store("charge_estimates.mzData",charge_map);
					file.store("intensity_scores.mzData",int_map);
					file.store("pattern_scores.mzData",pattern_map);
					file.store("selected_peaks.mzData",high_score_map_);
					
					FeatureXMLFile file2;
					file2.store("seeds.featureXML", seed_map);
				}
			}
			
			static FeatureFinderAlgorithm<PeakType,FeatureType>* create()
			{
				return new FeatureFinderAlgorithmPicked();
			}

			static const String getProductName()
			{
				return "picked_peak";
			}
	
		protected:
			/// The map of possible feature peaks
			FilteredMapType high_score_map_;
			
			///@name Members for parameters needed often in methods
			//@{
			DoubleReal pattern_tolerance_; ///< Stores mass_trace:mz_tolerance
			DoubleReal trace_tolerance_; ///< Stores isotopic_pattern:mz_tolerance
			UInt max_missing_trace_peaks_; ///< Stores mass_trace:max_missing
			DoubleReal slope_bound_; ///< Max slope of mass trace intensities
			DoubleReal overall_quality_; ///< Overall quality threshold
			//@}

			//Docu in base class
			virtual void updateMembers_()
			{
				pattern_tolerance_ = param_.getValue("mass_trace:mz_tolerance");
				trace_tolerance_ = param_.getValue("isotopic_pattern:mz_tolerance");
				max_missing_trace_peaks_ = param_.getValue("mass_trace:max_missing");
				slope_bound_ = param_.getValue("mass_trace:slope_bound");
				overall_quality_ = param_.getValue("quality:overall_threshold");
			}


			///Vector of precalculated isotope distributions for several mass winows
			std::vector< std::vector<DoubleReal> > isotope_distributions_;
			
			///Returns the isotope distribution for a certain mass window
			const std::vector<DoubleReal>& getIsotopeDistribution_(DoubleReal mass)
			{
				//calculate index in the vector
				UInt index = (UInt)std::floor(mass/100.0);
				
				//enlarge vector if necessary
				if (index>=isotope_distributions_.size())
				{
					isotope_distributions_.resize(index+1);
				}
				
				//calculate distribution if necessary
				if (isotope_distributions_[index].size()==0)
				{
					std::cout << "Calculating iso dist for mass: " << 50.0 + index * 100.0 << std::endl;
					IsotopeDistribution d;
					d.setMaxIsotope(7);
					d.estimateFromPeptideWeight(50.0 + index * 100.0);
					d.trimLeft(0.015);
					d.trimRight(0.015);
					for (IsotopeDistribution::Iterator it=d.begin(); it!=d.end(); ++it)
					{
						isotope_distributions_[index].push_back(it->second);
						std::cout << " - " << it->second << std::endl;
					}
				}
				//Return distribution
				return isotope_distributions_[index];
			}
						
			/**
				@brief Finds the best fitting position of the isotopic pattern estimate defined by @p center
				
				@param center the maximum peak of the isotope distribution (contains charge as well)
				@param begin Defines the first peak for searching
				@param end Defines the last peak for searching
				@param isotopic_peaks Returns the indices of the isotopic peaks. If a isopopic peak is missing -1 is returned.
				@param charge The charge of the pattern 
			*/
			//TODO: Leave out border isotopes if this increases the fitting quality
			//      Core pattern (e.g. >10%) + optional isotopes (e.g. >1%)
			DoubleReal findBestIsotopeFit_(const Seed& center, UInt begin, UInt end, std::vector<std::pair<Int,DoubleReal> >& isotopic_peaks, UInt charge)
			{
				std::cout << "Testing isotope patterns for charge " << charge << ": " << std::endl;			
				const FilteredMapType::SpectrumType& spectrum = high_score_map_[center.spectrum];
				std::cout << " - Seed: " << center.peak << " (mz:" << spectrum[center.peak].getMZ()<< ")" << std::endl;
				std::cout << " - Begin: " << begin << " (mz:" << spectrum[begin].getMZ()<< ")" << std::endl;
				std::cout << " - End: " << end << " (mz:" << spectrum[end].getMZ()<< ")" << std::endl;

				const std::vector<DoubleReal>& isotopes = getIsotopeDistribution_(spectrum[center.peak].getMZ()*charge);

				//fit isotope distribution to peaks (at least once, even if peaks are missing)
				std::vector<DoubleReal> peaks(isotopes.size());
				DoubleReal max_score = 0.0;
				for (UInt start=begin; start<=(UInt)std::max((Int)end+1-(Int)isotopes.size(),(Int)begin); ++start)
				{
					//find isotope peaks for the current start peak
					UInt peak_index = start;
					peaks[0] = spectrum[start].getIntensity();
					std::cout << " - Fitting at " << start << " (mz:" << spectrum[start].getMZ() << ")" << std::endl;
					std::cout << "   - Isotope 0: " << spectrum[start].getIntensity() << " (" << isotopes[0] << ")" << std::endl;
					std::vector<std::pair<Int,DoubleReal> > isotopic_peaks_tmp;
					isotopic_peaks_tmp.push_back(std::make_pair(start,0.0));
					DoubleReal mz_score = 1.0;
					UInt missing_peaks = 0;
					for (UInt iso=1; iso<isotopes.size(); ++iso)
					{
						DoubleReal pos = spectrum[start].getMZ() + iso*1.0/(DoubleReal)charge;
						peak_index = nearest_(pos, spectrum, peak_index);
						DoubleReal single_isotope_score = positionScore_(pos, spectrum[peak_index].getMZ(), pattern_tolerance_);
						if (single_isotope_score!=0.0)
						{
							std::cout << "   - Isotope " << iso << ": " << spectrum[peak_index].getIntensity() << " (" << isotopes[iso] << ")" << std::endl;

							mz_score += single_isotope_score;
							peaks[iso] = spectrum[peak_index].getIntensity();
							isotopic_peaks_tmp.push_back(std::make_pair(peak_index,0.0));
						}
						else
						{
							std::cout << "   - Isotope " << iso << ": missing" << " (" << isotopes[iso] << ")" << std::endl;
							peaks[iso] = 0.0;
							isotopic_peaks_tmp.push_back(std::make_pair(-1,pos));
							++missing_peaks;
						}
					}
					//check if the seed is contained, other wise abort
					bool seed_contained = false;
					for (UInt iso=0; iso<isotopic_peaks_tmp.size(); ++iso)
					{
						if (isotopic_peaks_tmp[iso].first==(Int)center.peak)
						{
							seed_contained = true;
						}
					}
					if(!seed_contained)
					{
						std::cout << "   - aborting as seed is not contained!" << std::endl;
					}
					else
					{
						//normlize m/z score with the number of isotopes
						mz_score /= isotopes.size();
						//scale m/z score to 0, if all but one peak are missing
						if (missing_peaks==isotopes.size()-1)
						{
							mz_score=0.0;
						}
						std::cout << "   - m/z score: " << mz_score << std::endl;
						
						//TODO Optimierung: Wenn mz_score gleich 0, braucht man die Korrelation nicht mehr auszurechnen
						//calculate correlation
						DoubleReal int_score = isotopeScore_(peaks, isotopes);
						std::cout << "   - int score: " << int_score << std::endl;
						
						DoubleReal score = mz_score * int_score;
						std::cout << "   - final score: " << score << std::endl;
						if (score>max_score)
						{
							max_score = score;
							isotopic_peaks = isotopic_peaks_tmp;
						}
					}
				}
				std::cout << " - best score: " << max_score << std::endl;
				return max_score;
			}

			///Extends a seed to a mass trace
			void extendSeedToMassTrace_(const Seed& start, MassTrace& trace)
			{
				const FilteredMapType::PeakType& seed = high_score_map_[start.spectrum][start.peak];
				//initialize trace with seed data
				trace.peaks.push_back(std::make_pair(high_score_map_[start.spectrum].getRT(), &seed));
				trace.max_peak = &seed;
				trace.max_rt = high_score_map_[start.spectrum].getRT();

				//search downstream
				std::list<DoubleReal> deltas;
				UInt missing_peaks = 0;
				Int scan_index = start.spectrum-1;
				DoubleReal last_intensity = seed.getIntensity();
				deltas.push_back(0.0);
				deltas.push_back(0.0);
				deltas.push_back(0.0);
				UInt added_peaks = 0;
				while(scan_index>=0 && missing_peaks<max_missing_trace_peaks_ && (std::accumulate(deltas.begin(),deltas.end(),0.0)/3.0)<=slope_bound_)
				{
					//TODO Optimierung: find nearest peak faster!
					UInt peak_index = nearest_(seed.getMZ(),high_score_map_[scan_index], 0);
					const FilteredMapType::PeakType& peak = high_score_map_[scan_index][peak_index];
					if ( positionScore_( seed.getMZ(), peak.getMZ(), trace_tolerance_)>0 )
					{
						deltas.push_back((peak.getIntensity() - last_intensity) / last_intensity);
						deltas.pop_front();
						last_intensity = peak.getIntensity();
						missing_peaks = 0;
						trace.peaks.push_back(std::make_pair(high_score_map_[scan_index].getRT(), &peak));
						++added_peaks;
						//update maximum peak of trace
						if (peak.getIntensity()>trace.max_peak->getIntensity())
						{
							trace.max_peak = &peak;
							trace.max_rt = high_score_map_[scan_index].getRT();
						}
					}
					else
					{
						++missing_peaks;
					}
					--scan_index;
				}
				//Go back two peaks if the average went up (we extended too far)
				if ((std::accumulate(deltas.begin(),deltas.end(),0.0)/3.0)>slope_bound_)
				{
					//We need to update the maximum, if if is removed
					bool max_removed = false;
					if (added_peaks>=1)
					{
						if (trace.peaks.back().second==trace.max_peak) max_removed = true;
						trace.peaks.pop_back();
					}
					if (added_peaks>=2)
					{
						if (trace.peaks.back().second==trace.max_peak) max_removed = true;
						trace.peaks.pop_back();
					}
					
					//Update maximum if necessary
					if (max_removed)
					{
						trace.max_peak = &seed;
						trace.max_rt = high_score_map_[start.spectrum].getRT();
						for (UInt i=0; i<trace.peaks.size(); ++i)
						{
							if (trace.peaks[i].second->getIntensity()>trace.max_peak->getIntensity())
							{
								trace.max_peak = trace.peaks[i].second;
								trace.max_rt = trace.peaks[i].first;
							}
						}
					}
				}
				//invert peak array to bring peaks in the correct cronological order
				std::reverse(trace.peaks.begin(), trace.peaks.end());
				
				//search upstream
				missing_peaks = 0;
				scan_index = start.spectrum+1;
				last_intensity = seed.getIntensity();
				deltas.clear();
				deltas.push_back(0.0);
				deltas.push_back(0.0);
				deltas.push_back(0.0);
				//count how many peaks we added (to avoid removing peaks that were not added in this step)
				added_peaks = 0;
				while(scan_index<(Int)high_score_map_.size() && missing_peaks<max_missing_trace_peaks_ && (std::accumulate(deltas.begin(),deltas.end(),0.0)/3.0)<=slope_bound_)
				{
					//TODO Optimierung: find nearest peak faster!
					UInt peak_index = nearest_(seed.getMZ(),high_score_map_[scan_index], 0);
					const FilteredMapType::PeakType& peak = high_score_map_[scan_index][peak_index];
					if ( positionScore_( seed.getMZ(), peak.getMZ(), trace_tolerance_)>0 )
					{
						deltas.push_back((peak.getIntensity() - last_intensity) / last_intensity);
						deltas.pop_front();
						last_intensity = peak.getIntensity();
						missing_peaks = 0;
						trace.peaks.push_back(std::make_pair(high_score_map_[scan_index].getRT(), &peak));
						++added_peaks;
						//update maximum peak of trace
						if (peak.getIntensity()>trace.max_peak->getIntensity())
						{
							trace.max_peak = &peak;
							trace.max_rt = high_score_map_[scan_index].getRT();
						}
					}
					else
					{
						++missing_peaks;
					}
					++scan_index;
				}
				//Go back two peaks if the average went up (we extended too far)
				if ((std::accumulate(deltas.begin(),deltas.end(),0.0)/3.0)>slope_bound_)
				{
					//We need to update the maximum, if if is removed
					bool max_removed = false;
					if (added_peaks>=1)
					{
						if (trace.peaks.back().second==trace.max_peak) max_removed = true;
						trace.peaks.pop_back();
					}
					if (added_peaks>=2)
					{
						if (trace.peaks.back().second==trace.max_peak) max_removed = true;
						trace.peaks.pop_back();
					}
					//Update maximum if necessary
					if (max_removed)
					{
						trace.max_peak = &seed;
						trace.max_rt = high_score_map_[start.spectrum].getRT();
						for (UInt i=0; i<trace.peaks.size(); ++i)
						{
							if (trace.peaks[i].second->getIntensity()>trace.max_peak->getIntensity())
							{
								trace.max_peak = trace.peaks[i].second;
								trace.max_rt = trace.peaks[i].first;
							}
						}
					}
				}
			}
			
			/// Returns the index of the peak nearest to m/z @p pos in spectrum @p spec starting from index @p start
			template <typename SpectrumType>
			UInt nearest_(CoordinateType pos, const SpectrumType& spec, UInt start) const
			{
				UInt index = start;
				CoordinateType dist = std::fabs(pos-spec[index].getMZ());
				++index;
				while (index < spec.size() && std::fabs(pos-spec[index].getMZ()) < dist)
				{
					dist = std::fabs(pos-spec[index].getMZ());
					++index;
				}
				return --index; 
			}
			
			/// Calculates a score between 0 and 1 for the m/z deviation of two peaks. The score becomes 0 at two times the @p allowed_deviation
			DoubleReal positionScore_(CoordinateType pos1, CoordinateType pos2, DoubleReal allowed_deviation) const
			{
				DoubleReal diff = fabs(pos1 - pos2);
				if (diff <= allowed_deviation)
				{
					return 0.1*(allowed_deviation-diff)/allowed_deviation+0.9;
				}
				else if (diff <= 2.0*allowed_deviation)
				{
					return 0.9*(2*allowed_deviation-diff)/allowed_deviation;
				}
				return 0.0;
			}

			/// Calculates a score between 0 and 1 for isotope patterns
			DoubleReal isotopeScore_(const std::vector<DoubleReal>& peaks, const std::vector<DoubleReal>& isotopes) const
			{
				DoubleReal result = Math::BasicStatistics<DoubleReal>::pearsonCorrelationCoefficient(isotopes.begin(), isotopes.end(), peaks.begin(), peaks.end());
				//truncate correlation to [0,1] interval
				return std::max(0.0, result);
			}
		
		private:
			
			/// Not implemented
			FeatureFinderAlgorithmPicked& operator=(const FeatureFinderAlgorithmPicked&);
			/// Not implemented
			FeatureFinderAlgorithmPicked(const FeatureFinderAlgorithmPicked&);

	};

} // namespace OpenMS

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMPICKED_H
