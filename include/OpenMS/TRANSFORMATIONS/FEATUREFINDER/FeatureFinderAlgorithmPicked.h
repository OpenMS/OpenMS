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

#include <numeric>

namespace OpenMS
{
	/** 
		@brief FeatureFinderAlgorithm for picked peaks.
		
		@improvement Estimate intensity cutoff from histogram of each bin (Marc)
		
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
			//@}
			
			using FeatureFinderAlgorithm<PeakType, FeatureType>::map_;
			using FeatureFinderAlgorithm<PeakType, FeatureType>::param_;
				
		protected:
			/// Helper structure that stores precalculated information for each peak used in FeatureFinderAlgorithmPicked
			struct PeakInfo
			{
				Real trace_score; ///< precalculated mass trace score
				Real intensity_score; ///< precalculated intensity score
				Real pattern_score; ///< precalculated isotope pattern score (for the current charge)
				Real overall_score; ///< overall score of the peak (for the current charge)
				bool local_max; ///< local maximum flag (possible seed)

				PeakInfo()
					: trace_score(0.0),
						intensity_score(0.0),
						pattern_score(0.0),
						overall_score(0.0),
						local_max(false)
				{
				}
			};

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
				const PeakType* max_peak;
				///RT of maximum peak
				DoubleReal max_rt;
				
				///Contained peaks (pair of RT and pointer to peak)
				std::vector<std::pair<DoubleReal, const PeakType*> > peaks;
				
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
				
				///Sets the maximum to the highest contained peak of the trace
				void updateMaximum()
				{
					if (peaks.size()==0) return;

					max_rt = peaks.begin()->first;					
					max_peak = peaks.begin()->second;
					
					for (UInt i=1; i<peaks.size(); ++i)
					{
						if (peaks[i].second->getIntensity()>max_peak->getIntensity())
						{
							max_rt = peaks[i].first;					
							max_peak = peaks[i].second;
						}
					}
				}
			};

			///Helper structure for a theoretical isotope pattern used in FeatureFinderAlgorithmPicked
			struct TheoreticalIsotopePattern
			{
				///Vector of intensity contributions 
				std::vector<DoubleReal> intensity;
				///Number of optional peaks at the beginning of the pattern
				UInt optional_begin;
				///Number of optional peaks at the end of the pattern
				UInt optional_end;

				/// Returns the size
				UInt size() const
				{
					return intensity.size();
				}
			};

			///Helper structure for a found isotope pattern used in FeatureFinderAlgorithmPicked
			struct IsotopePattern
			{
				///Peak index (-1 if peak was not found, -2 if it was removed to improve the isotope fit)
				std::vector<Int> peak;
				///Spectrum index (undefined if peak index is -1 or -2)
				std::vector<UInt> spectrum;
				///Peak intensity (0 if peak index is -1 or -2)
				std::vector<DoubleReal> intensity;
				///m/z score of peak (0 if peak index is -1 or -2)
				std::vector<DoubleReal> mz_score;
				///Theoretical m/z value of the isotope peak
				std::vector<DoubleReal> theoretical_mz;
				
				/// Constructor that resizes the internal vectors
				IsotopePattern(UInt size)
					: peak(size),
						spectrum(size),
						intensity(size),
						mz_score(size),
						theoretical_mz(size)
				{
				}
			};
			

		public:			
			/// default constructor 
			FeatureFinderAlgorithmPicked() 
				: FeatureFinderAlgorithm<PeakType,FeatureType>(),
					log_("featurefinder.log")
			{
				//debugging
				this->defaults_.setValue("debug",0,"If not 0 debug mode is activated. Then several files with intermediate results are written.");
				//intensity
				this->defaults_.setValue("intensity:bins",10,"Number of bins per dimension (RT and m/z).");
				this->defaults_.setValue("intensity:percentage",35.0,"Percentage of most intense peaks per bin that might be part of a feature.");
				this->defaults_.setSectionDescription("intensity","Settings for the calculation of a score indicating if a peak's intensity is significant (between 0 and 1)");
				//mass trace search parameters
				this->defaults_.setValue("mass_trace:mz_tolerance",0.06,"m/z difference tolerance of peaks belonging to the same mass trace.");
				this->defaults_.setValue("mass_trace:min_spectra",14,"Number of spectra the have to show the same peak mass for a mass trace.");
				this->defaults_.setValue("mass_trace:max_missing",4,"Number of spectra where a high mass deviation or missing peak is acceptable.");
				this->defaults_.setValue("mass_trace:slope_bound",0.1,"The maximum slope of mass trace intensities when extending from the highest peak", true);
				this->defaults_.setSectionDescription("mass_trace","Settings for the calculation of a score indicating if a peak is part of a mass trace (between 0 and 1).");
				//Isotopic pattern search paramters
				this->defaults_.setValue("isotopic_pattern:charge_low",1,"Lowest charge to search for.");
				this->defaults_.setValue("isotopic_pattern:charge_high",4,"Highest charge to search for.");
				this->defaults_.setValue("isotopic_pattern:mz_tolerance",0.06,"Tolerated mass deviation from the theoretical isotopic pattern.");		
				this->defaults_.setValue("isotopic_pattern:intensity_percentage",10.0,"Isotopic peaks that contribute more than this percentage to the overall isotope pattern intensity must be present.", true);
				this->defaults_.setValue("isotopic_pattern:intensity_percentage_optional",0.1,"Isotopic peaks that contribute more than this percentage to the overall isotope pattern intensity can be missing.", true);
				this->defaults_.setValue("isotopic_pattern:optional_fit_improvement",3.0,"Minimal percental improvement of isotope fit to allow leaving out an optional peak.", true);
				this->defaults_.setValue("isotopic_pattern:mass_window_width",100.0,"Window width in Dalton for precalcuation of estimated isotope distribtions.", true);
				this->defaults_.setSectionDescription("isotopic_pattern","Settings for the calculation of a score indicating if a peak is part of a isotoipic pattern (between 0 and 1).");
				//Feature settings
				this->defaults_.setValue("feature:intensity_as_max","true","Determines if feature intensity is reported as the maximum of the feature peaks (true) or the sum of all intensities (false).");
				this->defaults_.setValue("feature:minimum_quality",0.75, "Overall quality threshold for a feature to be reported.");
				this->defaults_.setValue("feature:min_isotope_fit",0.65,"Minimum isotope fit quality.", true);
				this->defaults_.setValue("feature:mass_trace_max_border_intensity",0.7, "Factor how much intensity the border peaks of a mass trace are allowed to have in comarison to the maximum.", true);
				this->defaults_.setSectionDescription("feature","Settings for the features (intensity, quality assessment, ...)");
				
				this->defaultsToParam_();
			}
			
			/// Main method for actual FeatureFinder
			virtual void run()
			{
				//-------------------------------------------------------------------------
				//initialize debugging
				bool debug = ((UInt)(param_.getValue("debug"))!=0);				
				MapType int_map;
				MapType trace_map;
				MapType pattern_map;
				MapType selected_map;
				FeatureMap<> seed_map;

				if (debug)
				{
					//resize maps
					int_map.resize(map_->size());
					trace_map.resize(map_->size());
					pattern_map.resize(map_->size());
					selected_map.resize(map_->size());
					seed_map.reserve(1000);

					//set RTs of scans in debug info maps
					for (UInt s=0; s<map_->size(); ++s)
					{
						DoubleReal rt = map_->at(s).getRT();
						if (debug)
						{
							int_map[s].setRT(rt);
							trace_map[s].setRT(rt);
							pattern_map[s].setRT(rt);
							selected_map[s].setRT(rt);
						}
					}
				}
					
				//-------------------------------------------------------------------------
				//General initialization
				this->features_->reserve(1000);
				//quality estimation
				DoubleReal mass_trace_max_border_intensity = param_.getValue("feature:mass_trace_max_border_intensity");
				DoubleReal min_feature_quality = param_.getValue("feature:minimum_quality");
				//feature intensity
				bool max_intensity = String(param_.getValue("feature:intensity_as_max"))=="true";
				//initialize info_
				info_.resize(map_->size());
				for (UInt s=0; s<map_->size(); ++s)
				{
					info_[s].resize(map_->at(s).size());
				}
				
				//---------------------------------------------------------------------------
				//Step 1:
				//Precalculate intensity scores for peaks
				//---------------------------------------------------------------------------
				log_ << "Precalculating intensity thresholds ..." << std::endl;
				//new scope to make local variables disapear
				{
					DoubleReal percentage = param_.getValue("intensity:percentage");
					DoubleReal rt_start = map_->getMinRT();
					DoubleReal mz_start = map_->getMinMZ();
					intensity_rt_step_ = (map_->getMaxRT() - rt_start ) / (DoubleReal)intensity_bins_;
				 	intensity_mz_step_ = (map_->getMaxMZ() - mz_start ) / (DoubleReal)intensity_bins_;
					intensity_thresholds_.resize(intensity_bins_);
					for (UInt rt=0; rt<intensity_bins_; ++rt)
					{
						intensity_thresholds_[rt].resize(intensity_bins_);
						DoubleReal min_rt = rt_start + rt * intensity_rt_step_;
						DoubleReal max_rt = rt_start + ( rt + 1 ) * intensity_rt_step_;
						std::vector<DoubleReal> tmp;
						for (UInt mz=0; mz<intensity_bins_; ++mz)
						{
							DoubleReal min_mz = mz_start + mz * intensity_mz_step_;
							DoubleReal max_mz = mz_start + ( mz + 1 ) * intensity_mz_step_;
							//std::cout << "rt range: " << min_rt << " - " << max_rt << std::endl;
							//std::cout << "mz range: " << min_mz << " - " << max_mz << std::endl;
							tmp.clear();
							for (typename MapType::ConstAreaIterator it = map_->areaBeginConst(min_rt,max_rt,min_mz,max_mz); it!=map_->areaEndConst(); ++it)
							{
								tmp.push_back(it->getIntensity());
							}
							std::sort(tmp.begin(), tmp.end());
							//std::cout << "number of peaks: " << tmp.size() << std::endl;
							if (tmp.size()==0)
							{
								intensity_thresholds_[rt][mz] = std::make_pair(0.0,0.0);
							}
							else
							{
								//std::cout << "rt:" << rt << " mz:" << mz << " index="<< index << std::endl;
								//std::cout << (min_rt+max_rt)/2.0 << " " << (min_mz+max_mz)/2.0 << " " << tmp.size() << std::endl;
								UInt index = (UInt)std::ceil(tmp.size()*(100.0-percentage)/100.0);
								intensity_thresholds_[rt][mz] = std::make_pair(tmp[index], tmp.back());
							}
						}
					}
					//store intensity score in PeakInfo
					for (UInt s=0; s<map_->size(); ++s)
					{
						for (UInt p=0; p<map_->at(s).size(); ++p)
						{
							info_[s][p].intensity_score = intensityScore_(map_->at(s)[p].getIntensity(),s,p);
							if (debug && info_[s][p].intensity_score>0)
							{
								PeakType tmp;
								tmp.setPos(map_->at(s)[p].getMZ());
								tmp.setIntensity(info_[s][p].intensity_score);
								int_map[s].push_back(tmp);
							}
						}
					}
					//Store intensity score map
					if (debug)
					{
						MzDataFile().store("intensity_scores.mzData",int_map);
					}
				}

				//---------------------------------------------------------------------------
				//Step 2:
				//Prealculate mass trace scores and local trace maximum for each peak
				//---------------------------------------------------------------------------
				//new scope to make local variables disapear
				{
					this->ff_->startProgress(0, map_->size(), "Precalculating mass trace scores");
					for (UInt s=0; s<map_->size(); ++s)
					{
						this->ff_->setProgress(s);
						//do nothing for the first few and last few spectra as the scans required to search for traces are missing
						if (s<min_spectra_ || s>=map_->size()-min_spectra_)
						{
							continue;
						}
						
						const SpectrumType& spectrum = map_->at(s);
						std::vector<UInt> indices_after(min_spectra_+1, 0);
						std::vector<UInt> indices_before(min_spectra_+1, 0);
						//iterate over all peaks of the scan
						while( indices_after[0] < spectrum.size() )
						{
							std::vector<DoubleReal> scores;
							
							CoordinateType pos = spectrum[indices_after[0]].getMZ();
							IntensityType inte = spectrum[indices_after[0]].getIntensity();
							
							//log_ << std::endl << "Peak: " << pos << std::endl;
							bool is_max_peak = true; //checking the maximum intensity peaks -> use them later as feature seeds.
							for (UInt i=1; i<=min_spectra_; ++i)
							{
								const SpectrumType& spec = map_->at(s+i);
								indices_after[i] = nearest_(pos, spec, indices_after[i]);
								DoubleReal position_score = positionScore_( pos, spec[indices_after[i]].getMZ(), trace_tolerance_);
								if (position_score >0 && spec[indices_after[i]].getIntensity()>inte) is_max_peak = false;
								scores.push_back(position_score);
							}
							indices_before[0] = indices_after[0];
							for (UInt i=1; i<=min_spectra_; ++i)
							{
								const SpectrumType& spec = map_->at(s-i);
								indices_before[i] = nearest_(pos, spec, indices_before[i]);
								DoubleReal position_score = positionScore_( pos, spec[indices_before[i]].getMZ(), trace_tolerance_);
								if (position_score>0 && spec[indices_before[i]].getIntensity()>inte) is_max_peak = false;
								scores.push_back(position_score);
							}
							
							//Calculate a consensus score out of the scores calculated before
							std::sort(scores.begin(), scores.end());
							DoubleReal trace_score = std::accumulate(scores.begin()+max_missing_trace_peaks_, scores.end(),0.0);
							trace_score /= (2*min_spectra_-max_missing_trace_peaks_);
							
							//store final score for later use
							info_[s][indices_after[0]].trace_score = trace_score;
							info_[s][indices_after[0]].local_max = is_max_peak;

							if (debug && trace_score>0.0)
							{
								PeakType tmp;
								tmp.setMZ(pos);
								tmp.setIntensity(trace_score);
								trace_map[s].push_back(tmp);
							}
							indices_after[0]++;
						}
					}
					//store mass trace score map
					if (debug)
					{
						MzDataFile().store("trace_scores.mzData",trace_map);
					}
					this->ff_->endProgress();
				}

				//-------------------------------------------------------------------------
				//Step 3:
				//Charge loop (create seeds and features for each charge separately)
				//-------------------------------------------------------------------------
				for (UInt c=(UInt)param_.getValue("isotopic_pattern:charge_low"); c<=(UInt)param_.getValue("isotopic_pattern:charge_high"); ++c)
				{
					std::vector<Seed> seeds;
					//initialize pattern scores
					for (UInt s=0; s<map_->size(); ++s)
					{
						for (UInt p=0; p<map_->at(s).size(); ++p)
						{
							info_[s][p].pattern_score = 0.0;
							info_[s][p].overall_score = 0.0;
						}
					}
					
					//-----------------------------------------------------------
					//Step 3.1: Precalculate IsotopePattern score
					//-----------------------------------------------------------
					this->ff_->startProgress(0, map_->size(), String("Calculating isotope pattern scores for charge ")+c);
					for (UInt s=0; s<map_->size(); ++s)
					{
						this->ff_->setProgress(s);
						const SpectrumType& spectrum = map_->at(s);
						for (UInt p=0; p<spectrum.size(); ++p)
						{
							DoubleReal mz = spectrum[p].getMZ();
							//get isotope distribution for this mass
							const TheoreticalIsotopePattern& isotopes = getIsotopeDistribution_(mz*c);
							//determine highest peak in isopope distribution
							UInt max_isotope = std::max_element(isotopes.intensity.begin(), isotopes.intensity.end()) - isotopes.intensity.begin();
							//Look up expected isotopic peaks (in the current spectrum or adjacent spectra)
							UInt peak_index = (UInt)spectrum.findNearest(mz-((DoubleReal)(isotopes.size()+1)/c));
							IsotopePattern pattern(isotopes.size());
							for (UInt i=0; i<isotopes.size(); ++i)
							{
								CoordinateType isotope_pos = mz + ((DoubleReal)i-max_isotope)/c;
								findIsotope_(isotope_pos, s, pattern, i, false, peak_index);
							}
							DoubleReal pattern_score = isotopeScore_(isotopes, pattern, true, false);
							
							//update pattern scores of all contained peaks (if necessary)
							if (pattern_score > 0.0)
							{
								//store debug info
								if (debug)
								{
									PeakType tmp;
									tmp.setPos(mz);
									tmp.setIntensity(pattern_score);
									pattern_map[s].push_back(tmp);
								}
								
								for (UInt i=0; i<pattern.peak.size(); ++i)
								{
									if (pattern.peak[i]>=0 && pattern_score>info_[pattern.spectrum[i]][pattern.peak[i]].pattern_score)
									{
										info_[pattern.spectrum[i]][pattern.peak[i]].pattern_score = pattern_score;
									}
								}
							}
						}
					}
					//store mass trace score map
					if (debug)
					{
						MzDataFile().store(String("pattern_scores_")+c+".mzData",pattern_map);
					}
					this->ff_->endProgress();
					
					//-----------------------------------------------------------
					//Step 3.2:
					//Find seeds for this charge
					//-----------------------------------------------------------		
					this->ff_->startProgress(0, map_->size(), String("Finding seeds for charge ")+c);
					for (UInt s=0; s<map_->size(); ++s)
					{
						this->ff_->setProgress(s);
						//do nothing for the first few and last few spectra as the scans required to search for traces are missing
						if (s<min_spectra_ || s>=map_->size()-min_spectra_)
						{
							continue;
						}
						const SpectrumType& spectrum = map_->at(s);
						//iterate over peaks
						for (UInt p=0; p<spectrum.size(); ++p)
						{	
							PeakInfo& info = info_[s][p];
							info.overall_score = std::pow(info.intensity_score*info.trace_score*info.pattern_score, 1.0f/3.0f);
							if (debug && info.overall_score>0.0)
							{
								PeakType tmp;
								tmp.setPos(map_->at(s)[p].getMZ());
								tmp.setIntensity(info.overall_score);
								selected_map[s].push_back(tmp);
							}
							//add seed to vector if certain conditions are fullfilled
							if (info.local_max && info.overall_score>0.2)
							{
								Seed seed;
								seed.spectrum = s;
								seed.peak = p;
								seed.intensity = map_->at(s)[p].getIntensity();								
								seeds.push_back(seed);
							}
						}
					}
					//sort seeds according to intensity
					std::sort(seeds.rbegin(),seeds.rend());
					//create and store seeds map
					if (debug)
					{
						for (UInt i=0; i<seeds.size(); ++i)
						{
							UInt spectrum = seeds[i].spectrum;
							UInt peak = seeds[i].peak;
							Feature tmp;
							tmp.setIntensity(seeds[i].intensity);
							tmp.setRT(map_->at(spectrum).getRT());
							tmp.setMZ(map_->at(spectrum)[peak].getMZ());
							seed_map.push_back(tmp);
						}
						FeatureXMLFile().store(String("seeds_")+c+".featureXML", seed_map);
						MzDataFile().store(String("selected_peaks_")+c+".mzData", selected_map);
					}
					this->ff_->endProgress();
					std::cout << "Found " << seeds.size() << " seeds for charge " << c << "." << std::endl;
					
					//------------------------------------------------------------------
					//Step 3.3:
					//Extension of seeds
					//------------------------------------------------------------------
					this->ff_->startProgress(0,seeds.size(), String("Extending seeds for charge ")+c);
					UInt added_features = 0;
					for (UInt i=0; i<seeds.size(); ++i)
					{
						//------------------------------------------------------------------
						//Step 3.3.1:
						//Extend all mass traces
						//------------------------------------------------------------------
						this->ff_->setProgress(i);
						log_ << std::endl << "Seed " << i+1 << ":" << std::endl;
						//If the intensity is zero this seed is already uses in another feature
						const SpectrumType& spectrum = map_->at(seeds[i].spectrum);
						const PeakType peak = spectrum[seeds[i].peak];
						log_ << " - Int: " << peak.getIntensity() << std::endl;
						log_ << " - RT: " << spectrum.getRT() << std::endl;
						log_ << " - MZ: " << peak.getMZ() << std::endl;
						if (seeds[i].intensity == 0.0)
						{
							abort_("Seed was already used");
							continue;
						}
						
						//----------------------------------------------------------------
						//Find best fitting isotope pattern for this charge (using averagene)
						IsotopePattern best_pattern(0);
						DoubleReal isotope_fit_quality = findBestIsotopeFit_(seeds[i], c, best_pattern);
						if (isotope_fit_quality<min_isotope_fit_)
						{
							abort_("Isotope pattern score too low");
							continue;
						}
						
						//extend the convex hull in m/z dimension (starting from the trace peaks)
						log_ << "Collecting mass traces" << std::endl;
						std::vector<MassTrace> traces;
						traces.reserve(best_pattern.peak.size());
						extendMassTraces_(best_pattern, traces);
						
						//Abort if too few traces were found
						if (traces.size()<2)
						{
							abort_("Found less than two mass traces");
							continue;
						}
						
						//------------------------------------------------------------------
						//Step 3.3.2:
						//Quality estimation
						//------------------------------------------------------------------
						log_ << "Quality estimation" << std::endl;
						
						//------------------------------------------------------------------
						//(1) isotope fit: isotope fit of the collected mass trace maxima
						log_ << " - Isotope fit: " << isotope_fit_quality << std::endl;
						
						//------------------------------------------------------------------
						//(2) overall shape: RT-spread of mass traces decreases with smaller intensities
						std::vector<DoubleReal> rts(traces.size());
						std::vector<DoubleReal> ints(traces.size());
						for (UInt j=0; j<traces.size(); ++j)
						{
							rts[j] = traces[j].peaks.back().first - traces[j].peaks[0].first;
							ints[j] = traces[j].max_peak->getIntensity();
						}
						DoubleReal overall_shape_quality = (Math::BasicStatistics<DoubleReal>::pearsonCorrelationCoefficient(rts.begin(),rts.end(),ints.begin(), ints.end())+1.0)/2.0;
						if (isnan(overall_shape_quality))
						{
							if (traces.size()==2) //for two traces it's ok to have the same width/intensity
							{
								overall_shape_quality = 0.5;
							}
							else //for more than two traces it is not ok to have the same width/intensity
							{
								overall_shape_quality = 0.1;
							}
						}
						log_ << " - overall shape: " << overall_shape_quality << std::endl;
		
						//------------------------------------------------------------------					
						//(3) trace m/z distances
						std::vector<DoubleReal> positions(traces.size());
						for (UInt j=0; j<traces.size(); ++j)
						{
							for (UInt k=0; k<traces[j].peaks.size(); ++k)
							{
								positions[j]+= traces[j].peaks[k].second->getMZ();
							}
							positions[j] /= traces[j].peaks.size();
						}
						DoubleReal mz_distance_quality = 0.0;
						for (UInt j=0; j<positions.size()-1; ++j)
						{
							mz_distance_quality += positionScore_(positions[j+1]-positions[j], 1.0/c, pattern_tolerance_);
						}
						mz_distance_quality /= positions.size()-1;
						log_ << " - mz distances: " << mz_distance_quality << std::endl;
	
						//------------------------------------------------------------------
						//(4) trace shape: trace intensity goes down towards the border
						UInt error_count = 0;
						for (UInt j=0; j<traces.size(); ++j)
						{
							UInt size = traces[j].peaks.size();
							if (size>=5)
							{
								DoubleReal max = traces[j].max_peak->getIntensity();
								//log_ << "------- max: " << max << std::endl;
								DoubleReal low_int = (traces[j].peaks[0].second->getIntensity() + traces[j].peaks[1].second->getIntensity()) / 2.0;
								//log_ << "------- low_int: " << low_int << std::endl;
								if (low_int / max > mass_trace_max_border_intensity)
								{
									//log_ << "------- error " << std::endl;
									++error_count;
								}
								DoubleReal high_int = (traces[j].peaks[size-2].second->getIntensity() + traces[j].peaks[size-1].second->getIntensity()) / 2.0;
								//log_ << "------- high_int: " << high_int << std::endl;
								if (high_int / max > mass_trace_max_border_intensity)
								{
									//log_ << "------- error " << std::endl;
									++error_count;
								}
							}
							else //increase error count by one as this trace is too small
							{
								++error_count;
							}
						}
						//Score: fraction of mal-formed rt profiles
						DoubleReal rt_shape_quality = 1.0 - (DoubleReal)error_count / (2.0*traces.size());					
						log_ << " - trace shape: " << rt_shape_quality << std::endl;
						
						//------------------------------------------------------------------					
						//(5) quality measure: maxima on one line
						//determine max peak RT and RT spread of that trace in both directions
						DoubleReal max = 0.0;
						DoubleReal max_rt = 0.0;
						DoubleReal spread_low = 0.0;
						DoubleReal spread_high = 0.0;
						for (UInt j=0; j<traces.size(); ++j)
						{
							if (traces[j].max_peak->getIntensity() > max)
							{
								max = traces[j].max_peak->getIntensity();
								max_rt = traces[j].max_rt;
								spread_low = std::max(0.01, max_rt - traces[j].peaks[0].first);
								spread_high = std::max(0.01, traces[j].peaks.back().first - max_rt);
							}
						}
						
						//look at max peak shifts of different scans
						DoubleReal rel_max_deviation = 0.0;
						for (UInt j=0; j<traces.size(); ++j)
						{
							if (traces[j].max_rt>max_rt)
							{
								rel_max_deviation += std::min(1.0,(traces[j].max_rt - max_rt) / spread_high);
							}
							else
							{
								rel_max_deviation += std::min(1.0,(max_rt - traces[j].max_rt) / spread_low);
							}
						}
						rel_max_deviation /= traces.size()-1;
						DoubleReal maxima_quality = 1.0 - rel_max_deviation;
						log_ << " - maxima positions: " << maxima_quality << std::endl;
	
						//----------------------------------------------------------------
						//abort if quality too low
						DoubleReal overall_quality_mean = std::pow(isotope_fit_quality * overall_shape_quality * mz_distance_quality * rt_shape_quality * maxima_quality, 1.0/5.0);
						log_ << " => final score: " << overall_quality_mean << std::endl;
						if (overall_quality_mean<min_feature_quality)
						{
							abort_("Feature quality too low");
							continue;
						}
	
						//------------------------------------------------------------------
						//Step 3.3.3:
						//Feature creation
						//------------------------------------------------------------------
						Feature f;
						f.setCharge(c);
						f.setOverallQuality(overall_quality_mean);
						if (debug)
						{
							f.setMetaValue("rt_shape",rt_shape_quality);
							f.setMetaValue("mz_distance",mz_distance_quality);
							f.setMetaValue("isotope_fit",isotope_fit_quality);
							f.setMetaValue("overall_shape",overall_shape_quality);
							f.setMetaValue("maxima_positions",maxima_quality);
						}
						
						//set feature position and intensity
						//feature RT is the average RT of the mass trace maxima
						//feature intensity and m/z are taken from the highest peak
						DoubleReal rt = 0.0;
						for (UInt j=0; j<traces.size(); ++j)
						{
							if (traces[j].max_peak->getIntensity()>f.getIntensity())
							{
								f.setIntensity(traces[j].max_peak->getIntensity());
								f.setMZ(traces[j].max_peak->getMZ());
							}
							rt += traces[j].max_rt;
						}
						f.setRT(rt/traces.size());
						//feature intensity is the sum of all the peak intensities
						if (!max_intensity)
						{
							DoubleReal int_sum = 0.0;
							for (UInt j=0; j<traces.size(); ++j)
							{
								for (UInt k=0; k<traces[j].peaks.size(); ++k)
								{
									int_sum += traces[j].peaks[k].second->getIntensity();
								}
							}
							f.setIntensity(int_sum);
						} 
						
						//add convex hulls of mass traces
						for (UInt j=0; j<traces.size(); ++j)
						{
							f.getConvexHulls().push_back(traces[j].getConvexhull());
						}
						//add feature to feature list
						this->features_->push_back(f);
						++added_features;
						log_ << "Feature number: " << this->features_->size() << std::endl;
	
	
						//----------------------------------------------------------------
						//Remove all seeds that lie inside the convex hull of the new feature
						DBoundingBox<2> bb = f.getConvexHull().getBoundingBox();
						for (UInt j=i+1; j<seeds.size(); ++j)
						{
							DoubleReal rt = map_->at(seeds[j].spectrum).getRT();
							DoubleReal mz = map_->at(seeds[j].spectrum)[seeds[j].peak].getMZ();
							if (bb.encloses(rt,mz) && f.encloses(rt,mz))
							{
								//set intensity to zero => the peak will be skipped!
								seeds[j].intensity = 0.0;
							}
						}
					}
					this->ff_->endProgress();
					std::cout << "Found " << added_features << " features candidates for charge " << c << "." << std::endl;
				}
				
				std::cout << std::endl;
				std::cout << "Abort reasons during feature construction:" << std::endl;
				for (std::map<String,UInt>::const_iterator it=aborts_.begin(); it!=aborts_.end(); ++it)
				{
					std::cout << "- " << it->first << ": " << it->second << std::endl;
				}
					
				//------------------------------------------------------------------
				//Step 4:
				//TODO: Resolve contradicting and overlapping features
				//------------------------------------------------------------------
				
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
			/// Output stream for log/debug info
			std::ofstream log_; 
			/// Array of abort reasons
			std::map<String, UInt> aborts_;
						
			///@name Members for parameters often needed in methods
			//@{
			DoubleReal pattern_tolerance_; ///< Stores mass_trace:mz_tolerance
			DoubleReal trace_tolerance_; ///< Stores isotopic_pattern:mz_tolerance
			UInt min_spectra_; ///< Number of spectra that have to show the same mass (for finding a mass trace)
			UInt max_missing_trace_peaks_; ///< Stores mass_trace:max_missing
			DoubleReal slope_bound_; ///< Max slope of mass trace intensities
			DoubleReal intensity_percentage_; ///< Isotope pattern intensity contribution of required peaks
			DoubleReal intensity_percentage_optional_; ///< Isotope pattern intensity contribution of optional peaks
			DoubleReal optional_fit_improvement_; ///< Minimal imrovment for leaving out optional isotope
			DoubleReal mass_window_width_; ///< Width of the isotope pattern mass bins
			UInt intensity_bins_; ///< Number of bins (in RT and MZ) for intensity significance estimation
			DoubleReal min_isotope_fit_; ///< Mimimum isotope pattern fit for a feature
			//@}

			///@name Members for intensity significance estimation
			//@{			
			/// RT bin width
			DoubleReal intensity_rt_step_;
			/// m/z bin width
			DoubleReal intensity_mz_step_;
			/// Precalculated threshold and maximum stored for each bin (rt,mz)
			std::vector< std::vector< std::pair<Real,Real> > > intensity_thresholds_;
			//@}

			///Precalculated information for each peak
			std::vector< std::vector<PeakInfo> > info_;
			
			///Vector of precalculated isotope distributions for several mass winows
			std::vector< TheoreticalIsotopePattern > isotope_distributions_;

			//Docu in base class
			virtual void updateMembers_()
			{
				pattern_tolerance_ = param_.getValue("mass_trace:mz_tolerance");
				trace_tolerance_ = param_.getValue("isotopic_pattern:mz_tolerance");
				min_spectra_ = (UInt)std::floor((DoubleReal)param_.getValue("mass_trace:min_spectra")*0.5);
				max_missing_trace_peaks_ = param_.getValue("mass_trace:max_missing");
				slope_bound_ = param_.getValue("mass_trace:slope_bound");
				intensity_percentage_ = (DoubleReal)param_.getValue("isotopic_pattern:intensity_percentage")/100.0;
				intensity_percentage_optional_ = (DoubleReal)param_.getValue("isotopic_pattern:intensity_percentage_optional")/100.0;
				optional_fit_improvement_ = (DoubleReal)param_.getValue("isotopic_pattern:optional_fit_improvement")/100.0;
				mass_window_width_ = param_.getValue("isotopic_pattern:mass_window_width");
				intensity_bins_ =  param_.getValue("intensity:bins");
				min_isotope_fit_ = param_.getValue("feature:min_isotope_fit");
				
			}
			
			///Writes the abort reason to the log file and counts occurences for each reason
			void abort_(const String& reason)
			{
				log_ << "Abort: " << reason << std::endl;
				aborts_[reason]++;
			}
			
			///Returns the isotope distribution for a certain mass window
			const TheoreticalIsotopePattern& getIsotopeDistribution_(DoubleReal mass)
			{
				//calculate index in the vector
				UInt index = (UInt)std::floor(mass/mass_window_width_);
				
				//enlarge vector if necessary
				if (index>=isotope_distributions_.size())
				{
					isotope_distributions_.resize(index+1);
				}
				
				//calculate distribution if necessary
				if (isotope_distributions_[index].intensity.size()==0)
				{
					//log_ << "Calculating iso dist for mass: " << 0.5*mass_window_width_ + index * mass_window_width_ << std::endl;
					IsotopeDistribution d;
					d.setMaxIsotope(10);
					d.estimateFromPeptideWeight(0.5*mass_window_width_ + index * mass_window_width_);
					d.trimLeft(intensity_percentage_optional_);
					d.trimRight(intensity_percentage_optional_);
					for (IsotopeDistribution::Iterator it=d.begin(); it!=d.end(); ++it)
					{
						isotope_distributions_[index].intensity.push_back(it->second);
						//log_ << " - " << it->second << std::endl;
					}
					//determine the number of optional peaks at the beginning/end
					UInt begin = 0;
					UInt end = 0;
					bool is_begin = true;
					bool is_end = false;
					for (UInt i=0; i<isotope_distributions_[index].intensity.size(); ++i)
					{
						if (isotope_distributions_[index].intensity[i]<intensity_percentage_)
						{
							if (!is_end && !is_begin) is_end = true;
							if (is_begin) ++begin;
							else if (is_end) ++end;
						}
						else if (is_begin)
						{
							is_begin = false;
						}
					}
					isotope_distributions_[index].optional_begin = begin;
					isotope_distributions_[index].optional_end = end;
					//log_ << " - optinal begin/end:" << begin << " / " << end << std::endl;
				}
				//Return distribution
				return isotope_distributions_[index];
			}
						
			/**
				@brief Finds the best fitting position of the isotopic pattern estimate defined by @p center
				
				@param center the maximum peak of the isotope distribution (contains charge as well)
				@param charge The charge of the pattern 
				@param best_pattern Returns the indices of the isotopic peaks. If a isopopic peak is missing -1 is returned.
			*/
			DoubleReal findBestIsotopeFit_(const Seed& center, UInt charge, IsotopePattern& best_pattern)
			{
				log_ << "Testing isotope patterns for charge " << charge << ": " << std::endl;			
				const SpectrumType& spectrum = map_->at(center.spectrum);
				const TheoreticalIsotopePattern& isotopes = getIsotopeDistribution_(spectrum[center.peak].getMZ()*charge);	
				log_ << " - Seed: " << center.peak << " (mz:" << spectrum[center.peak].getMZ()<< ")" << std::endl;
				
				//Find m/z boundaries of search space (linear search as this is local and we have the center already)
				DoubleReal mass_window = (DoubleReal)(isotopes.size()+1) / (DoubleReal)charge;
				log_ << " - Mass window: " << mass_window << std::endl;
				UInt end = center.peak;
				while(end<spectrum.size() && spectrum[end].getMZ()<spectrum[center.peak].getMZ()+mass_window)
				{
					++end;
				}
				--end;
				//search begin
				Int begin = center.peak;
				while(begin>=0 && spectrum[begin].getMZ()>spectrum[center.peak].getMZ()-mass_window)
				{
					--begin;
				}
				++begin;
				log_ << " - Begin: " << begin << " (mz:" << spectrum[begin].getMZ()<< ")" << std::endl;
				log_ << " - End: " << end << " (mz:" << spectrum[end].getMZ()<< ")" << std::endl;

				//fit isotope distribution to peaks
				DoubleReal max_score = 0.0;
				for (UInt start=begin; start<=end; ++start)
				{
					//find isotope peaks for the current start peak
					UInt peak_index = start;
					IsotopePattern pattern(isotopes.size());
					pattern.intensity[0] = spectrum[start].getIntensity();
					pattern.peak[0] = start;
					pattern.spectrum[0] = center.spectrum;
					pattern.mz_score[0] = 1.0;
					pattern.theoretical_mz[0] = spectrum[start].getMZ();
					log_ << " - Fitting at " << start << " (mz:" << spectrum[start].getMZ() << ")" << std::endl;
					log_ << "   - Isotope 0: " << pattern.intensity[0] << std::endl;
					for (UInt iso=1; iso<isotopes.size(); ++iso)
					{
						DoubleReal pos = spectrum[start].getMZ() + iso/(DoubleReal)charge;
						findIsotope_(pos, center.spectrum, pattern, iso, true, peak_index);
					}
					
					//check if the seed is contained, otherwise abort
					bool seed_contained = false;
					for (UInt iso=0; iso<pattern.peak.size(); ++iso)
					{
						if (pattern.peak[iso]==(Int)center.peak && pattern.spectrum[iso]==center.spectrum)
						{
							seed_contained = true;
							break;
						}
					}
					if(!seed_contained)
					{
						log_ << "   - aborting: seed is not contained!" << std::endl;
						continue;
					}

					DoubleReal score = isotopeScore_(isotopes, pattern, false, true);
	
					//check if the seed is still contained, otherwise abort
					seed_contained = false;
					for (UInt iso=0; iso<pattern.peak.size(); ++iso)
					{
						if (pattern.peak[iso]==(Int)center.peak && pattern.spectrum[iso]==center.spectrum)
						{
							seed_contained = true;
							break;
						}
					}
					if(!seed_contained)
					{
						log_ << "   - aborting: seed was removed during isotope fit!" << std::endl;
						continue;
					}
					
					log_ << "   - final score: " << score << std::endl;
					if (score>max_score)
					{
						max_score = score;
						best_pattern = pattern;
					}
				}
				log_ << " - best score: " << max_score << std::endl;
				return max_score;
			}
			
			//TODO Extend simultaniously
			///Extends all mass traces of a isotope pattern in one step
			void extendMassTraces_(const IsotopePattern& pattern, std::vector<MassTrace>& traces)
			{
				for (UInt p=0; p<pattern.peak.size(); ++p)
				{
					log_ << " - Trace " << p << std::endl;
					Seed starting_peak;
					starting_peak.spectrum = pattern.spectrum[p];
					starting_peak.peak = pattern.peak[p];
					if (pattern.peak[p]==-2)
					{
						log_ << "   - removed during isotope fit" << std::endl;
						continue;
					}
					else if (pattern.peak[p]==-1)
					{
						log_ << "   - missing" << std::endl;
						continue;
					}
					starting_peak.intensity = map_->at(starting_peak.spectrum)[starting_peak.peak].getIntensity();
					log_ << "   - trace seed: " << map_->at(starting_peak.spectrum).getRT() << " / " << map_->at(starting_peak.spectrum)[starting_peak.peak].getMZ() << " (int: " << map_->at(starting_peak.spectrum)[starting_peak.peak].getIntensity() << ")" << std::endl;
					
					//search for nearby maximum of the mass trace as the extension assumes that it starts at the maximum
					UInt begin = std::max(0u,starting_peak.spectrum-min_spectra_);
					UInt end = std::min(starting_peak.spectrum+min_spectra_,(UInt)map_->size());
					DoubleReal mz = map_->at(starting_peak.spectrum)[starting_peak.peak].getMZ();
					DoubleReal inte = map_->at(starting_peak.spectrum)[starting_peak.peak].getIntensity();
					for (UInt spectrum_index=begin; spectrum_index<end; ++spectrum_index)
					{
						//find better seeds (no-empty scan/low mz diff/higher intensity)
						Int peak_index = map_->at(spectrum_index).findNearest(map_->at(starting_peak.spectrum)[starting_peak.peak].getMZ());
						if (peak_index==-1 ||
								map_->at(spectrum_index)[peak_index].getIntensity()<=inte ||
								std::fabs(mz-map_->at(spectrum_index)[peak_index].getMZ())>=pattern_tolerance_
							 ) continue;
						starting_peak.spectrum = spectrum_index;
						starting_peak.peak = peak_index;
						inte = map_->at(spectrum_index)[peak_index].getIntensity();
					}
					log_ << "   - extending from: " << map_->at(starting_peak.spectrum).getRT() << " / " << map_->at(starting_peak.spectrum)[starting_peak.peak].getMZ() << " (int: " << map_->at(starting_peak.spectrum)[starting_peak.peak].getIntensity() << ")" << std::endl;
					
					//------------------------------------------------------------------
					//Extend seed to a mass trace
					MassTrace trace;
					const PeakType& seed = map_->at(starting_peak.spectrum)[starting_peak.peak];
					//initialize trace with seed data
					trace.max_peak = &seed;
					trace.max_rt = map_->at(starting_peak.spectrum).getRT();
					//extend in downstream direction
					extendMassTrace_(trace, starting_peak.spectrum -1, seed.getMZ(), false);
					//invert peak array to bring peaks in the correct cronological order
					std::reverse(trace.peaks.begin(), trace.peaks.end());
					//extend in upstream direction
					extendMassTrace_(trace, starting_peak.spectrum +1, seed.getMZ(), true);
					
					//check if enough peaks were found
					if (trace.peaks.size()<3)
					{
						log_ << "   - could not extend trace " << std::endl;
						continue;
					}
					traces.push_back(trace);
				}
			}

			/**
				@brief Extends a single mass trace in one RT direction
				
				@note this method assumes that it extends from a local maximum. Otherwise it will not work!
			*/
			void extendMassTrace_(MassTrace& trace, Int spectrum_index, DoubleReal mz, bool inc_rt)
			{
				std::vector<DoubleReal> ints(2, trace.max_peak->getIntensity());
				std::vector<DoubleReal> smoothed(3, trace.max_peak->getIntensity());
				UInt missing_peaks = 0;
				UInt added_peaks = 0;
				bool remove_last_peaks = false;
				std::pair<DoubleReal, const PeakType*> last_peak = std::make_pair(trace.max_rt, trace.max_peak);
				while((!inc_rt && spectrum_index>=0) || (inc_rt && spectrum_index<(Int)map_->size()))
				{
					Int peak_index = map_->at(spectrum_index).findNearest(mz);
					if (peak_index<0 || info_[spectrum_index][peak_index].overall_score<0.1 || positionScore_( mz, map_->at(spectrum_index)[peak_index].getMZ(), trace_tolerance_)<=0.0)
					{
						++missing_peaks;
						if(missing_peaks>max_missing_trace_peaks_)
						{
							//log_ << " # Too many peaks missing!" << std::endl;
							//add last peak to trace
							trace.peaks.push_back(last_peak);
							break;
						}
					}
					else
					{
						missing_peaks = 0;
						ints.push_back(map_->at(spectrum_index)[peak_index].getIntensity());
						//add last peak to trace
						trace.peaks.push_back(last_peak);
						++added_peaks;
						//update smoothed data
						smoothed.push_back(std::accumulate(ints.end()-3, ints.end(),0.0)/3.0);
						//update maximum peak of trace (smoothed to avoid maxima caused by noise)
						if (smoothed.back()>trace.max_peak->getIntensity())
						{
							trace.max_peak = last_peak.second;
							trace.max_rt = last_peak.first;
							//log_ << " # Updated Maximum to " << smoothed.back() << " / " << trace.max_peak->getIntensity() << std::endl;
						}
						//update last peak
						last_peak = std::make_pair(map_->at(spectrum_index).getRT(), &(map_->at(spectrum_index)[peak_index]));
						//Calculate deltas
						UInt last = smoothed.size()-1;
						DoubleReal delta1 = (smoothed[last-2]-smoothed[last-3])/smoothed[last-3];
						DoubleReal delta2 = (smoothed[last-1]-smoothed[last-2])/smoothed[last-2];
						DoubleReal delta3 = (smoothed[last]-smoothed[last-1])/smoothed[last-1];
						//Abort if the last three deltas are positive
						if (delta1>0.0 && delta2>0.0 && delta3>0.0)
						{
							//log_ << " # Deltas all positive!" << std::endl;
							remove_last_peaks = true;
							break;
						}
						//Abort if the average delta is too big (as intensity increases then)
						if ((delta1+delta2+delta3)/3.0>slope_bound_)
						{
							//log_ << " # Deltas too positive!" << std::endl;
							remove_last_peaks = true;
							break;
						}
					}
					//increase/decrease scan index
					if (inc_rt) ++spectrum_index; else --spectrum_index;
				}
				//Go back several peaks if we extended too far
				//Update maximum, if we removed it
				if (remove_last_peaks)
				{
					bool max_removed = false;
					UInt remove = std::min(added_peaks,2u);
					for(UInt i=remove; i>0; --i)
					{
						if (trace.peaks.back().second==trace.max_peak) max_removed = true;
						trace.peaks.pop_back();
					}
					if (max_removed)
					{
						trace.updateMaximum();
					}
					added_peaks -= remove;
				}
				log_ << "   - Added " << added_peaks << " peaks" << std::endl;
			}

			/// Returns the index of the peak nearest to m/z @p pos in spectrum @p spec (linear search starting from index @p start)
			template <typename SpectrumType>
			UInt nearest_(CoordinateType pos, const SpectrumType& spec, UInt start) const
			{
				UInt index = start;
				CoordinateType dist = std::fabs(pos-spec[index].getMZ());
				++index;
				while (index < spec.size())
				{
					CoordinateType new_dist = std::fabs(pos-spec[index].getMZ());
					if (new_dist<dist)
					{
						dist = new_dist;
						++index;	
					}
					else
					{
						break;
					}
				}
				return --index; 
			}
			
			/**
				@brief Searches for an isotopic peak in the current spectrum and the adjacent spectra
				
				@param pos m/z position of the searched for peak
				@param spectrum_index index of the central spectrum
				@param pattern IsotopePattern to store found peaks
				@param pattern_index index of the isotope in the pattern
				@param debug Flag that turn on debug info
				@param peak_index starting index of the search (to avoid multiple binary searches)
			*/
			void findIsotope_(CoordinateType pos, UInt spectrum_index, IsotopePattern& pattern, UInt pattern_index, bool debug, UInt& peak_index)
			{
				//search in the center spectrum
				const SpectrumType& spectrum = map_->at(spectrum_index);
				peak_index = nearest_(pos, spectrum, peak_index);
				DoubleReal mz_score = positionScore_(pos, spectrum[peak_index].getMZ(), pattern_tolerance_);
				pattern.theoretical_mz[pattern_index] = pos;
				if (mz_score!=0.0)
				{
					if (debug) log_ << "   - Isotope " << pattern_index << ": " << spectrum[peak_index].getIntensity() << std::endl;
					pattern.peak[pattern_index] = peak_index;
					pattern.spectrum[pattern_index] = spectrum_index;
					pattern.mz_score[pattern_index] = mz_score;
					pattern.intensity[pattern_index] = spectrum[peak_index].getIntensity();
					return;
				}
				//try to find the mass in the previous spectrum
				if (spectrum_index!=0)
				{
					const SpectrumType& spectrum_before = map_->at(spectrum_index-1);
					Int index_before = spectrum_before.findNearest(pos);
					if (index_before!=-1 && positionScore_(pos, spectrum_before[index_before].getMZ(), pattern_tolerance_)!=0.0)
					{
						if (debug) log_ << "   - Isotope " << pattern_index << ": " << spectrum_before[index_before].getIntensity() << " - previous spectrum" << std::endl;
						pattern.peak[pattern_index] = index_before;
						pattern.spectrum[pattern_index] = spectrum_index-1;
						pattern.mz_score[pattern_index] = positionScore_(pos, spectrum_before[index_before].getMZ(), pattern_tolerance_);
						pattern.intensity[pattern_index] = spectrum_before[index_before].getIntensity();
						return;
					}
				}
				//try to find the mass in the next spectrum
				if (spectrum_index!=map_->size()-1)
				{
					const SpectrumType& spectrum_after = map_->at(spectrum_index+1);
					Int index_after = spectrum_after.findNearest(pos);
					if (index_after!=-1 && positionScore_(pos, spectrum_after[index_after].getMZ(), pattern_tolerance_)!=0.0)
					{
						if (debug) if (debug) log_ << "   - Isotope " << pattern_index << ": " << spectrum_after[index_after].getIntensity() << " - next spectrum" << std::endl;
						pattern.peak[pattern_index] = index_after;
						pattern.spectrum[pattern_index] = spectrum_index+1;
						pattern.mz_score[pattern_index] = positionScore_(pos, spectrum_after[index_after].getMZ(), pattern_tolerance_);
						pattern.intensity[pattern_index] = spectrum_after[index_after].getIntensity();
						return;
					}
				}
				//no isotope found
				if (debug) log_ << "   - Isotope " << pattern_index << ": missing" << std::endl;
				pattern.peak[pattern_index] = -1;
				pattern.mz_score[pattern_index] = 0.0;
				pattern.intensity[pattern_index] = 0.0;
			}

			/// Calculates a score between 0 and 1 for the m/z deviation of two peaks.
			DoubleReal positionScore_(CoordinateType pos1, CoordinateType pos2, DoubleReal allowed_deviation) const
			{
				DoubleReal diff = fabs(pos1 - pos2);
				if (diff <= 0.5*allowed_deviation)
				{
					return 0.1*(0.5*allowed_deviation-diff)/(0.5*allowed_deviation)+0.9;
				}
				else if (diff <= allowed_deviation)
				{
					return 0.9*(allowed_deviation-diff)/(0.5*allowed_deviation);
				}
				return 0.0;
			}

			/// Calculates a score between 0 and 1 for the correlation between theoretical and found isotope pattern
			DoubleReal isotopeScore_(const TheoreticalIsotopePattern& isotopes, IsotopePattern& pattern, bool consider_mz_distances, bool debug)
			{
				//Abort if a core peak is missing
				for (UInt iso=0+isotopes.optional_begin; iso<pattern.peak.size()-isotopes.optional_end; ++iso)
				{
					if (pattern.peak[iso]==-1)
					{
						if (debug) log_ << "   - aborting: core peak is missing" << std::endl;
						return 0.0;
					}
				}
				//Find best isotope fit
				// - try to leave out optional isotope peaks to improve the fit
				// - do not allow gaps inside the pattern
				DoubleReal best_int_score = 0.01; //Not 0 as this would result in problems when checking for the percental improvement
				UInt best_begin = 0;
				for (UInt i=isotopes.optional_begin; i>0; --i)
				{
					if (pattern.peak[i-1]==-1)
					{
						best_begin = i;
						break;
					}
				}
				UInt best_end = 0;
				for (UInt i=isotopes.optional_end; i>0; --i)
				{
					if (pattern.peak[pattern.peak.size()-i]==-1)
					{
						best_end = i;
						break;
					}
				}
				if (debug) log_ << "   - best_begin/end: " << best_begin << "/" << best_end << std::endl;
				for (UInt b=best_begin; b<=isotopes.optional_begin; ++b)
				{
					for (UInt e=best_end; e<=isotopes.optional_end; ++e)
					{
						//Make sure we have more than 2 peaks (unless in the first loop interation) 
						if (isotopes.size()-b-e>2 || (b==best_begin && e==best_end))
						{
							DoubleReal int_score = Math::BasicStatistics<DoubleReal>::pearsonCorrelationCoefficient(isotopes.intensity.begin()+b, isotopes.intensity.end()-e, pattern.intensity.begin()+b, pattern.intensity.end()-e);	
							if (isnan(int_score)) int_score = 0.0;
							if (isotopes.size()-b-e==2 && int_score>min_isotope_fit_) int_score = min_isotope_fit_; //special case for the first loop iteration (otherwise the score is 1)
							if (debug) log_ << "   - fit (" << b << "/" << e << "): " << int_score;
							if (int_score/best_int_score>=1.0+optional_fit_improvement_)
							{
								if (debug) log_ << " - new best fit ";
								best_int_score = int_score;
								best_begin = b;
								best_end = e;
							}
							if (debug) log_ << std::endl;
						}
					}
				}
				//remove left out peaks from the beginning
				for (UInt i=0; i<best_begin; ++i)
				{
					pattern.peak[i] = -2;
					pattern.intensity[i] = 0.0;
					pattern.mz_score[i] = 0.0;
				}
				//remove left out peaks from the end
				for (UInt i=0; i<best_end; ++i)
				{
					pattern.peak[isotopes.size()-1-i] = -2;
					pattern.intensity[isotopes.size()-1-i] = 0.0;
					pattern.mz_score[isotopes.size()-1-i] = 0.0;
				}
				//calculate m/z score (if required)
				if (consider_mz_distances)
				{
					best_int_score *= std::accumulate(pattern.mz_score.begin()+best_begin, pattern.mz_score.end()-best_end,0.0) / (pattern.mz_score.size()-best_begin-best_end);
				}
				//return final score
				return best_int_score;
			}
			
			DoubleReal intensityScore_(DoubleReal intensity, UInt spectrum, UInt peak)
			{
				UInt rt_bin = std::min(intensity_bins_-1,(UInt)std::floor((map_->at(spectrum).getRT() - map_->getMinRT()) / intensity_rt_step_));
				UInt mz_bin = std::min(intensity_bins_-1,(UInt)std::floor((map_->at(spectrum)[peak].getMZ() - map_->getMinMZ()) / intensity_mz_step_));
						
				DoubleReal threshold = intensity_thresholds_[rt_bin][mz_bin].first;
				DoubleReal maximum = intensity_thresholds_[rt_bin][mz_bin].second;
				if (intensity>threshold)
				{
					return (intensity-threshold)/(maximum-threshold);
				}
				return 0.0;
			}

		private:
			
			/// Not implemented
			FeatureFinderAlgorithmPicked& operator=(const FeatureFinderAlgorithmPicked&);
			/// Not implemented
			FeatureFinderAlgorithmPicked(const FeatureFinderAlgorithmPicked&);

	};

} // namespace OpenMS

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMPICKED_H
