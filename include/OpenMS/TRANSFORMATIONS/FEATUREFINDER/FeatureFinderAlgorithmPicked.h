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
// $Maintainer: Marc Sturm $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMPICKED_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMPICKED_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithm.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>

#include <boost/math/special_functions/fpclassify.hpp>

#include <numeric>
#include <fstream>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>

#include <QtCore/QDir>

#ifdef _OPENMP
#ifdef OPENMS_WINDOWSPLATFORM
#include <omp.h>
#endif
#endif

namespace OpenMS
{
	/** 
		@brief FeatureFinderAlgorithm for picked peaks.

    @htmlinclude OpenMS_FeatureFinderAlgorithmPicked.parameters

		@improvement RT model with tailing/fronting (Marc)
		@improvement More general MZ model - e.g. based on co-elution or with sulphur-averagines (Marc)
		
		@todo Fix output in parallel mode, change assignment of charges to threads, add parallel TOPP test (Marc)
		@todo Implement user-specified seed lists support (Marc)
		
	
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
			typedef typename FeatureFinderAlgorithm<PeakType, FeatureType>::FeatureMapType FeatureMapType;
			typedef typename MapType::SpectrumType SpectrumType;
			typedef typename SpectrumType::FloatDataArrays FloatDataArrays;
			//@}
			
			using FeatureFinderAlgorithm<PeakType, FeatureType>::param_;
			using FeatureFinderAlgorithm<PeakType, FeatureType>::features_;
			using FeatureFinderAlgorithm<PeakType, FeatureType>::ff_;
			using FeatureFinderAlgorithm<PeakType, FeatureType>::defaults_;
				
		protected:
			///Helper structure for seeds used in FeatureFinderAlgorithmPicked
			struct Seed
			{
				///Spectrum index
				Size spectrum;
				///Peak index
				Size peak;
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

				///Theoretical intensity value (scaled to [0,1])
				DoubleReal theoretical_int;
				
				///Contained peaks (pair of RT and pointer to peak)
				std::vector<std::pair<DoubleReal, const PeakType*> > peaks;
				
				///determindes the convex hull of the trace
				ConvexHull2D getConvexhull() const
				{
					ConvexHull2D::PointArrayType hull_points(peaks.size());
					for (Size i=0; i<peaks.size(); ++i)
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
					
					for (Size i=1; i<peaks.size(); ++i)
					{
						if (peaks[i].second->getIntensity()>max_peak->getIntensity())
						{
							max_rt = peaks[i].first;					
							max_peak = peaks[i].second;
						}
					}
				}

				///Returns the average m/z of all peaks in this trace (weighted by intensity)
				DoubleReal getAvgMZ() const
				{
					DoubleReal sum = 0.0;
					DoubleReal intensities = 0.0;
					for (Size i=0; i<peaks.size(); ++i)
					{
						sum += peaks[i].second->getMZ()*peaks[i].second->getIntensity();
						intensities += peaks[i].second->getIntensity();
					}
					return sum / intensities;
				}
				
				///Checks if this Trace is valid (has more than 2 points)
				bool isValid() const
				{
					return (peaks.size()>=3);
				}
				
			};
			
			///Helper struct for a collection of mass traces used in FeatureFinderAlgorithmPicked
			struct MassTraces
				: public std::vector<MassTrace>
			{
				/// Constructor
				MassTraces()
					: max_trace(0)
				{
				}
				
				/// Returns the peak count of all traces
				Size getPeakCount() const
				{
					Size sum = 0;
					for (Size i=0; i<this->size(); ++i)
					{
						sum += this->at(i).peaks.size();
					}
					return sum;
				}
				
				///Checks if still valid (seed still contained and enough traces)
				bool isValid(DoubleReal seed_mz, DoubleReal trace_tolerance)
				{
					//Abort if too few traces were found
					if (this->size()<2) return false;

					//Abort if the seed was removed
					for (Size j=0; j<this->size(); ++j)
					{
						if (std::fabs(seed_mz-this->at(j).getAvgMZ())<=trace_tolerance)
						{
							return true;
						}
					}
					return false;
				}
				
				/**
				  @brief Returns the theoretical maximum trace index

				  @exception Exception::Precondition is thrown if there are not mass traces (not only in debug mode)
				*/
				Size getTheoreticalMax() const
				{
					if (!this->size())
					{
						throw Exception::Precondition(__FILE__,__LINE__,__PRETTY_FUNCTION__,"There must be at least one trace to determine the theoretical maximum trace!");
					}
					
					Size max=0;
					DoubleReal max_int=this->at(0).theoretical_int;
					for (Size i=1; i<this->size(); ++i)
					{
						if (this->at(i).theoretical_int>max_int)
						{
							max_int = this->at(i).theoretical_int;
							max = i;
						}
					}
					return max;
				}

				///Sets the baseline to the lowest contained peak of the trace
				void updateBaseline()
				{
					if (this->size()==0)
					{
						baseline = 0.0;
						return;
					}
					bool first = true;					
					for (Size i=0; i<this->size(); ++i)
					{
						for (Size j=0; j<this->at(i).peaks.size(); ++j)
						{
							if (first)
							{
								baseline = baseline = this->at(i).peaks[j].second->getIntensity();
								first = false;
							}
							if (this->at(i).peaks[j].second->getIntensity()<baseline)
							{
								baseline = this->at(i).peaks[j].second->getIntensity();
							}
						}
					}
				}

				/**
				  @brief Returns the RT boundaries of the mass traces

				  @exception Exception::Precondition is thrown if there are no mass traces (not only in debug mode)
				*/
				std::pair<DoubleReal,DoubleReal> getRTBounds() const
				{
					if (!this->size())
					{
						throw Exception::Precondition(__FILE__,__LINE__,__PRETTY_FUNCTION__,"There must be at least one trace to determine the RT boundaries!");
					}
					
					DoubleReal min = std::numeric_limits<DoubleReal>::max();
					DoubleReal max = -std::numeric_limits<DoubleReal>::max();
					//Abort if the seed was removed
					for (Size i=0; i<this->size(); ++i)
					{
						for (Size j=0; j<this->at(i).peaks.size(); ++j)
						{
							DoubleReal rt = this->at(i).peaks[j].first;
							if (rt>max) max = rt;
							if (rt<min) min = rt;
						}
					}
					return std::make_pair(min,max);
				}

				/// Maximum intensity trace
				Size max_trace;
				/// Estimated baseline in the region of the feature (used for the fit)
				DoubleReal baseline;
			};
			
			///Helper structure for a theoretical isotope pattern used in FeatureFinderAlgorithmPicked
			struct TheoreticalIsotopePattern
			{
				///Vector of intensity contributions 
				std::vector<DoubleReal> intensity;
				///Number of optional peaks at the beginning of the pattern
				Size optional_begin;
				///Number of optional peaks at the end of the pattern
				Size optional_end;
				///The maximum intensity contribution before scaling the pattern to 1
				DoubleReal max;
				///The number of isotopes trimmed on the left side. This is needed to reconstruct the monoisotopic peak.
				Size trimmed_left;
				/// Returns the size
				Size size() const
				{
					return intensity.size();
				}
			};

			///Helper structure for a found isotope pattern used in FeatureFinderAlgorithmPicked
			struct IsotopePattern
			{
				///Peak index (-1 if peak was not found, -2 if it was removed to improve the isotope fit)
				std::vector<SignedSize> peak;
				///Spectrum index (undefined if peak index is -1 or -2)
				std::vector<Size> spectrum;
				///Peak intensity (0 if peak index is -1 or -2)
				std::vector<DoubleReal> intensity;
				///m/z score of peak (0 if peak index is -1 or -2)
				std::vector<DoubleReal> mz_score;
				///Theoretical m/z value of the isotope peak
				std::vector<DoubleReal> theoretical_mz;
				///Theoretical isotope pattern
				TheoreticalIsotopePattern theoretical_pattern;
				
				/// Constructor that resizes the internal vectors
				IsotopePattern(Size size)
					: peak(size,-1),
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
					map_(),
					log_()
			{
				//debugging
				defaults_.setValue("debug","false","When debug mode is activated, several files with intermediate results are written to the folder 'debug' (do not use in parallel mode).");
				defaults_.setValidStrings("debug",StringList::create("true,false"));
				//intensity
				defaults_.setValue("intensity:bins",10,"Number of bins per dimension (RT and m/z). The higher this value, the more local the intensity significance score is.\nThis parameter should be decreased, if the algorithm is used on small regions of a map.");
				defaults_.setMinInt("intensity:bins",1);
				defaults_.setSectionDescription("intensity","Settings for the calculation of a score indicating if a peak's intensity is significant in the local environment (between 0 and 1)");
				//mass trace search parameters
				defaults_.setValue("mass_trace:mz_tolerance",0.03,"Tolerated m/z deviation of peaks belonging to the same mass trace.\nIt should be larger than the m/z resolution of the instument.\nThis value must be smaller than that 1/charge_high!");
				defaults_.setMinFloat("mass_trace:mz_tolerance",0.0);
				defaults_.setValue("mass_trace:min_spectra",10,"Number of spectra that have to show a similar peak mass in a mass trace.");
				defaults_.setMinInt("mass_trace:min_spectra",1);
				defaults_.setValue("mass_trace:max_missing",1,"Number of spectra where a high mass deviation or missing peak is acceptable.\nThis parameter should be well below 'min_spectra'!");
				defaults_.setMinInt("mass_trace:max_missing",0);
				defaults_.setValue("mass_trace:slope_bound",0.1,"The maximum slope of mass trace intensities when extending from the highest peak.\nThis parameter is important to seperate overlapping elution peaks.\nIt should be increased if feature elution profiles fluctuate a lot.");
				defaults_.setMinFloat("mass_trace:slope_bound",0.0);
				defaults_.setSectionDescription("mass_trace","Settings for the calculation of a score indicating if a peak is part of a mass trace (between 0 and 1).");
				//Isotopic pattern search parameters
				defaults_.setValue("isotopic_pattern:charge_low",1,"Lowest charge to search for.");
				defaults_.setMinInt("isotopic_pattern:charge_low",1);
				defaults_.setValue("isotopic_pattern:charge_high",4,"Highest charge to search for.");
				defaults_.setMinInt("isotopic_pattern:charge_high",1);
				defaults_.setValue("isotopic_pattern:mz_tolerance",0.03,"Tolerated m/z deviation from the theoretical isotopic pattern.\nIt should be larger than the m/z resolution of the instument.\nThis value must be smaller than that 1/charge_high!");		
				defaults_.setMinFloat("isotopic_pattern:mz_tolerance",0.0);
				defaults_.setValue("isotopic_pattern:intensity_percentage",10.0,"Isotopic peaks that contribute more than this percentage to the overall isotope pattern intensity must be present.", StringList::create("advanced"));
				defaults_.setMinFloat("isotopic_pattern:intensity_percentage",0.0);
				defaults_.setMaxFloat("isotopic_pattern:intensity_percentage",100.0);
				defaults_.setValue("isotopic_pattern:intensity_percentage_optional",0.1,"Isotopic peaks that contribute more than this percentage to the overall isotope pattern intensity can be missing.", StringList::create("advanced"));
				defaults_.setMinFloat("isotopic_pattern:intensity_percentage_optional",0.0);
				defaults_.setMaxFloat("isotopic_pattern:intensity_percentage_optional",100.0);
				defaults_.setValue("isotopic_pattern:optional_fit_improvement",2.0,"Minimal percental improvement of isotope fit to allow leaving out an optional peak.", StringList::create("advanced"));
				defaults_.setMinFloat("isotopic_pattern:optional_fit_improvement",0.0);
				defaults_.setMaxFloat("isotopic_pattern:optional_fit_improvement",100.0);
				defaults_.setValue("isotopic_pattern:mass_window_width",25.0,"Window width in Dalton for precalcuation of estimated isotope distribtions.", StringList::create("advanced"));
				defaults_.setMinFloat("isotopic_pattern:mass_window_width",1.0);
				defaults_.setMaxFloat("isotopic_pattern:mass_window_width",200.0);
				defaults_.setSectionDescription("isotopic_pattern","Settings for the calculation of a score indicating if a peak is part of a isotoipic pattern (between 0 and 1).");
				//Seed settings
				defaults_.setValue("seed:min_score",0.8,"Minimum seed score a peak has to reach to be used as seed.\nThe seed score is the geometric mean of intensity score, mass trace score and isotope pattern score.\nIf your features show a large deviation from the averagene isotope distribution or from an gaussian elution profile, lower this score.");
				defaults_.setMinFloat("seed:min_score",0.0);
				defaults_.setMaxFloat("seed:min_score",1.0);
				defaults_.setSectionDescription("seed","Settings that determine which peaks are considered a seed");
				//Fitting settings
				defaults_.setValue("fit:epsilon_abs",0.0001,"Absolute epsilon used for convergence of the fit.", StringList::create("advanced"));
				defaults_.setMinFloat("fit:epsilon_abs",0.0);
				defaults_.setValue("fit:epsilon_rel",0.0001,"Relative epsilon used for convergence of the fit.", StringList::create("advanced"));
				defaults_.setMinFloat("fit:epsilon_rel",0.0);
				defaults_.setValue("fit:max_iterations",500,"Maximum number of iterations of the fit.", StringList::create("advanced"));
				defaults_.setMinInt("fit:max_iterations",1);
				defaults_.setSectionDescription("fit","Settings for the model fitting");
				//Feature settings
				defaults_.setValue("feature:min_score",0.7, "Feature score threshold for a feature to be reported.\nThe feature score is the geometric mean of the average relative deviation and the correlation between the model and the observed peaks.");
				defaults_.setMinFloat("feature:min_score",0.0);
				defaults_.setMaxFloat("feature:min_score",1.0);
				defaults_.setValue("feature:min_isotope_fit",0.8,"Minimum isotope fit of the feature before model fitting.", StringList::create("advanced"));
				defaults_.setMinFloat("feature:min_isotope_fit",0.0);
				defaults_.setMaxFloat("feature:min_isotope_fit",1.0);
				defaults_.setValue("feature:min_trace_score",0.5, "Trace score threshold.\nTraces below this threshold are removed after the model fitting.\nThis parameter is important for features that overlap in m/z dimension.", StringList::create("advanced"));
				defaults_.setMinFloat("feature:min_trace_score",0.0);
				defaults_.setMaxFloat("feature:min_trace_score",1.0);
				defaults_.setValue("feature:min_rt_span",0.333, "Minimum RT span in relation to extended area that has to remain after model fitting.", StringList::create("advanced"));
				defaults_.setMinFloat("feature:min_rt_span",0.0);
				defaults_.setMaxFloat("feature:min_rt_span",1.0);
				defaults_.setValue("feature:max_rt_span",2.5, "Maximum RT span in relation to extended area that the model is allowed to have.", StringList::create("advanced"));
				defaults_.setMinFloat("feature:max_rt_span",0.5);
				defaults_.setValue("feature:max_intersection",0.35, "Maximum allowed intersection of features.", StringList::create("advanced"));
				defaults_.setMinFloat("feature:max_intersection",0.0);
				defaults_.setMaxFloat("feature:max_intersection",1.0);
				defaults_.setValue("feature:reported_mz","monoisotopic", "The mass type that is reported for features.\n'maximum' returns the m/z value of the highest mass trace.\n'average' returns the intensity-weighted average m/z value of all contained peaks.\n'monoisotopic' returns the monoisotopic m/z value derived from the fitted isotope model.");
				defaults_.setValidStrings("feature:reported_mz",StringList::create("maximum,average,monoisotopic"));
				defaults_.setSectionDescription("feature","Settings for the features (intensity, quality assessment, ...)");
				//user-specified seed settings
				defaults_.setValue("user-seed:rt_tolerance",5.0, "Allowed RT deviation of seeds from the user-specified seed position.");
				defaults_.setMinFloat("user-seed:rt_tolerance",0.0);
				defaults_.setValue("user-seed:mz_tolerance",1.1, "Allowed m/z deviation of seeds from the user-specified seed position.");
				defaults_.setMinFloat("user-seed:mz_tolerance",0.0);				
				defaults_.setValue("user-seed:min_score",0.5, "Overwrites 'seed:min_score' for user-specified seeds. The cutoff is typically a bit lower in this case.");
				defaults_.setMinFloat("user-seed:min_score",0.0);
				defaults_.setMaxFloat("user-seed:min_score",1.0);
				defaults_.setSectionDescription("user-seed","Settings for user-specified seeds.");
				//debug settings
				defaults_.setValue("debug:pseudo_rt_shift",500.0,"Pseudo RT shift used when .", StringList::create("advanced"));
				defaults_.setMinFloat("debug:pseudo_rt_shift",1.0);
				this->defaultsToParam_();
			}

			// docu in base class
			virtual void setSeeds(const FeatureMapType& seeds)
			{
				seeds_ = seeds;
			}
			
			/// Main method for actual FeatureFinder
			virtual void run()
			{
				//-------------------------------------------------------------------------
				//General initialization
				//---------------------------------------------------------------------------
				
				//quality estimation
				DoubleReal min_feature_score = param_.getValue("feature:min_score");
				//charges to look at				
				SignedSize charge_low = (Int)param_.getValue("isotopic_pattern:charge_low");
				SignedSize charge_high = (Int)param_.getValue("isotopic_pattern:charge_high");
				//fitting settings
				UInt max_iterations = param_.getValue("fit:max_iterations");
				DoubleReal epsilon_abs = param_.getValue("fit:epsilon_abs");
				DoubleReal epsilon_rel = param_.getValue("fit:epsilon_rel");
				
				//copy the input map
				map_ = *(FeatureFinderAlgorithm<PeakType, FeatureType>::map_);
				
				//flag for user-specified seed mode
				bool user_seeds = (seeds_.size()>0);
				if (user_seeds)
				{
					seeds_.sortByMZ();
				}
				DoubleReal user_rt_tol = param_.getValue("user-seed:rt_tolerance");
				DoubleReal user_mz_tol = param_.getValue("user-seed:mz_tolerance");
				DoubleReal user_seed_score = param_.getValue("user-seed:min_score");
				
				//reserve space for calculated scores
				UInt charge_count = charge_high - charge_low + 1;
				for (Size s=0; s<map_.size(); ++s)
				{
					Size scan_size = map_[s].size();
					map_[s].getFloatDataArrays().resize(3 + 2*charge_count);
					map_[s].getFloatDataArrays()[0].setName("trace_score");
					map_[s].getFloatDataArrays()[0].assign(scan_size,0.0);
					map_[s].getFloatDataArrays()[1].setName("intensity_score");
					map_[s].getFloatDataArrays()[1].assign(scan_size,0.0);
					map_[s].getFloatDataArrays()[2].setName("local_max");
					map_[s].getFloatDataArrays()[2].assign(scan_size,0.0);
					//create isotope pattern score arrays
					UInt charge = charge_low;
					for (Size i = 3; i< 3+charge_count; ++i)
					{
						map_[s].getFloatDataArrays()[i].setName(String("pattern_score_")+charge);
						map_[s].getFloatDataArrays()[i].assign(scan_size,0.0);
						++charge;
					}
					//create overall score arrays
					charge = charge_low;
					for (Size i = 3+charge_count; i< 3+2*charge_count; ++i)
					{
						map_[s].getFloatDataArrays()[i].setName(String("overall_score_")+charge);
						map_[s].getFloatDataArrays()[i].assign(scan_size,0.0);
						++charge;
					}
				}

				bool debug = ( (String)(param_.getValue("debug"))=="true" );
				//clean up / create folders for debug information
				if (debug)
				{
					QDir dir(".");
					dir.mkpath("debug/features");
					log_.open("debug/log.txt");
				}
				
				//time
				StopWatch stop_watch;
				stop_watch.start();
				
				//---------------------------------------------------------------------------
				//Step 1:
				//Precalculate intensity scores for peaks
				//---------------------------------------------------------------------------
				log_ << "Precalculating intensity thresholds ..." << std::endl;
				//new scope to make local variables disappear
				{
					ff_->startProgress(0, intensity_bins_*intensity_bins_, "Precalculating intensity scores");
					DoubleReal rt_start = map_.getMinRT();
					DoubleReal mz_start = map_.getMinMZ();
					intensity_rt_step_ = (map_.getMaxRT() - rt_start ) / (DoubleReal)intensity_bins_;
				 	intensity_mz_step_ = (map_.getMaxMZ() - mz_start ) / (DoubleReal)intensity_bins_;
					intensity_thresholds_.resize(intensity_bins_);
					for (Size rt=0; rt<intensity_bins_; ++rt)
					{
						intensity_thresholds_[rt].resize(intensity_bins_);
						DoubleReal min_rt = rt_start + rt * intensity_rt_step_;
						DoubleReal max_rt = rt_start + ( rt + 1 ) * intensity_rt_step_;
						std::vector<DoubleReal> tmp;
						for (Size mz=0; mz<intensity_bins_; ++mz)
						{
							ff_->setProgress(rt*intensity_bins_+ mz);
							DoubleReal min_mz = mz_start + mz * intensity_mz_step_;
							DoubleReal max_mz = mz_start + ( mz + 1 ) * intensity_mz_step_;
							//std::cout << "rt range: " << min_rt << " - " << max_rt << std::endl;
							//std::cout << "mz range: " << min_mz << " - " << max_mz << std::endl;
							tmp.clear();
							for (typename MapType::ConstAreaIterator it = map_.areaBeginConst(min_rt,max_rt,min_mz,max_mz); it!=map_.areaEndConst(); ++it)
							{
								tmp.push_back(it->getIntensity());
							}
							//init vector
							intensity_thresholds_[rt][mz].assign(21, 0.0);
							//store quantiles (20)
							if (tmp.size()!=0)
							{
								std::sort(tmp.begin(), tmp.end());
								for (Size i=0;i<21;++i)
								{
									Size index = (Size)std::floor(0.05*i*(tmp.size()-1));
									intensity_thresholds_[rt][mz][i] = tmp[index];
								}
							}
						}
					}
					
					//store intensity score in PeakInfo
					for (Size s=0; s<map_.size(); ++s)
					{
						for (Size p=0; p<map_[s].size(); ++p)
						{
							map_[s].getFloatDataArrays()[1][p] = intensityScore_(s,p);
						}
					}
					ff_->endProgress();
				}
				
				//---------------------------------------------------------------------------
				//Step 2:
				//Prealculate mass trace scores and local trace maximum for each peak
				//---------------------------------------------------------------------------
				//new scope to make local variables disappear
				{
					ff_->startProgress(0, map_.size(), "Precalculating mass trace scores");
					for (Size s=0; s<map_.size(); ++s)
					{
						ff_->setProgress(s);
						//do nothing for the first few and last few spectra as the scans required to search for traces are missing
						if (s<min_spectra_ || s>=map_.size()-min_spectra_)
						{
							continue;
						}
						
						const SpectrumType& spectrum = map_[s];
						//iterate over all peaks of the scan
						for (Size p=0; p<spectrum.size(); ++p)
						{
							std::vector<DoubleReal> scores;
							scores.reserve(2*min_spectra_);
							
							DoubleReal pos = spectrum[p].getMZ();
							Real inte = spectrum[p].getIntensity();
							
							//log_ << std::endl << "Peak: " << pos << std::endl;
							bool is_max_peak = true; //checking the maximum intensity peaks -> use them later as feature seeds.
							for (Size i=1; i<=min_spectra_; ++i)
							{
								try
								{
									Size spec_index = map_[s+i].findNearest(pos);
									DoubleReal position_score = positionScore_(pos, map_[s+i][spec_index].getMZ(), trace_tolerance_);
									if (position_score >0 && map_[s+i][spec_index].getIntensity()>inte) is_max_peak = false;
									scores.push_back(position_score);
								}
								catch(...) //no peaks in the spectrum
								{
									scores.push_back(0.0);
								}
							}
							for (Size i=1; i<=min_spectra_; ++i)
							{
								try
								{
									Size spec_index = map_[s-i].findNearest(pos);
									DoubleReal position_score = positionScore_(pos, map_[s-i][spec_index].getMZ(), trace_tolerance_);
									if (position_score>0 && map_[s-i][spec_index].getIntensity()>inte) is_max_peak = false;
									scores.push_back(position_score);
								}
								catch(...) //no peaks in the spectrum
								{
									scores.push_back(0.0);
								}
							}
							//Calculate a consensus score out of the scores calculated before
							DoubleReal trace_score = std::accumulate(scores.begin(), scores.end(),0.0) / scores.size();
							
							//store final score for later use
							map_[s].getFloatDataArrays()[0][p] = trace_score;
							map_[s].getFloatDataArrays()[2][p] = is_max_peak;
						}
					}
					ff_->endProgress();
				}

				//---------------------------------------------------------------------------
				//Step 2.5:
				//Prealculate isotope distributions for interesting mass ranges
				//---------------------------------------------------------------------------
				//new scope to make local variables disappear
				{
					DoubleReal max_mass = map_.getMaxMZ()*charge_high;
					Size num_isotopes = std::ceil(max_mass/mass_window_width_);
					ff_->startProgress(0, num_isotopes, "Precalculating isotope distributions");
					
					//reserve enough space
					isotope_distributions_.resize(num_isotopes);

					//calculate distribution if necessary
					for (Size index=0; index<num_isotopes; ++index)
					{
						//log_ << "Calculating iso dist for mass: " << 0.5*mass_window_width_ + index * mass_window_width_ << std::endl;
						IsotopeDistribution d;
						d.setMaxIsotope(20);
						d.estimateFromPeptideWeight(0.5*mass_window_width_ + index * mass_window_width_);
						//trim left and right. And store the number of isotopes on the left, to reconstruct the monoisotopic peak
						Size size_before = d.size();
						d.trimLeft(intensity_percentage_optional_);
						isotope_distributions_[index].trimmed_left = size_before - d.size();
						d.trimRight(intensity_percentage_optional_);
						for (IsotopeDistribution::Iterator it=d.begin(); it!=d.end(); ++it)
						{
							isotope_distributions_[index].intensity.push_back(it->second);
							//log_ << " - " << it->second << std::endl;
						}
						//determine the number of optional peaks at the beginning/end
						Size begin = 0;
						Size end = 0;
						bool is_begin = true;
						bool is_end = false;
						for (Size i=0; i<isotope_distributions_[index].intensity.size(); ++i)
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
						//scale the distibution to a maximum of 1
						DoubleReal max = 0.0;
						for (Size i=0; i<isotope_distributions_[index].intensity.size(); ++i)
						{
							if (isotope_distributions_[index].intensity[i]>max) max = isotope_distributions_[index].intensity[i];
						}
						isotope_distributions_[index].max = max;
						for (Size i=0; i<isotope_distributions_[index].intensity.size(); ++i)
						{
							isotope_distributions_[index].intensity[i] /= max;
						}
						
						//log_ << " - optinal begin/end:" << begin << " / " << end << std::endl;
					}
					
					ff_->endProgress();
				}
				
				//-------------------------------------------------------------------------
				//Step 3:
				//Charge loop (create seeds and features for each charge separately)
				//-------------------------------------------------------------------------
				Int plot_nr_global = -1; //counter for the number of plots (debug info)
#ifdef _OPENMP
#pragma omp parallel for
#endif
				for (SignedSize c=charge_low; c<=charge_high; ++c)
				{
					UInt meta_index_isotope = 3 + c - charge_low;
					UInt meta_index_overall = 3 + charge_count + c - charge_low;
					
					Size feature_candidates = 0;
					std::vector<Seed> seeds;
					
					//-----------------------------------------------------------
					//Step 3.1: Precalculate IsotopePattern score
					//-----------------------------------------------------------
					ff_->startProgress(0, map_.size(), String("Calculating isotope pattern scores for charge ")+String(c));
					for (Size s=0; s<map_.size(); ++s)
					{
						ff_->setProgress(s);
						const SpectrumType& spectrum = map_[s];
						for (Size p=0; p<spectrum.size(); ++p)
						{
							DoubleReal mz = spectrum[p].getMZ();
							
							//get isotope distribution for this mass
							const TheoreticalIsotopePattern& isotopes = getIsotopeDistribution_(mz*c);
							//determine highest peak in isopope distribution
							Size max_isotope = std::max_element(isotopes.intensity.begin(), isotopes.intensity.end()) - isotopes.intensity.begin();
							//Look up expected isotopic peaks (in the current spectrum or adjacent spectra)
							Size peak_index = spectrum.findNearest(mz-((DoubleReal)(isotopes.size()+1)/c));
							IsotopePattern pattern(isotopes.size());
							for (Size i=0; i<isotopes.size(); ++i)
							{
								DoubleReal isotope_pos = mz + ((DoubleReal)i-max_isotope)/c;
								findIsotope_(isotope_pos, s, pattern, i, false, peak_index);
							}
							DoubleReal pattern_score = isotopeScore_(isotopes, pattern, true, false);
							
							//update pattern scores of all contained peaks (if necessary)
							if (pattern_score > 0.0)
							{
								for (Size i=0; i<pattern.peak.size(); ++i)
								{
									if (pattern.peak[i]>=0 && pattern_score>map_[pattern.spectrum[i]].getFloatDataArrays()[meta_index_isotope][pattern.peak[i]])
									{
										map_[pattern.spectrum[i]].getFloatDataArrays()[meta_index_isotope][pattern.peak[i]] = pattern_score;
									}
								}
							}
						}
					}
					ff_->endProgress();
					//-----------------------------------------------------------
					//Step 3.2:
					//Find seeds for this charge
					//-----------------------------------------------------------		
					ff_->startProgress(0, map_.size(), String("Finding seeds for charge ")+String(c));
					DoubleReal min_seed_score = param_.getValue("seed:min_score");
					for (Size s=0; s<map_.size(); ++s)
					{
						ff_->setProgress(s);
						//do nothing for the first few and last few spectra as the scans required to search for traces are missing
						if (s<min_spectra_ || s>=map_.size()-min_spectra_)
						{
							continue;
						}
						//iterate over peaks
						for (Size p=0; p<map_[s].size(); ++p)
						{	
							FloatDataArrays& meta = map_[s].getFloatDataArrays();
							DoubleReal overall_score = std::pow(meta[0][p]*meta[1][p]*meta[meta_index_isotope][p], 1.0f/3.0f);
							meta[meta_index_overall][p] = overall_score;
							//add seed to vector if certain conditions are fullfilled
							if (meta[2][p]!=0.0) // local maximum of mass trace is prerequisite for all features
							{
								//automatic seeds: overall score greater than the min seed score
								if (!user_seeds && overall_score>=min_seed_score)
								{
									Seed seed;
									seed.spectrum = s;
									seed.peak = p;
									seed.intensity = map_[s][p].getIntensity();								
									seeds.push_back(seed);
								}
								//user-specified seeds: overall score greater than USER min seed score
								else if (user_seeds && overall_score>=user_seed_score)
								{
									//only consider seeds, if they are near a user-specified seed
									FeatureType tmp;
									tmp.setMZ(map_[s][p].getMZ() - user_mz_tol);
									for (typename FeatureMapType::const_iterator it = std::lower_bound(seeds_.begin(), seeds_.end(), tmp, typename FeatureType::MZLess()); it<seeds_.end(); ++it)
									{
										if (it->getMZ()>map_[s][p].getMZ() + user_mz_tol) break;
										if (fabs(it->getMZ()-map_[s][p].getMZ())<user_mz_tol && fabs(it->getRT()-map_[s].getRT())<user_rt_tol)
										{
											Seed seed;
											seed.spectrum = s;
											seed.peak = p;
											seed.intensity = map_[s][p].getIntensity();								
											seeds.push_back(seed);
											break;
										}
									}
								} 
							}
						}
					}
					//sort seeds according to intensity
					std::sort(seeds.rbegin(),seeds.rend());
					//create and store seeds map and selected peak map
					if (debug)
					{
						//seeds
						FeatureMap<> seed_map;
						seed_map.reserve(seeds.size());
						for (Size i=0; i<seeds.size(); ++i)
						{
							Size spectrum = seeds[i].spectrum;
							Size peak = seeds[i].peak;
							const FloatDataArrays& meta = map_[spectrum].getFloatDataArrays();
							Feature tmp;
							tmp.setIntensity(seeds[i].intensity);
							tmp.setOverallQuality(meta[meta_index_overall][peak]);
							tmp.setRT(map_[spectrum].getRT());
							tmp.setMZ(map_[spectrum][peak].getMZ());
							tmp.setMetaValue("intensity_score",	meta[1][peak]);
							tmp.setMetaValue("pattern_score",	meta[meta_index_isotope][peak]);
							tmp.setMetaValue("trace_score",	meta[0][peak]);
							seed_map.push_back(tmp);
						}
						FeatureXMLFile().store(String("debug/seeds_")+String(c)+".featureXML", seed_map);
					}
					ff_->endProgress();
					std::cout << "Found " << seeds.size() << " seeds for charge " << c << "." << std::endl;
					
					//------------------------------------------------------------------
					//Step 3.3:
					//Extension of seeds
					//------------------------------------------------------------------
					ff_->startProgress(0,seeds.size(), String("Extending seeds for charge ")+String(c));
					for (Size i=0; i<seeds.size(); ++i)
					{
						//------------------------------------------------------------------
						//Step 3.3.1:
						//Extend all mass traces
						//------------------------------------------------------------------
						ff_->setProgress(i);
						log_ << std::endl << "Seed " << i << ":" << std::endl;
						//If the intensity is zero this seed is already uses in another feature
						const SpectrumType& spectrum = map_[seeds[i].spectrum];
						const PeakType& peak = spectrum[seeds[i].peak];
						log_ << " - Int: " << peak.getIntensity() << std::endl;
						log_ << " - RT: " << spectrum.getRT() << std::endl;
						log_ << " - MZ: " << peak.getMZ() << std::endl;
						if (seeds[i].intensity == 0.0)
						{
							log_ << "Already used in another feature => aborting!" << std::endl;
							continue;
						}
						
						//----------------------------------------------------------------
						//Find best fitting isotope pattern for this charge (using averagine)
						IsotopePattern best_pattern(0);
						DoubleReal isotope_fit_quality = findBestIsotopeFit_(seeds[i], c, best_pattern);
						if (isotope_fit_quality<min_isotope_fit_)
						{
							abort_(debug, seeds[i], "Could not find good enough isotope pattern containing the seed");
							continue;
						}
						
						//extend the convex hull in RT dimension (starting from the trace peaks)
						log_ << "Collecting mass traces" << std::endl;
						MassTraces traces;
						traces.reserve(best_pattern.peak.size());
						extendMassTraces_(best_pattern, traces, meta_index_overall);

						//check if the traces are still valid
						DoubleReal seed_mz = map_[seeds[i].spectrum][seeds[i].peak].getMZ();
						if (!traces.isValid(seed_mz, trace_tolerance_))
						{
							abort_(debug, seeds[i], "Could not extend seed");
							continue;
						}
						
						//------------------------------------------------------------------
						//Step 3.3.2:
						//Gauss fit (first fit to find the feature boundaries)
						//------------------------------------------------------------------
						Int plot_nr=-1;
#ifdef _OPENMP
#pragma omp critical (FeatureFinderAlgorithmPicked)
#endif
						plot_nr = ++plot_nr_global;
						log_ << "Fitting model" << std::endl;
					  const gsl_multifit_fdfsolver_type *T;
					  gsl_multifit_fdfsolver *s;
					  int status;
					  const size_t param_count = 3;
						const size_t data_count = traces.getPeakCount();
					  gsl_multifit_function_fdf func;

						//TODO try fit with baseline term once more
						//baseline estimate
						traces.updateBaseline();
						traces.baseline = 0.75 * traces.baseline;

					  //parameter estimates (height, x0, sigma)
						traces[traces.max_trace].updateMaximum();
						DoubleReal height = traces[traces.max_trace].max_peak->getIntensity() - traces.baseline;
						DoubleReal x0 = traces[traces.max_trace].max_rt;
						const DoubleReal region_rt_span = traces[traces.max_trace].peaks.back().first-traces[traces.max_trace].peaks[0].first;
						DoubleReal sigma = region_rt_span/20.0;
					  double x_init[param_count] = {height, x0, sigma};
						log_ << " - estimates - height: " << height << " x0: " << x0 <<  " sigma: " << sigma  << std::endl;
						
						//fit					  
					  gsl_vector_view x = gsl_vector_view_array(x_init, param_count);	
					  const gsl_rng_type * type;
					  gsl_rng * r;
					  gsl_rng_env_setup();
					  type = gsl_rng_default;
					  r = gsl_rng_alloc(type);
					  func.f = &gaussF_;
					  func.df = &gaussDF_;
					  func.fdf = &gaussFDF_;
					  func.n = data_count;
					  func.p = param_count;
					  func.params = &traces;
					  T = gsl_multifit_fdfsolver_lmsder;
					  s = gsl_multifit_fdfsolver_alloc(T, data_count, param_count);
					  gsl_multifit_fdfsolver_set(s, &func, &x.vector);
					  size_t iter = 0;					
					  do
					  {
					    iter++;
					    status = gsl_multifit_fdfsolver_iterate(s);
					    //log_ << "iter " << iter << ": " << gsl_vector_get(s->x, 0) << " " << gsl_vector_get(s->x, 1) << " " << gsl_vector_get(s->x, 2) << std::endl;
					    if (status) break;
					    status = gsl_multifit_test_delta(s->dx, s->x, epsilon_abs, epsilon_rel);
					  } 
					  while (status == GSL_CONTINUE && iter < max_iterations);
						height = gsl_vector_get(s->x, 0);
						x0 = gsl_vector_get(s->x, 1);
						sigma = std::fabs(gsl_vector_get(s->x, 2));						
						gsl_multifit_fdfsolver_free(s);
						log_ << " - fit - height: " << height  << " x0: " << x0 << " sigma: " << sigma << std::endl;
						
						//------------------------------------------------------------------
						//Step 3.3.3:
						//Crop feature according to RT fit (2.5*sigma) and remove badly fitting traces
						//------------------------------------------------------------------
						MassTraces new_traces;
						DoubleReal low_bound = x0 - 2.5 * sigma;
						DoubleReal high_bound = x0 + 2.5 * sigma;
						log_ << "    => RT bounds: " << low_bound << " - " << high_bound << std::endl;
						for (Size t=0; t< traces.size(); ++t)
						{
							MassTrace& trace = traces[t];
							log_ << "   - Trace " << t << ": (" << trace.theoretical_int << ")" << std::endl;

							MassTrace new_trace;
							//compute average relative deviation and correlation
							DoubleReal deviation = 0.0;
							std::vector<DoubleReal> v_theo, v_real;
							for (Size k=0; k<trace.peaks.size(); ++k)
							{
								//consider peaks when inside RT bounds only
								if (trace.peaks[k].first>=low_bound && trace.peaks[k].first<= high_bound)
								{
									new_trace.peaks.push_back(trace.peaks[k]);
									DoubleReal theo = traces.baseline + trace.theoretical_int *  height * exp(-0.5 * pow(trace.peaks[k].first - x0, 2) / pow(sigma, 2) );
									v_theo.push_back(theo);
									DoubleReal real = trace.peaks[k].second->getIntensity();
									v_real.push_back(real);
									deviation += std::fabs(real-theo)/theo;
								}
							}
							DoubleReal fit_score = 0.0;
							DoubleReal correlation = 0.0;							
							DoubleReal final_score = 0.0;
							if (new_trace.peaks.size()!=0)
							{
								fit_score = deviation / new_trace.peaks.size();
								correlation = std::max(0.0, Math::pearsonCorrelationCoefficient(v_theo.begin(),v_theo.end(),v_real.begin(), v_real.end()));
								final_score = std::sqrt(correlation * std::max(0.0, 1.0-fit_score));
							}
							log_ << "     - peaks: " << new_trace.peaks.size() << " / " << trace.peaks.size() << " - relative deviation: " << fit_score << " - correlation: " << correlation << " - final score: " << correlation << std::endl;
							//remove badly fitting traces
							if ( !new_trace.isValid() || final_score < min_trace_score_)
							{
								if (t<traces.max_trace)
								{
									new_traces = MassTraces();
									log_ << "     - removed this and previous traces due to bad fit" << std::endl;
									new_traces.clear(); //remove earlier traces
									continue;
								}
								else if (t==traces.max_trace)
								{
									new_traces = MassTraces();
									log_ << "     - aborting (max trace was removed)" << std::endl;
									break;
								}
								else if (t>traces.max_trace)
								{
									log_ << "     - removed due to bad fit => omitting the rest" << std::endl;
		            	break; //no more traces are possible
								}
							}
							//add new trace
							else
							{
								new_trace.theoretical_int = trace.theoretical_int;
								new_traces.push_back(new_trace);
								if (t==traces.max_trace)
								{
									new_traces.max_trace = new_traces.size()-1;
								}
							}
						}
						new_traces.baseline = traces.baseline;			

						//------------------------------------------------------------------
						//Step 3.3.4:
						//Check if feature is ok
						//------------------------------------------------------------------
					  bool feature_ok = true;
					  String error_msg = "";
						//check if the sigma fit was ok (if it is larger than 'max_rt_span')
						if (feature_ok)
						{
							if (5.0*sigma > max_rt_span_*region_rt_span )
							{
						  	feature_ok = false;
						  	error_msg = "Invalid fit: Fitted model is bigger than 'max_rt_span'";
							}
						}
					  //check if the feature is valid
					  if (!new_traces.isValid(seed_mz, trace_tolerance_))
					  {
					  	feature_ok = false;
					  	error_msg = "Invalid feature after fit - too few traces or peaks left";
					  }
						//check if x0 is inside feature bounds
						if (feature_ok)
						{
							std::pair<DoubleReal,DoubleReal> rt_bounds = new_traces.getRTBounds();
							if (x0<rt_bounds.first || x0>rt_bounds.second)
							{
						  	feature_ok = false;
						  	error_msg = "Invalid fit: Center outside of feature bounds";
							}
						}
						//check if the remaining traces fill out at least 'min_rt_span' of the RT span
						if (feature_ok)
						{
							std::pair<DoubleReal,DoubleReal> rt_bounds = new_traces.getRTBounds();
							if ((rt_bounds.second-rt_bounds.first) < min_rt_span_*5.0*sigma )
							{
						  	feature_ok = false;
						  	error_msg = "Invalid fit: Less than 'min_rt_span' left after fit";
							}
						}
					  //check if feature quality is high enough (average relative deviation and correlation of the whole feature)
						DoubleReal fit_score = 0.0;
						DoubleReal correlation = 0.0;
						DoubleReal final_score = 0.0;

						if(feature_ok)
						{
							std::vector<DoubleReal> v_theo, v_real;
							DoubleReal deviation = 0.0;
							for (Size t=0; t< new_traces.size(); ++t)
							{
								MassTrace& trace = new_traces[t];
								for (Size k=0; k<trace.peaks.size(); ++k)
								{
									DoubleReal theo = new_traces.baseline + trace.theoretical_int *  height * exp(-0.5 * pow(trace.peaks[k].first - x0, 2) / pow(sigma, 2) );
									v_theo.push_back(theo);
									DoubleReal real = trace.peaks[k].second->getIntensity();
									v_real.push_back(real);
									deviation += std::fabs(real-theo)/theo;
								}
							}
							fit_score = std::max(0.0, 1.0 - (deviation / new_traces.getPeakCount()));
							correlation = std::max(0.0,Math::pearsonCorrelationCoefficient(v_theo.begin(),v_theo.end(),v_real.begin(), v_real.end()));
							final_score = std::sqrt(correlation * fit_score);
						  if (final_score<min_feature_score)
						  {
						  	feature_ok = false;
						  	error_msg = "Feature quality too low after fit";
						  }
							//quality output
							log_ << "Quality estimation:" << std::endl;
							log_ << " - relative deviation: " << fit_score << std::endl;
							log_ << " - correlation: " << correlation << std::endl;
							log_ << " => final score: " << final_score << std::endl;
						}
					  				  
						//write debug output of feature
						if (debug)
						{
							DoubleReal pseudo_rt_shift = param_.getValue("debug:pseudo_rt_shift");
							TextFile tf;
							//gnuplot script	
							String script = String("plot \"debug/features/") + plot_nr + ".dta\" title 'before fit (RT: " +  String::number(x0,2) + " m/z: " +  String::number(peak.getMZ(),4) + ")' with points 1";
							//feature before fit
							for (Size k=0; k<traces.size(); ++k)
							{
								for (Size j=0; j<traces[k].peaks.size(); ++j)
								{
									tf.push_back(String(pseudo_rt_shift*k+traces[k].peaks[j].first) + "	" + traces[k].peaks[j].second->getIntensity());
								}
							}
							tf.store(String("debug/features/") + plot_nr + ".dta");
							//fitted feature
							if (new_traces.getPeakCount()!=0)
							{
								tf.clear();
								for (Size k=0; k<new_traces.size(); ++k)
								{
									for (Size j=0; j<new_traces[k].peaks.size(); ++j)
									{
										tf.push_back(String(pseudo_rt_shift*k+new_traces[k].peaks[j].first) + "	" + new_traces[k].peaks[j].second->getIntensity());
									}
								}
								tf.store(String("debug/features/") + plot_nr + "_cropped.dta");
								script = script + ", \"debug/features/" + plot_nr + "_cropped.dta\" title 'feature ";
								if (!feature_ok)
								{
									script = script + " - " + error_msg;
								}
								else
								{
									script = script + (features_->size()+1) + " (score: " +  String::number(final_score,3) + ")";
								}
								script = script + "' with points 3";
							}
							//fitted functions
							tf.clear();
							for (Size k=0; k<traces.size(); ++k)
							{
								char fun = 'f';
								fun += (char)k;
								tf.push_back(String(fun)+"(x)= " + traces.baseline + " + " + (traces[k].theoretical_int*height) + " * exp(-0.5*(x-" + (pseudo_rt_shift*k+x0) + ")**2/(" + sigma + ")**2)");
								script =  script + ", " + fun + "(x) title 'Trace " + k + " (m/z: " + String::number(traces[k].getAvgMZ(),4) + ")'";
							}
							//output
							tf.push_back("set xlabel \"pseudo RT (mass traces side-by-side)\"");
							tf.push_back("set ylabel \"intensity\"");
							tf.push_back("set samples 1000");
							tf.push_back(script);
							tf.push_back("pause -1");
							tf.store(String("debug/features/") + plot_nr + ".plot");
						}
						traces = new_traces;
						
						log_ << "Feature label: " << plot_nr << std::endl;
						
						//validity output
						if (!feature_ok)
						{
							abort_(debug, seeds[i], error_msg);
							continue;
						}
						
						//------------------------------------------------------------------
						//Step 3.3.5:
						//Feature creation
						//------------------------------------------------------------------
						Feature f;
						//set label
						f.setMetaValue(3,plot_nr);
						f.setCharge(c);
						f.setOverallQuality(final_score);
						if (debug)
						{
							f.setMetaValue("score_fit",fit_score);
							f.setMetaValue("score_correlation",correlation);
						}
						f.setRT(x0);
						
						//Calculate the mass of the feature: maximum, average, monoisotopic
						String reported_mz = param_.getValue("feature:reported_mz");
						if(reported_mz=="maximum")
						{
							f.setMZ(traces[traces.getTheoreticalMax()].getAvgMZ());
						}
						else if(reported_mz=="average")
						{
							DoubleReal total_intensity = 0.0;
							DoubleReal average_mz = 0.0;
	 						for (Size t=0; t<traces.size(); ++t)
							{
								for (Size p=0; p<traces[t].peaks.size(); ++p)
								{
									average_mz += traces[t].peaks[p].second->getMZ()*traces[t].peaks[p].second->getIntensity();
									total_intensity+=traces[t].peaks[p].second->getIntensity();
								}
							}
							average_mz /= total_intensity;
							f.setMZ(average_mz);
						}
						else if(reported_mz=="monoisotopic")
						{
							DoubleReal mono_mz = traces[traces.getTheoreticalMax()].getAvgMZ();
							mono_mz -= (1.005/c) * (traces.getTheoreticalMax() + best_pattern.theoretical_pattern.trimmed_left);
							f.setMZ(mono_mz);
						}
						
						//Calculate intensity based on model only
						// - the model does not include the baseline, so we ignore it here
						// - as we scaled the isotope distribution to 
						f.setIntensity(2.5 * height * sigma / getIsotopeDistribution_(f.getMZ()).max);
						//add convex hulls of mass traces
						for (Size j=0; j<traces.size(); ++j)
						{
							f.getConvexHulls().push_back(traces[j].getConvexhull());
						}
#ifdef _OPENMP
#pragma omp critical (FeatureFinderAlgorithmPicked)
#endif
						features_->push_back(f);
						feature_candidates++;
						
						//----------------------------------------------------------------
						//Remove all seeds that lie inside the convex hull of the new feature
						DBoundingBox<2> bb = f.getConvexHull().getBoundingBox();
						for (Size j=i+1; j<seeds.size(); ++j)
						{
							DoubleReal rt = map_[seeds[j].spectrum].getRT();
							DoubleReal mz = map_[seeds[j].spectrum][seeds[j].peak].getMZ();
							if (bb.encloses(rt,mz) && f.encloses(rt,mz))
							{
								//set intensity to zero => the peak will be skipped!
								seeds[j].intensity = 0.0;
							}
						}
					}
					ff_->endProgress();
					std::cout << "Found " << feature_candidates << " features candidates for charge " << c << "." << std::endl;
				}
					
				//------------------------------------------------------------------
				//Step 4:
				//Resolve contradicting and overlapping features
				//------------------------------------------------------------------
				ff_->startProgress(0, features_->size()*features_->size(), "Resolving overlapping features");
				log_ << "Resolving intersecting features (" << features_->size() << " candidates)" << std::endl;
				//sort features according to m/z in order to speed up the resolution
				features_->sortByMZ();
				//precalculate BBs and maximum mz span
				std::vector< DBoundingBox<2> > bbs(features_->size());
				DoubleReal max_mz_span = 0.0;
				for (Size i=0; i<features_->size(); ++i)
				{
					bbs[i] = features_->at(i).getConvexHull().getBoundingBox();
					if (bbs[i].height()>max_mz_span)
					{
						max_mz_span = bbs[i].height();
					}
				}
				//intersect
				for (Size i=0; i<features_->size(); ++i)
				{
					Feature& f1(features_->at(i));
					for (Size j=i+1; j<features_->size(); ++j)
					{
						ff_->setProgress(i*features_->size()+j);
						Feature& f2(features_->at(j));
						//features that are more than 2 times the maximum m/z span apart do not overlap => abort 
						if (f2.getMZ()-f1.getMZ()>2.0*max_mz_span) break;
						//do nothting if one of the features is alreay removed
						if (f1.getIntensity()==0.0 || f2.getIntensity()==0.0) continue;
						//do nothing if the overall convex hulls du not overlap
						if (!bbs[i].intersects(bbs[j])) continue;
						//act depending on the intersection
						DoubleReal intersection = intersection_(f1, f2);
						if (intersection>=max_feature_intersection_)
						{
							log_ << " - Intersection (" << (i+1) << "/" << (j+1) << "): " << intersection << std::endl;
							if (f1.getCharge()==f2.getCharge())
							{
								if (f1.getIntensity()*f1.getOverallQuality()>f2.getIntensity()*f2.getOverallQuality())
								{
									log_ << "   - same charge -> removing duplicate " << (j+1) << std::endl;
									f1.getSubordinates().push_back(f2);
									f2.setIntensity(0.0);
								}
								else
								{
									log_ << "   - same charge -> removing duplicate " << (i+1) << std::endl;
									f2.getSubordinates().push_back(f1);
									f1.setIntensity(0.0);
								}
							}
							else if (f2.getCharge()%f1.getCharge()==0)
							{
								log_ << "   - different charge (one is the multiple of the other) -> removing lower charge " << (i+1) << std::endl;
								f2.getSubordinates().push_back(f1);
								f1.setIntensity(0.0);
							}
							else if (f1.getCharge()%f2.getCharge()==0)
							{
								log_ << "   - different charge (one is the multiple of the other) -> removing lower charge " << (i+1) << std::endl;
								f1.getSubordinates().push_back(f2);
								f2.setIntensity(0.0);
							}
							else
							{
								if (f1.getOverallQuality()>f2.getOverallQuality())
								{
									log_ << "   - different charge -> removing lower score " << (j+1) << std::endl;
									f1.getSubordinates().push_back(f2);
									f2.setIntensity(0.0);
								}
								else
								{
									log_ << "   - different charge -> removing lower score " << (i+1) << std::endl;
									f2.getSubordinates().push_back(f1);
									f1.setIntensity(0.0);
								}
							}
						}
					}
				}
				//finally remove features with intensity 0
				FeatureMap<> tmp;
				tmp.reserve(features_->size());
				UInt removed = 0;
				for (Size i=0; i<features_->size(); ++i)
				{
					if (features_->operator[](i).getIntensity()!=0.0)
					{
						tmp.push_back(features_->operator[](i));
					}
					else
					{
						++removed;
					}
				}
				tmp.Base::swap(*features_);
				//sort features by intensty
				features_->sortByIntensity(true);
				ff_->endProgress();
				std::cout << "Removed " << removed << " overlapping features." << std::endl;
				std::cout << features_->size() << " features left." << std::endl;
				
				//Abort reasons 
				std::cout << std::endl;
				std::cout << "Abort reasons during feature construction:" << std::endl;
				for (std::map<String,UInt>::const_iterator it=aborts_.begin(); it!=aborts_.end(); ++it)
				{
					std::cout << "- " << it->first << ": " << it->second << std::endl;
				}
				if (debug)
				{
					//store map of abort reasons for failed seeds
					FeatureMap<> abort_map;
					abort_map.reserve(abort_reasons_.size());
					for (typename std::map<Seed, String>::iterator it2=abort_reasons_.begin(); it2!=abort_reasons_.end(); ++it2)
					{
						Feature f;
						f.setRT(map_[it2->first.spectrum].getRT());
						f.setMZ(map_[it2->first.spectrum][it2->first.peak].getMZ());
						f.setIntensity(map_[it2->first.spectrum][it2->first.peak].getIntensity());
						f.setMetaValue("abort_reason",it2->second);
						abort_map.push_back(f);
					}
					FeatureXMLFile().store("debug/abort_reasons.featureXML", abort_map);
					
					//store input map with calculated scores (without overall score)
					for (Size s=0; s<map_.size(); ++s)
					{
						map_[s].getFloatDataArrays().erase(map_[s].getFloatDataArrays().begin()+2);
					}					
					MzMLFile().store("debug/input.mzML", map_);
				}

				//Execution time
				stop_watch.stop();
				std::cout << std::endl;
				std::cout << "Execution time: " << stop_watch.getCPUTime() << " s" << std::endl;
			}
			
			static FeatureFinderAlgorithm<PeakType,FeatureType>* create()
			{
				return new FeatureFinderAlgorithmPicked();
			}

			static const String getProductName()
			{
				return "centroided";
			}
	
		protected:
			/// editable copy of the map
			MapType map_;
			/// Output stream for log/debug info
			mutable std::ofstream log_; 
			/// Array of abort reasons
			std::map<String, UInt> aborts_;
			/// Array of abort reasons
			std::map<Seed, String> abort_reasons_;
			/// User-specified seed list
			FeatureMapType seeds_;
		
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
			DoubleReal min_trace_score_; ///< Minimum quality of a traces
			DoubleReal min_rt_span_; ///< Minimum RT range that has to be left after the fit
			DoubleReal max_rt_span_; ///< Maximum RT range the model is allowed to span
			DoubleReal max_feature_intersection_; ///< Maximum allowed feature intersection (if larger, that one of the feature is removed)
			//@}

			///@name Members for intensity significance estimation
			//@{			
			/// RT bin width
			DoubleReal intensity_rt_step_;
			/// m/z bin width
			DoubleReal intensity_mz_step_;
			/// Precalculated intensity 20-quantiles (binned)
			std::vector< std::vector< std::vector<DoubleReal> > > intensity_thresholds_;
			//@}

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
				min_trace_score_ = param_.getValue("feature:min_trace_score");
				min_rt_span_ = param_.getValue("feature:min_rt_span");
				max_rt_span_ = param_.getValue("feature:max_rt_span");
				max_feature_intersection_ = param_.getValue("feature:max_intersection");
			}
			
			///Writes the abort reason to the log file and counts occurences for each reason
			void abort_(bool debug, const Seed& seed, const String& reason)
			{
				log_ << "Abort: " << reason << std::endl;
				aborts_[reason]++;
				if (debug) abort_reasons_[seed] = reason;
			}

			///Calculates the intersection between features.
			///The value is normalized by the size of the smaller feature, so it rages from 0 to 1.
			DoubleReal intersection_(const Feature& f1, const Feature& f2) const
			{
				//calculate the RT range sum of feature 1
				DoubleReal s1 = 0.0;
				const std::vector<ConvexHull2D>& hulls1 = f1.getConvexHulls();
				for (Size i=0; i<hulls1.size(); ++i)
				{
					s1 += hulls1[i].getBoundingBox().width();
				}
				
				//calculate the RT range sum of feature 2
				DoubleReal s2 = 0.0;
				const std::vector<ConvexHull2D>& hulls2 = f2.getConvexHulls();
				for (Size j=0; j<hulls2.size(); ++j)
				{
					s2 += hulls2[j].getBoundingBox().width();
				}
				
				//calculate overlap
				DoubleReal overlap = 0.0;
				for (Size i=0; i<hulls1.size(); ++i)
				{
					DBoundingBox<2> bb1 = hulls1[i].getBoundingBox();
					for (Size j=0; j<hulls2.size(); ++j)
					{
						DBoundingBox<2> bb2 = hulls2[j].getBoundingBox();
						if (bb1.intersects(bb2))
						{
							if (bb1.min()[0]<=bb2.min()[0] && bb1.max()[0]>=bb2.max()[0]) //bb1 contains bb2
							{
								overlap += bb2.width();
							}
							else if (bb2.min()[0]<=bb1.min()[0] && bb2.max()[0]>=bb1.max()[0]) //bb2 contains bb1
							{
								overlap += bb1.width();
							}
							else if (bb1.min()[0]<=bb2.min()[0] && bb1.max()[0]<=bb2.max()[0]) //the end of bb1 overlaps with bb2
							{
								overlap += bb1.max()[0]-bb2.min()[0];
							}
							else if (bb2.min()[0]<=bb1.min()[0] && bb2.max()[0]<=bb1.max()[0]) //the end of bb2 overlaps with bb1
							{
								overlap += bb2.max()[0]-bb1.min()[0];
							}
						}
					}
				}
				
				return overlap/std::min(s1,s2);
			}
			
			///Returns the isotope distribution for a certain mass window
			const TheoreticalIsotopePattern& getIsotopeDistribution_(DoubleReal mass) const
			{
				//calculate index in the vector
				Size index = (Size)std::floor(mass/mass_window_width_);
				
				if (index>=isotope_distributions_.size()) std::cout << "INDEX: " << index << " SIZE: " << isotope_distributions_.size() << std::endl;
				
				//Return distribution
				return isotope_distributions_[index];
			}
						
			/**
				@brief Finds the best fitting position of the isotopic pattern estimate defined by @p center
				
				@param center the maximum peak of the isotope distribution (contains charge as well)
				@param charge The charge of the pattern 
				@param best_pattern Returns the indices of the isotopic peaks. If a isopopic peak is missing -1 is returned.
			*/
			DoubleReal findBestIsotopeFit_(const Seed& center, UInt charge, IsotopePattern& best_pattern) const
			{
				log_ << "Testing isotope patterns for charge " << charge << ": " << std::endl;			
				const SpectrumType& spectrum = map_[center.spectrum];
				const TheoreticalIsotopePattern& isotopes = getIsotopeDistribution_(spectrum[center.peak].getMZ()*charge);	
				log_ << " - Seed: " << center.peak << " (mz:" << spectrum[center.peak].getMZ()<< ")" << std::endl;
				
				//Find m/z boundaries of search space (linear search as this is local and we have the center already)
				DoubleReal mass_window = (DoubleReal)(isotopes.size()+1) / (DoubleReal)charge;
				log_ << " - Mass window: " << mass_window << std::endl;
				Size end = center.peak;
				while(end<spectrum.size() && spectrum[end].getMZ()<spectrum[center.peak].getMZ()+mass_window)
				{
					++end;
				}
				--end;
				//search begin
				SignedSize begin = center.peak;
				while(begin>=0 && spectrum[begin].getMZ()>spectrum[center.peak].getMZ()-mass_window)
				{
					--begin;
				}
				++begin;
				log_ << " - Begin: " << begin << " (mz:" << spectrum[begin].getMZ()<< ")" << std::endl;
				log_ << " - End: " << end << " (mz:" << spectrum[end].getMZ()<< ")" << std::endl;

				//fit isotope distribution to peaks
				DoubleReal max_score = 0.0;
				for (Size start=begin; start<=end; ++start)
				{
					//find isotope peaks for the current start peak
					Size peak_index = start;
					IsotopePattern pattern(isotopes.size());
					log_ << " - Fitting at " << start << " (mz:" << spectrum[start].getMZ() << ")" << std::endl;
					for (Size iso=0; iso<isotopes.size(); ++iso)
					{
						DoubleReal pos = spectrum[start].getMZ() + iso/(DoubleReal)charge;
						findIsotope_(pos, center.spectrum, pattern, iso, true, peak_index);
					}
					
					//check if the seed is contained, otherwise abort
					bool seed_contained = false;
					for (Size iso=0; iso<pattern.peak.size(); ++iso)
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
					for (Size iso=0; iso<pattern.peak.size(); ++iso)
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
					if (score> max_score)
					{
						max_score = score;
						best_pattern = pattern;
					}
				}
				log_ << " - best score              : " << max_score << std::endl;
				best_pattern.theoretical_pattern = isotopes;
				return max_score;
			}
			
			///Extends all mass traces of a isotope pattern in one step
			void extendMassTraces_(const IsotopePattern& pattern, MassTraces& traces, Size meta_index_overall) const
			{
				//find index of the trace with the maximum intensity
				DoubleReal max_int =  0.0;
				Size max_trace_index = 0;
				for (Size p=0; p<pattern.peak.size(); ++p)
				{
					if (pattern.peak[p]<0) continue; //skip missing and removed traces
					if (map_[pattern.spectrum[p]][pattern.peak[p]].getIntensity()>max_int)
					{
						max_int = map_[pattern.spectrum[p]][pattern.peak[p]].getIntensity();
						max_trace_index = p;
					}
				}
				
				//extend the maximum intensity trace to determine the boundaries in RT dimension
				Size start_index = pattern.spectrum[max_trace_index];
				const PeakType* start_peak = &(map_[pattern.spectrum[max_trace_index]][pattern.peak[max_trace_index]]);
				DoubleReal start_mz = start_peak->getMZ();
				DoubleReal start_rt = map_[start_index].getRT();
				log_ << " - Trace " << max_trace_index << " (maximum intensity)" << std::endl;
				log_ << "   - extending from: " << map_[start_index].getRT() << " / " << start_mz << " (int: " << start_peak->getIntensity() << ")" << std::endl;
				//initialize the trace and extend
				MassTrace max_trace;
				max_trace.peaks.push_back(std::make_pair(start_rt,start_peak));
				extendMassTrace_(max_trace, start_index, start_mz, false, meta_index_overall);
				extendMassTrace_(max_trace, start_index, start_mz, true, meta_index_overall);

				DoubleReal rt_max = max_trace.peaks.back().first;
				DoubleReal rt_min = max_trace.peaks.begin()->first;
				log_ << "   - rt bounds: " << rt_min << "-" << rt_max << std::endl;
				//Abort if too few peak were found
				if (!max_trace.isValid() || max_trace.peaks.size()<2*min_spectra_-max_missing_trace_peaks_)
				{
					log_ << "   - could not extend trace with maximum intensity => abort" << std::endl;
					return;
				}
				for (Size p=0; p<pattern.peak.size(); ++p)
				{
					log_ << " - Trace " << p << std::endl;
					if (p==max_trace_index)
					{
						log_ << "   - previously extended maximum trace" << std::endl;
						traces.push_back(max_trace);
						traces.back().theoretical_int = pattern.theoretical_pattern.intensity[p];
						traces.max_trace = traces.size()-1;
						continue;
					}
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
					starting_peak.intensity = map_[starting_peak.spectrum][starting_peak.peak].getIntensity();
					log_ << "   - trace seed: " << map_[starting_peak.spectrum].getRT() << " / " << map_[starting_peak.spectrum][starting_peak.peak].getMZ() << " (int: " << map_[starting_peak.spectrum][starting_peak.peak].getIntensity() << ")" << std::endl;
					
					//search for nearby maximum of the mass trace as the extension assumes that it starts at the maximum
					Size begin = std::max((Size)0,starting_peak.spectrum-min_spectra_);
					Size end = std::min(starting_peak.spectrum+min_spectra_,(Size)map_.size());
					DoubleReal mz = map_[starting_peak.spectrum][starting_peak.peak].getMZ();
					DoubleReal inte = map_[starting_peak.spectrum][starting_peak.peak].getIntensity();
					for (Size spectrum_index=begin; spectrum_index<end; ++spectrum_index)
					{
						//find better seeds (no-empty scan/low mz diff/higher intensity)
						SignedSize peak_index = -1;
						try
						{
							peak_index = map_[spectrum_index].findNearest(map_[starting_peak.spectrum][starting_peak.peak].getMZ());
						}
						catch(...) //no peaks in the spectrum
						{
							peak_index=-1;
						}
						if (peak_index<0 ||
								map_[spectrum_index][peak_index].getIntensity()<=inte ||
								std::fabs(mz-map_[spectrum_index][peak_index].getMZ())>=pattern_tolerance_
							 ) continue;
						starting_peak.spectrum = spectrum_index;
						starting_peak.peak = peak_index;
						inte = map_[spectrum_index][peak_index].getIntensity();
					}
					log_ << "   - extending from: " << map_[starting_peak.spectrum].getRT() << " / " << map_[starting_peak.spectrum][starting_peak.peak].getMZ() << " (int: " << map_[starting_peak.spectrum][starting_peak.peak].getIntensity() << ")" << std::endl;
					
					//------------------------------------------------------------------
					//Extend seed to a mass trace
					MassTrace trace;
					const PeakType* seed = &(map_[starting_peak.spectrum][starting_peak.peak]);
					//initialize trace with seed data and extend
					trace.peaks.push_back(std::make_pair(map_[starting_peak.spectrum].getRT(),seed));
					extendMassTrace_(trace, starting_peak.spectrum, seed->getMZ(), false, meta_index_overall, rt_min, rt_max);
					extendMassTrace_(trace, starting_peak.spectrum, seed->getMZ(), true, meta_index_overall, rt_min, rt_max);
					
					//check if enough peaks were found
					if (!trace.isValid())
					{
						log_ << "   - could not extend trace " << std::endl;
						//Missing traces in the middle of a pattern are not acceptable => fix this
						if (p<traces.max_trace)
						{
							traces.clear(); //remove earlier traces
							continue;
						}
						else if (p>traces.max_trace)
						{
            	break; //no more traces are possible
						}
					}
					traces.push_back(trace);
					traces.back().theoretical_int = pattern.theoretical_pattern.intensity[p];
				}
			}

			/**
				@brief Extends a single mass trace in one RT direction
				
				How to use this method:
				- Add the starting peak to the @p trace
				- Extend in downstream direction
				- extend in upstream direction
				
				@note this method assumes that it extends from a local maximum.
			*/
			void extendMassTrace_(MassTrace& trace, SignedSize spectrum_index, DoubleReal mz, bool inc_rt, Size meta_index_overall, DoubleReal min_rt=0.0, DoubleReal max_rt = 0.0) const
			{
				//Reverse peaks if we run the method for the second time (to keep them in chronological order)
				if (inc_rt)
				{
					++spectrum_index;
					std::reverse(trace.peaks.begin(), trace.peaks.end());
				}
				else
				{
					--spectrum_index;
				}
				//check if boundaries are set
				bool boundaries = false;
				if (max_rt != min_rt)
				{
					boundaries = true;
				}
				//Relax slope theshold if there is a hard boundary for the extension
				DoubleReal current_slope_bound = (1.0 + (DoubleReal)boundaries) * slope_bound_;
				Size delta_count = min_spectra_;
				DoubleReal last_int = trace.peaks.back().second->getIntensity();
				std::vector<DoubleReal> deltas(delta_count-1, 0);
				//deltas.reserve(2*delta_count);
				UInt missing_peaks = 0;
				Size peaks_before = trace.peaks.size();
				String abort_reason = "";
				while((!inc_rt && spectrum_index>=0) || (inc_rt && spectrum_index<(SignedSize)map_.size()))
				{
					if(boundaries && ((!inc_rt && map_[spectrum_index].getRT()<min_rt) || (inc_rt && map_[spectrum_index].getRT()>max_rt)) )
					{
						abort_reason = "Hit upper/lower boundary";
						break;
					}
					SignedSize peak_index = -1;
					try
					{
						peak_index = map_[spectrum_index].findNearest(mz);
					}
					catch(...) //no peaks in the spectrum
					{
						peak_index=-1;
					}
					if (peak_index<0 || map_[spectrum_index].getFloatDataArrays()[meta_index_overall][peak_index]<0.01 || positionScore_( mz, map_[spectrum_index][peak_index].getMZ(), trace_tolerance_)==0.0)
					{
						++missing_peaks;
						if(missing_peaks>max_missing_trace_peaks_)
						{
							abort_reason = "too many peaks missing";
							break;
						}
					}
					else
					{
						missing_peaks = 0;
						//add last peak to trace
						trace.peaks.push_back(std::make_pair(map_[spectrum_index].getRT(),&(map_[spectrum_index][peak_index])));
						//update ints and deltas 
						deltas.push_back( (map_[spectrum_index][peak_index].getIntensity() - last_int) / last_int);
						last_int = map_[spectrum_index][peak_index].getIntensity();

						//Abort if the average delta is too big (as intensity increases then)
						DoubleReal average_delta = std::accumulate(deltas.end()-delta_count,deltas.end(),0.0) / (DoubleReal)delta_count;
						if ( average_delta > current_slope_bound)
						{
							abort_reason = String("Average delta above threshold: ")+average_delta+"/"+current_slope_bound;
							//remove last peaks as we extended too far
							Size remove = std::min((Size)(trace.peaks.size()-peaks_before),delta_count-1);
							trace.peaks.erase(trace.peaks.end()-remove,trace.peaks.end());
							break;
						}
					}
					//increase/decrease scan index
					if (inc_rt) ++spectrum_index; else --spectrum_index;
				}
				log_ << "   - Added " << (trace.peaks.size()-peaks_before) << " peaks (abort: " << abort_reason << ")" << std::endl;
			}

			/// Returns the index of the peak nearest to m/z @p pos in spectrum @p spec (linear search starting from index @p start)
			template <typename SpectrumType>
			Size nearest_(DoubleReal pos, const SpectrumType& spec, Size start) const
			{
				Size index = start;
				DoubleReal dist = std::fabs(pos-spec[index].getMZ());
				++index;
				while (index < spec.size())
				{
					DoubleReal new_dist = std::fabs(pos-spec[index].getMZ());
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
			void findIsotope_(DoubleReal pos, Size spectrum_index, IsotopePattern& pattern, Size pattern_index, bool debug, Size& peak_index) const
			{
				if (debug) log_ << "   - Isotope " << pattern_index << ": "; 				
				
				DoubleReal intensity = 0.0;
				DoubleReal pos_score = 0.0;
				UInt matches = 0;
				
				//search in the center spectrum
				const SpectrumType& spectrum = map_[spectrum_index];
				peak_index = nearest_(pos, spectrum, peak_index);
				DoubleReal mz_score = positionScore_(pos, spectrum[peak_index].getMZ(), pattern_tolerance_);
				pattern.theoretical_mz[pattern_index] = pos;
				if (mz_score!=0.0)
				{
					if (debug) log_ << String::number(spectrum[peak_index].getIntensity(),1) << " ";
					pattern.peak[pattern_index] = peak_index;
					pattern.spectrum[pattern_index] = spectrum_index;
					intensity += spectrum[peak_index].getIntensity();
					pos_score += mz_score;
					++matches;
				}
				
				//previous spectrum
				if (spectrum_index!=0 && map_[spectrum_index-1].size()>0)
				{
					const SpectrumType& spectrum_before = map_[spectrum_index-1];
					Size index_before = spectrum_before.findNearest(pos);
					DoubleReal mz_score = positionScore_(pos, spectrum_before[index_before].getMZ(), pattern_tolerance_);
					if (mz_score!=0.0)
					{
						if (debug) log_ << String::number(spectrum_before[index_before].getIntensity(),1) << "b ";
						intensity += spectrum_before[index_before].getIntensity();
						pos_score += mz_score;
						++matches;

						if (pattern.peak[pattern_index]==-1)
						{
							pattern.peak[pattern_index] = index_before;
							pattern.spectrum[pattern_index] = spectrum_index-1;
						}
					}
				}
				//next spectrum
				if (spectrum_index!=map_.size()-1 && map_[spectrum_index+1].size()>0)
				{
					const SpectrumType& spectrum_after = map_[spectrum_index+1];
					Size index_after = spectrum_after.findNearest(pos);
					DoubleReal mz_score = positionScore_(pos, spectrum_after[index_after].getMZ(), pattern_tolerance_);
					if (mz_score!=0.0)
					{
						if (debug) log_ << String::number(spectrum_after[index_after].getIntensity(),1) << "a ";
						intensity += spectrum_after[index_after].getIntensity();
						pos_score += mz_score;
						++matches;
						
						if (pattern.peak[pattern_index]==-1)
						{
							pattern.peak[pattern_index] = index_after;
							pattern.spectrum[pattern_index] = spectrum_index+1;
						}
					}
				}
				//no isotope found
				if (matches==0)
				{
					if (debug) log_ << " missing" << std::endl;
					pattern.peak[pattern_index] = -1;
					pattern.mz_score[pattern_index] = 0.0;
					pattern.intensity[pattern_index] = 0.0;
				}
				else
				{
					if (debug) log_ << "=> " << intensity/matches << std::endl;
					pattern.mz_score[pattern_index] = pos_score/matches;
					pattern.intensity[pattern_index] = intensity/matches;
				}
			}

			/// Calculates a score between 0 and 1 for the m/z deviation of two peaks.
			DoubleReal positionScore_(DoubleReal pos1, DoubleReal pos2, DoubleReal allowed_deviation) const
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
			DoubleReal isotopeScore_(const TheoreticalIsotopePattern& isotopes, IsotopePattern& pattern, bool consider_mz_distances, bool debug) const
			{
				if (debug) log_ << "   - fitting " << pattern.intensity.size() << " peaks" << std::endl;
				//Abort if a core peak is missing
				for (Size iso=0+isotopes.optional_begin; iso<pattern.peak.size()-isotopes.optional_end; ++iso)
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
				Size best_begin = 0;
				for (Size i=isotopes.optional_begin; i>0; --i)
				{
					if (pattern.peak[i-1]==-1)
					{
						best_begin = i;
						break;
					}
				}
				Size best_end = 0;
				for (Size i=isotopes.optional_end; i>0; --i)
				{
					if (pattern.peak[pattern.peak.size()-i]==-1)
					{
						best_end = i;
						break;
					}
				}
				if (debug) log_ << "   - best_begin/end: " << best_begin << "/" << best_end << std::endl;
				for (Size b=best_begin; b<=isotopes.optional_begin; ++b)
				{
					for (Size e=best_end; e<=isotopes.optional_end; ++e)
					{
						//Make sure we have more than 2 peaks (unless in the first loop interation, there we allow two points) 
						if (isotopes.size()-b-e>2 || (b==best_begin && e==best_end && isotopes.size()-b-e>1))
						{
							DoubleReal int_score = Math::pearsonCorrelationCoefficient(isotopes.intensity.begin()+b, isotopes.intensity.end()-e, pattern.intensity.begin()+b, pattern.intensity.end()-e);	
							if (boost::math::isnan(int_score)) int_score = 0.0;
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
				
				//if the best fit is empty, abort
				if (pattern.mz_score.size()-best_begin-best_end==0)
				{
					return 0.0;
				}
				
				//remove left out peaks from the beginning
				for (Size i=0; i<best_begin; ++i)
				{
					pattern.peak[i] = -2;
					pattern.intensity[i] = 0.0;
					pattern.mz_score[i] = 0.0;
				}
				//remove left out peaks from the end
				for (Size i=0; i<best_end; ++i)
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
				OPENMS_POSTCONDITION(best_int_score>=0.0,  (String("Internal error: Isotope score (") + best_int_score + ") should be >=0.0").c_str())
				OPENMS_POSTCONDITION(best_int_score<=1.0,  (String("Internal error: Isotope score (") + best_int_score + ") should be <=1.0").c_str())
				return best_int_score;
			}
			
			DoubleReal intensityScore_(Size spectrum, Size peak) const
			{
				//calculate (half) bin numbers
				DoubleReal intensity  = map_[spectrum][peak].getIntensity();
				DoubleReal rt = map_[spectrum].getRT();
				DoubleReal mz = map_[spectrum][peak].getMZ();
				DoubleReal rt_min = map_.getMinRT();
				DoubleReal mz_min = map_.getMinMZ();
				UInt rt_bin = std::min(2*intensity_bins_-1,(UInt)std::floor((rt - rt_min) / intensity_rt_step_ * 2.0));
				UInt mz_bin = std::min(2*intensity_bins_-1,(UInt)std::floor((mz - mz_min) / intensity_mz_step_ * 2.0));
				//determine mz bins
				UInt ml,mh;
				if (mz_bin==0 || mz_bin==2*intensity_bins_-1)
				{
					ml = mz_bin/2;
					mh = mz_bin/2;
				}
				else if (Math::isOdd(mz_bin))
				{
					ml = mz_bin/2;
					mh = mz_bin/2+1;
				}
				else
				{
					ml = mz_bin/2-1; 
					mh = mz_bin/2;
				}
				//determine rt bins
				UInt rl,rh;
				if (rt_bin==0 || rt_bin==2*intensity_bins_-1)
				{
					rl = rt_bin/2;
					rh = rt_bin/2;
				}
				else if (Math::isOdd(rt_bin))
				{
					rl = rt_bin/2;
					rh = rt_bin/2+1;
				}
				else
				{
					rl = rt_bin/2-1; 
					rh = rt_bin/2;
				}
				//calculate distances to surrounding points (normalized to [0,1])
				DoubleReal drl = std::fabs(rt_min+(0.5+rl)*intensity_rt_step_-rt)/intensity_rt_step_;
				DoubleReal drh = std::fabs(rt_min+(0.5+rh)*intensity_rt_step_-rt)/intensity_rt_step_;
				DoubleReal dml = std::fabs(mz_min+(0.5+ml)*intensity_mz_step_-mz)/intensity_mz_step_;
				DoubleReal dmh = std::fabs(mz_min+(0.5+mh)*intensity_mz_step_-mz)/intensity_mz_step_;
				//Calculate weights for the intensity scores (the nearer to better)
				DoubleReal d1 = std::sqrt(std::pow(1.0-drl,2.0)+std::pow(1.0-dml,2.0));
				DoubleReal d2 = std::sqrt(std::pow(1.0-drh,2.0)+std::pow(1.0-dml,2.0));
				DoubleReal d3 = std::sqrt(std::pow(1.0-drl,2.0)+std::pow(1.0-dmh,2.0));
				DoubleReal d4 = std::sqrt(std::pow(1.0-drh,2.0)+std::pow(1.0-dmh,2.0));
				DoubleReal d_sum = d1 + d2 + d3 + d4;				
				//Final score
				DoubleReal final = intensityScore_(rl, ml, intensity)*(d1/d_sum)
													+ intensityScore_(rh, ml, intensity)*(d2/d_sum)
													+ intensityScore_(rl, mh, intensity)*(d3/d_sum)
													+ intensityScore_(rh, mh, intensity)*(d4/d_sum);
				
				OPENMS_POSTCONDITION(final>=0.0, (String("Internal error: Intensity score (") + final + ") should be >=0.0").c_str())
				OPENMS_POSTCONDITION(final<=1.0001, (String("Internal error: Intensity score (") + final + ") should be <=1.0").c_str())
				return final;
			}

			DoubleReal intensityScore_(Size rt_bin, Size mz_bin, DoubleReal intensity) const
			{
				//interpolate score value according to quantiles(20)
				const std::vector<DoubleReal>& quantiles20 = intensity_thresholds_[rt_bin][mz_bin];
				std::vector<DoubleReal>::const_iterator it = std::lower_bound(quantiles20.begin(),quantiles20.end(),intensity);
				//bigger than the biggest value => return 1.0
				if (it==quantiles20.end())
				{
					return 1.0;
				}
				//interpolate inside the bin
				DoubleReal bin_score = 0.0;
				if (it==quantiles20.begin())
				{
					bin_score = 0.05 * intensity / *it;
				}
				else
				{
					bin_score = 0.05 * (intensity-*(it-1)) / (*it-*(it-1));
				}
				
				DoubleReal final = bin_score + 0.05*((it - quantiles20.begin()) -1.0);
				
				//fix numerical problems
				if (final<0.0) final = 0.0;
				if (final>1.0) final = 1.0;				
				
				return final;
			}

			static int gaussF_(const gsl_vector* param, void* data, gsl_vector* f)
			{
				MassTraces* traces = static_cast<MassTraces*>(data);
				double height = gsl_vector_get (param, 0);
				double x0 = gsl_vector_get (param, 1);
				double sig = gsl_vector_get (param, 2);
				
				UInt count = 0;
				for (Size t=0; t< traces->size(); ++t)
				{
					MassTrace& trace = traces->at(t);			
					for (Size i=0; i<trace.peaks.size(); ++i)
					{
						gsl_vector_set(f, count, traces->baseline + trace.theoretical_int * height * exp(-0.5 * pow(trace.peaks[i].first - x0, 2)  / pow(sig, 2)) - trace.peaks[i].second->getIntensity() );
						++count;
					}
				}
				return GSL_SUCCESS;
			}
		
			static int gaussDF_(const gsl_vector* param, void* data, gsl_matrix* J)
			{
				MassTraces* traces = static_cast<MassTraces*>(data);
				double height = gsl_vector_get (param, 0);
				double x0 = gsl_vector_get (param, 1);
				double sig = gsl_vector_get (param, 2);
				
				UInt count = 0;
				for (Size t=0; t<traces->size(); ++t)
				{
					MassTrace& trace = traces->at(t);
					for (Size i=0; i< trace.peaks.size(); ++i)
					{
						DoubleReal rt = trace.peaks[i].first;
						gsl_matrix_set (J, count, 0, trace.theoretical_int * exp(-0.5 * pow(rt-x0,2) / pow(sig,2)));
						gsl_matrix_set (J, count, 1, trace.theoretical_int * height * exp(-0.5 * pow(rt-x0,2) / pow(sig,2)) * (rt-x0) / pow(sig,2));
						gsl_matrix_set (J, count, 2, 0.125 * trace.theoretical_int * height * exp(-0.5 * pow(rt-x0,2) / pow(sig,2)) * pow(rt-x0,2) / pow(sig,3));				
						++count;
					}
				}
			  return GSL_SUCCESS;
			}
		
			static int gaussFDF_(const gsl_vector* param, void* data, gsl_vector* f, gsl_matrix* J)
			{
			  gaussF_(param, data, f);
			  gaussDF_(param, data, J);
			  return GSL_SUCCESS;
			}

		private:
			
			/// Not implemented
			FeatureFinderAlgorithmPicked& operator=(const FeatureFinderAlgorithmPicked&);
			/// Not implemented
			FeatureFinderAlgorithmPicked(const FeatureFinderAlgorithmPicked&);
	};

} // namespace OpenMS

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMPICKED_H
