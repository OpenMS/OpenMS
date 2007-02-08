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
// $Maintainer: Alexandra Zerck $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_TWODOPTIMIZATION_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_TWODOPTIMIZATION_H

//#define DEBUG_2D
#undef DEBUG_2D

#ifdef DEBUG_2D
#include<iostream>
#include<fstream>
#endif

#include <vector>
#include <utility>
#include <cmath>

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakShape.h>
#include <OpenMS/KERNEL/DPickedPeak.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/DPeak.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/Param.h>

#ifndef OPENMS_SYSTEM_STOPWATCH_H
# include <OpenMS/SYSTEM/StopWatch.h>
#endif


#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePeakDeconvolution.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePick.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakShape.h>

namespace OpenMS
{   
	/// stores information about an isotopic cluser (i.e. potential peptide charge variant)
	struct IsotopeCluster
	{
		IsotopeCluster()
			: charge_(0), peaks_(), scans_()
		{}
      
		// predicted charge state of this peptide
		UnsignedInt charge_;
		// peaks in this cluster
		std::vector< std::pair<UnsignedInt,UnsignedInt> > peaks_;
		// the scans of this cluster
		std::vector<double> scans_;
	};
	
  /**
		 @brief Namespace for all functions and classes needed for 2D optimization.
	
	*/
  namespace TwoDOptimizationFunctions
  {
  
    /// Raw data point type
    typedef DRawDataPoint<1> RawDataPointType;
    extern std::vector<std::pair<int,int> > signal2D; 
    extern std::multimap<double,IsotopeCluster>::iterator iso_map_iter;
    extern unsigned int total_nr_peaks;
    extern std::map<int, std::vector<MSSpectrum<DPickedPeak<1> >::Iterator > > matching_peaks;
    extern MSExperiment<DPickedPeak<1> >::Iterator picked_peaks_iter;
    extern MSExperiment<RawDataPointType>::ConstIterator raw_data_first;


    /**
       @name Functions provided to the gsl Levenberg-Marquardt
		*/
    //@{
    /// Function computing estimated signal and its deviation to the experimental signal*/
    int residual2D(const gsl_vector* x, void* params , gsl_vector* f);
    /// Function computing the Jacobian */
    int jacobian2D(const gsl_vector* x, void* params, gsl_matrix* J);
    /// Function that calls residual2D and jacobian2D*/
    int evaluate2D(const gsl_vector* x, void* params, gsl_vector* f, gsl_matrix* J);

  }
 
	/**
		 @brief This class provides the two-dimensional optimization of the picked peak parameters.
			
		 Given the picked peaks, this class optimizes the peak parameters of each isotope pattern using
		 a non-linear optimization. The peaks of adjacent scans are adjusted to achieve that a peak occuring in
		 several scans has always the same m/z position. For the optimization the Levenberg-Marquardt algorithm
		 provided from the GSL is used. The optimized parameters are the m/z values,
		 the left and right width, which shall be equal for a peak in all scans,
		 and the peaks' heights.

		 @ingroup PeakPicking
			
		 @todo use DefaultParamHandler (Alexandra)
	*/
	class TwoDOptimization
	{
	public:


      
		/// Constructor
		TwoDOptimization()
		{
			// 2D optimization parameters
			param_.setValue("2D_optimization:penalties:position",0.0);
			param_.setValue("2D_optimization:penalties:height",1.0);
			param_.setValue("2D_optimization:penalties:left_width",0.0);
			param_.setValue("2D_optimization:penalties:right_width",0.0);
			param_.setValue("2D_optimization:thresholds:tolerance_mz",0.2);
			param_.setValue("2D_optimization:thresholds:max_peak_distance",1.0);
			param_.setValue("2D_optimization:skip_optimization","yes");
			param_.setValue("2D_optimization:delta_abs_error",1e-05f);
			param_.setValue("2D_optimization:delta_rel_error",1e-05f);
			param_.setValue("2D_optimization:iterations",10);
			param_.setValue("2D_optimization:real_2D","yes");
		}

		/// Constructor using a Param-object
		TwoDOptimization(const Param& param) 
		{
			param_ = param;
	
	
			DataValue dv = param_.getValue("2D_optimization:thresholds:tolerance_mz");
			if (dv.isEmpty() || dv.toString() == "") tolerance_mz_ = 0.1;
			else tolerance_mz_ = (float)dv;
	
			dv = param_.getValue("2D_optimization:thresholds:max_peak_distance");
			if (dv.isEmpty() || dv.toString() == "") max_peak_distance_ = 1;
			else max_peak_distance_ = (float)dv;
	
			dv = param_.getValue("2D_optimization:delta_abs_error");
			if (dv.isEmpty() || dv.toString() == "") eps_abs_ = 0.00000001;
			else eps_abs_ = (float)dv;

			dv = param_.getValue("2D_optimization:delta_rel_error");
			if (dv.isEmpty() || dv.toString() == "") eps_rel_ = 0.00000001;
			else eps_rel_ = (float)dv;
	
			dv = param_.getValue("2D_optimization:iterations");
			if (dv.isEmpty() || dv.toString() == "") max_iteration_ = 15;
			else max_iteration_ = (unsigned int)dv;



			penalties_ = OptimizationFunctions::PenaltyFactorsInt();
			dv = param_.getValue("2D_optimization:penalties:height");
			if (dv.isEmpty() || dv.toString() == "") penalties_.height = 2;
			else penalties_.height = (float)dv;

			dv = param_.getValue("2D_optimization:penalties:left_width");
			if (dv.isEmpty() || dv.toString() == "") penalties_.lWidth = 1;
			else penalties_.lWidth = (float)dv;
	
			dv = param_.getValue("2D_optimization:penalties:right_width");
			if (dv.isEmpty() || dv.toString() == "") penalties_.rWidth = 1;
			else penalties_.rWidth = (float)dv;

			dv = param_.getValue("2D_optimization:penalties:position");
			if (dv.isEmpty() || dv.toString() == "") penalties_.pos = 2;
			else penalties_.pos = (float)dv;
		}

		/// Copy constructor
		TwoDOptimization(const TwoDOptimization& opt)
			:tolerance_mz_(opt.tolerance_mz_),
			 eps_abs_(opt.eps_abs_),
			 eps_rel_(opt.eps_rel_),
			 max_iteration_(opt.max_iteration_)
		{
			param_ = opt.param_;
			penalties_ = opt.penalties_;
		}

		/// Destructor
		virtual ~TwoDOptimization(){}

		/// Assignment operator
		TwoDOptimization& operator=(const TwoDOptimization& opt)
		{
			if(&opt == this) return *this; 
	
			tolerance_mz_ = opt.tolerance_mz_;
			eps_abs_ = opt.eps_abs_;
			eps_rel_ = opt.eps_rel_;
			max_iteration_ = opt.max_iteration_;
			param_ = opt.param_;
			penalties_ = opt.penalties_;
			return *this;
		}
			
		///Non-mutable access to the matching epsilon
		inline const double& getMZTolerance() const {return tolerance_mz_;}
		///Mutable access to the matching epsilon
		inline double& getMZTolerance() {return tolerance_mz_;}
		///Mutable access to the matching epsilon
		inline void setMZTolerance(double tolerance_mz) { tolerance_mz_ = tolerance_mz;}

		///Non-mutable access to the maximal peak distance in a cluster
		inline const double& getMaxPeakDistance() const {return max_peak_distance_;}
		///Mutable access to the maximal peak distance in a cluster
		inline double& getMaxPeakDistance() {return max_peak_distance_;}
		///Mutable access to the maximal peak distance in a cluster
		inline void setMaxPeakDistance(double max_peak_distance) { max_peak_distance_ = max_peak_distance;}

		///Non-mutable access to the maximal absolute error
		inline const double& getMaxAbsError() const {return eps_abs_;}
		///Mutable access to the  maximal absolute error
		inline double& getMaxAbsError() {return eps_abs_;}
		///Mutable access to the  maximal absolute error
		inline void setMaxAbsError(double eps_abs) { eps_abs_ = eps_abs;}
      
		///Non-mutable access to the maximal relative error
		inline const double& getMaxRelError() const {return eps_rel_;}
		///Mutable access to the maximal relative error
		inline double& getMaxRelError() {return eps_rel_;}
		///Mutable access to the maximal relative error
		inline void setMaxRelError(double eps_rel) { eps_rel_ = eps_rel;}

		///Non-mutable access to the maximal number of iterations
		inline const int& getMaxIterations() const {return max_iteration_;}
		///Mutable access to the  maximal number of iterations
		inline int& getMaxIterations() {return max_iteration_;}
		///Mutable access to the  maximal number of iterations
		inline void setMaxIterations(int max_iteration) { max_iteration_ = max_iteration;}

		///Non-mutable access to the minimal number of adjacent scans
		inline const OptimizationFunctions::PenaltyFactorsInt& getPenalties() const {return penalties_;}
		///Mutable access to the minimal number of adjacent scans
		inline OptimizationFunctions::PenaltyFactorsInt& getPenalties() {return penalties_;}
		///Mutable access to the minimal number of adjacent scans
		inline void setPenalties(OptimizationFunctions::PenaltyFactorsInt& penalties) { penalties_ = penalties;}



		/** Find two dimensional peak clusters and optimize their peak parameters */
		template <typename InputSpectrumIterator,typename OutputPeakType>
		void twoDOptimize(InputSpectrumIterator& first,
											InputSpectrumIterator& last,
											MSExperiment< OutputPeakType >& ms_exp,bool real2D=true);


	protected:

		/// stores the retention time of each isotopic cluster
		std::multimap<double,IsotopeCluster> iso_map_;
      
		/// Pointer to the current region
		std::multimap<double,IsotopeCluster>::const_iterator curr_region_;

		/// upper bound for distance between two peaks belonging to the same region
		double max_peak_distance_;

		/// threshold for the difference in the peak position of two matching peaks
		double tolerance_mz_;
      
		/// Indices of peaks in the adjacent scans matching peaks in the scan with no. ref_scan
		std::map<int, std::vector<MSSpectrum<DPickedPeak<1> >::Iterator > > matching_peaks_;
      
		/// Convergence Parameter: Maximal absolute error
		double eps_abs_;
      
		/// Convergence Parameter: Maximal relative error 
		double eps_rel_;

		/// Convergence Parameter: Maximal number of iterations 
		int max_iteration_;

		/// Optimization considering all scans of a cluster or optimization of each scan separately
		bool real_2D_;

		/// Parameter object 
		Param param_;
      
		/// Penalty factors for some parameters in the optimization
		OptimizationFunctions::PenaltyFactorsInt penalties_;
      
		/**
			 @name Auxiliary Functions for the search of matching regions
		*/
		//@{
		std::vector<double>::iterator searchInScan_(std::vector<double>::iterator scan_begin,
																								std::vector<double>::iterator scan_end ,
																								double current_mz);

		/** Performs 2D optimization of all regions */
		template <typename InputSpectrumIterator,typename OutputPeakType>
		void optimizeRegions_(InputSpectrumIterator& first,
													InputSpectrumIterator& last,
													MSExperiment<OutputPeakType>& ms_exp);

		/** Performs an optimization of all regions by calling OptimizePick */
		template <typename InputSpectrumIterator,typename OutputPeakType>
		void optimizeRegionsScanwise_(InputSpectrumIterator& first,
																	InputSpectrumIterator& last,
																	MSExperiment<OutputPeakType>& ms_exp);
    

		/// Get the indices of the first and last raw data point of this region
		template<typename InputSpectrumIterator,typename OutputPeakType>
		void getRegionEndpoints_(MSExperiment<OutputPeakType>& exp,
														 InputSpectrumIterator& first,
														 InputSpectrumIterator& last,
														 unsigned int iso_map_idx,
														 double noise_level);

		/// Identify matching peak in a peak cluster
		void findMatchingPeaks_(std::multimap<double, IsotopeCluster>::iterator& it,
														MSExperiment< DPickedPeak<1> >& ms_exp);
      
		//@}
	};


	template <typename InputSpectrumIterator,typename OutputPeakType>
	void TwoDOptimization::twoDOptimize(InputSpectrumIterator& first,
																			InputSpectrumIterator& last,
																			MSExperiment< OutputPeakType >& ms_exp,bool real2D)
	{
		real_2D_ = real2D;
			
		typedef typename InputSpectrumIterator::value_type InputSpectrumType;
		typedef typename InputSpectrumType::value_type RawDataPointType;
		typedef MSSpectrum<RawDataPointType> SpectrumType;

		typename MSExperiment<OutputPeakType>::Iterator ms_exp_it = ms_exp.begin();
		typename MSExperiment<OutputPeakType>::Iterator ms_exp_it_end = ms_exp.end();
		// stores the monoisotopic peaks of isotopic clusters
		std::vector<double> iso_last_scan;
		std::vector<double> iso_curr_scan;
		std::vector<std::multimap<double,IsotopeCluster>::iterator> clusters_last_scan;
		std::vector<std::multimap<double,IsotopeCluster>::iterator> clusters_curr_scan;
		std::multimap<double,IsotopeCluster>::iterator cluster_iter;
		double current_rt=ms_exp_it->getRetentionTime(),last_rt  = 0;

		// retrieve values for accepted peaks distances
		max_peak_distance_ = param_.getValue("2D_optimization:thresholds:max_peak_distance");
		double tolerance_mz = param_.getValue("2D_optimization:thresholds:tolerance_mz");
	
		UnsignedInt current_charge     = 0;			// charge state of the current isotopic cluster
		double mz_in_hash   = 0;			// used as reference to the current isotopic peak			
	
		// sweep through scans
		for (unsigned int curr_scan =0; ms_exp_it+curr_scan != ms_exp_it_end;++curr_scan)
			{
				unsigned int nr_peaks_in_scan = (ms_exp_it +curr_scan)->size();
				//last_rt = current_rt;
				current_rt = (ms_exp_it+curr_scan)->getRetentionTime();
				typename MSSpectrum<OutputPeakType>::Iterator peak_it  = (ms_exp_it+curr_scan)->begin();
				typename MSSpectrum<OutputPeakType>::Iterator peak_it_last  = (ms_exp_it+curr_scan)->end();

				// copy cluster information of least scan
				iso_last_scan = iso_curr_scan;
				iso_curr_scan.clear();
				clusters_last_scan = clusters_curr_scan;
				clusters_curr_scan.clear();

#ifdef DEBUG_2D
				std::cout << "Next scan with rt: " << current_rt << std::endl;
				std::cout << "Next scan, rt = "<<current_rt<<" last_rt: "<<last_rt<< std::endl;
				std::cout << "---------------------------------------------------------------------------" << std::endl;
#endif	  
				MSSpectrum<RawDataPointType> s;
				s.setRetentionTime(current_rt);
				// check if there were scans in between
				if ( last_rt == 0|| // are we in the first scan
						 ( (lower_bound(first,last,s, typename SpectrumType::RTLess())-1)->getRetentionTime() == last_rt))
					{
	      
	  
						for(unsigned int curr_peak=0; peak_it+curr_peak < peak_it_last-1;++curr_peak)
							{
		  
								// store the m/z of the current peak
								double curr_mz         = (peak_it+curr_peak)->getPos();
								double dist2nextpeak = (peak_it+curr_peak+1)->getPos() - curr_mz;
		  
								if (dist2nextpeak <= max_peak_distance_) // one single peak without neighbors isn't optimized
									{
#ifdef DEBUG_2D	      
										std::cout << "Isotopic pattern found ! " << std::endl;
										std::cout << "We are at: " << (peak_it+curr_peak)->getPos()  << " " << curr_mz << std::endl;
#endif      
										if (iso_last_scan.size() > 0)  // Did we find any isotopic cluster in the last scan?
											{
												// there were some isotopic clustures in the last scan...
												std::vector<double>::iterator it =
													searchInScan_(iso_last_scan.begin(),iso_last_scan.end(),curr_mz);
			  
												double delta_mz = fabs(*it - curr_mz);
												std::vector<double>::iterator itneu = iso_last_scan.begin();
			  
												if ( delta_mz > tolerance_mz) // check if first peak of last cluster is close enough
													{
														mz_in_hash = curr_mz; // update current hash key
			      
														// create new isotopic cluster
#ifdef DEBUG_2D
														std::cout << "Last peak cluster too far, creating new cluster at "<<curr_mz << std::endl;
#endif
														IsotopeCluster new_cluster;
														new_cluster.charge_  = current_charge;
														new_cluster.scans_.push_back( current_rt );					
														cluster_iter = iso_map_.insert(std::pair<double,IsotopeCluster>(mz_in_hash,new_cluster));
			      
													}
												else
													{
#ifdef DEBUG_2D
														std::cout << "Found neighbouring peak with distance (m/z) " << delta_mz << std::endl;
#endif
														cluster_iter = clusters_last_scan[distance(iso_last_scan.begin(),it)];

														// check whether this scan is already contained
														if(find(cluster_iter->second.scans_.begin(),cluster_iter->second.scans_.end(),current_rt)
															 == cluster_iter->second.scans_.end())
															{
																cluster_iter->second.scans_.push_back( current_rt );
															}
			      
#ifdef DEBUG_2D
														std::cout << "Cluster with " << cluster_iter->second.peaks_.size()
																			<< " peaks retrieved." << std::endl;
#endif
													}
			  
											}
										else // last scan did not contain any isotopic cluster
											{	
#ifdef DEBUG_2D
												std::cout << "Last scan was empty => creating new cluster." << std::endl;
												std::cout << "Creating new cluster at m/z: " << curr_mz << std::endl;
#endif
			  
												mz_in_hash = curr_mz; // update current hash key
			  
												// create new isotopic cluster
												IsotopeCluster new_cluster;
												new_cluster.charge_  = current_charge;
												new_cluster.scans_.push_back( current_rt );					
												cluster_iter = iso_map_.insert(std::pair<double,IsotopeCluster>(mz_in_hash,new_cluster));

											}
		  
#ifdef DEBUG_2D
										std::cout << "Storing found peak in current isotopic cluster" << std::endl;
#endif


		      
										cluster_iter->second.peaks_.push_back(std::pair<UnsignedInt,UnsignedInt>(curr_scan,curr_peak));
		      
										iso_curr_scan.push_back(  mz_in_hash );
										clusters_curr_scan.push_back(cluster_iter);
										++curr_peak;
		      
										cluster_iter->second.peaks_.push_back(std::pair<UnsignedInt,UnsignedInt>(curr_scan,curr_peak));
										iso_curr_scan.push_back((peak_it+curr_peak)->getPos());
										clusters_curr_scan.push_back(cluster_iter);
		      
										// check distance to next peak
										if ( (curr_peak+1) >= nr_peaks_in_scan ) break;
										dist2nextpeak = (peak_it+curr_peak+1)->getPos() -  (peak_it+curr_peak)->getPos();
		      
		
										// loop until end of isotopic pattern in this scan
										while (dist2nextpeak <= max_peak_distance_
													 &&  curr_peak < (nr_peaks_in_scan-1) )
											{
												cluster_iter->second.peaks_.push_back(std::pair<UnsignedInt,UnsignedInt>(curr_scan,curr_peak+1));				// save peak in cluster
												iso_curr_scan.push_back((peak_it+curr_peak+1)->getPos());
												clusters_curr_scan.push_back(cluster_iter);
												// std::cout << "new enter'd: "<<(peak_it+curr_peak+1)->getPos()<<" im while"<<std::endl;
												++curr_peak;			
												if(curr_peak >= nr_peaks_in_scan-1) break;
												dist2nextpeak = (peak_it+curr_peak+1)->getPos() -  (peak_it+curr_peak)->getPos(); // get distance to next peak

			  
											} // end while(...)
		      
		    
		      
									} // end of if (charge > 0)
	      
								current_charge = 0; // reset charge
							} // end for (...)
					}
				last_rt = current_rt;
			}
		curr_region_ = iso_map_.begin();
#ifdef DEBUG_2D
		std::cout << iso_map_.size() << " isotopic clusters were found ! " << std::endl;
#endif

		if(real_2D_)    optimizeRegions_(first,last,ms_exp);
		else           optimizeRegionsScanwise_(first,last,ms_exp);
	}



    
	template < >
	void TwoDOptimization::
	optimizeRegions_<MSExperiment<DRawDataPoint<1> >::const_iterator, DPickedPeak<1> >
	(MSExperiment<DRawDataPoint<1> >::const_iterator& first,
	 MSExperiment<DRawDataPoint<1> >::const_iterator& last,
	 MSExperiment<DPickedPeak<1> >& ms_exp);

	template < >
	void TwoDOptimization::
	optimizeRegionsScanwise_<MSExperiment<DRawDataPoint<1> >::const_iterator, DPickedPeak<1> >
	(MSExperiment<DRawDataPoint<1> >::const_iterator& first,
	 MSExperiment<DRawDataPoint<1> >::const_iterator& last,
	 MSExperiment<DPickedPeak<1> >& ms_exp);
    
	template <typename InputSpectrumIterator,typename OutputPeakType>
	void TwoDOptimization::optimizeRegions_(InputSpectrumIterator& /*first*/,
																					InputSpectrumIterator& /*last*/,
																					MSExperiment< OutputPeakType >& /*ms_exp*/)
	{
		throw Exception::IllegalArgument(__FILE__,
																		 __LINE__,
																		 __PRETTY_FUNCTION__,
																		 "wrong input peak type, must be DPickedPeak<1>,");
      
      
	}

  template <typename InputSpectrumIterator,typename OutputPeakType>
	void TwoDOptimization::optimizeRegionsScanwise_(InputSpectrumIterator& /*first*/,
																									InputSpectrumIterator& /*last*/,
																									MSExperiment< OutputPeakType >& /*ms_exp*/)
	{
		throw Exception::IllegalArgument(__FILE__,
																		 __LINE__,
																		 __PRETTY_FUNCTION__,
																		 "wrong input peak type, must be DPickedPeak<1>,");
      
      
	}


	template<typename InputSpectrumIterator,typename OutputPeakType>
	void TwoDOptimization::getRegionEndpoints_(MSExperiment<OutputPeakType>& exp,
																						 InputSpectrumIterator& first,
																						 InputSpectrumIterator& last,
																						 unsigned int iso_map_idx,
																						 double noise_level)
	{
		TwoDOptimizationFunctions::signal2D.clear();
		typedef typename InputSpectrumIterator::value_type InputExperimentType;
		typedef typename InputExperimentType::value_type InputPeakType;
		typedef std::multimap<double,IsotopeCluster> MapType;

		double rt, first_peak_mz, last_peak_mz;
		MSSpectrum<InputPeakType> spec;
		InputPeakType peak;

		MapType::iterator iso_map_iter = iso_map_.begin();
		for(unsigned int i=0;i<iso_map_idx;++i)  ++iso_map_iter;
	
		// get left and right endpoint for all scans in the current cluster
		for(unsigned int i=0; i< iso_map_iter->second.scans_.size();++i)
			{
				typename MSExperiment<OutputPeakType>::iterator exp_it;
	  
				// first the right scan through binary search
				rt = iso_map_iter->second.scans_[i];
				spec.setRetentionTime(rt);
				InputSpectrumIterator iter = lower_bound(first, last, spec, typename MSSpectrum<InputPeakType>::RTLess());

				exp_it = exp.RTBegin(rt);
				// now the right mz
				unsigned int j=0;
				while(j < (iso_map_iter->second.peaks_.size()-1) &&
							iso_map_iter->second.peaks_[j].first != iso_map_iter->second.peaks_[0].first + i)
					{
						++j;
					}
				// consider a bit more of the signal to the left
				first_peak_mz = (exp_it->begin() + iso_map_iter->second.peaks_[j].second)->getPos() - 1;
				if(j == iso_map_iter->second.peaks_.size()-1)
					{
						last_peak_mz = first_peak_mz;
					}
				else{
					while(j < iso_map_iter->second.peaks_.size()-1  &&
								iso_map_iter->second.peaks_[j].first == iso_map_iter->second.peaks_[0].first + i) 
						{
							++j;
						}
					if(j >= iso_map_iter->second.peaks_.size() || 
						 iso_map_iter->second.peaks_[j].first != iso_map_iter->second.peaks_[0].first + i)  --j;
					last_peak_mz = (exp_it->begin() + iso_map_iter->second.peaks_[j].second)->getPos();
				}

				peak.setPos(first_peak_mz);
				typename MSSpectrum<InputPeakType>::const_iterator raw_data_iter
					= lower_bound(iter->begin(), iter->end(), peak, typename InputPeakType::PositionLess());
				if(raw_data_iter != iter->begin())
					{
						--raw_data_iter;
					}
				double intensity = raw_data_iter->getIntensity();
				// while the intensity is falling go to the left
				while(raw_data_iter != iter->begin() && (raw_data_iter-1)->getIntensity() < intensity &&
							(raw_data_iter-1)->getIntensity() > noise_level)
					{
						--raw_data_iter;
						intensity = raw_data_iter->getIntensity();
					}
				++raw_data_iter;
				std::pair<unsigned int,unsigned int> left,right;
				left.first = distance(first,iter);
				left.second = distance(iter->begin(),raw_data_iter);

				// consider a bit more of the signal to the right
				peak.setPos(last_peak_mz + 1);
				raw_data_iter
					= upper_bound(iter->begin(), iter->end(), peak, typename InputPeakType::PositionLess());
				if(raw_data_iter == iter->end()) --raw_data_iter;
				intensity = raw_data_iter->getIntensity();
				// while the intensity is falling go to the right
				while(raw_data_iter+1 != iter->end() && (raw_data_iter+1)->getIntensity() < intensity)
					{
						++raw_data_iter;
						intensity = raw_data_iter->getIntensity();
						if( (raw_data_iter+1!= iter->end()) && (raw_data_iter+1)->getIntensity() > noise_level ) break;
					}
				right.first = left.first;
				right.second = distance(iter->begin(),raw_data_iter);
				// region endpoints are stored in global vector
				TwoDOptimizationFunctions::signal2D.push_back(left);
				TwoDOptimizationFunctions::signal2D.push_back(right);

			}


	}
    

    
  
}

#endif //OPENMS_TRANSFORMATIONS_RAW2PEAK_TWODOPTIMIZATION_H
