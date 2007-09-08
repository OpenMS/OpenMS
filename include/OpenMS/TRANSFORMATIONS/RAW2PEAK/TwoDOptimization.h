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
#include <set>

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakShape.h>
#include <OpenMS/KERNEL/PickedPeak1D.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/DPeak.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/IsotopeCluster.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

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
	
	typedef std::pair<unsigned int,unsigned int> Idx  ;
	typedef std::set<Idx> IndexSet;

  namespace OptimizationFunctions
  {
  
    /// Raw data point type
    typedef RawDataPoint1D RawDataPointType;
    extern std::vector<std::pair<int,int> > signal2D; 
    extern std::multimap<double,IsotopeCluster>::iterator iso_map_iter;
    extern unsigned int total_nr_peaks;
    extern std::map<int, std::vector<MSSpectrum<PickedPeak1D>::Iterator > > matching_peaks;
    extern MSExperiment<PickedPeak1D>::Iterator picked_peaks_iter;
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
		 
		 @ref TwoDOptimization_Parameters are explained on a separate page.

		 @ingroup PeakPicking
	*/
	class TwoDOptimization : public DefaultParamHandler
	{
	public:

		///Comparator for the retention time.
		struct IndexLess
			: public std::binary_function <Idx,Idx, bool>
		{
			inline bool operator () (const Idx& a, const Idx& b) const
			{
				return (a.first < b.first);
			}
		};
      
		/// Constructor
		TwoDOptimization();

		/// Copy constructor
		TwoDOptimization(const TwoDOptimization& opt);

		/// Destructor
		virtual ~TwoDOptimization(){}

		/// Assignment operator
		TwoDOptimization& operator=(const TwoDOptimization& opt);

			
		///Non-mutable access to the matching epsilon
		inline DoubleReal getMZTolerance() const {return tolerance_mz_;}
		///Mutable access to the matching epsilon
		inline void setMZTolerance(double tolerance_mz)
		{
			tolerance_mz_ = tolerance_mz;
			param_.setValue("thresholds:tolerance_mz",tolerance_mz);
		}

		///Non-mutable access to the maximal peak distance in a cluster
		inline DoubleReal getMaxPeakDistance() const {return max_peak_distance_;}
		///Mutable access to the maximal peak distance in a cluster
		inline void setMaxPeakDistance(double max_peak_distance)
		{
			max_peak_distance_ = max_peak_distance;
			param_.setValue("thresholds:max_peak_distance",max_peak_distance);
		}

		///Non-mutable access to the maximal absolute error
		inline DoubleReal getMaxAbsError() const {return eps_abs_;}
		///Mutable access to the  maximal absolute error
		inline void setMaxAbsError(double eps_abs)
		{
			eps_abs_ = eps_abs;
			param_.setValue("delta_abs_error",eps_abs);
		}
      
		///Non-mutable access to the maximal relative error
		inline DoubleReal getMaxRelError() const {return eps_rel_;}
		///Mutable access to the maximal relative error
		inline void setMaxRelError(double eps_rel)
		{
			eps_rel_ = eps_rel;
			param_.setValue("delta_rel_error",eps_rel);
		}

		///Non-mutable access to the maximal number of iterations
		inline Int getMaxIterations() const {return max_iteration_;}
		///Mutable access to the  maximal number of iterations
		inline void setMaxIterations(int max_iteration)
		{
			max_iteration_ = max_iteration;
			param_.setValue("iterations",max_iteration);
		}

		///Non-mutable access to the minimal number of adjacent scans
		inline const OptimizationFunctions::PenaltyFactorsIntensity& getPenalties() const {return penalties_;}
		///Mutable access to the minimal number of adjacent scans
		inline void setPenalties(OptimizationFunctions::PenaltyFactorsIntensity& penalties)
		{
			penalties_ = penalties;
			param_.setValue("penalties:position",penalties.pos);
			param_.setValue("penalties:height",penalties.height);
			param_.setValue("penalties:left_width",penalties.lWidth);
			param_.setValue("penalties:right_width",penalties.rWidth);
		}



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
		std::map<int, std::vector<MSSpectrum<PickedPeak1D>::Iterator > > matching_peaks_;
      
		/// Convergence Parameter: Maximal absolute error
		double eps_abs_;
      
		/// Convergence Parameter: Maximal relative error 
		double eps_rel_;

		/// Convergence Parameter: Maximal number of iterations 
		int max_iteration_;

		/// Optimization considering all scans of a cluster or optimization of each scan separately
		bool real_2D_;

    
		/// Penalty factors for some parameters in the optimization
		OptimizationFunctions::PenaltyFactorsIntensity penalties_;
      
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
														MSExperiment<PickedPeak1D>& ms_exp);
      
		//@}

		/// update members method from DefaultParamHandler to update the members 
		void updateMembers_();
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
		double current_rt=ms_exp_it->getRT(),last_rt  = 0;

		// retrieve values for accepted peaks distances
		max_peak_distance_ = param_.getValue("thresholds:max_peak_distance");
		double tolerance_mz = param_.getValue("thresholds:tolerance_mz");
	
		UInt current_charge     = 0;			// charge state of the current isotopic cluster
		double mz_in_hash   = 0;			// used as reference to the current isotopic peak			
	
		// sweep through scans
		for (unsigned int curr_scan =0; ms_exp_it+curr_scan != ms_exp_it_end;++curr_scan)
			{
				unsigned int nr_peaks_in_scan = (ms_exp_it +curr_scan)->size();
				//last_rt = current_rt;
				current_rt = (ms_exp_it+curr_scan)->getRT();
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
				s.setRT(current_rt);
				// check if there were scans in between
				if ( last_rt == 0|| // are we in the first scan
						 ( (lower_bound(first,last,s, typename SpectrumType::RTLess())-1)->getRT() == last_rt))
					{
	      
	  
						for(unsigned int curr_peak=0; peak_it+curr_peak < peak_it_last-1;++curr_peak)
							{
		  
								// store the m/z of the current peak
								double curr_mz         = (peak_it+curr_peak)->getMZ();
								double dist2nextpeak = (peak_it+curr_peak+1)->getMZ() - curr_mz;
		  
								if (dist2nextpeak <= max_peak_distance_) // one single peak without neighbors isn't optimized
									{
#ifdef DEBUG_2D	      
										std::cout << "Isotopic pattern found ! " << std::endl;
										std::cout << "We are at: " << (peak_it+curr_peak)->getMZ()  << " " << curr_mz << std::endl;
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
														new_cluster.peaks_.charge_  = current_charge;
														new_cluster.scans_.push_back( curr_scan );					
														cluster_iter = iso_map_.insert(std::pair<double,IsotopeCluster>(mz_in_hash,new_cluster));
			      
													}
												else
													{
#ifdef DEBUG_2D
														std::cout << "Found neighbouring peak with distance (m/z) " << delta_mz << std::endl;
#endif
														cluster_iter = clusters_last_scan[distance(iso_last_scan.begin(),it)];

														// check whether this scan is already contained
														if(find(cluster_iter->second.scans_.begin(),cluster_iter->second.scans_.end(),curr_scan)
															 == cluster_iter->second.scans_.end())
															{
																cluster_iter->second.scans_.push_back( curr_scan );
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
												new_cluster.peaks_.charge_  = current_charge;
												new_cluster.scans_.push_back( curr_scan );					
												cluster_iter = iso_map_.insert(std::pair<double,IsotopeCluster>(mz_in_hash,new_cluster));

											}
		  
#ifdef DEBUG_2D
										std::cout << "Storing found peak in current isotopic cluster" << std::endl;
#endif


		      
										cluster_iter->second.peaks_.insert(std::pair<UInt,UInt>(curr_scan,curr_peak));
		      
										iso_curr_scan.push_back(  mz_in_hash );
										clusters_curr_scan.push_back(cluster_iter);
										++curr_peak;
		      
										cluster_iter->second.peaks_.insert(std::pair<UInt,UInt>(curr_scan,curr_peak));
										iso_curr_scan.push_back((peak_it+curr_peak)->getMZ());
										clusters_curr_scan.push_back(cluster_iter);
		      
										// check distance to next peak
										if ( (curr_peak+1) >= nr_peaks_in_scan ) break;
										dist2nextpeak = (peak_it+curr_peak+1)->getMZ() -  (peak_it+curr_peak)->getMZ();
		      
		
										// loop until end of isotopic pattern in this scan
										while (dist2nextpeak <= max_peak_distance_
													 &&  curr_peak < (nr_peaks_in_scan-1) )
											{
												cluster_iter->second.peaks_.insert(std::pair<UInt,UInt>(curr_scan,curr_peak+1));				// save peak in cluster
												iso_curr_scan.push_back((peak_it+curr_peak+1)->getMZ());
												clusters_curr_scan.push_back(cluster_iter);
												// std::cout << "new enter'd: "<<(peak_it+curr_peak+1)->getMZ()<<" im while"<<std::endl;
												++curr_peak;			
												if(curr_peak >= nr_peaks_in_scan-1) break;
												dist2nextpeak = (peak_it+curr_peak+1)->getMZ() -  (peak_it+curr_peak)->getMZ(); // get distance to next peak

			  
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
	void TwoDOptimization::optimizeRegions_<MSExperiment<RawDataPoint1D >::const_iterator, PickedPeak1D >
	(MSExperiment<RawDataPoint1D >::const_iterator& first,
	 MSExperiment<RawDataPoint1D >::const_iterator& last,
	 MSExperiment<PickedPeak1D >& ms_exp);

	template < >
	void TwoDOptimization::optimizeRegionsScanwise_<MSExperiment<RawDataPoint1D >::const_iterator, PickedPeak1D >
	(MSExperiment<RawDataPoint1D >::const_iterator& first,
	 MSExperiment<RawDataPoint1D >::const_iterator& last,
	 MSExperiment<PickedPeak1D >& ms_exp);
    
	template <typename InputSpectrumIterator,typename OutputPeakType>
	void TwoDOptimization::optimizeRegions_(InputSpectrumIterator& /*first*/,
																					InputSpectrumIterator& /*last*/,
																					MSExperiment< OutputPeakType >& /*ms_exp*/)
	{
		throw Exception::IllegalArgument(__FILE__,
																		 __LINE__,
																		 __PRETTY_FUNCTION__,
																		 "wrong input peak type, must be PickedPeak1D,");
      
      
	}

  template <typename InputSpectrumIterator,typename OutputPeakType>
	void TwoDOptimization::optimizeRegionsScanwise_(InputSpectrumIterator& /*first*/,
																									InputSpectrumIterator& /*last*/,
																									MSExperiment< OutputPeakType >& /*ms_exp*/)
	{
		throw Exception::IllegalArgument(__FILE__,
																		 __LINE__,
																		 __PRETTY_FUNCTION__,
																		 "wrong input peak type, must be PickedPeak1D,");
      
      
	}


	template<typename InputSpectrumIterator,typename OutputPeakType>
	void TwoDOptimization::getRegionEndpoints_(MSExperiment<OutputPeakType>& exp,
																						 InputSpectrumIterator& first,
																						 InputSpectrumIterator& last,
																						 unsigned int iso_map_idx,
																						 double noise_level)
	{
		OptimizationFunctions::signal2D.clear();
		typedef typename InputSpectrumIterator::value_type InputExperimentType;
		typedef typename InputExperimentType::value_type InputPeakType;
		typedef std::multimap<double,IsotopeCluster> MapType;

		double rt, first_peak_mz, last_peak_mz;
		MSSpectrum<InputPeakType> spec;
		InputPeakType peak;

		MapType::iterator iso_map_iter = iso_map_.begin();
		for(unsigned int i=0;i<iso_map_idx;++i)  ++iso_map_iter;

#ifdef DEBUG2D
		std::cout << "rt begin: "<<exp[iso_map_iter->second.scans_[0]].getRT()
							<< "\trt end: "<<exp[iso_map_iter->second.scans_[iso_map_iter->second.scans_.size()-1]].getRT()
							<< " \t"<<iso_map_iter->second.scans_.size()<<" scans"
							<< std::endl;
#endif
		
		// get left and right endpoint for all scans in the current cluster
		for(unsigned int i=0; i< iso_map_iter->second.scans_.size();++i)
			{
				typename MSExperiment<OutputPeakType>::iterator exp_it;
	  
				// first the right scan through binary search
				rt = exp[iso_map_iter->second.scans_[i]].getRT();
				spec.setRT(rt);
				InputSpectrumIterator iter = lower_bound(first, last, spec, typename MSSpectrum<InputPeakType>::RTLess());
				//				if(iter->getRT() != rt) --iter;
				exp_it = exp.RTBegin(rt);
#ifdef DEBUG2D
				std::cout << exp_it->getRT() << " vs "<< iter->getRT()<<std::endl;
#endif
				// now the right mz
				IndexSet::const_iterator j=(iso_map_iter->second.peaks_.begin());
																		
				Idx pair;
				pair.first =  iso_map_iter->second.peaks_.begin()->first + i;
				IndexSet::const_iterator set_iter = lower_bound(iso_map_iter->second.peaks_.begin(),
																												iso_map_iter->second.peaks_.end(),
																												pair,IndexLess());
				
				// consider a bit more of the signal to the left
				first_peak_mz = (exp_it->begin() + set_iter->second)->getMZ() - 1;
				
				// find the last entry with this rt-value
				++pair.first;
				IndexSet::const_iterator set_iter2 = lower_bound(iso_map_iter->second.peaks_.begin(),
																												 iso_map_iter->second.peaks_.end(),
																												 pair,IndexLess());
				--set_iter2;
				last_peak_mz = (exp_it->begin() + set_iter2->second)->getMZ() + 1;
				
				//std::cout << rt<<": first peak mz "<<first_peak_mz << "\tlast peak mz "<<last_peak_mz <<std::endl;
				peak.setPosition(first_peak_mz);
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
				Idx left,right;
				left.first = distance(first,iter);
				left.second = distance(iter->begin(),raw_data_iter);
#ifdef DEBUG2D
				std::cout << "left: "<<iter->getRT()<<"\t"<<raw_data_iter->getMZ()<<std::endl;
#endif
				// consider a bit more of the signal to the right
				peak.setPosition(last_peak_mz + 1);
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
#ifdef DEBUG2D
				std::cout << "rightt: "<<iter->getRT()<<"\t"<<raw_data_iter->getMZ()<<std::endl;
#endif
				// region endpoints are stored in global vector
				OptimizationFunctions::signal2D.push_back(left);
				OptimizationFunctions::signal2D.push_back(right);
			}
#ifdef DEBUG2D
		//std::cout << "fertig"<< std::endl;
		std::cout << first_peak_mz <<"\t"<<last_peak_mz<<std::endl;
#endif
	}
    

    
  
}

#endif //OPENMS_TRANSFORMATIONS_RAW2PEAK_TWODOPTIMIZATION_H
