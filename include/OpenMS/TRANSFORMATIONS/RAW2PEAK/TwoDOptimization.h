// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
#include <iostream>
#include <fstream>
#endif

#include <vector>
#include <utility>
#include <cmath>
#include <set>

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakShape.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/PeakIndex.h>
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

 
	/**
		 @brief This class provides the two-dimensional optimization of the picked peak parameters.
			
		 Given the picked peaks, this class optimizes the peak parameters of each isotope pattern using
		 a non-linear optimization. The peaks of adjacent scans are adjusted to achieve that a peak occuring in
		 several scans has always the same m/z position. For the optimization the Levenberg-Marquardt algorithm
		 provided from the GSL is used. The optimized parameters are the m/z values,
		 the left and right width, which shall be equal for a peak in all scans,
		 and the peaks' heights.
		 
		 @todo Works only with defined types due to pointers to the data in the optimization namespace! Change that or remove templates (Alexandra)
		 
		 @htmlinclude OpenMS_TwoDOptimization.parameters
	*/
	class OPENMS_DLLAPI TwoDOptimization : public DefaultParamHandler
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
			param_.setValue("2d:tolerance_mz",tolerance_mz);
		}

		///Non-mutable access to the maximal peak distance in a cluster
		inline DoubleReal getMaxPeakDistance() const {return max_peak_distance_;}
		///Mutable access to the maximal peak distance in a cluster
		inline void setMaxPeakDistance(double max_peak_distance)
		{
			max_peak_distance_ = max_peak_distance;
			param_.setValue("2d:max_peak_distance",max_peak_distance);
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



		/**
			@brief Find two dimensional peak clusters and optimize their peak parameters
						
			@note For the peak spectra, the following meta data arrays (see DSpectrum) have to be present and have to be named just as listed here:
			- area (index:1)
			- leftWidth (index:3)
			- rightWidth (index:4)
			- peakShape (index:5) 

			@param first begin of the raw data spectra iterator range
			@param last end of the raw data spectra interator range
			@param ms_exp peak map corresponding to the raw data in the range from @p first to @p last
			@param real2D flag if the optimization should be two dimensional or on each scan separately
			@exception Exception::IllegalArgument is thrown if required meta information from peak picking is missing (area, shape, left width, right width) or if the input data is invalid in some other way
			
		*/
		template <typename InputSpectrumIterator,typename OutputPeakType>
		void optimize(InputSpectrumIterator& first,
									InputSpectrumIterator& last,
									MSExperiment< OutputPeakType >& ms_exp,bool real2D=true);


	protected:
		/// Helper struct (contains the size of an area and a raw data container)
		struct Data
		{
			std::vector<std::pair<int,int> > signal2D; 
			std::multimap<double,IsotopeCluster>::iterator iso_map_iter;
			unsigned int total_nr_peaks;
			std::map<int, std::vector<PeakIndex> > matching_peaks;
			MSExperiment<> picked_peaks;
			MSExperiment<Peak1D>::ConstIterator raw_data_first;
			OptimizationFunctions::PenaltyFactorsIntensity penalties;
			std::vector<double> positions;
			std::vector<double> signal;
		};
		
		/// stores the retention time of each isotopic cluster
		std::multimap<double,IsotopeCluster> iso_map_;
      
		/// Pointer to the current region
		std::multimap<double,IsotopeCluster>::const_iterator curr_region_;

		/// upper bound for distance between two peaks belonging to the same region
		double max_peak_distance_;

		/// threshold for the difference in the peak position of two matching peaks
		double tolerance_mz_;
      
		/// Indices of peaks in the adjacent scans matching peaks in the scan with no. ref_scan
		//		std::map<int, std::vector<MSExperiment<>::SpectrumType::Iterator > > matching_peaks_;
		std::map<int,std::vector<PeakIndex> > matching_peaks_;

		
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
       @name Functions provided to the gsl Levenberg-Marquardt
		*/
    //@{
    /// Function computing estimated signal and its deviation to the experimental signal*/
    static int residual2D_(const gsl_vector* x, void* params , gsl_vector* f);
    /// Function computing the Jacobian */
    static int jacobian2D_(const gsl_vector* x, void* params, gsl_matrix* J);
    /// Function that calls residual2D and jacobian2D*/
    static int evaluate2D_(const gsl_vector* x, void* params, gsl_vector* f, gsl_matrix* J);

		
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
														 double noise_level,
														 TwoDOptimization::Data& d);

		/// Identify matching peak in a peak cluster
		void findMatchingPeaks_(std::multimap<double, IsotopeCluster>::iterator& it,
														MSExperiment<>& ms_exp);
      
		//@}

		/// update members method from DefaultParamHandler to update the members 
		void updateMembers_();
	};


	template <typename InputSpectrumIterator,typename OutputPeakType>
	void TwoDOptimization::optimize(InputSpectrumIterator& first, InputSpectrumIterator& last, MSExperiment< OutputPeakType >& ms_exp,bool real2D)
	{
		//#define DEBUG_2D
		//check if the input maps have the same number of spectra
		if ((UInt)distance(first,last)!=ms_exp.size())
		{
			throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Error in Two2Optimization: Raw and peak map do not have the same number of spectra");
		}
		//do nothing if there are no scans
		if(ms_exp.size() == 0)
		{
			return;
		}
		//check if required meta data arrays are present (for each scan)
		for (UInt i=0; i<ms_exp.size(); ++i)
		{
			//check if enough meta data arrays are present
			if (ms_exp[i].getMetaDataArrays().size()<6)
			{
				throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Error in Two2Optimization: Not enough meta data arrays present (1:area, 5:shape, 3:left width, 4:right width)");
			}
			bool area = ms_exp[i].getMetaDataArrays()[1].getName() == "maximumIntensity";
			bool wleft = ms_exp[i].getMetaDataArrays()[3].getName() == "leftWidth";
			bool wright = ms_exp[i].getMetaDataArrays()[4].getName() == "rightWidth";
			bool shape = ms_exp[i].getMetaDataArrays()[5].getName() == "peakShape";

			if (!area || !wleft || !wright || !shape)
			{
				throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Error in Two2Optimization: One or several meta data arrays missing (1:area, 5:shape, 3:left width, 4:right width)");
			}
		}
		real_2D_ = real2D;
		typedef typename InputSpectrumIterator::value_type InputSpectrumType;
		typedef typename InputSpectrumType::value_type PeakType;
		typedef MSSpectrum<PeakType> SpectrumType;

		typename MSExperiment<OutputPeakType>::Iterator ms_exp_it = ms_exp.begin();
		typename MSExperiment<OutputPeakType>::Iterator ms_exp_it_end = ms_exp.end();
		if(ms_exp.size() == 0)
			{
				std::cout << "empty experiment"<<std::endl;
				return;
			}
		// stores the monoisotopic peaks of isotopic clusters
		std::vector<double> iso_last_scan;
		std::vector<double> iso_curr_scan;
		std::vector<std::multimap<double,IsotopeCluster>::iterator> clusters_last_scan;
		std::vector<std::multimap<double,IsotopeCluster>::iterator> clusters_curr_scan;
		std::multimap<double,IsotopeCluster>::iterator cluster_iter;
		double current_rt=ms_exp_it->getRT(),last_rt  = 0;

		// retrieve values for accepted peaks distances
		max_peak_distance_ = param_.getValue("2d:max_peak_distance");
		double tolerance_mz = param_.getValue("2d:tolerance_mz");
	
		UInt current_charge     = 0;			// charge state of the current isotopic cluster
		double mz_in_hash   = 0;			// used as reference to the current isotopic peak			
	
		// sweep through scans
		for (unsigned int curr_scan =0; ms_exp_it+curr_scan != ms_exp_it_end;++curr_scan)
			{
				unsigned int nr_peaks_in_scan = (ms_exp_it +curr_scan)->size();
				//last_rt = current_rt;
				current_rt = (ms_exp_it+curr_scan)->getRT();
				typename MSExperiment<OutputPeakType>::SpectrumType::Iterator peak_it  = (ms_exp_it+curr_scan)->begin();
				typename MSExperiment<OutputPeakType>::SpectrumType::Iterator peak_it_last  = (ms_exp_it+curr_scan)->end();

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
				MSSpectrum<PeakType> s;
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
												//std::cout << delta_mz << " "<< tolerance_mz << std::endl;
												if ( delta_mz > tolerance_mz) // check if first peak of last cluster is close enough
													{
														mz_in_hash = curr_mz; // update current hash key
			      
														// create new isotopic cluster
// #ifdef DEBUG_2D
// 														std::cout << "Last peak cluster too far, creating new cluster at "<<curr_mz << std::endl;
// #endif
														IsotopeCluster new_cluster;
														new_cluster.peaks_.charge_  = current_charge;
														new_cluster.scans_.push_back( curr_scan );					
														cluster_iter = iso_map_.insert(std::pair<double,IsotopeCluster>(mz_in_hash,new_cluster));
			      
													}
												else
													{
// //#ifdef DEBUG_2D
// 														std::cout << "Found neighbouring peak with distance (m/z) " << delta_mz << std::endl;
// //#endif
														cluster_iter = clusters_last_scan[distance(iso_last_scan.begin(),it)];

														// check whether this scan is already contained
														if(find(cluster_iter->second.scans_.begin(),cluster_iter->second.scans_.end(),curr_scan)
															 == cluster_iter->second.scans_.end())
															{
																cluster_iter->second.scans_.push_back( curr_scan );
															}
			      
// 														//#ifdef DEBUG_2D
// 														std::cout << "Cluster with " << cluster_iter->second.peaks_.size()
// 																			<< " peaks retrieved." << std::endl;
// 														//#endif
													}
			  
											}
										else // last scan did not contain any isotopic cluster
											{	
// 												//#ifdef DEBUG_2D
// 												std::cout << "Last scan was empty => creating new cluster." << std::endl;
// 												std::cout << "Creating new cluster at m/z: " << curr_mz << std::endl;
// 												//#endif
			  
												mz_in_hash = curr_mz; // update current hash key
			  
												// create new isotopic cluster
												IsotopeCluster new_cluster;
												new_cluster.peaks_.charge_  = current_charge;
												new_cluster.scans_.push_back( curr_scan );					
												cluster_iter = iso_map_.insert(std::pair<double,IsotopeCluster>(mz_in_hash,new_cluster));

											}
		  
// 										//#ifdef DEBUG_2D
// 										std::cout << "Storing found peak in current isotopic cluster" << std::endl;
// 										//#endif


		      
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
												//	std::cout << "new enter'd: "<<(peak_it+curr_peak+1)->getMZ()<<" im while"<<std::endl;
												++curr_peak;			
												if(curr_peak >= nr_peaks_in_scan-1) break;
												dist2nextpeak = (peak_it+curr_peak+1)->getMZ() -  (peak_it+curr_peak)->getMZ(); // get distance to next peak

			  
											} // end while(...)
		      
		    
		      
									} // end of if (dist2nextpeak <= max_peak_distance_)
								else
									{
										if (iso_last_scan.size() > 0)  // Did we find any isotopic cluster in the last scan?
											{
												// there were some isotopic clusters in the last scan...
												std::vector<double>::iterator it =
													searchInScan_(iso_last_scan.begin(),iso_last_scan.end(),curr_mz);
			  
												double delta_mz = fabs(*it - curr_mz);
												std::vector<double>::iterator itneu = iso_last_scan.begin();
												//												std::cout << delta_mz << " "<< tolerance_mz << std::endl;
												if ( delta_mz > tolerance_mz) // check if first peak of last cluster is close enough
													{
														mz_in_hash = curr_mz; // update current hash key
			      
														// create new isotopic cluster
// 														//#ifdef DEBUG_2D
// 														std::cout << "Last peak cluster too far, creating new cluster at "<<curr_mz << std::endl;
// 														//#endif
														IsotopeCluster new_cluster;
														new_cluster.peaks_.charge_  = current_charge;
														new_cluster.scans_.push_back( curr_scan );					
														cluster_iter = iso_map_.insert(std::pair<double,IsotopeCluster>(mz_in_hash,new_cluster));
			      
													}
												else
													{
// 														//#ifdef DEBUG_2D
// 														std::cout << "Found neighbouring peak with distance (m/z) " << delta_mz << std::endl;
// 														//#endif
														cluster_iter = clusters_last_scan[distance(iso_last_scan.begin(),it)];

														// check whether this scan is already contained
														if(find(cluster_iter->second.scans_.begin(),cluster_iter->second.scans_.end(),curr_scan)
															 == cluster_iter->second.scans_.end())
															{
																cluster_iter->second.scans_.push_back( curr_scan );
															}
			      
// 														//#ifdef DEBUG_2D
// 														std::cout << "Cluster with " << cluster_iter->second.peaks_.size()
// 																			<< " peaks retrieved." << std::endl;
// 														//#endif
													}
			  
											}
										else // last scan did not contain any isotopic cluster
											{	
// 												//#ifdef DEBUG_2D
// 												std::cout << "Last scan was empty => creating new cluster." << std::endl;
// 												std::cout << "Creating new cluster at m/z: " << curr_mz << std::endl;
// 												//#endif
			  
												mz_in_hash = curr_mz; // update current hash key
			  
												// create new isotopic cluster
												IsotopeCluster new_cluster;
												new_cluster.peaks_.charge_  = current_charge;
												new_cluster.scans_.push_back( curr_scan );					
												cluster_iter = iso_map_.insert(std::pair<double,IsotopeCluster>(mz_in_hash,new_cluster));

											}
		  
// 										//#ifdef DEBUG_2D
// 										std::cout << "Storing found peak in current isotopic cluster" << std::endl;
// 										//#endif


		      
										cluster_iter->second.peaks_.insert(std::pair<UInt,UInt>(curr_scan,curr_peak));
		      
										iso_curr_scan.push_back(  mz_in_hash );
										clusters_curr_scan.push_back(cluster_iter);
								
		    						
									}
	      
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
		//#undef DEBUG_2D
	}

	template <typename InputSpectrumIterator,typename OutputPeakType>
	void TwoDOptimization::optimizeRegions_(InputSpectrumIterator& first,
																					InputSpectrumIterator& last,
																					MSExperiment< OutputPeakType >& ms_exp)
  {

		//#define DEBUG_2D
    int counter =0;
    //#endif
    // go through the clusters
    for (std::multimap<double, IsotopeCluster>::iterator it = iso_map_.begin();
         it != iso_map_.end();
         ++it)
			{
#ifdef DEBUG_2D
				std::cout << "element: " << counter<< std::endl;
				std::cout << "mz: "<< it->first << std::endl<<"rts: ";
// 				for(unsigned int i=0;i<it->second.scans_.size();++i) std::cout << it->second.scans_[i] << "\n";
				std::cout<<std::endl<<"peaks: ";
				IndexSet::const_iterator iter = it->second.peaks_.begin();
				for(; iter != it->second.peaks_.end(); ++iter)std::cout << ms_exp[iter->first].getRT() << " "<<(ms_exp[iter->first][iter->second]).getMZ()<<std::endl;
				
//for(unsigned int i=0;i<it->second.peaks_.size();++i) std::cout << ms_exp[it->first].getRT() << " "<<(ms_exp[it->first][it->second]).getMZ()<<std::endl;
				std::cout << std::endl << std::endl;

#endif

				// prepare for optimization:
				// determine the matching peaks
				matching_peaks_.clear();
				findMatchingPeaks_(it,ms_exp);
				TwoDOptimization::Data d;
				d.penalties = penalties_;
				d.matching_peaks = matching_peaks_;
				// and the endpoints of each isotope pattern in the cluster
				getRegionEndpoints_(ms_exp,first,last,counter,400,d);

				// peaks have to be stored globally
				d.iso_map_iter = it;

				d.picked_peaks = ms_exp;
				d.raw_data_first =  first;

				int nr_diff_peaks = matching_peaks_.size();
				d.total_nr_peaks = it->second.peaks_.size();

				int nr_parameters = nr_diff_peaks*3 + d.total_nr_peaks;

				gsl_vector *start_value=gsl_vector_alloc(nr_parameters);
				gsl_vector_set_zero(start_value);

				// initialize parameters for optimization
				std::map<int, std::vector<PeakIndex> >::iterator m_peaks_it
					= d.matching_peaks.begin();
				double av_mz=0,av_lw=0,av_rw=0,avr_height=0,height;
				int peak_counter = 0;
				int diff_peak_counter =0;
				// go through the matching peaks
				for(;m_peaks_it != d.matching_peaks.end();++m_peaks_it)
					{
						av_mz=0,av_lw=0,av_rw=0,avr_height=0;
						std::vector<PeakIndex>::iterator iter_iter = (m_peaks_it)->second.begin();
						for(;iter_iter != m_peaks_it->second.end();++iter_iter)
							{
								height = (iter_iter)->getPeak(ms_exp).getIntensity();
								avr_height += height;
								av_mz += (iter_iter)->getPeak(ms_exp).getMZ() * height;
								av_lw += ms_exp[(iter_iter)->spectrum].getMetaDataArrays()[3][(iter_iter)->peak]* height; //left width
								av_rw +=	ms_exp[(iter_iter)->spectrum].getMetaDataArrays()[4][(iter_iter)->peak]* height;//right width
								gsl_vector_set(start_value,peak_counter,height);
								++peak_counter;
							}
						gsl_vector_set(start_value,d.total_nr_peaks+3*diff_peak_counter, av_mz/avr_height);
						gsl_vector_set(start_value,d.total_nr_peaks+3*diff_peak_counter+1, av_lw/avr_height);
						gsl_vector_set(start_value,d.total_nr_peaks+3*diff_peak_counter+2, av_rw/avr_height);
						++diff_peak_counter;
					}

#ifdef DEBUG_2D
				std::cout << "----------------------------\n\nstart_value: "<<std::endl;
				for(unsigned int k=0;k<start_value->size;++k)
					{
						std::cout << gsl_vector_get(start_value,k)<<std::endl;
					}
#endif
				int num_positions = 0;
				for(unsigned int i=0; i<d.signal2D.size(); i+=2)
					{
						num_positions += (d.signal2D[i+1].second - d.signal2D[i].second +1);
#ifdef DEBUG_2D
						std::cout << d.signal2D[i+1].second << " - "<< d.signal2D[i].second <<" +1 "<< std::endl;
#endif

					}
#ifdef DEBUG_2D
				std::cout << "num_positions : "<< num_positions << std::endl;
#endif
				// The gsl algorithms require us to provide function pointers for the evaluation of
				// the target function.
				gsl_multifit_function_fdf fit_function;
				fit_function.f   = (int (*)(const gsl_vector * x, void * params, gsl_vector * f))&OpenMS::TwoDOptimization::residual2D_;
				fit_function.df  = (int (*)(const gsl_vector * x, void * params, gsl_matrix * J))&OpenMS::TwoDOptimization::jacobian2D_;
				fit_function.fdf = (int (*)(const gsl_vector * x, void * params, gsl_vector * f, gsl_matrix * J))&OpenMS::TwoDOptimization::evaluate2D_;
				// dunno why, but gsl crashes when n is smaller than p!!!!!!!! ????
				fit_function.n   = std::max(num_positions+1,(int)(nr_parameters));
				fit_function.p   = nr_parameters;
				fit_function.params = &d;
#ifdef DEBUG_2D
				std::cout << "fit_function.n "<<fit_function.n
									<< "\tfit_function.p "<<fit_function.p<<std::endl;
#endif
				const gsl_multifit_fdfsolver_type *type = gsl_multifit_fdfsolver_lmsder;

				gsl_multifit_fdfsolver *fit=gsl_multifit_fdfsolver_alloc(type,
																																 std::max(num_positions+1,(int)(nr_parameters)),
																																 nr_parameters);

				gsl_multifit_fdfsolver_set(fit, &fit_function, start_value);



				// initial norm
#ifdef DEBUG_2D
				std::cout << "Before optimization: ||f|| = " << gsl_blas_dnrm2(fit->f) << std::endl;
#endif
				// Iteration
				int iteration = 0;
				int status;

				do
					{
						iteration++;
						status = gsl_multifit_fdfsolver_iterate(fit);
#ifdef DEBUG_2D
						std::cout << "Iteration " << iteration << "; Status " << gsl_strerror(status) << "; " << std::endl;
						std::cout << "||f|| = " << gsl_blas_dnrm2(fit->f) << std::endl;
						std::cout << "Number of parms: " << nr_parameters << std::endl;
						std::cout << "Delta: " << gsl_blas_dnrm2(fit->dx) << std::endl;
#endif

						status = gsl_multifit_test_delta(fit->dx, fit->x, eps_abs_, eps_rel_);
						if (status != GSL_CONTINUE)
							break;

					}
				while (status == GSL_CONTINUE && iteration < max_iteration_);

#ifdef DEBUG_2D
				std::cout << "Finished! No. of iterations" << iteration << std::endl;
				std::cout << "Delta: " << gsl_blas_dnrm2(fit->dx) << std::endl;
				double chi = gsl_blas_dnrm2(fit->f);
				std::cout << "After optimization: || f || = " << gsl_blas_dnrm2(fit->f) << std::endl;
				std::cout << "chisq/dof = " << pow(chi, 2.0) / (num_positions - nr_parameters);


				std::cout << "----------------------------------------------\n\nnachher"<<std::endl;
				for(unsigned int k=0;k<fit->x->size;++k)
					{
						std::cout << gsl_vector_get(fit->x,k)<<std::endl;
					}
#endif
				int peak_idx =0;
				std::map<int, std::vector<PeakIndex> >::iterator it
					=d.matching_peaks.begin();
				for(; it != d.matching_peaks.end(); ++it)
					{
						int i = distance(d.matching_peaks.begin(),it);
						for(unsigned int j=0;j<it->second.size();++j)
							{

#ifdef DEBUG_2D
								std::cout << "pos: "<<it->second[j].getPeak(ms_exp).getMZ()<<"\nint: "<<it->second[j].getPeak(ms_exp).getIntensity()
													<<"\nlw: "<<it->second[j].getSpectrum(ms_exp).getMetaDataArrays()[3][it->second[j].peak]
													<<"\nrw: "<<it->second[j].getSpectrum(ms_exp).getMetaDataArrays()[4][it->second[j].peak] << "\n";

#endif

								ms_exp[it->second[j].spectrum][it->second[j].peak].setMZ(gsl_vector_get(fit->x,d.total_nr_peaks+3*i));

								ms_exp[it->second[j].spectrum][it->second[j].peak].setIntensity(gsl_vector_get(fit->x,peak_idx));
								ms_exp[it->second[j].spectrum].getMetaDataArrays()[3][it->second[j].peak] =
									gsl_vector_get(fit->x,d.total_nr_peaks+3*i+1);
								ms_exp[it->second[j].spectrum].getMetaDataArrays()[4][it->second[j].peak] =
									gsl_vector_get(fit->x,d.total_nr_peaks+3*i+2);


#ifdef DEBUG_2D
								std::cout << "pos: "<<it->second[j].getPeak(ms_exp).getMZ()<<"\nint: "<<it->second[j].getPeak(ms_exp).getIntensity()
													<<"\nlw: "<<it->second[j].getSpectrum(ms_exp).getMetaDataArrays()[3][it->second[j].peak]
													<<"\nrw: "<<it->second[j].getSpectrum(ms_exp).getMetaDataArrays()[4][it->second[j].peak] << "\n";

// 								std::cout << "pos: "<<it->second[j]->getMZ()<<"\nint: "<<it->second[j]->getIntensity()
// 													<<"\nlw: "<<it->second[j]->getLeftWidthParameter()
// 													<<"\nrw: "<<it->second[j]->getRightWidthParameter() << "\n";

#endif

								++peak_idx;


							}
					}

				gsl_multifit_fdfsolver_free(fit);
				gsl_vector_free(start_value);
				++counter;
			}// end for
		//#undef DEBUG_2D
  }


  template <typename InputSpectrumIterator,typename OutputPeakType>
	void TwoDOptimization::optimizeRegionsScanwise_(InputSpectrumIterator& first,
																									InputSpectrumIterator& last,
																									MSExperiment< OutputPeakType >& ms_exp)
  {
		typedef typename InputSpectrumIterator::value_type InputSpectrumType;
		typedef typename InputSpectrumType::value_type PeakType;
		
		
    int counter =0;
		TwoDOptimization::Data d;
    d.picked_peaks = ms_exp;
    d.raw_data_first =  first;

		std::cout << "richtig hier" << std::endl;
    struct OpenMS::OptimizationFunctions::PenaltyFactors penalties;


    DataValue dv = param_.getValue("penalties:position");
    if (dv.isEmpty() || dv.toString() == "")
      penalties.pos = 0.;
    else
      penalties.pos = (float)dv;

    dv = param_.getValue("penalties:left_width");
    if (dv.isEmpty() || dv.toString() == "")
      penalties.lWidth = 1.;
    else
      penalties.lWidth = (float)dv;

    dv = param_.getValue("penalties:right_width");
    if (dv.isEmpty() || dv.toString() == "")
      penalties.rWidth = 1.;
    else
      penalties.rWidth = (float)dv;
#ifdef DEBUG_2D
		std::cout << penalties.pos << " "
							<< penalties.rWidth << " "
							<< penalties.lWidth << std::endl;
#endif
// 		MSExperiment<Peak1D >::const_iterator help = first;
// 		//		std::cout << "\n\n\n\n---------------------------------------------------------------";
// 		while(help!=last)
// 			{
// 				//	std::cout<<help->getRT()<<std::endl;
// 				++help;
// 			}
		//		std::cout << "---------------------------------------------------------------\n\n\n\n";
		
    unsigned int max_iteration;
    dv = param_.getValue("iterations");
    if (dv.isEmpty() || dv.toString() == "")
      max_iteration = 15;
    else
      max_iteration = (unsigned int)dv;

    double eps_abs;
    dv = param_.getValue("delta_abs_error");
    if (dv.isEmpty() || dv.toString() == "")
      eps_abs = 1e-04f;
    else
      eps_abs = (double)dv;

    double eps_rel;
    dv = param_.getValue("delta_rel_error");
    if (dv.isEmpty() || dv.toString() == "")
      eps_rel = 1e-04f;
    else
      eps_rel = (double)dv;

    std::vector<PeakShape> peak_shapes;


    // go through the clusters
    for (std::multimap<double, IsotopeCluster>::iterator it = iso_map_.begin();
         it != iso_map_.end();
         ++it)
			{
				d.iso_map_iter = it;
				//#ifdef DEBUG_2D
				std::cerr << "element: " << counter<< std::endl;
				std::cerr << "mz: "<< it->first <<std::endl<<"rts: ";
				for(unsigned int i=0;i<it->second.scans_.size();++i) std::cerr << it->second.scans_[i] << "\n";
				std::cerr<<std::endl<<"peaks: ";
				IndexSet::const_iterator iter = it->second.peaks_.begin();
				for(; iter != it->second.peaks_.end(); ++iter)std::cerr << ms_exp[iter->first].getRT() << " "<<(ms_exp[iter->first][iter->second]).getMZ()<<std::endl;
				
//for(unsigned int i=0;i<it->second.peaks_.size();++i) std::cout << ms_exp[it->first].getRT() << " "<<(ms_exp[it->first][it->second]).getMZ()<<std::endl;
				std::cerr << std::endl << std::endl;

				//#endif
				// prepare for optimization:
				// determine the matching peaks
				// and the endpoints of each isotope pattern in the cluster
			
				getRegionEndpoints_(ms_exp,first,last,counter,400,d);
				
				
				
				unsigned int idx = 0;
				for(unsigned int i=0; i < d.signal2D.size()/2; ++i)
					{
						OptimizationFunctions::positions_.clear();
						OptimizationFunctions::signal_.clear();

						MSExperiment<Peak1D >::SpectrumType::const_iterator	ms_it =
							(d.raw_data_first + d.signal2D[2*i].first)->begin()+d.signal2D[2*i].second;
						int size = distance(ms_it,(d.raw_data_first + d.signal2D[2*i].first)->begin()+d.signal2D[2*i+1].second);
						OptimizationFunctions::positions_.reserve(size);
						OptimizationFunctions::signal_.reserve(size);

						while(ms_it != (d.raw_data_first + d.signal2D[2*i].first)->begin()+d.signal2D[2*i+1].second)
							{
								OptimizationFunctions::positions_.push_back(ms_it->getMZ());
								OptimizationFunctions::signal_.push_back(ms_it->getIntensity());
								++ms_it;
							}


						Idx pair;
						pair.first =  d.iso_map_iter->second.peaks_.begin()->first + idx;

						IndexSet::const_iterator set_iter = lower_bound(d.iso_map_iter->second.peaks_.begin(),
																														d.iso_map_iter->second.peaks_.end(),
																														pair,IndexLess());


						// find the last entry with this rt-value
						++pair.first;
						IndexSet::const_iterator set_iter2 = lower_bound(d.iso_map_iter->second.peaks_.begin(),
																														 d.iso_map_iter->second.peaks_.end(),
																														 pair,IndexLess());

						while(set_iter != set_iter2)
							{
								const UInt peak_index = set_iter->second;
								const MSSpectrum<>& spec = ms_exp[set_iter->first];
								PeakShape shape(spec[peak_index].getIntensity(),
																spec[peak_index].getMZ(),
																spec.getMetaDataArrays()[3][peak_index], //left width
																spec.getMetaDataArrays()[4][peak_index], //right width
																spec.getMetaDataArrays()[1][peak_index], //area
																std::vector<Peak1D>::iterator(),
																std::vector<Peak1D>::iterator(),
																PeakShape::Type(Int(spec.getMetaDataArrays()[5][peak_index]))); //shape
								peak_shapes.push_back(shape);
								++set_iter;
							}
						//#ifdef DEBUG_2D
						std::cout << "rt "
											<<(d.raw_data_first + d.signal2D[2*i].first)->getRT()
											<<"\n";
						//#endif
						OptimizePick opt(penalties,max_iteration,eps_abs,eps_rel);
						//#ifdef DEBUG_2D
						std::cout << "vorher\n";

						for(unsigned int p =0; p < peak_shapes.size();++p)
							{
								std::cout << peak_shapes[p].mz_position <<"\t" << peak_shapes[p].height
													<<"\t" << peak_shapes[p].left_width <<"\t" << peak_shapes[p].right_width  <<std::endl;
							}
						//#endif
						opt.optimize(peak_shapes);
						//#ifdef DEBUG_2D
						std::cout << "nachher\n";
						for(unsigned int p =0; p < peak_shapes.size();++p)
							{
								std::cout << peak_shapes[p].mz_position <<"\t" << peak_shapes[p].height
													<<"\t" << peak_shapes[p].left_width <<"\t" << peak_shapes[p].right_width  <<std::endl;
							}
						//#endif
						std::sort(peak_shapes.begin(),peak_shapes.end(),PeakShape::PositionLess());
						pair.first =  d.iso_map_iter->second.peaks_.begin()->first + idx;

						set_iter = lower_bound(d.iso_map_iter->second.peaks_.begin(),
																	 d.iso_map_iter->second.peaks_.end(),
																	 pair,IndexLess());
						unsigned int i=0;
						while(i < peak_shapes.size())
							{
								MSSpectrum<>& spec = ms_exp[set_iter->first];
								spec[set_iter->second].setMZ(peak_shapes[i].mz_position);
								spec[set_iter->second].setIntensity(peak_shapes[i].height);
								spec.getMetaDataArrays()[3][set_iter->second] = peak_shapes[i].left_width;
								spec.getMetaDataArrays()[4][set_iter->second] = peak_shapes[i].right_width;

								++set_iter;
								++i;
							}
						++idx;
						peak_shapes.clear();
					}

				++counter;
			}
  }

	template<typename InputSpectrumIterator,typename OutputPeakType>
	void TwoDOptimization::getRegionEndpoints_(MSExperiment<OutputPeakType>& exp,
																						 InputSpectrumIterator& first,
																						 InputSpectrumIterator& last,
																						 unsigned int iso_map_idx,
																						 double noise_level,
																						 TwoDOptimization::Data& d)
	{
		d.signal2D.clear();
		typedef typename InputSpectrumIterator::value_type InputExperimentType;
		typedef typename InputExperimentType::value_type InputPeakType;
		typedef std::multimap<double,IsotopeCluster> MapType;

		double rt, first_peak_mz, last_peak_mz;
		
		//MSSpectrum<InputPeakType> spec;
		typename MSExperiment<InputPeakType>::SpectrumType spec;
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
				typename MSExperiment<InputPeakType>::SpectrumType::const_iterator raw_data_iter
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
				left.second = raw_data_iter-iter->begin();
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
				right.second = raw_data_iter-iter->begin();
#ifdef DEBUG2D
				std::cout << "rightt: "<<iter->getRT()<<"\t"<<raw_data_iter->getMZ()<<std::endl;
#endif
				// region endpoints are stored in global vector
				d.signal2D.push_back(left);
				d.signal2D.push_back(right);
			}
#ifdef DEBUG2D
		//std::cout << "fertig"<< std::endl;
		std::cout << first_peak_mz <<"\t"<<last_peak_mz<<std::endl;
#endif
	}
    

    
  
}

#endif //OPENMS_TRANSFORMATIONS_RAW2PEAK_TWODOPTIMIZATION_H
