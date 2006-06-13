
// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: DPeakPickerCWT.h,v 1.49 2006/06/02 14:42:37 elange Exp $
// $Author: elange $
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_DPEAKPICKERCWT_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_DPEAKPICKERCWT_H

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/DPeakPicker.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakShape.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/DExtractSignalRegions.h>
#include <OpenMS/FILTERING/NOISEESTIMATION/DSignalToNoiseEstimatorWindowing.h>

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/ContinuousWaveletTransformNumIntegration.h>

#ifdef GSL_DEF
# include <OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePick.h>
#endif

#include <math.h>
#include <vector>
#include <algorithm>

#ifdef DEBUG_PEAK_PICKING
#include<iostream>
#include<fstream>
#endif

#define OPENMS_TRANSFORMATIONS_RAW2PEAK_DEFAULT_FILE "TRANSFORMATIONS/RAW2PEAK/PeakPicking.xml"

//#define DEBUG_PEAK_PICKING
#undef DEBUG_PEAK_PICKING

namespace OpenMS
{
  /**
     @defgroup PeakPickingCWT PeakPickingCWT

     @brief Classes for the peak picking method as described by Lange et al. (2006) Proc. PSB-06.

     This module contains all important classes that are involved in the peak picking method.

     @ingroup Raw2Peak
  */

  /**
     @brief This class implements a peak picking algorithm using wavelet techniques.

     This peak picking algorithm uses the continuous wavelet transform of a raw data signal to detect mass peaks.
     Afterwards a given asymmetric peak function is fitted to the raw data and important peak parameters (e.g. fwhm)
     are extracted.
     In an optional step these parameters can be optimized using a non-linear opimization method.

     @ingroup PeakPickingCWT

     @todo write test

  */

  template <Size D = 1, 
	    typename MapType = MSExperiment<DRawDataPoint<1> >, 
	    typename MapTypeOut = MSExperiment<DPickedPeak<1> > >
  class DPeakPickerCWT : public DPeakPicker<D, MapType, MapTypeOut>
  {
  public:
    /** @name Type definitions
     */
    //@{
    typedef DPeakArrayNonPolymorphic< D,DRawDataPoint<D> > RawData;
    ///
    typedef typename RawData::const_iterator RawConstIterator;
    ///
    typedef typename RawData::iterator RawIterator;
    ///
    typedef DPickedPeak<D> OutputPeak;
    ///
    typedef DPeakArray<D, DPickedPeak<D> > PeakData;
    ///
    typedef typename PeakData::iterator PeakIterator;
    ///
    typedef typename PeakData::iterator PeakBackInserter;
    ///
    typedef std::vector<RawConstIterator> IteratorVector;
    ///
    typedef typename IteratorVector::const_iterator IteratorVectorConstIterator;
    ///
    typedef typename MapType::const_iterator SpectrumConstIterator;
    ///
    typedef typename MapType::SpectrumType SpectrumType;
    ///
    typedef typename MapTypeOut::SpectrumType SpectrumTypeOut;
    ///
    typedef typename SpectrumType::const_iterator PointConstIterator;
    ///
    //@}

    using DPeakPicker<D, MapType, MapTypeOut>::mz_dim_;
    using DPeakPicker<D, MapType, MapTypeOut>::rt_dim_;
    using DPeakPicker<D, MapType, MapTypeOut>::peak_bound_;
    using DPeakPicker<D, MapType, MapTypeOut>::peak_bound_ms2_level_;
    using DPeakPicker<D, MapType, MapTypeOut>::signal_to_noise_;
    using DPeakPicker<D, MapType, MapTypeOut>::peaks_;
    using DPeakPicker<D, MapType, MapTypeOut>::ms_exp_peaks_;
    using DPeakPicker<D, MapType, MapTypeOut>::param_;


    /** @name Constructors and Destructor
     */
    //@{
    DPeakPickerCWT();
    ///
    DPeakPickerCWT(const Param& parameters);
    ///
    DPeakPickerCWT(const String& filename);
    ///
    void init();
    ///
    DPeakPickerCWT(const DPeakPickerCWT& pp)
      : DPeakPicker<D, MapType, MapTypeOut>(pp),
        scale_(pp.scale_),
        peak_bound_cwt_(pp.peak_bound_cwt_),
        peak_bound_ms2_level_cwt_(pp.peak_bound_ms2_level_cwt_),
        radius_(pp.radius_),
        peak_asymm_bound_(pp.peak_asymm_bound_),
        peak_corr_bound_(pp.peak_corr_bound_),
        peak_fwhm_bound_(pp.peak_fwhm_bound_),
        noise_level_(pp.noise_level_),
        optimization_(pp.optimization_),
        num_integration_(pp.num_integration_)
    {}
    ///
    virtual ~DPeakPickerCWT();
    ///
    //@}

    /** @name Assignment
     */
    //@{
    inline DPeakPickerCWT& operator=(const DPeakPickerCWT& pp)
    {
      mz_dim_ = pp.mz_dim_;
      rt_dim_ = pp.rt_dim_;
      peak_bound_ = pp.peak_bound_;
      signal_to_noise_ = pp.signal_to_noise_;
      peaks_ = pp.peaks_; // TODO copy the peaks
      param_ = pp.param_;
      scale_ = pp.scale_;
      peak_bound_cwt_ = pp.peak_bound_cwt_;
      peak_bound_ms2_level_cwt_ = pp.peak_bound_ms2_level_cwt_;
      radius_ = pp.radius_;
      peak_asymm_bound_ = pp.peak_asymm_bound_;
      peak_corr_bound_ = pp.peak_corr_bound_;
      peak_fwhm_bound_ = pp.peak_fwhm_bound_;
      noise_level_ = pp.noise_level_;
      optimization_ = pp.optimization_;
      num_integration_ = pp.num_integration_;

      return *this;
    }
    //@}


    /** Accessors
     */
    //@{
    /// Mutable access to the noise level
    virtual void setPeakBound(const float peak_bound) 
    { 
      peak_bound_ = peak_bound; 
      calculatePeakBoundCWT_();
    }

    /// Non-mutable access to the vector with peak shapes
    inline const std::vector<PeakShape>& getPeakShapes() const { return peak_shapes_; }
    /// Mutable access to the vector with peak shapes
    inline std::vector<PeakShape>& getPeakShapes() { return peak_shapes_;  }
    /// Mutable access to the vector with peak shapes
    inline void setPeakShapes(const std::vector<PeakShape>&  peak_shapes) { peak_shapes_ = peak_shapes; }

    /// Non-mutable access to the wavelet transform
    inline const ContinuousWaveletTransform<D>& getWaveletTransform() const { return *wt_; }
    /// Mutable access to the wavelet transform
    inline ContinuousWaveletTransform<D>& getWaveletTransform() { return *wt_;  }
    /// Mutable access to the wavelet transform
    inline void setWaveletTransform(const ContinuousWaveletTransform<D>& wt) { wt_ = &wt; }

    /// Non-mutable access to the search radius for the peak maximum
    inline const int getSearchRadius() const { return radius_; }
    /// Mutable access to the search radius for the peak maximum
    inline int getSearchRadius() { return radius_;  }
    /// Mutable access to the search radius for the peak maximum
    inline void setSearchRadius(const int radius) { radius_ = radius; }

    /// Non-mutable access to the scale of the wavelet transform
    inline const float getWaveletScale() const { return scale_; }
    /// Mutable access to the scale of the wavelet transform
    inline float getWaveletScale() { return scale_; }
    /// Mutable access to the scale of the wavelet transform
    inline void setWaveletScale(const float scale) { scale_ = scale; }

    /// Non-mutable access to the peak bound in the wavelet transform
    inline const float getPeakBoundCWT() const { return peak_bound_cwt_; }
    /// Mutable access to the peak bound in the wavelet transform
    inline float getPeakBoundCWT() { return peak_bound_cwt_; }
    /// Mutable access to the peak bound in the wavelet transform
    inline void setPeakBoundCWT(const float peak_bound_cwt) { peak_bound_cwt_ = peak_bound_cwt; }

    /// Non-mutable access to the peak bound in the wavelet transform
    inline const float getPeakBoundMs2LevelCWT() const { return peak_bound_ms2_level_cwt_; }
    /// Mutable access to the peak bound in the wavelet transform
    inline float getPeakBoundMs2LevelCWT() { return peak_bound_ms2_level_cwt_; }
    /// Mutable access to the peak bound in the wavelet transform
    inline void setPeakBoundMs2LevelCWT(const float peak_bound_ms2_level_cwt) { peak_bound_ms2_level_cwt_ = peak_bound_ms2_level_cwt; }

    /// Non-mutable access to the value for the maximum peak asymmetry
    inline const float getPeakAsymBound() const { return peak_asymm_bound_; }
    /// Mutable access to the value for the maximum peak asymmetry
    inline float getPeakAsymBound() { return peak_asymm_bound_; }
    /// Mutable access to the maximum peak asymmetry
    inline void setAsymBound(const float peak_asymm_bound) { peak_asymm_bound_ = peak_asymm_bound; }

    /// Non-mutable access to the minimum peak correlation coefficient
    inline const float getPeakCorrBound() const { return peak_corr_bound_; }
    /// Mutable access to the minimum peak correlation coefficient
    inline float getPeakCorrBound() { return peak_corr_bound_; }
    /// Mutable access to the minimum peak correlation coefficient
    inline void setPeakCorrBound(const float peak_corr_bound) { peak_corr_bound_ = peak_corr_bound; }

    /// Non-mutable access to the value of minimum peak full-width-at-half-max
    inline const float getPeakFwhmBound() const { return peak_fwhm_bound_; }
    /// Mutable access to the value of minimum peak full-width-at-half-max
    inline float getPeakFwhmBound() { return peak_corr_bound_; }
    /// Mutable access to the value of minimum peak full-width-at-half-max
    inline void setPeakFwhmBound(const float peak_fwhm_bound) { peak_fwhm_bound_ = peak_fwhm_bound; }

    /// Non-mutable access to the noise level
    inline const float getNoiseLevel() const { return noise_level_; }
    /// Mutable access to the noise level
    inline float getNoiseLevel() { return noise_level_; }
    /// Mutable access to the noise level
    inline void setNoiseLevel(const float noise_level) { noise_level_ = noise_level; }

    /// Non-mutable access to the optimization switch
    inline const bool getOptimizationValue() const { return optimization_; }
    /// Mutable access to the optimization switch
    inline float getOptimizationValue() { return optimization_; }
    /// Mutable access to the optimization switch
    inline void setOptimizationValue(const bool optimization) { optimization_ = optimization; }

    /// Non-mutable access to the numerical integration switch
    inline const bool getNumOptValue() const { return num_integration_; }
    /// Mutable access to the numerical integration switch
    inline float getNumOptValue() { return num_integration_; }
    /// Mutable access to the numerical integration switch
    inline void setNumOptValue(const bool num_integration) { num_integration_ = num_integration; }
    //@}

    /** @name Functions to manage the workflow of peak picking
     */
    //@{
    virtual void pick(const MapType& ms_exp_raw);
	    
    void pick(const MapType& ms_exp_raw, DPeakPickerCWT<1, MapType, MapTypeOut> const*);
		
    void pick(const MapType& ms_exp_raw, DPeakPickerCWT<2, MapType, MapTypeOut> const*);

    /// To enhance the runtime of the peak picking algorithm every mass spectrum is decomposed in boxes
    /// of a typical length of 10 Da. The start and endpoints of each box are stored in an iterator vector.
    virtual void pick(RawConstIterator first, RawConstIterator last, std::back_insert_iterator<PeakData >& output);
    ///
    /// For each box the peak picking algorithm is performed
    template<typename ContainerType>
    void pick(typename std::vector<typename ContainerType::const_iterator>::const_iterator first,
              typename std::vector<typename ContainerType::const_iterator>::const_iterator last,
              DSignalToNoiseEstimatorWindowing<D,ContainerType>& sne,
              bool ms_experiment = false, double current_rt = 0, int ms_level = 1)
    {

#ifdef DEBUG_PEAK_PICKING
      std::cout << "DPeakPickerCWT<D>::pick" << std::endl;
#endif
      std::cout << "****************** PICK ******************" << std::endl;
      RawData raw_peak_array(0);
      // The begin and end iterators for the split regions
      RawIterator it_pick_begin;
      RawIterator it_pick_end;
      // The maximum position in the data (peak)
      RawIterator it_max_pos;
      IteratorVectorConstIterator split_first=first;

      /** We need the spacing, so we take the first and last point of the first split **/
      RawConstIterator first_first, first_last;
      first_first = *first;
      first_last  = *(first+1)-1;

      // temporary peakshape vector
      std::vector<PeakShape> peak_shapes;
      std::vector<double> peak_endpoints;

      //???
      unsigned int counter_cwt=0;

      // For every splitted part of the array apply the cwt and detect peak positions
      while ((split_first+1) < last)
	{

	  // try to find peaks only if the split contains more than 3 data points
	  if(distance(*split_first, *(split_first+1)) > 3 )
	    {
#ifdef GSL_DEF
	      if (optimization_)
		{
		  OptimizationFunctions::positions_.clear();
		  OptimizationFunctions::signal_.clear();
		}
#endif

	      // copy the raw data of the next split intervall into a DPeakArrayNonPolymorphic<DRawDataPoint<D> >
	      raw_peak_array.assign(*split_first, *(split_first+1));

	      //   std::cout << "*************** Raw_Peak_Array *****************" << std::endl;
	      // 	  for (int i = 0; i < raw_peak_array.size(); ++i)
	      // 	    {
	      // 	      std::cout << raw_peak_array[i].getPosition()[mz_dim_] << ' ' << raw_peak_array[i].getIntensity() << std::endl;
	      // 	    } 
	      // 	  std::cout << "*************** Raw_Peak_Array *****************" << std::endl;

	      it_pick_begin = raw_peak_array.begin();
	      it_pick_end   = raw_peak_array.end();

#ifdef GSL_DEF
	      if (optimization_)
		{
		  unsigned int l=raw_peak_array.size();

		  OptimizationFunctions::positions_.resize(l);
		  OptimizationFunctions::signal_.resize(l);

		  for (unsigned int i = 0; i < l ;++i)
		    {
		      OptimizationFunctions::positions_[i] = raw_peak_array[i].getPosition()[mz_dim_];
		      OptimizationFunctions::signal_[i] = raw_peak_array[i].getIntensity();
		    }
		}
#endif

#ifdef DEBUG_PEAK_PICKING
	      if (D==2)
		{
		  std::cout << "SPLIT " << it_pick_begin->getPosition()[rt_dim_] << " " ;
		}
	      std::cout << it_pick_begin->getPosition()[mz_dim_] << " UNTIL " << (it_pick_end-1)->getPosition()[mz_dim_] << std::endl;
#endif

	      //???
	      std::cout << "Counter: " << counter_cwt << std::endl;
	      std::cout << "SPLIT :\n" 
			<< "MZ intervall " <<  fabs(it_pick_begin->getPosition()[mz_dim_] - (it_pick_end-1)->getPosition()[mz_dim_]) << std::endl;
	      counter_cwt = 0; 
	  
	      unsigned int number_of_peaks = 0;
	      do 
		{  
		  number_of_peaks = 0;
		  int peak_left_index, peak_right_index;
		  // compute the cwt of the split
		  wt_->transform(it_pick_begin, it_pick_end,1.);
		  ++counter_cwt;  
		  std::cout << "COMPUTE THE CWT  " << counter_cwt << " TIME " << std::endl;
		  
		  PeakArea_ area;
		  bool centroid_fit=false;
		  bool regular_endpoints=true;

		  // search for maxima in the cwt
		  int direction=1;
		  int distance_from_scan_border = 0;
		  while ((distance(it_pick_begin, it_pick_end) > 3) && getMaxPosition_(it_pick_begin, it_pick_end, wt_, area, distance_from_scan_border, ms_level, direction))
		    {
		      //search for the endpoints of the peak
		      regular_endpoints=getPeakEndPoints_(it_pick_begin,   it_pick_end,area,
							  peak_left_index, peak_right_index);
		      getPeakCentroid_(area);

		      //if the endpoints meet the claim of minimal width
		      if (regular_endpoints)
			{
#ifdef DEBUG_PEAK_PICKING
			  std::cout << "The endpoints are "
				    << area.left->getPosition()[mz_dim_]
				    << " and "
				    << area.right->getPosition()[mz_dim_]
				    << std::endl;
#endif
			  // determine the best fitting lorezian or sech2 function
			  PeakShape shape = fitPeakShape_(area,centroid_fit,ms_experiment,current_rt);

			  // Use the centroid for Optimization
			  shape.mz_position=area.centroid_position[mz_dim_];

			  // TEST!!!!!
			  if ( (shape.r_value > peak_corr_bound_) && ((sne.getSignalToNoise(area.max)) >= signal_to_noise_))
			    {
			      //  shape.getSymmetricMeasure();
			      shape.signal_to_noise = sne.getSignalToNoise(area.max);
			      peak_shapes.push_back(shape);
			      peak_endpoints.push_back(area.left->getPosition()[mz_dim_]);
			      peak_endpoints.push_back(area.right->getPosition()[mz_dim_]);
			      ++number_of_peaks;
			    }

			  else
			    {
#ifdef DEBUG_PEAK_PICKING
			      std::cout << "Bad fitting peak "<< std::endl;
#endif

			    }
			}

		      // remove the peak from the signal
		      // TODO: does this work as expected???
		      for (RawIterator pi=area.left; pi!=area.right+1; pi++)
			{
			  pi->setIntensity(0.);
			}

		      // and transform again
		      //   wt_->transform(it_pick_begin, it_pick_end,1.);

		      // search for the next peak
		      it_pick_begin = area.right;
		      distance_from_scan_border = distance(raw_peak_array.begin(),it_pick_begin);

		      // direction*=-1;
		    } //end while (getMaxPosition_(it_pick_begin, it_pick_end, wt_, area, distance_from_scan_border, ms_level, direction))
		  it_pick_begin = raw_peak_array.begin();
		  std::cout << "FOUND " << number_of_peaks << " PEAKS " << std::endl;
		}
	      while (number_of_peaks != 0);
		
	      // start the nonlinear optimization for all peaks in split
#ifdef DEBUG_PEAK_PICKING
		std::cout << "Try the optimization run... with " << peak_shapes.size() << std::endl;
#endif
	      
	      if (peak_shapes.size() > 0)
		{
#ifdef GSL_DEF
		  
		  if (optimization_)
		    {
		      struct OpenMS::OptimizationFunctions::PenaltyFactors penalties;
		      
		      DataValue dv = param_.getValue("Optimization:Penalties:Position");
		      if (dv.isEmpty() || dv.toString() == "") penalties.pos = 0.;
		      else penalties.pos = (float)dv;
		      
		      dv = param_.getValue("Optimization:Penalties:LeftWidth");
		      if (dv.isEmpty() || dv.toString() == "") penalties.lWidth = 1.;
		      else penalties.lWidth = (float)dv;
		      
		      dv = param_.getValue("Optimization:Penalties:RightWidth");
		      if (dv.isEmpty() || dv.toString() == "") penalties.rWidth = 1.;
		      else penalties.rWidth = (float)dv;
		      
		      unsigned int max_iteration;
		      dv = param_.getValue("Optimization:Iterations");
		      if (dv.isEmpty() || dv.toString() == "") max_iteration = 15;
		      else max_iteration = (unsigned int)dv;

		      double eps_abs;
		      dv = param_.getValue("Optimization:DeltaAbsError");
		      if (dv.isEmpty() || dv.toString() == "") eps_abs = 1e-04f;
		      else eps_abs = (double)dv;

		      double eps_rel;
		      dv = param_.getValue("Optimization:DeltaRelError");
		      if (dv.isEmpty() || dv.toString() == "") eps_rel = 1e-04f;
		      else eps_rel = (double)dv;

		      OptimizePick opt(penalties,max_iteration,eps_abs,eps_rel);

		      opt.optimize(peak_shapes);

		      // compute the new correlation coefficients
		      for (unsigned int i=0, j=0; i < peak_shapes.size(); ++i, j+=2)
			{
			  peak_shapes[i].r_value=opt.correlate(peak_shapes[i],peak_endpoints[j], peak_endpoints[j+1]);
			}
		    }
#endif

		  int last = peak_shapes_.size();
		  number_of_peaks += last;
		  peak_shapes_.resize(peak_shapes_.size() + peak_shapes.size());
		  std::copy(peak_shapes.begin(),peak_shapes.end(),peak_shapes_.begin()+last);

		  peak_shapes.clear();
		  peak_endpoints.clear();
		} // end if (optimization_)
	    } // end if (peak_shapes.size() > 0)
	  split_first+=2;
	} // end while ((split_first+1) < last)
    }
    //@}

    // compute the resolution of the data (spacing in Th)
    inline double computeResolution(const MapType& ms_raw)
    {
      int max_index = 0;
      for (unsigned int i=0; i<ms_raw.size(); ++i)
      {
        if (ms_raw[i].size() > ms_raw[max_index].size())
	  {
	    max_index = i;
	  }
      }

      float spacing = 1000;
      MSSpectrum<DRawDataPoint<1> >::const_iterator first = ms_raw[max_index].begin();
      MSSpectrum<DRawDataPoint<1> >::const_iterator last = ms_raw[max_index].end();

      for (; first != (last-2) ;++first)
      {
        float act_spacing = fabs((first+1)->getPosition()[0] - first->getPosition()[0]);
        if (act_spacing < spacing)
	  {
	    spacing = act_spacing;
	  }
      }

      return spacing;
    }

  protected:
    /** @name Attributes
     */
    //@{
    // The computed peak shapes for all peaks we have found
    std::vector<PeakShape> peak_shapes_;
    ///
    // The class providing the continuous wavelet transform
    ContinuousWaveletTransform<D>* wt_;
    ///
    // The search radius for the determination of the maximum
    int radius_;
    ///
    // The dilation of the wavelet
    float scale_;
    ///
    // the minimal height which defines a peak in the CWT
    float peak_bound_cwt_;
    ///
    // the minimal height which defines a peak in the CWT
    float peak_bound_ms2_level_cwt_;
    ///
    // The threshold for asymmetry
    ///
    float peak_asymm_bound_;
    ///
    // The threshold for correlation
    float peak_corr_bound_;
    ///
    // The threshold for the minimal fwhm
    float peak_fwhm_bound_;
    // The threshold for the noise level (TODO: Use the information of the signal to noise estimator)
    float noise_level_;
    // Use the optimization of peak parameters
    bool optimization_;
    // Compute the integration with numerical methods
    bool num_integration_;
    //@}

    /** A regularData-Object which contents some additional useful informations
	for analysing peaks and their properties
    */
    class PeakArea_
    {
      typedef typename std::vector<DRawDataPoint<D> >::iterator RawIterator;

    public:
      PeakArea_(){}

      RawIterator left, max, right, left_behind_centroid;
      DPosition<D> centroid_position;
    }
    ;

    /** @name Functions for Stick Conversion
     */
    //@{
    void getPeakArea_(const PeakArea_& area, double &area_left, double &area_right);

    /** Fit a peak form.
     */
    PeakShape fitPeakShape_(const PeakArea_& area, bool enable_centroid_fit, bool ms_experiment, double current_rt);

    /** Returns a measure for the similarity of s and phi. For that, the increase of the
	linear regression of the points (s(i),phi(i)) for each grid-point i in s is computed
	If the value is near 1, s and phi are expected to be very similar.
    */
    double correlate_(const PeakShape& peak, const PeakArea_& area, int direction=0) const;
    //@}

    /** @name Functions for the Separation of Overlapping Peaks
     */
    //@{
    // fit a symmetric sech peak in the data at the centroid position
    PeakShape fitSymmetricSechPeakShape_(const PeakArea_& area, bool enable_centroid_fit, int direction);
    //@}

    /** @name Auxiliary Functions
     */
    //@{
    // calculate the weighted mean (weights are the square root of the peak areas)
    void changeCentroidPositions_(PeakData& sort_wrt_mz, std::vector<int>& index_of_peaks, double monoisotopic_position);
    void subtractPeakShapeVector_(RawData& raw_array, const PeakShape& peak);
    // Given a m/z value pos search for the corresponding rawdata point in the intervall [first,last).
    // The function returns -1 if pos doesn't lie in the intervall
    int searchPosInRawData(RawIterator first, RawIterator last,
                           unsigned int start_index,
                           double mz_value);
    //@}

    /** @name Functions for Peak Detection
     */
    //@{
    /** Finds the next peak-position in f, using the cwt Wf. Only peaks with values greater
	than max are relvant. If direction is +1, the method starts at the left end and runs
	to the right, if direction is -1, from right to the left.
	Either take the position out of the CWT (when the resolution of the CWT is bigger than
	one) or assign it in the signal.
	The return value is a pair of an iterator which either points to the maximum in the signal
	(resolution of the CWT == 1) or to the point in the signal which corresponds the
	maximum position in the CWT.
    */
    bool getMaxPosition_(RawIterator first, 
			 RawIterator last,
                         ContinuousWaveletTransform<D>* wt, 
			 PeakArea_& area, 
			 int distance_from_scan_border,
			 int ms_level, 
			 int direction=1);


    /** This method extracts the Endpoints of the Peak given the position of a local maximum.
	'Walk' simultaniously left and right down of the maximum until a value is smaller than
	the peak bound. With the difference of this value to the maximum position jump to the
	other side of the maximum. If the signal next to this points left and right beside the maximum
	are still falling extend the intervall until real minima on both sides of the maximum are found.
    */
    bool getPeakEndPoints_(RawIterator first, RawIterator last,
                           PeakArea_ &area, int& peak_left_index, int& peak_right_index);


    /**
     *  Compute the transformation of a peak area with a higher resolution (more supporting points) and define the maximum of
     *  the transformed peak as the centroid.
     *  If the calculated centroid lies between two data point positions the function returns true otherwise return false.
     */
    void getPeakCentroid_(PeakArea_& area);
    //@}


    /** @name Functions for the Separation of Overlapping Peaks
     */
    //@{
    void getSymmetricPeakEndPoints_(RawIterator first, RawIterator last,
                                    PeakArea_& area, int& left_index, int& right_index,
                                    int& direction);

    /** This function separates overlapping peaks.
     *  First the transformation of the assumed overlapped peak with a smaller
     *  dilation is computed. Afterwards the maxima, the endpoints and the
     *  centroid of the located peaks are iteratively determined
     *  and to every peak maxima is correlated with a symmetric sech peak.
     */
    bool separateOverlappingPeaks_(PeakArea_ area, std::vector<PeakShape>& overlapped_peaks,
                                   bool enable_centroid_fit);
    //@}

    /** @name Auxiliary Functions
     */
    //@{
    /// Compute a theroretical lorentz peak with height peak_bound
    double lorentz_(double height, double lambda, double pos, double x);

    /// Compute the minimal Intensity in the cwt for a data point to be considered as a peak
    void calculatePeakBoundCWT_();

    /** Yields a peak iterator which points to peak[i], whereby: peak[i] < value < peak[i+1]
     */
    RawIterator getIteratorLeftDataPoint_(RawIterator first,
                                          RawIterator last,
                                          double value);
    //@}
  }
  ; // end of internal class PeakArea_


  template <Size D,  typename MapType, typename MapTypeOut>
  DPeakPickerCWT<D, MapType, MapTypeOut>::DPeakPickerCWT()
    : DPeakPicker<D, MapType, MapTypeOut>(),
      radius_(3),
      scale_(0.15),
      peak_asymm_bound_(0),
      peak_corr_bound_(0),
      peak_fwhm_bound_(0.2),
      noise_level_(10),
      optimization_(false),
      num_integration_(true)
  {

    wt_= new ContinuousWaveletTransformNumIntegration<D>();

    // estimate the peak bound in the wavelet transform concerning the peak bound in the original signal
    calculatePeakBoundCWT_();
  }

  template <Size D,  typename MapType, typename MapTypeOut>
  DPeakPickerCWT<D, MapType, MapTypeOut>::DPeakPickerCWT(const String& filename) : DPeakPicker<D, MapType, MapTypeOut>(filename)
  {
    init();
  }

  template <Size D,  typename MapType, typename MapTypeOut>
  DPeakPickerCWT<D, MapType, MapTypeOut>::DPeakPickerCWT(const Param& parameters) : DPeakPicker<D, MapType, MapTypeOut>(parameters)
  {
    init();
  }

  template <Size D,  typename MapType, typename MapTypeOut>
  void DPeakPickerCWT<D, MapType, MapTypeOut>::init()
  {

    //std::cout << param_ << std::endl;
    // if a peak picking parameter is missed in the param object the value should be substituted by a default value
    DataValue dv = param_.getValue("Thresholds:Asymmetry");
    if (dv.isEmpty() || dv.toString() == "") peak_asymm_bound_ = 0;
    else peak_asymm_bound_ = (float)dv;

    dv =  param_.getValue("Thresholds:Correlation");
    if (dv.isEmpty() || dv.toString() == "") peak_corr_bound_ = 0;
    else peak_corr_bound_ = (float)dv;

    dv = param_.getValue("Thresholds:Fwhm");
    if (dv.isEmpty() || dv.toString() == "") peak_fwhm_bound_ = 0.2;
    else peak_fwhm_bound_ = (float)dv;

    dv = (param_.getValue("Optimization:SkipOptimization"));
    if (dv.isEmpty() || dv.toString() == "") optimization_ = false;
    else optimization_ = (dv.toString() == "no");

    dv = param_.getValue("WaveletTransform:Scale");
    if (dv.isEmpty() || dv.toString() == "") scale_ = 0.15;
    else scale_ = (float)dv;

    dv = param_.getValue("Thresholds:NoiseLevel");
    if (dv.isEmpty() || dv.toString() == "") noise_level_ = 10.;
    else noise_level_ = (float)dv;

    //    std::cout << "Noise Level " << noise_level_ << " scale " << scale_ << std::endl;

    dv =param_.getValue("Thresholds:SearchRadius");
    if (dv.isEmpty() || dv.toString() == "") radius_ = 3;
    else radius_ = (int)dv;

    dv = param_.getValue("WaveletTransform:NumIntegration");
    if (dv.isEmpty() || dv.toString() == "") num_integration_ = true;
    else num_integration_ = (dv.toString() == "yes");

    wt_= new ContinuousWaveletTransformNumIntegration<D>();

    // estimate the peak bound in the wavelet transform concerning the peak bound in the original signal
    calculatePeakBoundCWT_();
  }

  template <Size D,  typename MapType, typename MapTypeOut>
  DPeakPickerCWT<D, MapType, MapTypeOut>::~DPeakPickerCWT()
  {
    if (wt_)
      {
	delete wt_;
#ifdef DEBUG_PEAK_PICKING
	std::cout << "delete wt" << std::endl;
#endif
	wt_=0;
      }
  }

  template <Size D,  typename MapType, typename MapTypeOut>
  void DPeakPickerCWT<D, MapType, MapTypeOut>::pick(const  MapType& ms_exp_raw)
  {
    pick(ms_exp_raw,this);
  }

  template <Size D,  typename MapType, typename MapTypeOut>
  void DPeakPickerCWT<D, MapType, MapTypeOut>::pick(const MapType&, DPeakPickerCWT<2, MapType, MapTypeOut> const*)
  {
    OPENMS_PRECONDITION(D==1, "Use the one dimensional peak picker for instances of MSExperiment.");
  }

  template <Size D,  typename MapType, typename MapTypeOut>
  void DPeakPickerCWT<D, MapType, MapTypeOut>::pick(const MapType& ms_exp_raw, DPeakPickerCWT<1, MapType, MapTypeOut> const*)
  {
    IteratorVector vec(2);
    //    DExtractSignalRegions<D,RawData >  split(param_);
    DSignalToNoiseEstimatorWindowing<D,RawData> sne;

    //    split.setParam(param_);

    // iterate over all scans
    for (unsigned int i=0; i < ms_exp_raw.size(); ++i)
      {
	PointConstIterator first_data_point = ms_exp_raw[i].begin();
	PointConstIterator last_data_point = ms_exp_raw[i].end();

	wt_->init(scale_, 0.001, mz_dim_);


	sne.init(first_data_point,last_data_point); 
	  
	//std::cout << "scan " << ms_exp_raw[i].getRetentionTime() << std::endl;
      
	//split.splitScan(first_data_point,last_data_point,noise_level_,vec);
	vec[0]=first_data_point;
	vec[1]=last_data_point;
      
	if (vec.size() > 0)
	  {

	    pick(vec.begin(),vec.end(),sne,true,
		 ms_exp_raw[i].getRetentionTime(),
		 ms_exp_raw[i].getMSLevel());
			   
	    std::sort(peak_shapes_.begin(),peak_shapes_.end(),PeakShape::PositionLess(2));

	    // write the peak shapes into the MSExperiment
	    typename std::vector<PeakShape>::const_iterator iter = peak_shapes_.begin();
	    typedef typename SpectrumTypeOut::PeakType Peak;

	    SpectrumTypeOut spec;

	    spec.setRetentionTime(ms_exp_raw[i].getRetentionTime(), ms_exp_raw[i].getRetentionTimeStart(), ms_exp_raw[i].getRetentionTimeStop());
	    spec.setMSLevel(ms_exp_raw[i].getMSLevel());
	    spec.setName(ms_exp_raw[i].getName());

	    for (; iter != peak_shapes_.end(); iter++)
	      {
		// create temporary peak and insert it into spectrum
		Peak input_peak;
		input_peak.getIntensity()   = iter->height;
		input_peak.getPosition()[0] = iter->mz_position;
		input_peak.getArea() = iter->area;
		input_peak.getRValue() = iter->r_value;
		input_peak.getLeftWidthParameter() = iter->left_width;
		input_peak.getRightWidthParameter() = iter->right_width;
		input_peak.getFWHM() = iter->getFWHM();
		input_peak.getPeakShape() = iter->type;
		input_peak.getCharge() = -1;
		input_peak.getSN() = iter->signal_to_noise;

		spec.getContainer().push_back(input_peak);
	      }
	    //clear the member peak_shapes_
	    peak_shapes_.clear();
       
	    ms_exp_peaks_->push_back(spec);
	  }

      }
  }

  template <Size D,  typename MapType, typename MapTypeOut>
  void DPeakPickerCWT<D,  MapType, MapTypeOut>::pick(RawConstIterator first, RawConstIterator last, std::back_insert_iterator<PeakData >& output)
  {
    IteratorVector vec;
    DExtractSignalRegions<D>  split;
    DSignalToNoiseEstimatorWindowing<D,DPeakArrayNonPolymorphic<D, DRawDataPoint<D> > > sne;

    double precision;
    DataValue dv = param_.getValue("Thresholds:Precision");
    if (dv.isEmpty() || dv.toString() == "") precision = 1e-5;
    else precision = (float)dv;
    split.setParam(param_);

    /** Initialize the Wavelet Transform **/
    float wavelet_spacing = 0.001; // ??? TODO: should be parameter in the Param object
    wt_->init(scale_, wavelet_spacing, mz_dim_);

    RawConstIterator scan_first=first;
    RawConstIterator scan_last;

    if (D==1)
      {
	// 1D raw data (one mass spectrum)
	scan_last=last;
      }
    else
      {
	scan_last=first+1;
	// search for the scans in the 2D raw data and pick peaks on each scan gradually
	while (scan_last != last)
	  {
	    double f1 = scan_first->getPosition()[rt_dim_];
	    double f2 = scan_last->getPosition()[rt_dim_];
	    // new scan
	    if (fabs(f1 - f2) > precision)
	      {
		sne.init(scan_first, scan_last);
		split.splitScan(scan_first,scan_last,noise_level_,vec);

		if (vec.size() == 0)
		  {
		    std::cerr << "No Peaks detected" << std::endl;
		  }
		else
		  {
		    pick(vec.begin(),vec.end(),sne,false,-1);
		  }
		scan_first=scan_last;
	      }
	    ++scan_last;
	  }
      }

    // pick the peaks on the last scan or rather on the 1D mass spectrum
    sne.init(scan_first, scan_last);
    // TODO: Split depending on the DSignalToNoiseEstimator !!!!!!
    split.splitScan(scan_first,scan_last,noise_level_,vec);

    if (vec.size() == 0)
      {
	std::cerr << "No Peaks detected " << std::endl;
      }
    else
      {
	pick(vec.begin(),vec.end(),sne,false,-1);
      }

    std::sort(peak_shapes_.begin(),peak_shapes_.end(),PeakShape::PositionLess(D));

    //
    for (unsigned int i=0; i < peak_shapes_.size(); ++i)
      {
	OutputPeak p;
	p.setIntensity(peak_shapes_[i].height);
	DPosition<D> pos;
	if (D==2)
	  {
	    pos[rt_dim_]=peak_shapes_[i].rt_position;
	  }
	pos[mz_dim_]=peak_shapes_[i].mz_position;
	p.setPosition(pos);
	p.setArea(peak_shapes_[i].area);
	*(output++)=p;
      }
  }

  template <Size D,  typename MapType, typename MapTypeOut>
  int DPeakPickerCWT<D, MapType, MapTypeOut>::searchPosInRawData(RawIterator first, RawIterator last,
								 unsigned int start_index,
								 double mz_value)
  {
    if ((mz_value < first->getPosition()[mz_dim_]) || (mz_value > (last-1)->getPosition()[mz_dim_]))
      {
	return -1;
      }
    else
      {
	unsigned int index = start_index;
	RawIterator help = first+index;
	double left_pos = help->getPosition()[mz_dim_];
	double right_pos = (help+1)->getPosition()[mz_dim_];
	int direction = (mz_value < left_pos) ? -1 : 1;

	while ((help >= first) || (help < (last-1)))
	  {
	    left_pos  = help->getPosition()[mz_dim_];
	    right_pos = (help+1)->getPosition()[mz_dim_];

	    // if mz_value lies between the help and help+1 raw data points return the index to the closest raw data point
	    if ((left_pos < mz_value) && (mz_value < right_pos))
	      {
		if (fabs(left_pos-mz_value) < fabs(right_pos-mz_value))
		  {
		    return index;
		  }
		else
		  {
		    return index + 1;
		  }
	      }
	    help+=direction;
	    index+=direction;
	  }
      }
    return -1;
  }

  template <Size D,  typename MapType, typename MapTypeOut>
  bool DPeakPickerCWT<D,  MapType, MapTypeOut>::getMaxPosition_(RawIterator first, 
								RawIterator last,
								ContinuousWaveletTransform<D>* wt, 
								PeakArea_& area,
								int distance_from_scan_border,
								int ms_level, 
								int direction)
  {
    // ATTENTION! It is assumed that the resolution==1 (no resolution higher than 1).
    double noise_level=0.;
    double noise_level_cwt=0.;
    if (ms_level==1)
      {
	noise_level = peak_bound_;
	noise_level_cwt = peak_bound_cwt_;
      }
    else
      {
	noise_level = peak_bound_ms2_level_;
	noise_level_cwt = peak_bound_ms2_level_cwt_;
      }

    int zeros_left_index  = wt->getLeftPaddingIndex();
    int zeros_right_index = wt->getRightPaddingIndex();

    // Points to most intensive data point in the signal
    RawIterator it_max_pos;
    double max_value;
    
    // Given direction, start the search from left or right
    int start = (direction > 0) ? ((zeros_left_index + 2) + distance_from_scan_border) : ((zeros_right_index - 2) - distance_from_scan_border) ;
    int end   = (direction > 0) ? (zeros_right_index - 1)  : zeros_left_index+1;

    int i=0, j=0, k, max_pos;

    for(i=start, k=0; i!=end; i+=direction, ++k)
      {
	// Check for maximum in cwt at position i
	if((((*wt)[i-1] - (*wt)[i]  ) < 0)
	   && (((*wt)[i]   - (*wt)[i+1]) > 0)
	   && ( (*wt)[i]   >  noise_level_cwt))
	  {
	    max_pos = (direction > 0) ? (i - distance_from_scan_border)  : i;
#ifdef DEBUG_PEAK_PICKING
	    std::cout << "MAX in CWT at " << (first + max_pos)->getPosition()[mz_dim_]<< " with " << (*wt)[i]
		      << std::endl;
#endif
	    max_value=(first + max_pos)->getIntensity();


	    // search for the corresponding maximum in the signal (consider the radius left and right adjacent points)
	    int start_intervall = ((max_pos - radius_) < 0 ) ? 0 : (max_pos - radius_);
	    int end_intervall= ((max_pos + radius_) >= distance(first,last)) ? 0 : (max_pos + radius_);

	    for(j = start_intervall; j <= end_intervall; ++j)
	      {
		if((first + j)->getIntensity() > max_value)
		  {
		    max_pos = j;
		    max_value = (first + j)->getIntensity();
		  }

	      }

	    // if the maximum position is high enough and isn't one of the border points, we return it
	    if(((first + max_pos)->getIntensity() >= noise_level)
	       && (((first + max_pos) != first)
		   && (first + max_pos)   != (last-1)))
	      {
		area.max = first + max_pos;
#ifdef DEBUG_PEAK_PICKING
		std::cout << "_________Max in data at__________ " << area.max->getPosition()[mz_dim_];
		if (D==2)
		  std::cout << "   in Scan " << area.max->getPosition()[rt_dim_] ;

		std::cout << std::endl;
#endif
		return true;
	      }

#ifdef DEBUG_PEAK_PICKING
	    std::cout << "NO REAL MAX!!!!!!!" << std::endl;
#endif

	  }
      }

    // No relevant peak was found
    return false;
  }



  template <Size D,  typename MapType, typename MapTypeOut>
  bool DPeakPickerCWT<D, MapType, MapTypeOut>::getPeakEndPoints_(RawIterator first,
								 RawIterator last,
								 PeakArea_& area,
								 int& peak_left_index,
								 int& peak_right_index)
  {
    // the Maximum may neither be the first or last point in the signal
    if ((area.max <= first) || (area.max >= last-1))
      {
	return false;
      }

    /*  The algorithm does the following:
     *    - let x_m be the position of the maximum in the data and let (x_l, x_r) be
     *      the left and right neighbours
     *
     *
     *  (1) starting from x_l', walk left until one of the following happens
     *         - the new point is lower than the original bound
     *              => we found our left endpoint
     *
     *         - the new point is larger than the last, but the point left from
     *           the new point is smaller. In that case, we either ran into another
     *           peak, or we encounter some noise. Therefore we now look in the cwt
     *           at the position corresponding to this value. If the cwt here is
     *           monotonous, we consider the point as noise and continue further to the
     *           left. Otherwise, we probably found the beginning of a new peak and
     *           therefore stop here.
     *
     *  (2) analogous procedure to the right of x_r
     */

    RawIterator it_help=area.max-1;
    double vec_pos;
    int cwt_pos;
    int ep_radius=2;
    int start;
    int stop;
    bool monoton;

    int zeros_left_index  = wt_->getLeftPaddingIndex();

    // search for the left endpoint
    while (((it_help-1) > first) && (it_help->getIntensity() > noise_level_))
      {
	// if the values are still falling to the left, everything is ok.
	if ((it_help-1)->getIntensity() < it_help->getIntensity())
	  {
	    --it_help;
	  }
	// if the values are _rising_, we have to check the cwt
	else
	  {
	    if ((it_help-2) <= first)
	      {
		break;
	      }
	    // now check the value to the left of the problematic value
	    if ((it_help-2)->getIntensity() > (it_help-1)->getIntensity()) // we probably ran into another peak
	      {
		break;
	      }


	    // to the left, the values are falling again => let the cwt decide if we
	    // are seeing a new peak or just noise

	    // compute the position of the corresponding point in the cwt
	    cwt_pos = distance(first, it_help);
	    vec_pos=it_help->getPosition()[mz_dim_];

	    // since the cwt is pretty smooth usually, we consider the point as noise
	    // if the cwt is monotonous in this region
	    // TODO: better monotonicity test... say two or three points more
	    monoton=true;
	    start   =   ((cwt_pos-ep_radius) < 0)
	      ? zeros_left_index+1
	      : cwt_pos-ep_radius+zeros_left_index+1;
	    stop    =   ((cwt_pos+ep_radius) > wt_->getSignalLength())
	      ?  (wt_->getSignalLength() + zeros_left_index)
	      : (cwt_pos+ep_radius+zeros_left_index);

	    for (; start < stop; ++start)
	      {
		if (   ((*wt_)[start-1] - (*wt_)[start]  )
		       * ((*wt_)[start]   - (*wt_)[start+1]) < 0 )
		  {
		    // different slopes at the sides => stop here
		    monoton=false;
		    break;
		  }
	      }

	    if (!monoton)
	      {
		break;
	      }
	    --it_help;
	  }
      }
    area.left=it_help;

    it_help=area.max+1;
    // search for the right endpoint ???
    while (((it_help+1) < last) && (it_help->getIntensity() > noise_level_))
      {
	// if the values are still falling to the right, everything is ok.
	if (it_help->getIntensity() > (it_help+1)->getIntensity())
	  {
	    ++it_help;
	  }
	// if the values are _rising_, we have to check the cwt
	else
	  {
	    if ((it_help+2) >= last)
	      {
		break;
	      }
	    // now check the value to the right of the problematic value
	    if ((it_help+2)->getIntensity() > (it_help+1)->getIntensity()) // we probably ran into another peak
	      {
		break;
	      }

	    // to the left, the values are falling again => let the cwt decide if we
	    // are seeing a new peak or just noise
	    // compute the position of the corresponding point in the cwt
	    cwt_pos = distance(first, it_help);
	    //cwt_pos = distance(first, it_help);
	    vec_pos=it_help->getPosition()[mz_dim_];

	    // since the cwt is pretty smooth usually, we consider the point as noise
	    // if the cwt is monotonous in this region
	    // TODO: better monotonicity test... say two or three points more
	    monoton = true;

	    start   =   ((cwt_pos-ep_radius) < 0)
	      ? zeros_left_index+1
	      : cwt_pos-ep_radius+zeros_left_index+1;
	    stop    =   ((cwt_pos+ep_radius) > wt_->getSignalLength())
	      ?  (wt_->getSignalLength() + zeros_left_index)
	      : (cwt_pos+ep_radius+zeros_left_index);

	    for (; start < stop; ++start)
	      {
		if (   ((*wt_)[start-1] - (*wt_)[start])
		       * ((*wt_)[start]   - (*wt_)[start+1]) < 0 )
		  {
		    // different slopes at the sides => stop here
		    monoton=false;
		    break;
		  }
	      }

	    if (!monoton)
	      {
		break;
	      }
	    ++it_help;
	  }
      }
    area.right=it_help;

    peak_left_index=distance(first,area.left);
    peak_right_index=distance(first,area.right);

    // The minimal raw data points per peak should be 2
    if ((distance(area.left,area.max) > 0) && (distance(area.max,area.right) > 0))
      {
	return true;
      }

    return false;
  }

  template <Size D,  typename MapType, typename MapTypeOut>
  void DPeakPickerCWT<D, MapType, MapTypeOut>::getPeakCentroid_(PeakArea_& area)
  {
    RawIterator left_it=area.max-1, right_it=area.max;
    double max_intensity=area.max->getIntensity();
    double rel_peak_height=max_intensity*0.6;
    double sum=0., w=0.;
    area.centroid_position[mz_dim_]=area.max->getPosition()[mz_dim_];

    // compute the centroid position (use weighted mean)
    while ((left_it >= area.left) && (left_it->getIntensity() >=rel_peak_height) )
      {
	if (left_it->getIntensity() >=rel_peak_height)
	  {
	    w+=left_it->getIntensity()*left_it->getPosition()[mz_dim_];
	    sum+=left_it->getIntensity();
	    --left_it;
	  }
      }

    while ((right_it < area.right) && (right_it->getIntensity() >=rel_peak_height) )
      {
	if (right_it->getIntensity() >=rel_peak_height)
	  {
	    w+=right_it->getIntensity()*right_it->getPosition()[mz_dim_];
	    sum+=right_it->getIntensity();
	    ++right_it;
	  }
      }

    area.centroid_position[mz_dim_]=w / sum;


    if (D==2)
      {
	area.centroid_position[rt_dim_]=area.max->getPosition()[rt_dim_];
      }


    // test if the centroid is an existing data point
    /*  if (fabs(area.centroid_position[mz_dim_]-area.max->getPosition()[mz_dim_]) < 0.001)
	{
	area.left_behind_centroid=area.max;
	return false;
	}

	area.left_behind_centroid = getIteratorLeftDataPoint_(area.left,
	area.right+1,
	area.centroid_position[mz_dim_]);
    */
#ifdef DEBUG_PEAK_PICKING
    std::cout << "________Centroid is___________ " << area.centroid_position[mz_dim_] << std::endl;
#endif

  }


  template <Size D,  typename MapType, typename MapTypeOut>
  void DPeakPickerCWT<D, MapType, MapTypeOut>::getSymmetricPeakEndPoints_(RawIterator first,
									  RawIterator last,
									  PeakArea_& area,
									  int& left_index,
									  int& right_index,
									  int& direction)
  {
#ifdef DEBUG_PEAK_PICKING
    std::cout << "*********getSymmetricPeakEndPoints_***********" << std::endl;
#endif

    // search the endpoints of the peak simultaneous and assume a symmetric peak form
    int pos, left=0, right=0, step=1, n=distance(first,last);
    RawIterator it_max_pos = area.max;
    pos = distance(it_max_pos, first);//(first,it_max_pos);

    while(true)
      {
	if ((pos-step > 0) && (it_max_pos-step)->getIntensity() < peak_bound_)
	  {
	    direction=1;
	    break;
	  }
	if ((pos+step < n) && (it_max_pos+step)->getIntensity() < peak_bound_)
	  {
	    direction=-1;
	    break;
	  }
	++step;
      }

    left = (pos-step < 0) ? pos : step;
    right = (pos+step > n) ? n-1-pos : step;

    area.left=it_max_pos-left;
    area.right=it_max_pos+right;

    left_index=distance(first,area.left);
    right_index=distance(first,area.right);
  }



  template <Size D,  typename MapType, typename MapTypeOut>
  double DPeakPickerCWT<D, MapType, MapTypeOut>::lorentz_(double height, double lambda, double pos, double x)
  {
    return height/(1+pow(lambda*(x-pos),2));
  }

  template <Size D, typename MapType, typename MapTypeOut>
  void DPeakPickerCWT<D, MapType, MapTypeOut>::calculatePeakBoundCWT_()
  {
#ifdef DEBUG_PEAK_PICKING
    std::cout << "DPeakPickerCWT<D>::calculatePeakBoundCWT_ peak_bound_" << peak_bound_ <<  std::endl;
#endif

    // build a lorentz peak of height peak_bound_
    // compute its cwt, and compute the resulting height
    // of the transformed peak

    //compute the peak in the intervall [-2*scale,2*scale]
    double spacing=0.001;
    int n = (int)((4*scale_)/spacing)+1;

    // compute the width parameter using height=peak_bound_ and the peak endpoints should be -scale and +scale, so at
    // positions -scale and +scale the peak value should correspond to the noise_level_
    double lambda = sqrt((-noise_level_*(-peak_bound_+noise_level_)))/(noise_level_*scale_);

    RawData lorentz_peak(n);
    RawData lorentz_peak2(n);

    /** TODO: switch the type of the transform **/

    ContinuousWaveletTransform<D>* lorentz_cwt;
    ContinuousWaveletTransform<D>* lorentz_ms2_cwt;

    lorentz_cwt = new ContinuousWaveletTransformNumIntegration<D>();
    lorentz_ms2_cwt = new ContinuousWaveletTransformNumIntegration<D>();

    lorentz_cwt->init(scale_, spacing, mz_dim_);
    lorentz_ms2_cwt->init(scale_, spacing, mz_dim_);
    double start = -2*scale_;
    for (int i=0; i < n; ++i)
      {
	DPosition<D> p;
	p[mz_dim_]=i*spacing + start;
	lorentz_peak[i].setPosition(p);
	lorentz_peak[i].setIntensity(lorentz_(peak_bound_,lambda,0,i*spacing + start));
	lorentz_peak2[i].setPosition(p);
	lorentz_peak2[i].setIntensity(lorentz_(peak_bound_ms2_level_,lambda,0,i*spacing + start));
      }

    float resolution = 1.;
    lorentz_cwt->transform(lorentz_peak.begin(), lorentz_peak.end(),resolution);
    lorentz_ms2_cwt->transform(lorentz_peak2.begin(), lorentz_peak2.end(),resolution);

    float peak_max=0;
    float peak_max2=0;

    for (int i=0; i<lorentz_cwt->getSignalLength(); i++)
      {
	if ((*lorentz_cwt)[i] > peak_max)
	  {
	    peak_max = (*lorentz_cwt)[i];
	  }
	if ((*lorentz_ms2_cwt)[i] > peak_max2)
	  {
	    peak_max2 = (*lorentz_ms2_cwt)[i];
	  }
      }

    peak_bound_cwt_ = peak_max;
    peak_bound_ms2_level_cwt_ = peak_max2;
#ifdef DEBUG_PEAK_PICKING
    std::cout << "PEAK BOUND IN CWT " << peak_bound_cwt_ << std::endl;
    std::cout << "PEAK BOUND IN CWT (MS 2 Level)" << peak_bound_ms2_level_cwt_ << std::endl;
#endif

    delete lorentz_cwt;
    delete lorentz_ms2_cwt;
  }

  template <Size D,  typename MapType, typename MapTypeOut>
  typename DPeakPickerCWT<D,  MapType, MapTypeOut>::RawIterator
  DPeakPickerCWT<D, MapType, MapTypeOut>::getIteratorLeftDataPoint_(RawIterator first,
								    RawIterator last,
								    double value)
  {
    int length = distance(first,last);
    double origin = first->getPosition()[mz_dim_];

    OPENMS_PRECONDITION(((origin < value) && (value < (last-1)->getPosition()[mz_dim_])),
                        "The position can't be found in this peak array.");

    double spacing = ((last-1)->getPosition()[mz_dim_]-origin)/(length-1);
    double distance=value-origin;

    int value_index=(int)(distance/spacing);

    RawIterator it_pos=first+value_index;
    while (true)
      {
	if (it_pos->getPosition()[mz_dim_] < value)
	  {
	    if ((it_pos+1)->getPosition()[mz_dim_] < value)
	      {
		++it_pos;
	      }
	    else
	      {
		return it_pos;
	      }
	  }
	else
	  {
	    --it_pos;
	  }
      }
  }

  template <Size D,  typename MapType, typename MapTypeOut>
  void DPeakPickerCWT<D, MapType, MapTypeOut>::getPeakArea_(const DPeakPickerCWT<D, MapType, MapTypeOut>::PeakArea_& area, double& area_left, double& area_right)
  {
    area_left += area.left->getIntensity() * ((area.left+1)->getPosition()[mz_dim_] - area.left->getPosition()[mz_dim_]) * 0.5;
    area_left += area.max->getIntensity() *  (area.max->getPosition()[mz_dim_] - (area.max-1)->getPosition()[mz_dim_]) * 0.5;

    for (RawIterator pi=area.left+1; pi<area.max; pi++)
      {
	double step = ((pi)->getPosition()[mz_dim_] - (pi-1)->getPosition()[mz_dim_]);
	area_left += step * pi->getIntensity();
      }

    area_right += area.right->getIntensity() * ((area.right)->getPosition()[mz_dim_] - (area.right-1)->getPosition()[mz_dim_]) * 0.5;
    area_right += (area.max+1)->getIntensity() *  ((area.max+2)->getPosition()[mz_dim_] - (area.max+1)->getPosition()[mz_dim_]) * 0.5;

    for (RawIterator pi=area.max+2; pi<area.right; pi++)
      {
	double step = ((pi)->getPosition()[mz_dim_] - (pi-1)->getPosition()[mz_dim_]);
	area_right += step * pi->getIntensity();
      }
  }

  template <Size D,  typename MapType, typename MapTypeOut>
  PeakShape DPeakPickerCWT<D, MapType, MapTypeOut>::fitPeakShape_
  (const DPeakPickerCWT<D, MapType, MapTypeOut>::PeakArea_& area,
   bool enable_centroid_fit,
   bool ms_experiment,
   double current_rt)
  {

#ifdef DEBUG_PEAK_PICKING
    std::cout << "Left end point: " << area.left->getPosition()[mz_dim_] << std::endl;
#endif
    double max_intensity   =   area.max->getIntensity();
    double left_intensity  =  area.left->getIntensity();
    double right_intensity = area.right->getIntensity();

    //avoid zero width
    float minimal_endpoint_centroid_distance=0.01;
    if (  (fabs( area.left->getPosition()[mz_dim_]-area.centroid_position[mz_dim_]) < minimal_endpoint_centroid_distance)
          ||(fabs(area.right->getPosition()[mz_dim_]-area.centroid_position[mz_dim_]) < minimal_endpoint_centroid_distance) )
      {
#ifdef DEBUG_PEAK_PICKING
	std::cout << "The distance between centroid and the endpoints is too small!" << std::endl;
#endif
	return PeakShape();
      }

    if (enable_centroid_fit)
      {
#ifdef DEBUG_PEAK_PICKING
	std::cout << "Fit at the peak centroid" << std::endl;
#endif
	// the maximal position was taken directly from the cwt.
	// first we do a "regular" fit of the left half
	// TODO: avoid zero width!

#ifdef DEBUG_PEAK_PICKING
	std::cout << "Left end point: "         << area.left->getPosition()[mz_dim_]
		  << " centroid: "              << area.centroid_position[mz_dim_]
		  << " right end point: "       << area.right->getPosition()[mz_dim_]
		  << std::endl;
	std::cout << " point left of centroid:" << (area.left_behind_centroid)->getPosition()[mz_dim_]
		  << std::endl;
	if (D==2)
	  std::cout << "Peak in Scan " << area.centroid_position[rt_dim_] << std::endl;
#endif

	// lorentzian fit

	// estimate the width parameter of the left peak side
	RawIterator left_it=area.left_behind_centroid;
	double x0=area.centroid_position[mz_dim_];
	double l_sqrd=0.;
	int n=0;
	while(left_it-1 >= area.left)
	  {
	    double x1=left_it->getPosition()[mz_dim_];
	    double x2=(left_it-1)->getPosition()[mz_dim_];
	    double c=(left_it-1)->getIntensity()/left_it->getIntensity();
	    l_sqrd+=(1-c)/(c*(pow((x2-x0),2))-pow((x1-x0),2));
	    --left_it;
	    ++n;
	  }
	double left_heigth=area.left_behind_centroid->getIntensity()/(1+l_sqrd*pow(area.left_behind_centroid->getPosition()[mz_dim_]-area.centroid_position[mz_dim_],2));

	// estimate the width parameter of the right peak side
	RawIterator right_it=area.left_behind_centroid+1;
	l_sqrd=0.;
	n=0;
	while(right_it+1 <= area.right)
	  {
	    double x1=right_it->getPosition()[mz_dim_];
	    double x2=(right_it+1)->getPosition()[mz_dim_];
	    double c=(right_it+1)->getIntensity()/right_it->getIntensity();
	    l_sqrd+=(1-c)/(c*(pow((x1-x0),2))-pow((x2-x0),2));
	    ++right_it;
	    ++n;
	  }

	//estimate the heigth
	double right_heigth=(area.left_behind_centroid+1)->getIntensity()/(1+l_sqrd*pow((area.left_behind_centroid+1)->getPosition()[mz_dim_]-area.centroid_position[mz_dim_],2));

	double height=std::min(left_heigth,right_heigth);

	// compute the left and right area
	double peak_area_left = 0.;
	peak_area_left += area.left->getIntensity() * (  (area.left+1)->getPosition()[mz_dim_]
							 -    area.left->getPosition()[mz_dim_]  ) * 0.5;
	peak_area_left += height * (area.centroid_position[mz_dim_]-area.left_behind_centroid->getPosition()[mz_dim_]) * 0.5;

	for (RawIterator pi=area.left+1; pi <= area.left_behind_centroid; pi++)
	  {
	    double step = ((pi)->getPosition()[mz_dim_] - (pi-1)->getPosition()[mz_dim_]);
	    peak_area_left += step * pi->getIntensity();
	  }

	double peak_area_right = 0.;
	peak_area_right += area.right->getIntensity() * ((area.right)->getPosition()[mz_dim_]
							 - (area.right-1)->getPosition()[mz_dim_]  ) * 0.5;
	peak_area_right += height * ( (area.left_behind_centroid+1)->getPosition()[mz_dim_]-area.centroid_position[mz_dim_]) * 0.5;

	for (RawIterator pi=area.left_behind_centroid+1; pi < area.right; pi++)
	  {
	    double step = ((pi)->getPosition()[mz_dim_] - (pi-1)->getPosition()[mz_dim_]);
	    peak_area_right += step * pi->getIntensity();
	  }

	double left_width =    height/peak_area_left
	  * atan( sqrt( height/area.left->getIntensity() - 1. ) );
	double right_width =  height/peak_area_right
	  * atan( sqrt( height/area.right->getIntensity() - 1. ) );


	// TODO: test different heights; recompute widths; compute area
	PeakShape lorentz(height, area.centroid_position[mz_dim_], -1, left_width, right_width,
			  peak_area_left + peak_area_right, PeakShapeType::LORENTZ_PEAK);
	if (D==2)
	  {
	    lorentz.rt_position = area.centroid_position[rt_dim_];
	  }

	if (ms_experiment)
	  {
	    lorentz.rt_position = current_rt;
	  }

	lorentz.r_value = correlate_(lorentz, area);

#ifdef DEBUG_PEAK_PICKING
	std::cout << "lorentz r: " << lorentz.r_value << " " << "pos: " << lorentz.mz_position << std::endl;
	std::cout << "w1, w2: " << lorentz.left_width << " " << lorentz.right_width << std::endl;
	std::cout << "h: " << lorentz.height << std::endl;
#endif
	return lorentz;
      }

    else // no fitting on centroids
      {
#ifdef DEBUG_PEAK_PICKING
	std::cout << "fit at the peak maximum " << std::endl;
#endif
	// determine the left half of the area of the PeakArea_...
	double peak_area_left = 0.;
	peak_area_left += area.left->getIntensity() * ((area.left+1)->getPosition()[mz_dim_] - area.left->getPosition()[mz_dim_]) * 0.5;
	peak_area_left += area.max->getIntensity() *  (area.max->getPosition()[mz_dim_] - (area.max-1)->getPosition()[mz_dim_]) * 0.5;

	for (RawIterator pi=area.left+1; pi<area.max; pi++)
	  {
	    double step = ((pi)->getPosition()[mz_dim_] - (pi-1)->getPosition()[mz_dim_]);
	    peak_area_left += step * pi->getIntensity();
	  }

	double peak_area_right = 0.;
	peak_area_right += area.right->getIntensity() * ((area.right)->getPosition()[mz_dim_] - (area.right-1)->getPosition()[mz_dim_]) * 0.5;
	peak_area_right += area.max->getIntensity() *  ((area.max+1)->getPosition()[mz_dim_] - (area.max)->getPosition()[mz_dim_]) * 0.5;

	for (RawIterator pi=area.max+1; pi<area.right; pi++)
	  {
	    double step = ((pi)->getPosition()[mz_dim_] - (pi-1)->getPosition()[mz_dim_]);
	    peak_area_right += step * pi->getIntensity();
	  }

	// first the lorentz-peak...

	double left_width = max_intensity / peak_area_left * atan(sqrt(max_intensity / left_intensity - 1.));
	double right_width = max_intensity / peak_area_right * atan(sqrt(max_intensity / right_intensity - 1.));



	PeakShape lorentz(max_intensity, area.max->getPosition()[mz_dim_], -1,
			  left_width, right_width, peak_area_left + peak_area_right,
			  PeakShapeType::LORENTZ_PEAK);
	if (D==2)
	  {
	    lorentz.rt_position=area.max->getPosition()[rt_dim_];
	  }

	if (ms_experiment)
	  {
	    lorentz.rt_position = current_rt;
	  }

	lorentz.r_value = correlate_(lorentz, area);

	// now the sech-peak...
	left_width  = max_intensity /peak_area_left * sqrt(1. - left_intensity / max_intensity);
	right_width  = max_intensity /peak_area_right * sqrt(1. - right_intensity / max_intensity);


	PeakShape sech(max_intensity, area.max->getPosition()[mz_dim_], -1,
		       left_width, right_width,
		       peak_area_left + peak_area_right,
		       PeakShapeType::SECH_PEAK);

	if (D==2)
	  {
	    sech.rt_position=area.max->getPosition()[rt_dim_];
	  }

	if (ms_experiment)
	  {
	    sech.rt_position = current_rt;
	  }

	sech.r_value = correlate_(sech, area);

#ifdef DEBUG_PEAK_PICKING
	std::cout << "r: " << lorentz.r_value << " " << sech.r_value << std::endl;
	std::cout << "pos: " << lorentz.mz_position <<  " " << sech.mz_position << std::endl;
	std::cout << "w1, w2: " << lorentz.left_width << " " << lorentz.right_width << " "
		  << sech.left_width << " " << sech.right_width << std::endl;
	std::cout << "h: " << lorentz.height << std::endl;
#endif

	if ((lorentz.r_value > sech.r_value) && isnan(sech.r_value))
	  {
	    return lorentz;
	  }
	else
	  {
	    return sech;
	  }
      }
  }

  template <Size D,  typename MapType, typename MapTypeOut>
  double DPeakPickerCWT<D, MapType, MapTypeOut>::correlate_(const PeakShape& peak,
							    const DPeakPickerCWT<D, MapType, MapTypeOut>::PeakArea_& area,
							    int direction) const
  {
    double SSxx = 0., SSyy = 0., SSxy = 0.;

    // compute the averages
    double data_average=0., fit_average=0.;
    double data_sqr=0., fit_sqr=0.;
    double cross=0.;

    int number_of_points = 0;
    RawIterator corr_begin=area.left;
    RawIterator corr_end=area.right;

    // for separate overlapping peak correlate until the max position...
    if (direction > 0)
      corr_end=area.max;
    else
      if (direction < 0)
        corr_begin=area.max;

    for (RawIterator pi = corr_begin; pi<=corr_end; pi++)
      {
	double data_val = pi->getIntensity();
	double peak_val = peak(pi->getPosition()[mz_dim_]);

	data_average += data_val;
	fit_average  += peak_val;

	data_sqr += data_val * data_val;
	fit_sqr  += peak_val * peak_val;

	cross += data_val * peak_val;

	number_of_points++;
      }

    if (number_of_points == 0)
      return 0.;

    data_average /= number_of_points;
    fit_average  /= number_of_points;

    SSxx = data_sqr - number_of_points * (data_average * data_average);
    SSyy = fit_sqr - number_of_points * (fit_average * fit_average);
    SSxy = cross - number_of_points * (data_average * fit_average);

    return (SSxy * SSxy) / (SSxx * SSyy);
  }

  template <Size D,  typename MapType, typename MapTypeOut>
  PeakShape DPeakPickerCWT<D, MapType, MapTypeOut>::fitSymmetricSechPeakShape_(const DPeakPickerCWT<D, MapType, MapTypeOut>::PeakArea_& area,
									       bool enable_centroid_fit,
									       int direction)
  {
    double peak_area = 0;

#ifdef PEAK_PICKING_DEBUG
    std::cout << "Left endpoint: "  << area.left->getPosition()[mz_dim_]
	      << " maximum "        << area.max->getPosition()[mz_dim_]
	      << " right endpoint " << area.right->getPosition()[mz_dim_]
	      << std::endl;
#endif

    // compute the area to the left endpoint
    if (direction > 0)
      {
	peak_area +=      area.left->getIntensity()
	  * ( (area.left+1)->getPosition()[mz_dim_]
	      -area.left->getPosition()[mz_dim_]
	      ) * 0.5;

	peak_area +=     area.max->getIntensity()
	  *  (   area.max->getPosition()[mz_dim_]
		 -(area.max-1)->getPosition()[mz_dim_]
		 ) * 0.5;

	for (RawIterator pi=area.left+1; pi<area.max; pi++)
	  {
	    double step = ((pi)->getPosition()[mz_dim_] - (pi-1)->getPosition()[mz_dim_]);
	    peak_area += step * pi->getIntensity();
	  }
      }
    else // direction = 0
      {
	peak_area +=      area.right->getIntensity()
	  * (  (area.right)->getPosition()[mz_dim_]
	       -(area.right-1)->getPosition()[mz_dim_]
	       ) * 0.5;
	peak_area +=      (area.max+1)->getIntensity()
	  * ( (area.max+2)->getPosition()[mz_dim_]
	      -(area.max+1)->getPosition()[mz_dim_]
	      ) * 0.5;

	for (RawIterator pi=area.max+2; pi<area.right; pi++)
	  {
	    double step = ((pi)->getPosition()[mz_dim_] - (pi-1)->getPosition()[mz_dim_]);
	    peak_area += step * pi->getIntensity();
	  }
      }

    double h=area.max->getIntensity();
    double width=2*h/peak_area;

    PeakShape sech(h,area.max->getPosition()[mz_dim_],-1,width/2,width/2,peak_area,PeakShapeType::SECH_PEAK);
    sech.r_value=correlate_(sech,area,direction);


#ifdef DEBUG_PEAK_PICKING
    std::cout << "sech r: " << sech.r_value << " " << "pos: " << sech.mz_position << std::endl;
    std::cout << "w: " << width << std::endl;
    std::cout << "h: " << sech.height << std::endl;
#endif

    return sech;
  }

  template <Size D,  typename MapType, typename MapTypeOut>
  void DPeakPickerCWT<D, MapType, MapTypeOut>::changeCentroidPositions_(PeakData& sort_wrt_mz,
									std::vector<int>& index_of_peaks,
									double monoisotopic_position)
  {
    const double const_distance=1.0005;
    double sum=0,k=0,sqr_area;
    int m=index_of_peaks.size();
    int i,l;

    PeakIterator it=sort_wrt_mz.begin();
    for (i=0; i < m-1; ++i)
      {
	for (l=(index_of_peaks[i]+1); l <= index_of_peaks[i+1]; l++)
	  {
#ifdef DEBUG_PEAK_PICKING
	    std::cout << "old position " << (it+l)->getPosition()[mz_dim_];
#endif

	    (it+l)->setPosition(monoisotopic_position+i*const_distance)[mz_dim_];

#ifdef DEBUG_PEAK_PICKING
	    std::cout << " new position " << (it+l)->getPosition()[mz_dim_] << std::endl;
#endif

	  }
      }
  }

  template <Size D,  typename MapType, typename MapTypeOut>
  void DPeakPickerCWT<D, MapType, MapTypeOut>::subtractPeakShapeVector_(RawData& raw_array, const PeakShape& peak)
  {
    for (RawIterator pi=raw_array.begin(); pi<raw_array.end(); ++pi)
      pi->setIntensity(pi->getIntensity()-peak(pi->getPosition()[mz_dim_]));
  }


} // namespace OpenMS

#endif
