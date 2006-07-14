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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKPICKERCWT_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKPICKERCWT_H

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPicker.h>
#include <OpenMS/KERNEL/DPickedPeak.h>
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

//#define DEBUG_PEAK_PICKING
namespace OpenMS
{
  /**
     @defgroup PeakPicking PeakPicking

     @brief Classes for the transformation of raw ms data into peak data.

     This module contains all important classes that are involved in the peak picking as described by Lange et al. (2006) Proc. PSB-06.

     @ingroup Transformations
  */

  /**
     @brief This class implements a peak picking algorithm using wavelet techniques.

     This peak picking algorithm uses the continuous wavelet transform of a raw data signal to detect mass peaks.
     Afterwards a given asymmetric peak function is fitted to the raw data and important peak parameters (e.g. fwhm)
     are extracted.
     In an optional step these parameters can be optimized using a non-linear opimization method.

     @ingroup PeakPicking

     @todo write test (Eva)

  */

  class PeakPickerCWT : public PeakPicker
  {
  public:
    /** @name Type definitions
      */
    //@{
    typedef PeakPicker BaseType;
    typedef DRawDataPoint<1> RawDataPointType;
    typedef DPeakArrayNonPolymorphic<1, RawDataPointType > RawDataArrayType;
    typedef RawDataArrayType::iterator RawDataPointIterator;
    typedef DPosition<1> PositionType;
    ///
    //@}

    using PeakPicker::param_;
    using PeakPicker::peak_bound_;
    using PeakPicker::peak_bound_ms2_level_;
    using PeakPicker::signal_to_noise_;

    /** @name Constructors and Destructor
        */
    //@{
    PeakPickerCWT();
    ///
    PeakPickerCWT(const Param& parameters);
    ///
    PeakPickerCWT(const String& filename);
    ///
    PeakPickerCWT(const PeakPickerCWT& pp)
        : PeakPicker(pp),
        radius_(pp.radius_),
        scale_(pp.scale_),
        peak_bound_cwt_(pp.peak_bound_cwt_),
        peak_bound_ms2_level_cwt_(pp.peak_bound_ms2_level_cwt_),
        peak_corr_bound_(pp.peak_corr_bound_),
        noise_level_(pp.noise_level_),
        optimization_(pp.optimization_)
    {}
    ///
    virtual ~PeakPickerCWT();
    ///
    void init();
    //@}

    /** @name Assignment
     */
    //@{
    inline PeakPickerCWT& operator=(const PeakPickerCWT& pp)
    {
      param_ = pp.param_;
      peak_bound_ = pp.peak_bound_;
      signal_to_noise_ = pp.signal_to_noise_;
      scale_ = pp.scale_;
      peak_bound_cwt_ = pp.peak_bound_cwt_;
      peak_bound_ms2_level_cwt_ = pp.peak_bound_ms2_level_cwt_;
      radius_ = pp.radius_;
      peak_corr_bound_ = pp.peak_corr_bound_;
      noise_level_ = pp.noise_level_;
      optimization_ = pp.optimization_;

      return *this;
    }
    //@}


    /** Accessors
     */
    //@{


    /// Non-mutable access to the vector with peak shapes
    inline const std::vector<PeakShape>& getPeakShapes() const { return peak_shapes_; }
    /// Mutable access to the vector with peak shapes
    inline std::vector<PeakShape>& getPeakShapes() { return peak_shapes_;  }
    /// Mutable access to the vector with peak shapes
    inline void setPeakShapes(const std::vector<PeakShape>&  peak_shapes) { peak_shapes_ = peak_shapes; }

    /// Non-mutable access to the wavelet transform
    inline const ContinuousWaveletTransform<1>& getWaveletTransform() const { return *wt_; }
    /// Mutable access to the wavelet transform
    inline ContinuousWaveletTransform<1>& getWaveletTransform() { return *wt_;  }
    /// Mutable access to the wavelet transform
    inline void setWaveletTransform(const ContinuousWaveletTransform<1>& wt)
    {
      if (wt_)
      {
        *wt_ = wt;
      }
      else
      {
        wt_= new ContinuousWaveletTransformNumIntegration<1>();
        *wt_ = wt;
      }
    }

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
    inline void setPeakBoundCWT(const float peak_bound_cwt) { peak_bound_cwt_ = peak_bound_cwt; }

    /// Non-mutable access to the peak bound in the wavelet transform
    inline const float getPeakBoundMs2LevelCWT() const { return peak_bound_ms2_level_cwt_; }
    /// Mutable access to the peak bound in the wavelet transform
    inline float getPeakBoundMs2LevelCWT() { return peak_bound_ms2_level_cwt_; }
    /// Mutable access to the peak bound in the wavelet transform
    inline void setPeakBoundMs2LevelCWT(const float peak_bound_ms2_level_cwt) { peak_bound_ms2_level_cwt_ = peak_bound_ms2_level_cwt; }

    /// Non-mutable access to the minimum peak correlation coefficient
    inline const float getPeakCorrBound() const { return peak_corr_bound_; }
    /// Mutable access to the minimum peak correlation coefficient
    inline float getPeakCorrBound() { return peak_corr_bound_; }
    /// Mutable access to the minimum peak correlation coefficient
    inline void setPeakCorrBound(const float peak_corr_bound) { peak_corr_bound_ = peak_corr_bound; }

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
    //@}


    template <typename InputPeakIterator, typename OutputPeakContainer>
    void pick(InputPeakIterator first, InputPeakIterator last, OutputPeakContainer& picked_peak_container, int ms_level = 1)
    {
      //DSignalToNoiseEstimatorWindowing<InputPeakType> sne;
      typedef typename OutputPeakContainer::value_type OutputPeakType;

      /// Initialize the wavelet transform
      double wavelet_spacing;
      DataValue dv = param_.getValue("WaveletTransform:Spacing");
      if (dv.isEmpty() || dv.toString() == "") wavelet_spacing= 0.001;
      else wavelet_spacing = (double)dv;

      wt_->init(scale_, wavelet_spacing, 0);

#ifdef DEBUG_PEAK_PICKING
      std::cout << "****************** PICK ******************" << std::endl;
#endif

      // vector of peak endpoint positions
      std::vector<double> peak_endpoints;

#ifdef GSL_DEF
      if (optimization_)
      {
        OptimizationFunctions::positions_.clear();
        OptimizationFunctions::signal_.clear();
      }
#endif

      // copy the raw data into a DPeakArrayNonPolymorphic<DRawDataPoint<D> >
      RawDataArrayType raw_peak_array;
      unsigned int n = distance(first, last);
      raw_peak_array.resize(n);

      for (unsigned int i = 0; i < n; ++i)
      {
        RawDataPointType raw_data_point;
        raw_data_point.getIntensity() = (first + i)->getIntensity();
        raw_data_point.getPosition() = (first + i)->getPosition();
        raw_peak_array[i] = raw_data_point;
      }

      RawDataPointIterator it_pick_begin = raw_peak_array.begin();
      RawDataPointIterator it_pick_end   = raw_peak_array.end();

#ifdef GSL_DEF
      if (optimization_)
      {
        unsigned int l=raw_peak_array.size();

        OptimizationFunctions::positions_.resize(l);
        OptimizationFunctions::signal_.resize(l);

        for (unsigned int i = 0; i < l ;++i)
        {
          OptimizationFunctions::positions_[i] = raw_peak_array[i].getPos();
          OptimizationFunctions::signal_[i] = raw_peak_array[i].getIntensity();
        }
      }
#endif

      // Points to the actual maximum position in the raw data
      RawDataPointIterator it_max_pos;

      // start the peak picking until no more maxima can be found in the wavelet transform
      unsigned int number_of_peaks = 0;
      do
      {
        number_of_peaks = 0;
        int peak_left_index, peak_right_index;

        // compute the continious wavelet transform with resolution 1
        wt_->transform(it_pick_begin, it_pick_end,1.);

        PeakArea_ area;
        bool centroid_fit=false;
        bool regular_endpoints=true;

        // search for maximum positions in the cwt and extract potential peaks
        int direction=1;
        int distance_from_scan_border = 0;
        while ((distance(it_pick_begin, it_pick_end) > 3) && getMaxPosition_(it_pick_begin, it_pick_end, wt_, area, distance_from_scan_border, ms_level, direction))
        {
          //search for the endpoints of the peak
          regular_endpoints = getPeakEndPoints_(it_pick_begin,
                                                it_pick_end,
                                                area,
                                                peak_left_index,
                                                peak_right_index);

          // compute the centroid position
          getPeakCentroid_(area);

          // if the peak achieves a minimal width, start the peak fitting
          if (regular_endpoints)
          {
#ifdef DEBUG_PEAK_PICKING
            std::cout << "The endpoints are "
            << area.left->getPosition()
            << " and "
            << area.right->getPosition()
            << std::endl;
#endif
            // determine the best fitting lorezian or sech2 function
            PeakShape shape = fitPeakShape_(area,centroid_fit);

            // Use the centroid for Optimization
            shape.mz_position=area.centroid_position[0];

            // TEST!!!!!
            if ( (shape.r_value > peak_corr_bound_) )//&& ((sne.getSignalToNoise(area.max)) >= signal_to_noise_))
            {
              //  shape.getSymmetricMeasure();
              // shape.signal_to_noise = sne.getSignalToNoise(area.max);
              peak_shapes_.push_back(shape);
              peak_endpoints.push_back(area.left->getPos());
              peak_endpoints.push_back(area.right->getPos());
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
          for (RawDataPointIterator pi=area.left; pi!=area.right+1; pi++)
          {
            pi->setIntensity(0.);
          }

          // search for the next peak
          it_pick_begin = area.right;
          distance_from_scan_border = distance(raw_peak_array.begin(),it_pick_begin);

        } //end while (getMaxPosition_(it_pick_begin, it_pick_end, wt_, area, distance_from_scan_border, ms_level, direction))
        it_pick_begin = raw_peak_array.begin();
      }
      while (number_of_peaks != 0);

      // start the nonlinear optimization for all peaks in split
#ifdef DEBUG_PEAK_PICKING
      std::cout << "Try the optimization run... with " << peak_shapes_.size() << std::endl;
#endif

      if (peak_shapes_.size() > 0)
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
          opt.optimize(peak_shapes_);

          // compute the new correlation coefficients
          for (unsigned int i=0, j=0; i < peak_shapes_.size(); ++i, j+=2)
          {
            peak_shapes_[i].r_value=opt.correlate(peak_shapes_[i],peak_endpoints[j], peak_endpoints[j+1]);
          }
        } // if optimization
#endif

        // write the picked peaks to the outputcontainer
        for (unsigned int i = 0; i < peak_shapes_.size(); ++i)
        {
          OutputPeakType picked_peak;

          picked_peak.getIntensity() = peak_shapes_[i].height;
          picked_peak.getPos() = peak_shapes_[i].mz_position;

          fillPeak_(peak_shapes_[i],picked_peak);
	  picked_peak_container.push_back(picked_peak);
        }

      } // if (peak_shapes_.size() > 0)
    }

    template <typename InputPeakContainer, typename OutputPeakContainer>
    void pick(const InputPeakContainer& input_peak_container, OutputPeakContainer& picked_peaks_container)
    {
      pick(input_peak_container.begin(), input_peak_container.end(), picked_peaks_container);
    }

    template <typename InputSpectrumIterator, typename OutputSpectrumContainer>
    void pickExperiment(InputSpectrumIterator first, InputSpectrumIterator last, OutputSpectrumContainer& spectrum_container)
    {
      typedef typename OutputSpectrumContainer::value_type SpectrumType;

      InputSpectrumIterator help = first;
      while (help != last)
      {
        SpectrumType spectrum;
        pick(*help,spectrum);
        spectrum_container.push_back(spectrum);
      }
    }

    template <typename InputSpectrumContainer, typename OutputSpectrumContainer>
    void pickExperiment(const InputSpectrumContainer& input_peak_container, OutputSpectrumContainer& spectrum_container)
    {
      pickExperiment(input_peak_container.begin(),input_peak_container.end(),spectrum_container);
    }

    template <typename OutputPeakType>
    void fillPeak_(const PeakShape& /* peak_shape */, OutputPeakType& /* picked_peak */)
    {}

  protected:
    /** @name Attributes
     */
    //@{
    // Vector of the detected peak shapes
    std::vector<PeakShape> peak_shapes_;
    ///
    // This member perfomrs the continuous wavelet transform
    ContinuousWaveletTransform<1>* wt_;
    ///
    // The search radius for the determination of a peak's maximum position
    unsigned int radius_;
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
    // The threshold for correlation
    float peak_corr_bound_;
    ///
    // The threshold for the noise level (TODO: Use the information of the signal to noise estimator)
    float noise_level_;
    ///
    // Use the optimization of peak parameters
    bool optimization_;
    //@}



    /** A regularData-Object which contents some additional useful informations
       for analysing peaks and their properties
       */
    class PeakArea_
    {
      typedef std::vector<RawDataPointType>::iterator RawDataPointIterator;

    public:
      PeakArea_(){}

      RawDataPointIterator left, max, right, left_behind_centroid;
      DPosition<1> centroid_position;
    };

    /** @name Functions for Stick Conversion
     */
    //@{
    void getPeakArea_(const PeakArea_& area, double &area_left, double &area_right);

    /** Fit a peak form.
     */
    PeakShape fitPeakShape_(const PeakArea_& area, bool enable_centroid_fit);

    /** Returns a measure for the similarity of s and phi. For that, the increase of the
    linear regression of the points (s(i),phi(i)) for each grid-point i in s is computed
    If the value is near 1, s and phi are expected to be very similar.
    */
    double correlate_(const PeakShape& peak, const PeakArea_& area, int direction=0) const;
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
    bool getMaxPosition_(RawDataPointIterator first,
                         RawDataPointIterator last,
                         ContinuousWaveletTransform<1>* wt,
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
    bool getPeakEndPoints_(RawDataPointIterator first, RawDataPointIterator last,
                           PeakArea_ &area, int& peak_left_index, int& peak_right_index);


    /**
     *  Compute the transformation of a peak area with a higher resolution (more supporting points) and define the maximum of
     *  the transformed peak as the centroid.
     *  If the calculated centroid lies between two data point positions the function returns true otherwise return false.
     */
    void getPeakCentroid_(PeakArea_& area);
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
    RawDataPointIterator getIteratorLeftDataPoint_(RawDataPointIterator first,
        RawDataPointIterator last,
        double value);
    //@}
  }
  ; // end PeakPickerCWT


  template <>
  void PeakPickerCWT::fillPeak_< DPickedPeak<1> >(const PeakShape& peak_shape, DPickedPeak<1>& picked_peak);


}// namespace OpenMS

#endif
