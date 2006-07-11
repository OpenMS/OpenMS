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
//

#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_DPEAKPICKER_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_DPEAKPICKER_H

#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/FORMAT/Param.h>
#include <OpenMS/KERNEL/DimensionDescription.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>

//#define DEBUG_PEAK_PICKING

namespace OpenMS
{
  /**
    @defgroup Transformations Transformations

    @brief Classes for the transformation of ms data.

    This module contains all classes that are involved in a data-reduction
    method (e.g. the transformation of raw data into peak data).
  */

  /**
     @brief This class is the base class for every peak picker.

     @ingroup PeakPicking
     
     @todo Make it work on all classes derived from DRawDataPoint (Eva)

  */
  template <Size D, 
                   typename MapType = MSExperiment<DRawDataPoint<1> >,
  				   typename MapTypeOut = MSExperiment<DPickedPeak<1> > > 
  class DPeakPicker
  {

  public:
    /** @name Type definitions
     */
    //@{
    typedef DimensionDescription < DimensionDescriptionTagLCMS > DimensionDescription;
    ///
    typedef DPeakArrayNonPolymorphic<D, DRawDataPoint<D> > RawData;
    ///
    typedef typename RawData::const_iterator RawConstIterator;
    ///
    typedef typename RawData::iterator RawIterator;
    ///
    typedef DPickedPeak<1> OutputPeak;
    ///
    typedef DSpectrum<1,DPeakArrayNonPolymorphic<1,OutputPeak> > Spectrum;
    ///
    typedef DPeakArray<D, DPickedPeak<D> > PeakData;
    /// NOTE: Output goes into an "ordinary" MSExperiment. We will have to see if that goes well (memory issues).
    typedef MapTypeOut MSExperimentPeakData;
	///
	typedef typename PeakData::iterator PeakIterator;
    ///
    typedef std::vector<RawIterator> IteratorVector;
    ///
    typedef typename IteratorVector::iterator IteratorIterator;
    ///

    //@}

    /** @name Constructors and Destructor
     */
    //@{
    DPeakPicker()
        : peak_bound_(200),
        peak_bound_ms2_level_(50),
        signal_to_noise_(3),
        peaks_(0),
        ms_exp_peaks_(0)
    {
      if (D == 1)
      {
        rt_dim_ = -1;
        mz_dim_ = 0;
      }
      else
        if (D == 2)
        {
          rt_dim_ = DimensionDescription::RT;
          mz_dim_ = DimensionDescription::MZ;
        }
    }
    ///
    DPeakPicker(const String& filename);
    ///
    DPeakPicker(const Param& parameters);
    ///
    DPeakPicker(const DPeakPicker& pp);
    ///
    virtual ~DPeakPicker()
    {   }
    //@}

    /** @name
     */
    //@{
    // The member variable peaks_ points to the output data
    DPeakPicker& operator()(PeakData& peaks);
    DPeakPicker& operator()(MSExperimentPeakData& ms_exp_peaks);
    ///
    virtual void pick(const MapType& ms_exp_raw)=0;
	///
    // The pick function picks 2D peaks in the intervall [first,last)
    virtual void pick(RawConstIterator first, RawConstIterator last, std::back_insert_iterator<PeakData >& output)=0;
    //@}

    /** @name Assignment
     */
    //@{
    DPeakPicker& operator=(const DPeakPicker& pp);
    //@}

    /** Accessors
     */
    //@{
    /// Non-mutable access to he mz dimension
    inline const int getMZdim() const { return mz_dim_; }
    /// Mutable access to the mz dimensin
    inline int getMZdim() { return mz_dim_; }
    /// Mutable access to the mz dimensin
    inline void setMZdim(const int mz_dim) { mz_dim_ = mz_dim; }

    /// Non-mutable access to the rt dimension
    inline const int getRTdim() const { return rt_dim_; }
    /// Mutable access to the rt dimensin
    inline int getRTdim() { return rt_dim_; }
    /// Mutable access to the rt dimensin
    inline void setRTdim(const int rt_dim) { rt_dim_ = rt_dim; }

    /// Non-mutable access to the noise level
    inline const float getPeakBound() const { return peak_bound_; }
		/// Mutable access to the noise level
    virtual void setPeakBound(const float peak_bound) 
		{ 
			peak_bound_ = peak_bound; 
		}

    /// Non-mutable access to the noise level
    inline const float getPeakBoundMs2Level() const { return peak_bound_ms2_level_; }
    /// Mutable access to the noise level
    inline float getPeakBoundMs2Level() { return peak_bound_ms2_level_; }
    /// Mutable access to the noise level
    inline void setPeakBoundMs2Level(const float peak_bound_ms2_level) { peak_bound_ms2_level_ = peak_bound_ms2_level; }

    /// Non-mutable access to the signal to noise level
    inline const float getSignalToNoiseLevel() const { return signal_to_noise_; }
    /// Mutable access to the signal to noise level
    inline float getSignalToNoiseLevel() { return signal_to_noise_; }
    /// Mutable access to the signal to noise value
    inline void setSignalToNoiseLevel(const float signal_to_noise) { signal_to_noise_ = signal_to_noise; }

    /// Non-mutable access to the picked peaks
    inline const PeakData& getPeakData() const { return *peaks_; }
    /// Mutable access to the picked peaks
    inline PeakData& getPeakData() { return *peaks_; }
    /// Mutable access to the picked peaks
    inline void setPeakData(const PeakData& peaks) { peaks_ = &peaks; }

    /// Non-mutable access to the picked peaks
    inline const MSExperimentPeakData& getMSExperimentPeakData() const { return *ms_exp_peaks_; }
    /// Mutable access to the picked peaks
    inline MSExperimentPeakData& getMSExperimentPeakData() { return *ms_exp_peaks_; }
    /// Mutable access to the picked peaks
    inline void setMSExperimentPeakData(const MSExperimentPeakData& peaks) { peaks_ = &peaks; }

    /// Non-mutable access to the parameter object
    inline const Param& getParam() const { return param_; }
    /// Mutable access to the parameter object
    inline Param& getParam() { return param_; }
    /// Mutable access to the parameter object
    inline void setParam(const Param& param) { param_ = param; }
    //@}

  protected:
    ///
    // MZ dimension
    int mz_dim_;
    ///
    // RT dimension
    int rt_dim_;
    ///
    // The noise level (threshold for peaks in the MS 1 level)
    float peak_bound_;

    // The noise level (threshold for peaks in the MS 2 level)
    float peak_bound_ms2_level_;
    ///
    // Signal to noise threshold
    float signal_to_noise_;
    ///
    // peaks_ points to the picked peaks in PeakData
    PeakData* peaks_;
    // ms_exp_peaks_ points to MSExperiment containing the picked peaks
    MSExperimentPeakData* ms_exp_peaks_;
    ///
    // Parameter object
    Param param_;
  };

  // Usage: "ms_experiment_raw >> picker(ms_experiment_peaks)"
  template <typename MapType, typename MapTypeOut>
  const MapTypeOut&
  operator>>(const  MapType& ms_exp_raw, DPeakPicker<1,MapType, MapTypeOut>& p)
  {
// 		(p.getMSExperimentPeakData()).setSample(ms_exp_raw.getSample());
// 		(p.getMSExperimentPeakData()).setSourceFile(ms_exp_raw.getSourceFile());
// 		(p.getMSExperimentPeakData()).setContacts(ms_exp_raw.getContacts());
// 		(p.getMSExperimentPeakData()).setInstrument(ms_exp_raw.getInstrument());
// 		(p.getMSExperimentPeakData()).setSoftware(ms_exp_raw.getSoftware());
// 		(p.getMSExperimentPeakData()).setProcessingMethod(ms_exp_raw.getProcessingMethod());
// 		(p.getMSExperimentPeakData()).setHPLC(ms_exp_raw.getHPLC());
// 		(p.getMSExperimentPeakData()).setType(ms_exp_raw.getType());
// 		(p.getMSExperimentPeakData()).setDate(ms_exp_raw.getDate());

		p.pick(ms_exp_raw);
		return p.getMSExperimentPeakData();  
  }

   // Usage: "input_raw_data >> picker(output_peak_data)"
  template <Size D, typename MapType, typename MapTypeOut>
  typename DPeakPicker<D,MapType,MapTypeOut>::PeakData&
  operator>>(const typename DPeakPicker<D,MapType,MapTypeOut>::RawData& raw, DPeakPicker<D,MapType,MapTypeOut>& p)
  {
    std::back_insert_iterator< typename DPeakPicker<D,MapType,MapTypeOut>::PeakData > bi(p.getPeakData());
    p.pick(raw.begin(),raw.end(),bi);
    return p.getPeakData();
  }
  
  template <Size D, typename MapType, typename MapTypeOut>
  DPeakPicker<D, MapType, MapTypeOut>::DPeakPicker(const String& filename)
  {
    param_.load(filename);

    if (D == 1)
    {
      rt_dim_ = -1;
      mz_dim_ = 0;
    }
    else
      if (D == 2)
      {
        rt_dim_ = DimensionDescription::RT;
        mz_dim_ = DimensionDescription::MZ;
      }

    // if a peak picking parameter is missed in the param object the value should be substituted by a default value
    DataValue dv;
    dv = param_.getValue("Thresholds:SignalToNoise");
    if (dv.isEmpty() || dv.toString() == "") signal_to_noise_ = 5;
    else signal_to_noise_ = (float)dv;

    dv = param_.getValue("Thresholds:PeakBound");
    if (dv.isEmpty() || dv.toString() == "") peak_bound_ = 200;
    else peak_bound_ = (float)dv;

    dv = param_.getValue("Thresholds:PeakBoundMs2Level");
    if (dv.isEmpty() || dv.toString() == "") peak_bound_ms2_level_ = 30;
    else peak_bound_ms2_level_ = (float)dv;

#ifdef DEBUG_PEAK_PICKING
    std::cout << "RT " << rt_dim_ << "\nMZ " << mz_dim_ << "\nPeak Level " << peak_bound_<<std::endl;
#endif

  }


  template <Size D, typename MapType, typename MapTypeOut>
  DPeakPicker<D, MapType, MapTypeOut>::DPeakPicker(const Param& parameters)
  {
    param_ = parameters;

    if (D == 1)
    {
      rt_dim_ = -1;
      mz_dim_ = 0;
    }
    else
      if (D == 2)
      {
        rt_dim_ = DimensionDescription::RT;
        mz_dim_ = DimensionDescription::MZ;
      }

    // if a peak picking parameter is missed in the param object the value should be substituted by a dv value
    DataValue dv;
    dv = param_.getValue("Thresholds:SignalToNoise");
    if (dv.isEmpty() || dv.toString() == "") signal_to_noise_ = 3;
    else signal_to_noise_ = (float)dv;

    dv = param_.getValue("Thresholds:PeakBound");
    if (dv.isEmpty() || dv.toString() == "") peak_bound_ = 200;
    else peak_bound_ = (float)dv;

    dv = param_.getValue("Thresholds:PeakBoundMs2Level");
    if (dv.isEmpty() || dv.toString() == "") peak_bound_ms2_level_ = 30;
    else peak_bound_ms2_level_ = (float)dv;
  }

  template <Size D, typename MapType, typename MapTypeOut>
  DPeakPicker<D, MapType, MapTypeOut>::DPeakPicker(const DPeakPicker<D, MapType, MapTypeOut>& pp)
      : mz_dim_(pp.mz_dim_),
      rt_dim_(pp.rt_dim_),
      peak_bound_(pp.peak_bound_),
      peak_bound_ms2_level_(pp.peak_bound_ms2_level_),
      signal_to_noise_(pp.signal_to_noise_)
  {
    // ???? copy of *peaks_ ????
    peaks_= pp.peaks_;
    ms_exp_peaks_=pp.ms_exp_peaks_;
  }

  template <Size D, typename MapType, typename MapTypeOut>
  DPeakPicker<D, MapType, MapTypeOut>& DPeakPicker<D, MapType, MapTypeOut>::operator= (const DPeakPicker<D, MapType, MapTypeOut>& pp)
  {
    // take care of self assignments
    if (this == &pp)
    {
      return *this;
    }
    mz_dim_=pp.mz_dim_;
    rt_dim_=pp.rt_dim_;
    peak_bound_=pp.peak_bound_;
    peak_bound_ms2_level_=pp.peak_bound_ms2_level_;
    signal_to_noise_=pp.signal_to_noise_;
    peaks_= pp.peaks_;
    ms_exp_peaks_=pp.ms_exp_peaks_;


    return *this;
  }

  template <Size D, typename MapType, typename MapTypeOut>
  DPeakPicker<D, MapType, MapTypeOut>& DPeakPicker<D, MapType, MapTypeOut>::operator()(typename DPeakPicker<D, MapType, MapTypeOut>::PeakData& peaks)
  {
    peaks_=&peaks;

    return *this;
  }

  template <Size D, typename MapType, typename MapTypeOut>
  DPeakPicker<D, MapType, MapTypeOut>& DPeakPicker<D, MapType, MapTypeOut>::operator()(MapTypeOut& ms_exp_peaks)
  {
    OPENMS_PRECONDITION(D == 1, "Use a one-dimensional Peak Picker for MSExperiment");
    ms_exp_peaks_=&ms_exp_peaks;

    return *this;
  }
  
  

}// namespace OpenMS

#endif
