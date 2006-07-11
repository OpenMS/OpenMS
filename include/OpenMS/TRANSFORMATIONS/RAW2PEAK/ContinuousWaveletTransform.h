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

#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_CONTINUOUSWAVELETTRANSFORM_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_CONTINUOUSWAVELETTRANSFORM_H

#include <OpenMS/KERNEL/DPeakArrayNonPolymorphic.h>
#include <OpenMS/KERNEL/DRawDataPoint.h>

#include <vector>
#include <iostream>

//#define DEBUG_PEAK_PICKING

namespace OpenMS
{
  /**
  	@brief This class is the base class of the continuous wavelet transformation.
   
		@todo move the typenames into their corresponding header files! (Eva)
		@todo remove depencies on DPeakPicker (Eva)
  */
  template <Size D>
  class ContinuousWaveletTransform
  {
    /** @name Type definitions
    	
    */
    //@{
    typedef std::vector<DRawDataPoint<D> > RawDataVector;
    ///
    typedef typename RawDataVector::iterator RawDataPointIterator;
    ///
    typedef typename RawDataVector::const_iterator RawDataPointConstIterator;
    //@}

  public:

    /** @name Constructors and Destructor
       */
    //@{
    /** Default constructor. */
    ContinuousWaveletTransform();

    /** Destructor. */
    virtual ~ContinuousWaveletTransform();
    //@}


    /** @name Assignement
     */
    //@{
    ///
    inline ContinuousWaveletTransform& operator=(const ContinuousWaveletTransform& cwt)
    {
      signal_=cwt.signal_;
      wavelet_=cwt.wavelet_;
      scale_=cwt.scale_;
      spacing_=cwt.spacing_;
      signal_length_=cwt.signal_length_;
      end_left_padding_=cwt.end_left_padding_;
      begin_right_padding_=cwt.begin_right_padding_;
    }
    //@}

    /** Accessors
     */
    //@{
    /// Non-mutable access to the wavelet transform of the signal
    inline const DPeakArrayNonPolymorphic<2, DRawDataPoint<2> >& getSignal() const { return signal_; }
    /// Mutable access to the wavelet transform of the signal
    inline DPeakArrayNonPolymorphic<2, DRawDataPoint<2> >& getSignal() { return signal_; }
    /// Mutable access to the wavelet transform of the signal
    inline void setSignal(const DPeakArrayNonPolymorphic<2, DRawDataPoint<2> >& signal) { signal_ = signal; }

    /// Non-mutable access to the wavelet
    inline const std::vector<double>& getWavelet() const { return wavelet_; }
    /// Mutable access to the wavelet
    inline std::vector<double>& getWavelet() { return wavelet_; }
    /// Mutable access to the signal
    inline void setWavelet(const std::vector<double>& wavelet) { wavelet_ = wavelet; }

    /// Non-mutable access to the resolution
    inline const float getResolution() const { return 1.; }

    // Non-mutable access to the spacing of raw data
    inline const double getSpacing() const { return spacing_; }
    /// Mutable access to the spacing of raw data
    inline double getSpacing() { return spacing_; }
    /// Mutable access to the spacing of raw data
    inline void setSpacing(const double spacing) { spacing_ = spacing; }

    /// Non-mutable access to the position where the signal starts (in the intervall [0,end_left_padding_) are the padded zeros)
    inline const int getLeftPaddingIndex() const { return end_left_padding_; }
    /// Mutable access to the position where the signal starts
    inline int getLeftPaddingIndex() { return end_left_padding_; }
    /// Mutable access to position where the signal starts
    inline void setLeftPaddingIndex(const int end_left_padding) { end_left_padding_ = end_left_padding; }

    /// Non-mutable access to the position where the signal ends (in the intervall (begin_right_padding_,end] are the padded zeros)
    inline const int getRightPaddingIndex() const { return begin_right_padding_; }
    /// Mutable access to the position where the signal starts
    inline int getRightPaddingIndex() { return begin_right_padding_; }
    /// Mutable access to position where the signal starts
    inline void setRightPaddingIndex(const int begin_right_padding) { begin_right_padding_ = begin_right_padding; }

    /// Non-mutable access to signal length [end_left_padding,begin_right_padding]
    inline const int getSignalLength() const { return signal_length_; }
    /// Mutable access to signal length [end_left_padding,begin_right_padding]
    inline int getSignalLength() { return signal_length_; }
    /// Mutable access to signal length [end_left_padding,begin_right_padding]
    inline void setSignalLength(const int signal_length) { signal_length_ = signal_length; }

    /// Non-mutable access to signal length including padded zeros [0,end]
    inline const int getSize() const { return signal_.size(); }
    //@}

    /** Compute the continuous wavelet transformation (using a marr wavelet)
     * of the signal intervall [it_begin, it_end).
     */
    virtual void transform(RawDataPointIterator begin_input,
                           RawDataPointIterator end_input,
                           float resolution) = 0;

    /**
     *  Perform possibly necessary preprocessing steps, like tabulating the Wavelet.
     */
    virtual void init(double scale, double spacing, unsigned int MZ);

    /** Yields the signal at position i **/
    double& operator [] (const unsigned int i);
    const double& operator [] (const unsigned int i) const;

    /** Return the interpolated value out of the input array **/
    double getInterpolatedValue(double x, RawDataPointIterator it_left);

  protected:
    /** The transformed signal. **/
    DPeakArrayNonPolymorphic<2, DRawDataPoint<2> > signal_;

    /** The wavelet used for the transform **/
    std::vector<double> wavelet_;

    /** Spacing and scale of the wavelet and length of the signal. **/
    double scale_;
    double spacing_;
    int signal_length_;

    /** We often have to pad the transform at the left and right with
     *  zeros. Since we don't want to iterate over those as well, we
     *  have to store their positions.
     */
    int end_left_padding_;
    int begin_right_padding_;

    /** MZ dimension */
    unsigned int mz_dim_;
  };

  template <Size D>
  ContinuousWaveletTransform<D>::ContinuousWaveletTransform()
  {}

  template <Size D>
  ContinuousWaveletTransform<D>::~ContinuousWaveletTransform() {}

  template <Size D>
  double& ContinuousWaveletTransform<D>::operator [] (const unsigned int i)
  {
    return signal_[i].getIntensity();
  }

  template <Size D>
  const double& ContinuousWaveletTransform<D>::operator [] (const unsigned int i) const
  {
    return signal_[i].getIntensity();
  }

  template <Size D>
  double ContinuousWaveletTransform<D>::getInterpolatedValue(double x,
      RawDataPointIterator it_left)
  {
    // Interpolate between the point to the left and the point to the right.
    double left_position = it_left->getPosition()[mz_dim_];
    double right_position = (it_left+1)->getPosition()[mz_dim_];
    double d=(x-left_position)/(right_position-left_position);

    return ((it_left+1)->getIntensity()*d+it_left->getIntensity()*(1-d));
  }

  template <Size D>
  void ContinuousWaveletTransform<D>::init(double scale, double spacing, unsigned int MZ)
  {
    scale_ = scale;
    spacing_=spacing;
    mz_dim_ = MZ;
  }
}

#endif

