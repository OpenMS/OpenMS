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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_CONTINUOUSWAVELETTRANSFORM_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_CONTINUOUSWAVELETTRANSFORM_H

#include <OpenMS/KERNEL/DPeakArray.h>
#include <OpenMS/KERNEL/DRawDataPoint.h>

#include <vector>
#include <iostream>

namespace OpenMS
{
/**
	@brief This class is the base class of the continuous wavelet transformation. 
*/
class ContinuousWaveletTransform
{
public:
    /// Raw data const iterator type
    typedef std::vector<DRawDataPoint<1> >::const_iterator RawDataPointConstIterator;


    /// Constructor
    ContinuousWaveletTransform()
            : scale_(0),
            spacing_(0),
            signal_length_(0),
            end_left_padding_(0),
            begin_right_padding_(0),
            mz_dim_(0)
    {}

    /// Copy constructor
    ContinuousWaveletTransform(const ContinuousWaveletTransform& cwt)
            : signal_(cwt.signal_),
            wavelet_(cwt.wavelet_),
            scale_(cwt.scale_),
            spacing_(cwt.spacing_),
            signal_length_(cwt.signal_length_),
            end_left_padding_(cwt.end_left_padding_),
            begin_right_padding_(cwt.begin_right_padding_),
            mz_dim_(cwt.mz_dim_)
    {}

    /// Destructor.
    virtual ~ContinuousWaveletTransform()
    {}

    /// Assignment operator
    inline ContinuousWaveletTransform& operator=(const ContinuousWaveletTransform& cwt)
    {
        // take care of self assignments
        if (this == &cwt)
        {
            return *this;
        }

        signal_=cwt.signal_;
        wavelet_=cwt.wavelet_;
        scale_=cwt.scale_;
        spacing_=cwt.spacing_;
        signal_length_=cwt.signal_length_;
        end_left_padding_=cwt.end_left_padding_;
        begin_right_padding_=cwt.begin_right_padding_;
        mz_dim_ = cwt.mz_dim_;

        return *this;
    }

    /// Non-mutable access to the wavelet transform of the signal
    inline const DPeakArray<1, DRawDataPoint<1> >& getSignal() const
    {
        return signal_;
    }
    /// Mutable access to the wavelet transform of the signal
    inline DPeakArray<1, DRawDataPoint<1> >& getSignal()
    {
        return signal_;
    }
    /// Mutable access to the wavelet transform of the signal
    inline void setSignal(const DPeakArray<1, DRawDataPoint<1> >& signal)
    {
        signal_ = signal;
    }

    /// Non-mutable access to the wavelet
    inline const std::vector<double>& getWavelet() const
    {
        return wavelet_;
    }
    /// Mutable access to the wavelet
    inline std::vector<double>& getWavelet()
    {
        return wavelet_;
    }
    /// Mutable access to the signal
    inline void setWavelet(const std::vector<double>& wavelet)
    {
        wavelet_ = wavelet;
    }

    // Non-mutable access to the scale of the wavelet
    inline const double& getScale() const
    {
        return scale_;
    }
    /// Mutable access to the spacing of raw data
    inline double& getScale()
    {
        return scale_;
    }
    /// Mutable access to the spacing of raw data
    inline void setScale(const double& scale)
    {
        scale_ = scale;
    }

    // Non-mutable access to the spacing of raw data
    inline const double& getSpacing() const
    {
        return spacing_;
    }
    /// Mutable access to the spacing of raw data
    inline double& getSpacing()
    {
        return spacing_;
    }
    /// Mutable access to the spacing of raw data
    inline void setSpacing(const double spacing)
    {
        spacing_ = spacing;
    }

    /// Non-mutable access to the position where the signal starts (in the intervall [0,end_left_padding_) are the padded zeros)
    inline const int& getLeftPaddingIndex() const
    {
        return end_left_padding_;
    }
    /// Mutable access to the position where the signal starts
    inline int& getLeftPaddingIndex()
    {
        return end_left_padding_;
    }
    /// Mutable access to position where the signal starts
    inline void setLeftPaddingIndex(const int end_left_padding)
    {
        end_left_padding_ = end_left_padding;
    }

    /// Non-mutable access to the position where the signal ends (in the intervall (begin_right_padding_,end] are the padded zeros)
    inline const int& getRightPaddingIndex() const
    {
        return begin_right_padding_;
    }
    /// Mutable access to the position where the signal starts
    inline int& getRightPaddingIndex()
    {
        return begin_right_padding_;
    }
    /// Mutable access to position where the signal starts
    inline void setRightPaddingIndex(const int begin_right_padding)
    {
        begin_right_padding_ = begin_right_padding;
    }

    /// Non-mutable access to signal length [end_left_padding,begin_right_padding]
    inline const int& getSignalLength() const
    {
        return signal_length_;
    }
    /// Mutable access to signal length [end_left_padding,begin_right_padding]
    inline int& getSignalLength()
    {
        return signal_length_;
    }
    /// Mutable access to signal length [end_left_padding,begin_right_padding]
    inline void setSignalLength(const int signal_length)
    {
        signal_length_ = signal_length;
    }

    /// Non-mutable access to signal length including padded zeros [0,end]
    inline int getSize() const
    {
        return signal_.size();
    }


    /**
     @brief Perform possibly necessary preprocessing steps, like tabulating the Wavelet.
    */
    virtual void init(double scale, double spacing);


    /// Yields the signal (intensity) at position i
    inline double& operator [] (const unsigned int i)
    {
        return signal_[i].getIntensity();
    }

    inline const double& operator [] (const unsigned int i) const
    {
        return signal_[i].getIntensity();
    }
		
		template < typename InputPeakIterator >
		double getInterpolatedValue_(double x, InputPeakIterator it_left)
		{
    	// Interpolate between the point to the left and the point to the right.
    	double left_position = it_left->getPosition()[mz_dim_];
    	double right_position = (it_left+1)->getPosition()[mz_dim_];
    	double d=(x-left_position)/(right_position-left_position);

    	return ((it_left+1)->getIntensity()*d+it_left->getIntensity()*(1-d));
		}



protected:
    /// The transformed signal
    DPeakArray<1, DRawDataPoint<1> > signal_;

    /// The pretabulated wavelet used for the transform
    std::vector<double> wavelet_;

    /// Spacing and scale of the wavelet and length of the signal.
    double scale_;
    double spacing_;
    int signal_length_;

    /// We often have to pad the transform at the left and right with
    /// zeros. Since we don't want to iterate over those as well, we
    /// have to store their positions.
    int end_left_padding_;
    int begin_right_padding_;

    /// The index of the mass to charge dimension
    unsigned int mz_dim_;

    /// Computes the interpolated value at position x (mz) given the iterator it_left, which points
    /// to the left neighbour raw data point of x in the original data
    double getInterpolatedValue_(double x, RawDataPointConstIterator it_left);
};

} //namespace OpenMS

#endif

