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

#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_CONTINUOUSWAVELETTRANSFORMNUMINTEGRATION_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_CONTINUOUSWAVELETTRANSFORMNUMINTEGRATION_H

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/ContinuousWaveletTransform.h>

#include<math.h>

#ifdef DEBUG_PEAK_PICKING
#include<iostream>
#include<fstream>
#endif


//#define DEBUG_PEAK_PICKING
#undef DEBUG_PEAK_PICKING

namespace OpenMS
{
  /**
     @brief This class computes the continuous wavelet transformation using a marr wavelet.

     The convolution of the signal and the wavelet is computed by numerical integration.
     
     @ingroup PeakPicking
  */
  template <Size D>
  class ContinuousWaveletTransformNumIntegration : public ContinuousWaveletTransform<D>
  {
  public:
    /// Base type
    typedef ContinuousWaveletTransform<D> Base;
    /// Raw data const iterator type
    typedef typename Base::RawDataPointConstIterator RawDataPointConstIterator;

    using ContinuousWaveletTransform<D>::signal_;
    using ContinuousWaveletTransform<D>::wavelet_;
    using ContinuousWaveletTransform<D>::scale_;
    using ContinuousWaveletTransform<D>::spacing_;
    using ContinuousWaveletTransform<D>::end_left_padding_;
    using ContinuousWaveletTransform<D>::begin_right_padding_;
    using ContinuousWaveletTransform<D>::signal_length_;
    using ContinuousWaveletTransform<D>::mz_dim_;


    /// Constructor
    ContinuousWaveletTransformNumIntegration();

    /// Copy constructor
    ContinuousWaveletTransformNumIntegration(const ContinuousWaveletTransformNumIntegration& cwt)
        : ContinuousWaveletTransform<D>(cwt)
    {}
    /// Destructor.
    virtual ~ContinuousWaveletTransformNumIntegration();
    //@}


    /// Assignment operator
    inline ContinuousWaveletTransformNumIntegration& operator=(const ContinuousWaveletTransformNumIntegration& cwt)
    {
      // take care of self assignments
      if (this == &cwt)
      {
        return *this;
      }
      signal_=cwt.signal_;
      wavelet_=cwt.wavelet_;
      spacing_=cwt.spacing_;
      scale_=cwt.scale_;
      signal_length_=cwt.signal_length_;
      end_left_padding_=cwt.end_left_padding_;
      begin_right_padding_=cwt.begin_right_padding_;
      mz_dim_=cwt.mz_dim_;

      return *this;
    }
    /**
     @brief Computes the wavelet transform of a given raw data intervall [begin_input,end_input)
    */
    /// resolution = 1: the wavelet transform will be computed at every position of the raw data,
    /// resolution = 2: the wavelet transform will be computed at 2x(number of raw data positions) positions
    /// 								(the raw data are interpolated to get the intensity for missing positions)
    /// ...
    /// @note before starting the transformation you have to call the init function
    virtual void transform(RawDataPointConstIterator begin_input,
                           RawDataPointConstIterator end_input,
                           float resolution);

    /**
     @brief Perform necessary preprocessing steps like tabulating the Wavelet.
     
      Build a Marr-Wavelet for the current spacing and scale.
      	We store the wavelet in the vector<double> wavelet_;

      	We only need a finite amount of points since the Marr function
      	decays fast. We take 5*scale, since at that point the wavelet
      	has dropped to ~ -10^-4
    */
    virtual void init(double scale, double spacing, unsigned int mz_dim_);

  protected:

    /// Computes the convolution of the wavelet and the raw data at position x with resolution = 1
    double integrate_(RawDataPointConstIterator x, RawDataPointConstIterator first, RawDataPointConstIterator last);

    /// Computes the convolution of the wavelet and the raw data at position x with resolution > 1
    double integrate_(const std::vector<double>& processed_input, double spacing_data, int index);

    /// Computes the marr wavelet at position x
    double marr_(double x);
  };

  template <Size D>
  ContinuousWaveletTransformNumIntegration<D>::ContinuousWaveletTransformNumIntegration()
      : ContinuousWaveletTransform<D>()
  {
#ifdef DEBUG_PEAK_PICKING
    std::cout << "ContinuousWaveletTransformNumIntegration constructor." << std::endl;
#endif

  }

  template <Size D>
  ContinuousWaveletTransformNumIntegration<D>::~ContinuousWaveletTransformNumIntegration()
  {}

  template <Size D>
  void ContinuousWaveletTransformNumIntegration<D>::transform(typename ContinuousWaveletTransformNumIntegration<D>::RawDataPointConstIterator begin_input,
      typename ContinuousWaveletTransformNumIntegration<D>::RawDataPointConstIterator end_input,
      float resolution)
  {

#ifdef DEBUG_PEAK_PICKING
    std::cout << "ContinuousWaveletTransformNumIntegration::transform in dimension " << mz_dim_ << " from " << begin_input->getPosition()[mz_dim_] << " until " << (end_input-1)->getPosition()[mz_dim_] << std::endl;
#endif
    if (fabs(resolution-1) < 0.0001)
    {
      // resolution = 1 corresponds to the cwt at supporting points which have a distance corresponding to the minimal spacing in [begin_input,end_input)
      unsigned int n = distance(begin_input,end_input);
      signal_length_ = n;

      unsigned int i;

      signal_.clear();
      signal_.resize(n);

      // TODO avoid to compute the cwt for the zeros in signal
#ifdef DEBUG_PEAK_PICKING
      std::cout << "---------START TRANSFORM---------- \n";
#endif
      RawDataPointConstIterator help = begin_input;
      for (i=0; i < n; ++i)
      {
        signal_[i].getPos() = help->getPos();
        signal_[i].getIntensity()=integrate_(help,begin_input,end_input);
        ++help;
      }
#ifdef DEBUG_PEAK_PICKING
      std::cout << "---------END TRANSFORM----------" << std::endl;
#endif
      // no zeropadding
      begin_right_padding_=n;
      end_left_padding_=-1;
    }
    else
    {
      unsigned int n = (unsigned int) resolution * distance(begin_input, end_input);
      double origin  = begin_input->getPosition()[mz_dim_];
      double spacing = ((end_input-1)->getPosition()[mz_dim_]-origin)/(n-1);

      std::vector<double> processed_input(n);
      signal_.clear();
      signal_.resize(n);

      RawDataPointConstIterator it_help = begin_input;
      processed_input[0]=it_help->getPosition()[mz_dim_];

      double x;
      for (unsigned int k=1; k < n; ++k)
      {
        x = origin + k*spacing;
        // go to the real data point next to x
        while (((it_help+1) < end_input) && ((it_help+1)->getPosition()[mz_dim_] < x))
        {
          ++it_help;
        }
        processed_input[k] = getInterpolatedValue_(x,it_help);
      }

      // TODO avoid to compute the cwt for the zeros in signal
      for (unsigned int i=0; i < n; ++i)
      {
        signal_[i].getPos() = origin + i*spacing;
        signal_[i].getIntensity() = integrate_(processed_input,spacing,i);
      }

      // no zeropadding
      begin_right_padding_=n;
      end_left_padding_=-1;
    }
  }


  template <Size D>
  double ContinuousWaveletTransformNumIntegration<D>::integrate_(RawDataPointConstIterator x, RawDataPointConstIterator first, RawDataPointConstIterator last)
  {
#ifdef DEBUG_PEAK_PICKING
    std::cout << "integrate_" << std::endl;
#endif

    double v=0.;
    int middle = wavelet_.size();

    double start_pos = ((x->getPosition()[mz_dim_]-(middle*spacing_)) > first->getPosition()[mz_dim_]) ? (x->getPosition()[mz_dim_]-(middle*spacing_))
                       : first->getPosition()[mz_dim_];
    double end_pos = ((x->getPosition()[mz_dim_]+(middle*spacing_)) < (last-1)->getPosition()[mz_dim_]) ? (x->getPosition()[mz_dim_]+(middle*spacing_))
                     : (last-1)->getPosition()[mz_dim_];

    RawDataPointConstIterator help = x;

#ifdef DEBUG_PEAK_PICKING
    std::cout << "integrate from middle to start_pos "<< help->getPosition()[mz_dim_] << " until " << start_pos << std::endl;
#endif

    //integrate from middle to start_pos
    while ((help != first) && ((help-1)->getPosition()[mz_dim_] > start_pos))
    {
      // search for the corresponding datapoint of help in the wavelet (take the left most adjacent point)
      double distance = fabs(x->getPosition()[mz_dim_] - help->getPosition()[mz_dim_]);
      int index_w_r = (int)round(distance / spacing_);
      double wavelet_right =  wavelet_[index_w_r];

#ifdef DEBUG_PEAK_PICKING
      std::cout << "distance x help " << distance<< std::endl;
      std::cout << "distance in wavelet_ " << index_w_r*spacing_ << std::endl;
      std::cout << "wavelet_right "  <<  wavelet_right << std::endl;
#endif

      // search for the corresponding datapoint for (help-1) in the wavelet (take the left most adjacent point)
      distance = fabs(x->getPosition()[mz_dim_] - (help-1)->getPosition()[mz_dim_]);
      int index_w_l = (int)round(distance / spacing_);
      double wavelet_left =  wavelet_[index_w_l];

      // start the interpolation for the true value in the wavelet

#ifdef DEBUG_PEAK_PICKING
      std::cout << " help-1 " << (help-1)->getPosition()[mz_dim_] << " distance x, help-1" << distance << std::endl;
      std::cout << "distance in wavelet_ " << index_w_l*spacing_ << std::endl;
      std::cout << "wavelet_ at left " <<   wavelet_left << std::endl;

      std::cout << " intensity " << fabs((help-1)->getPosition()[mz_dim_]-help->getPosition()[mz_dim_]) / 2. << " * " << (help-1)->getIntensity() << " * " << wavelet_left <<" + " << (help)->getIntensity()<< "* " << wavelet_right
      << std::endl;
#endif

      v+= fabs((help-1)->getPosition()[mz_dim_]-help->getPosition()[mz_dim_]) / 2. * ((help-1)->getIntensity()*wavelet_left + help->getIntensity()*wavelet_right);
      --help;
    }


    //integrate from middle to end_pos
    help = x;
#ifdef DEBUG_PEAK_PICKING
    std::cout << "integrate from middle to endpos "<< (help)->getPosition()[mz_dim_] << " until " << end_pos << std::endl;
#endif
    while ((help != (last-1)) && ((help+1)->getPosition()[mz_dim_] < end_pos))
    {
      // search for the corresponding datapoint for help in the wavelet (take the left most adjacent point)
      double distance = fabs(x->getPosition()[mz_dim_] - help->getPosition()[mz_dim_]);
      int index_w_l = (int)round(distance / spacing_);
      double wavelet_left =  wavelet_[index_w_l];

#ifdef DEBUG_PEAK_PICKING
      std::cout << " help " << (help)->getPosition()[mz_dim_] << " distance x, help" << distance << std::endl;
      std::cout << "distance in wavelet_ " << index_w_l*spacing_ << std::endl;
      std::cout << "wavelet_ at left " <<   wavelet_left << std::endl;
#endif

      // search for the corresponding datapoint for (help+1) in the wavelet (take the left most adjacent point)
      distance = fabs(x->getPosition()[mz_dim_] - (help+1)->getPosition()[mz_dim_]);
      int index_w_r = (int)round(distance / spacing_);
      double wavelet_right =  wavelet_[index_w_r];

#ifdef DEBUG_PEAK_PICKING
      std::cout << " help+1 " << (help+1)->getPosition()[mz_dim_] << " distance x, help+1" << distance << std::endl;
      std::cout << "distance in wavelet_ " << index_w_r*spacing_ << std::endl;
      std::cout << "wavelet_ at right " <<   wavelet_right << std::endl;
#endif

      v+= fabs(help->getPosition()[mz_dim_] - (help+1)->getPosition()[mz_dim_]) / 2. * (help->getIntensity()*wavelet_left + (help+1)->getIntensity()*wavelet_right);
      ++help;
    }


#ifdef DEBUG_PEAK_PICKING
    std::cout << "return" << (v/sqrt(scale_)) << std::endl;
#endif
    return v / sqrt(scale_);
  }

  template <Size D>
  double ContinuousWaveletTransformNumIntegration<D>::marr_(double x)
  {
    return (1-x*x)*exp(-x*x/2);
  }


  template <Size D>
  double ContinuousWaveletTransformNumIntegration<D>::integrate_(const std::vector<double>& processed_input,
      double spacing_data,
      int index)
  {
    double v = 0.;
    int half_width = wavelet_.size();
    int index_in_data = (int)floor((half_width*spacing_) / spacing_data);
    int offset_data_left = ((index - index_in_data) < 0) ? 0 : (index-index_in_data);
    int offset_data_right = ((index + index_in_data) > (int)processed_input.size()) ? processed_input.size()-1 : (index+index_in_data);

    // integrate from i until offset_data_left
    for (int i = index; i > offset_data_left; --i)
    {
      int index_w_r = (int)round(((index-i)*spacing_data)/spacing_);
      int index_w_l = (int)round(((index-(i-1))*spacing_data)/spacing_);

      v += spacing_data / 2.*( processed_input[i]*wavelet_[index_w_r] + processed_input[i-1]*wavelet_[index_w_l] );
    }

    // integrate from i+1 until offset_data_right
    for (int i = index; i < offset_data_right; ++i)
    {
      int index_w_r = (int)round((((i+1)-index)*spacing_data)/spacing_);
      int index_w_l = (int)round(((i-index)*spacing_data)/spacing_);

      v += spacing_data / 2.*( processed_input[i+1]*wavelet_[index_w_r] + processed_input[i]*wavelet_[index_w_l]);
    }

    return v / sqrt(scale_);
  }


  template <Size D>
  void ContinuousWaveletTransformNumIntegration<D>::init(double scale, double spacing, unsigned int mz_dim)
  {
    ContinuousWaveletTransform<D>::init(scale, spacing, mz_dim);
    int number_of_points_right = (int)(ceil(5*scale_/spacing_));
    int number_of_points = number_of_points_right + 1;
    wavelet_.resize(number_of_points);
    wavelet_[0] = 1.;

    for (int i=1; i<number_of_points; i++)
    {
      wavelet_[i] = marr_(i*spacing_/scale_ );
    }

#ifdef DEBUG_PEAK_PICKING
    std::cout << "WAVELET" << std::endl;
    for (int i=0; i<number_of_points; i++)
    {
      std::cout <<  i*spacing_ << " " << wavelet_[i] << std::endl;
    }
#endif

  }
} //namespace OpenMS
#endif
