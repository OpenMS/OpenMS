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
// $Id: ContinuousWaveletTransformNumIntegration.h,v 1.22 2006/06/02 14:04:17 elange Exp $
// $Author: elange $
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

     @ingroup PeakPickingCWT

     @todo write test (Eva)
  */
  template <Size D>
  class ContinuousWaveletTransformNumIntegration : public ContinuousWaveletTransform<D>
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

    using ContinuousWaveletTransform<D>::signal_;
    using ContinuousWaveletTransform<D>::wavelet_;
    using ContinuousWaveletTransform<D>::scale_;
    using ContinuousWaveletTransform<D>::spacing_;
    using ContinuousWaveletTransform<D>::end_left_padding_;
    using ContinuousWaveletTransform<D>::begin_right_padding_;
    using ContinuousWaveletTransform<D>::signal_length_;
    using ContinuousWaveletTransform<D>::mz_dim_;

  public:
    /** @name Constructors and Destructor
     */
    //@{
    /// Default constructor.
    ContinuousWaveletTransformNumIntegration();

    /// Destructor.
    virtual ~ContinuousWaveletTransformNumIntegration();
    //@}


    /** @name Assignement
     */
    //@{
    ///
    inline ContinuousWaveletTransformNumIntegration& operator=(const ContinuousWaveletTransformNumIntegration& cwt)
    {
      signal_=cwt.signal_;
      wavelet_=cwt.wavelet_;
      spacing_=cwt.spacing_;
      scale_=cwt.scale_;
      signal_length_=cwt.signal_length_;
      end_left_padding_=cwt.end_left_padding_;
      begin_right_padding_=cwt.begin_right_padding_;
      mz_dim_=cwt.mz_dim_;
    }
    //@}


    // convolution with the marr wavelet (resolution = 1)
    double integrate(RawDataPointConstIterator x, RawDataPointConstIterator first, RawDataPointConstIterator last);

    // convolution with the marr wavelet (resolution > 1)
    double integrate(const std::vector<double>& processed_input, double spacing_data, int index);

    double marr(double x);

    /** Compute the continuous wavelet transformation (use a marr wavelet)
     * of the signal intervall [it_begin, it_end).
     * Using the numerical integration (trapezian rule) to compute the convolution of the
     * signal and the marr wavelet.
     */
    virtual
    void transform(RawDataPointIterator begin_input,
                   RawDataPointIterator end_input,
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

	//std::cout << "---------START TRANSFORM---------- \n";
        RawDataPointIterator help = begin_input;
        for (i=0; i < n; ++i)
        {
          signal_[i].getPos() = help->getPos();
          signal_[i].getIntensity()=integrate(help,begin_input,end_input);
	  //std::cout << signal_[i].getPos() << ' ' << signal_[i].getIntensity() << '\n';
          ++help;
        }
	//std::cout << "---------END TRANSFORM----------" << std::endl;

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

        RawDataPointIterator it_help = begin_input;
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
          processed_input[k]=getInterpolatedValue(x,it_help);
        }

        // TODO avoid to compute the cwt for the zeros in signal
        for (unsigned int i=0; i < n; ++i)
        {
          signal_[i].getPos() = origin + i*spacing;
          signal_[i].getIntensity() = integrate(processed_input,spacing,i);
					//  std::cout << origin + i*spacing << " " << signal_[i].getIntensity() << std::endl;
        }

        // no zeropadding
        begin_right_padding_=n;
        end_left_padding_=-1;
      }
    }

    /// Initialize the Wavelet
    virtual void init(double scale, double spacing, unsigned int mz_dim_);

  protected:
  };

  /// Default constructor.
  template <Size D>
  ContinuousWaveletTransformNumIntegration<D>::ContinuousWaveletTransformNumIntegration()
      : ContinuousWaveletTransform<D>()
  {
#ifdef DEBUG_PEAK_PICKING
    std::cout << "ContinuousWaveletTransformNumIntegration constructor." << std::endl;
#endif

  }

  /** Destructor. */
  template <Size D>
  ContinuousWaveletTransformNumIntegration<D>::~ContinuousWaveletTransformNumIntegration()
  {}


  template <Size D>
  double ContinuousWaveletTransformNumIntegration<D>::integrate(RawDataPointConstIterator x, RawDataPointConstIterator first, RawDataPointConstIterator last)
  {
#ifdef DEBUG_PEAK_PICKING
    std::cout << "integrate" << std::endl;
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
  double ContinuousWaveletTransformNumIntegration<D>::marr(double x)
  {
    return (1-x*x)*exp(-x*x/2);
  }


  template <Size D>
  double ContinuousWaveletTransformNumIntegration<D>::integrate(const std::vector<double>& processed_input,
      double spacing_data,
      int index)
  {
    double v = 0.;
    int half_width = wavelet_.size();
    int index_in_data = (int)floor((half_width*spacing_) / spacing_data);
    //   std::cout << "wavelet breite " << half_width*spacing_ << " spacing data " << spacing_data  << "index_in_data " << index_in_data << std::endl;

    int offset_data_left = ((index - index_in_data) < 0) ? 0 : (index-index_in_data);
    int offset_data_right = ((index + index_in_data) > (int)processed_input.size()) ? processed_input.size()-1 : (index+index_in_data);

    //   std::cout << "offset left " << offset_data_left << " " << offset_data_left*spacing_data << std::endl;
    //   std::cout << "offset right " << offset_data_right << " " << offset_data_right*spacing_data << std::endl;

    // integrate from i until offset_data_left
    for (int i = index; i > offset_data_left; --i)
    {
      int index_w_r = (int)round(((index-i)*spacing_data)/spacing_);
      int index_w_l = (int)round(((index-(i-1))*spacing_data)/spacing_);
      //  std::cout << "abstand in daten i " << (index-i)*spacing_data << " und im wavelet " << index_w_r*spacing_ << std::endl;
      //  std::cout << "abstand in daten i-1 " << (index-(i-1))*spacing_data << " und im wavelet " << (index_w_l)*spacing_ << std::endl;
      v += spacing_data / 2.*( processed_input[i]*wavelet_[index_w_r] + processed_input[i-1]*wavelet_[index_w_l] );
    }

    // integrate from i+1 until offset_data_right
    for (int i = index; i < offset_data_right; ++i)
    {
      int index_w_r = (int)round((((i+1)-index)*spacing_data)/spacing_);
      int index_w_l = (int)round(((i-index)*spacing_data)/spacing_);
      //  std::cout << "abstand in daten i " << (i-index)*spacing_data << " und im wavelet " << index_w_l*spacing_ << std::endl;
      //  std::cout << "abstand in daten i+1 " << ((i+1)-index)*spacing_data << " und im wavelet " << (index_w_r)*spacing_ << std::endl;
      v += spacing_data / 2.*( processed_input[i+1]*wavelet_[index_w_r] + processed_input[i]*wavelet_[index_w_l]);
    }

    return v / sqrt(scale_);
  }


  template <Size D>
  void ContinuousWaveletTransformNumIntegration<D>::init(double scale, double spacing, unsigned int mz_dim)
  {
    ContinuousWaveletTransform<D>::init(scale, spacing, mz_dim);

    /** Build a Marr-Wavelet for the current spacing and scale.
     *  We store the wavelet in the vector<double> wavelet_;
     *
     *  We only need a finite amount of points since the Marr function
     *  decays fast. We take 5*scale, since at that point the wavelet
     *  has dropped to ~ -10^-4
     */
    int number_of_points_right = (int)(ceil(5*scale_/spacing_));
    int number_of_points = number_of_points_right + 1;
    wavelet_.resize(number_of_points);
    wavelet_[0] = 1.;

    for (int i=1; i<number_of_points; i++)
    {
      wavelet_[i] = marr(i*spacing_/scale_ );
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
