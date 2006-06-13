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
// $Id: ContinuousWaveletTransformFFT.h,v 1.17 2006/05/29 16:18:55 elange Exp $
// $Author: elange $
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_CONTINUOUSWAVELETTRANSFORMFFT_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_CONTINUOUSWAVELETTRANSFORMFFT_H

# include <OpenMS/TRANSFORMATIONS/RAW2PEAK/ContinuousWaveletTransform.h>


#ifdef DEBUG_PEAK_PICKING
#include<iostream>
#include<fstream>
#endif


#include <math.h>
#include <complex>

#ifdef FFTW_DEF
#include <fftw3.h>
#endif

namespace OpenMS
{
  /**
     @brief This class computes the continuous wavelet transformation using a marr wavelet.

     The convolution of the signal and the wavelet is computed in the fourier space.

     @todo write test (Eva)
  */
  template <Size D>
  class ContinuousWaveletTransformFFT : public ContinuousWaveletTransform<D>
  {
    /** @name Type definitions
     */
    //@{
    typedef std::vector<DRawDataPoint<D> > RawDataVector;
    ///
    typedef typename RawDataVector::iterator RawDataPointIterator;
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
    //@{
    /// Constructors and Destructors.

    /** Default constructor. */
    ContinuousWaveletTransformFFT();

    /** Destructor. */
    virtual ~ContinuousWaveletTransformFFT();
    //@}

     /** @name Assignement
     */
    //@{
    ///
    inline ContinuousWaveletTransformFFT& operator=(const ContinuousWaveletTransformFFT& cwt)
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

    /** Compute the continuous wavelet transformation (use a marr wavelet)
     * of the signal intervall [it_begin, it_end).
     * Using the FFTW and IFFTW you have to imagine that x=1/N*fft(fft(x)).
     * Because we are only interested in finding the maxima positions in the
     * cwt (potential peak positions) we compute the cwt without scaling,
     * so we compute N*cwt(x).
     */
    virtual
    void transform(RawDataPointIterator begin_input, 
									 RawDataPointIterator end_input, 
									 float resolution);

    virtual
    void init(double scale, double spacing, unsigned int MZ);

    /** Stores a Marr wavelet of the desired parameters in the wavelet vector **/
    void buildMarrWavelet(double scale);

  protected:
    /** The processed input data used for the transform **/
    std::vector<double> processed_input_;

  };

  /** Default constructor. */
  template <Size D>
  ContinuousWaveletTransformFFT<D>::ContinuousWaveletTransformFFT()
    : ContinuousWaveletTransform<D>()
  {
#ifdef DEBUG_PEAK_PICKING
    std::cout << "ContinuousWaveletTransformFFT constructor." << std::endl;
#endif
  }

  /** Destructor. */
  template <Size D>
  ContinuousWaveletTransformFFT<D>::~ContinuousWaveletTransformFFT()
  {
#ifdef DEBUG_PEAK_PICKING
    std::cout << "ContinuousWaveletTransformFFT destructor." << std::endl;
#endif
  }

  template <Size D>
  void ContinuousWaveletTransformFFT<D>::transform(RawDataPointIterator begin_input,
                                                   RawDataPointIterator end_input,
																									 float resolution)
  {
    // make sure there is no signal and no input data left from a prior transform
    ContinuousWaveletTransform<D>::signal_.clear();
    processed_input_.clear();

    // first compute the length of the signal
    unsigned int n = distance(begin_input, end_input);

    signal_length_ = n;

    // compute the next power of two, so that we can perform an efficient transform
    unsigned int transform_length = (unsigned int) pow(2,(floor(log2(n)+1)));

    // ??? TODO The length of the array should be smaller than the max integer
    signal_.resize(transform_length);
    processed_input_.resize(transform_length);

    int m=transform_length/2+1;
    wavelet_.resize(m);//transform_length);

    float origin  = begin_input->getPosition()[mz_dim_];
    spacing_ = ((end_input-1)->getPosition()[mz_dim_]-origin)/(n-1);

    // now process the input data for the transform
    int number_of_zeros= (int) (transform_length-n)/2;

#ifdef DEBUG_PEAK_PICKING
    std::cout << "VECTOR " ;
    std::cout << begin_input->getPosition()[mz_dim_] << " UNTIL  " ;
    std::cout << (end_input-1)->getPosition()[mz_dim_] << std::endl;
    std::cout << "Number of zeros: " << number_of_zeros << " n: "
              << n << " transform_length: " << transform_length << std::endl;
#endif

    RawDataPointIterator it_help=begin_input;

    unsigned int i;
    for (i = 0; i < (unsigned int) number_of_zeros; ++i)
      processed_input_[i]=0;

    end_left_padding_   = i;
    processed_input_[i] = it_help->getIntensity();
    ++i;

    double x;
    for (unsigned int k=1; k < n; ++k,++i)
    {
      x = origin + k*spacing_;
      while (((it_help+1) < end_input) && ((it_help+1)->getPosition()[mz_dim_] < x))
      {
        ++it_help;
      }
      processed_input_[i]=getInterpolatedValue(x,it_help);
    }
    begin_right_padding_=i-1;

    // zeropadding
    for (;i<processed_input_.size();++i)
    {
      processed_input_[i]=0;
    }

#ifdef DEBUG_PEAK_PICKING
    std::cout << "ARRAY " ;
    std::cout << "spacing: " << spacing_ << " ";
    std::cout << origin << " UNTIL  " ;
    std::cout << origin+(begin_right_padding_-end_left_padding_)*spacing_ << std::endl;
#endif

    // fftw_complex *input_ft   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m);
    // fftw_complex *product_ft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m);
    std::vector<double> input_ft(transform_length);
    std::vector< std::complex<double>  > product_ft(transform_length);
    fftw_plan p_fft,p_ifft;

    fftw_complex* target = reinterpret_cast<fftw_complex*>(&(product_ft[0]));
    p_fft = fftw_plan_dft_r2c_1d(transform_length,&(processed_input_[0]),target,FFTW_ESTIMATE);

    // coumpute the fourier transform of the split region
    fftw_execute(p_fft);

    // Multiply the fourier transform W of the Marr wavelet with the fourier transform Y of the signal (Y_k=Y'_n-k (Y' is the conjugate of Y))
    // P_0=W_0*Y_0
    // P_1=W_1*Y_1
    // ...
    // P_n/2=W_n/2*Y_n/n
    // P_n/2+1=W_n/2-1*Y'_n/2-1
    // ...
    // P_n-1=W_1*Y'_1

    // The real and imaginery part of P are stored : r_0 r_1 r_2 ... r_n/2 i_(n+1)/2-1 ... i_2 i_1

    // Y_0 is real
    product_ft[0]=wavelet_[0]*product_ft[0];
    int middle = (int) transform_length / 2;

    int j;
    for (j=1; j < middle;++j)
    {
      // real part [0..n/2]
      product_ft[j]= wavelet_[j]*product_ft[j];
      // imaginery part [n/2+1..n-1]
      product_ft[transform_length-j]= wavelet_[j]*product_ft[transform_length-j];
    }
    // P_n/2=W_n/2*Y_n/n
    product_ft[j]= wavelet_[j]*product_ft[j];


    /*
    for (int i=0; i<product_ft.size(); i++)
      std::cerr << i << " " << product_ft[i].real()  << " " << product_ft[i].imag() << std::endl;
    */
    p_ifft = fftw_plan_dft_c2r_1d(transform_length,target,&(signal_[0]),FFTW_ESTIMATE);

    // ignore the first and last data points in the cwt signal (come from zeropadding in input)
    fftw_execute(p_ifft);

    /*
    for (int i=0; i<signal_.size(); i++)
          std::cerr << i << " " << signal_[i] << std::endl;
    */

    fftw_destroy_plan(p_fft);
    fftw_destroy_plan(p_ifft);
  }

  template <Size D>
  void ContinuousWaveletTransformFFT<D>::init(double scale, double spacing, unsigned int MZ)
  {
    ContinuousWaveletTransform<D>::init(scale, spacing, MZ);

    int transform_length  = wavelet_.size();
    int middle            = (int) transform_length/2+1;
    double fourierspacing =  (2.*M_PI) / (double)transform_length;
    double prefac         = sqrt(M_PI / transform_length) * pow(scale / 2., 3);
    double scale2_2       = scale*scale / 2.;

    // For the real and even wavelet w is W(f)=W(-f). Therefore we only store the data points of W in [0..n/2]
    double omega  = 0.;
    double omega2 = 0.;

#ifdef DEBUG_PEAK_PICKING
    std::ofstream wavelet("Wavelet.dta",std::ios::out);
#endif

    // (ost) made i signed to avoid compiler warning
    for (int i=0; i < middle; ++i)
    {
      omega += fourierspacing;
      omega2 = omega*omega;

      wavelet_[i]=prefac*omega2 * exp(-omega2*scale2_2);
#ifdef DEBUG_PEAK_PICKING
      wavelet << i << " " << wavelet_[i] << std::endl;
#endif
    }
  }

}

#endif

