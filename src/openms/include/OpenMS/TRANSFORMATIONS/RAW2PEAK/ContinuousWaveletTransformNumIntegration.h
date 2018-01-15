// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Eva Lange $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_CONTINUOUSWAVELETTRANSFORMNUMINTEGRATION_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_CONTINUOUSWAVELETTRANSFORMNUMINTEGRATION_H

#include <cmath>

#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/ContinuousWaveletTransform.h>


#ifdef DEBUG_PEAK_PICKING
#include <iostream>
#include <fstream>
#endif

namespace OpenMS
{
  /**
    @brief This class computes the continuous wavelet transformation using a marr wavelet.

    The convolution of the signal and the wavelet is computed by numerical integration.
  */
  class OPENMS_DLLAPI ContinuousWaveletTransformNumIntegration :
    public ContinuousWaveletTransform
  {
public:
    /// Profile data const iterator type
    typedef ContinuousWaveletTransform::PeakConstIterator PeakConstIterator;

    using ContinuousWaveletTransform::signal_;
    using ContinuousWaveletTransform::wavelet_;
    using ContinuousWaveletTransform::scale_;
    using ContinuousWaveletTransform::spacing_;
    using ContinuousWaveletTransform::end_left_padding_;
    using ContinuousWaveletTransform::begin_right_padding_;
    using ContinuousWaveletTransform::signal_length_;


    /// Constructor
    ContinuousWaveletTransformNumIntegration() :
      ContinuousWaveletTransform()
    {}

    /// Destructor.
    ~ContinuousWaveletTransformNumIntegration() override {}

    /**
        @brief Computes the wavelet transform of a given profile data interval [begin_input,end_input)

        - Resolution = 1: the wavelet transform will be computed at every position of the Profile data,
        - Resolution = 2: the wavelet transform will be computed at 2x(number of Profile data positions) positions
            (the profile data are interpolated to get the intensity for missing positions)
        .

        @note The InputPeakIterator should point to a Peak1D or a class derived from Peak1D.

        @note Before starting the transformation you have to call the init function

    */
    template <typename InputPeakIterator>
    void transform(InputPeakIterator begin_input,
                   InputPeakIterator end_input,
                   float resolution)
    {
#ifdef DEBUG_PEAK_PICKING
      std::cout << "ContinuousWaveletTransformNumIntegration::transform: start " << begin_input->getMZ() << " until " << (end_input - 1)->getMZ() << std::endl;
#endif
      if (fabs(resolution - 1) < 0.0001)
      {
        // resolution = 1 corresponds to the cwt at supporting points which have a distance corresponding to the minimal spacing in [begin_input,end_input)
        SignedSize n = distance(begin_input, end_input);
        signal_length_ = n;

        signal_.clear();
        signal_.resize(n);

        // TODO avoid to compute the cwt for the zeros in signal
#ifdef DEBUG_PEAK_PICKING
        std::cout << "---------START TRANSFORM---------- \n";
#endif
        InputPeakIterator help = begin_input;
        for (int i = 0; i < n; ++i)
        {
          signal_[i].setMZ(help->getMZ());
          signal_[i].setIntensity((Peak1D::IntensityType)integrate_(help, begin_input, end_input));
          ++help;
        }
#ifdef DEBUG_PEAK_PICKING
        std::cout << "---------END TRANSFORM----------" << std::endl;
#endif
        // no zeropadding
        begin_right_padding_ = n;
        end_left_padding_ = -1;
      }
      else
      {
        SignedSize n = SignedSize(resolution * distance(begin_input, end_input));
        double origin  = begin_input->getMZ();
        double spacing = ((end_input - 1)->getMZ() - origin) / (n - 1);

        std::vector<double> processed_input(n);
        signal_.clear();
        signal_.resize(n);

        InputPeakIterator it_help = begin_input;
        processed_input[0] = it_help->getIntensity();

        double x;
        for (SignedSize k = 1; k < n; ++k)
        {
          x = origin + k * spacing;
          // go to the real data point next to x
          while (((it_help + 1) < end_input) && ((it_help + 1)->getMZ() < x))
          {
            ++it_help;
          }
          processed_input[k] = getInterpolatedValue_(x, it_help);
        }
        
        // TODO avoid to compute the cwt for the zeros in signal
        for (Int i = 0; i < n; ++i)
        {
          signal_[i].setMZ(origin + i * spacing);
          signal_[i].setIntensity((Peak1D::IntensityType)integrate_(processed_input, spacing, i));
        }

        begin_right_padding_ = n;
        end_left_padding_ = -1;
      }
    }

    /**
        @brief Perform necessary preprocessing steps like tabulating the Wavelet.

        Build a Marr-Wavelet for the current spacing and scale.
        We store the wavelet in the vector<double> wavelet_;

        We only need a finite amount of points since the Marr function
        decays fast. We take 5*scale, since at that point the wavelet
        has dropped to ~ -10^-4
    */
    void init(double scale, double spacing) override;

protected:

    /// Computes the convolution of the wavelet and the profile data at position x with resolution = 1
    template <typename InputPeakIterator>
    double integrate_(InputPeakIterator x, InputPeakIterator first, InputPeakIterator last)
    {
#ifdef DEBUG_PEAK_PICKING
      std::cout << "integrate_" << std::endl;
#endif

      double v = 0.;
      double middle_spacing = wavelet_.size() * spacing_;

      double start_pos = ((x->getMZ() - middle_spacing) > first->getMZ()) ? (x->getMZ() - middle_spacing)
                         : first->getMZ();
      double end_pos = ((x->getMZ() + middle_spacing) < (last - 1)->getMZ()) ? (x->getMZ() + middle_spacing)
                       : (last - 1)->getMZ();

      InputPeakIterator help = x;

#ifdef DEBUG_PEAK_PICKING
      std::cout << "integrate from middle to start_pos " << help->getMZ() << " until " << start_pos << std::endl;
#endif

      //integrate from middle to start_pos
      while ((help != first) && ((help - 1)->getMZ() > start_pos))
      {
        // search for the corresponding data point of help in the wavelet (take the left most adjacent point)
        double distance = fabs(x->getMZ() - help->getMZ());
        Size index_w_r = (Size) Math::round(distance / spacing_);
        if (index_w_r >= wavelet_.size())
        {
          index_w_r = wavelet_.size() - 1;
        }
        double wavelet_right =  wavelet_[index_w_r];

#ifdef DEBUG_PEAK_PICKING
        std::cout << "distance x help " << distance << std::endl;
        std::cout << "distance in wavelet_ " << index_w_r * spacing_ << std::endl;
        std::cout << "wavelet_right "  <<  wavelet_right << std::endl;
#endif

        // search for the corresponding datapoint for (help-1) in the wavelet (take the left most adjacent point)
        distance = fabs(x->getMZ() - (help - 1)->getMZ());
        Size index_w_l = (Size) Math::round(distance / spacing_);
        if (index_w_l >= wavelet_.size())
        {
          index_w_l = wavelet_.size() - 1;
        }
        double wavelet_left =  wavelet_[index_w_l];

        // start the interpolation for the true value in the wavelet

#ifdef DEBUG_PEAK_PICKING
        std::cout << " help-1 " << (help - 1)->getMZ() << " distance x, help-1" << distance << std::endl;
        std::cout << "distance in wavelet_ " << index_w_l * spacing_ << std::endl;
        std::cout << "wavelet_ at left " <<   wavelet_left << std::endl;

        std::cout << " intensity " << fabs((help - 1)->getMZ() - help->getMZ()) / 2. << " * " << (help - 1)->getIntensity() << " * " << wavelet_left << " + " << (help)->getIntensity() << "* " << wavelet_right
                  << std::endl;
#endif

        v += fabs((help - 1)->getMZ() - help->getMZ()) / 2. * ((help - 1)->getIntensity() * wavelet_left + help->getIntensity() * wavelet_right);
        --help;
      }


      //integrate from middle to end_pos
      help = x;
#ifdef DEBUG_PEAK_PICKING
      std::cout << "integrate from middle to endpos " << (help)->getMZ() << " until " << end_pos << std::endl;
#endif
      while ((help != (last - 1)) && ((help + 1)->getMZ() < end_pos))
      {
        // search for the corresponding datapoint for help in the wavelet (take the left most adjacent point)
        double distance = fabs(x->getMZ() - help->getMZ());
        Size index_w_l = (Size) Math::round(distance / spacing_);
        if (index_w_l >= wavelet_.size())
        {
          index_w_l = wavelet_.size() - 1;
        }
        double wavelet_left =  wavelet_[index_w_l];

#ifdef DEBUG_PEAK_PICKING
        std::cout << " help " << (help)->getMZ() << " distance x, help" << distance << std::endl;
        std::cout << "distance in wavelet_ " << index_w_l * spacing_ << std::endl;
        std::cout << "wavelet_ at left " <<   wavelet_left << std::endl;
#endif

        // search for the corresponding datapoint for (help+1) in the wavelet (take the left most adjacent point)
        distance = fabs(x->getMZ() - (help + 1)->getMZ());
        Size index_w_r = (Size) Math::round(distance / spacing_);
        if (index_w_r >= wavelet_.size())
        {
          index_w_r = wavelet_.size() - 1;
        }
        double wavelet_right =  wavelet_[index_w_r];

#ifdef DEBUG_PEAK_PICKING
        std::cout << " help+1 " << (help + 1)->getMZ() << " distance x, help+1" << distance << std::endl;
        std::cout << "distance in wavelet_ " << index_w_r * spacing_ << std::endl;
        std::cout << "wavelet_ at right " <<   wavelet_right << std::endl;
#endif

        v += fabs(help->getMZ() - (help + 1)->getMZ()) / 2. * (help->getIntensity() * wavelet_left + (help + 1)->getIntensity() * wavelet_right);
        ++help;
      }


#ifdef DEBUG_PEAK_PICKING
      std::cout << "return" << (v / sqrt(scale_)) << std::endl;
#endif
      return v / sqrt(scale_);
    }

    /// Computes the convolution of the wavelet and the profile data at position x with resolution > 1
    double integrate_(const std::vector<double> & processed_input, double spacing_data, int index);

    /// Computes the Marr wavelet at position x
    inline double marr_(const double x) const
    {
      return (1 - x * x) * exp(-x * x / 2);
    }

  };
} //namespace OpenMS
#endif
