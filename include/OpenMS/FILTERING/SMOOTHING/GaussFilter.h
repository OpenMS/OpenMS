// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Alexandra Zerck $
// $Authors: Eva Lange $
// --------------------------------------------------------------------------

#ifndef OPENMS_FILTERING_SMOOTHING_GAUSSFILTER_H
#define OPENMS_FILTERING_SMOOTHING_GAUSSFILTER_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <cmath>

namespace OpenMS
{
  /**
    @brief This class represents a Gaussian lowpass-filter which works on uniform as well as on non-uniform profile data.

    Gaussian filters are important in many signal processing,
    image processing, and communication applications. These filters are characterized by narrow bandwidths,
    sharp cutoffs, and low passband ripple. A key feature of Gaussian filters is that the Fourier transform of a
    Gaussian is also a Gaussian, so the filter has the same response shape in both the time and frequency domains.
    The coefficients \f$ \emph{coeffs} \f$ of the Gaussian-window with length \f$ \emph{frameSize} \f$ are calculated
    from the gaussian distribution
    \f[ \emph{coeff}(x) = \frac{1}{\sigma \sqrt{2\pi}} e^{\frac{-x^2}{2\sigma^2}} \f]
    where \f$ x=[-\frac{frameSize}{2},...,\frac{frameSize}{2}] \f$ represents the window area and \f$ \sigma \f$
    is the standard derivation.

    @note The wider the kernel width the smoother the signal (the more detail information get lost!).
          Use a gaussian filter kernel which has approximately the same width as your mass peaks,
          whereas the gaussian peak width corresponds approximately to 8*sigma.

        @note The data must be sorted according to ascending m/z!

        @htmlinclude OpenMS_GaussFilter.parameters

    @ingroup SignalProcessing
  */
//#define DEBUG_FILTERING

  class OPENMS_DLLAPI GaussFilter :
    public ProgressLogger,
    public DefaultParamHandler
  {
public:
    /// Constructor
    GaussFilter();

    /// Destructor
    virtual ~GaussFilter();

    /**
      @brief Smoothes an MSSpectrum containing profile data.

      Convolutes the filter and the profile data and writes the result back to the spectrum.

        @exception Exception::IllegalArgument is thrown, if the @em gaussian_width parameter is too small.
    */
    template <typename PeakType>
    void filter(MSSpectrum<PeakType> & spectrum)
    {
      typedef std::vector<double> ContainerT;

      // make sure the right data type is set
      spectrum.setType(SpectrumSettings::RAWDATA);
      bool found_signal = false;
      Size data_size = spectrum.size();
      ContainerT mz_in(data_size), int_in(data_size), mz_out(data_size), int_out(data_size);

      // copy spectrum to container
      for (Size p = 0; p < spectrum.size(); ++p)
      {
        mz_in[p] = spectrum[p].getMZ();
        int_in[p] = spectrum[p].getIntensity();
      }

      // apply filter
      ContainerT::iterator mz_out_it = mz_out.begin();
      ContainerT::iterator int_out_it = int_out.begin();
      found_signal = filter(mz_in.begin(), mz_in.end(), int_in.begin(), mz_out_it, int_out_it);

      // If all intensities are zero in the scan and the scan has a reasonable size, throw an exception.
      // This is the case if the gaussian filter is smaller than the spacing of raw data
      if (!found_signal && spectrum.size() >= 3)
      {
        String error_message = "Found no signal. The gaussian width is probably smaller than the spacing in your profile data. Try to use a bigger width.";
        if (spectrum.getRT() > 0.0)
        {
          error_message += String(" The error occured in the spectrum with retention time ") + spectrum.getRT() + ".\n";
        }
        std::cerr << error_message;
      }
      else
      {
        // copy the new data into the spectrum
        ContainerT::iterator mz_it = mz_out.begin();
        ContainerT::iterator int_it = int_out.begin();
        for (Size p = 0; mz_it != mz_out.end(); mz_it++, int_it++, p++)
        {
          spectrum[p].setIntensity(*int_it);
          spectrum[p].setMZ(*mz_it);
        }
      }
    }

    /**
      @brief Smoothes an MSSpectrum containing profile data.

      Convolutes the filter and the profile data and writes the results into the output iterators mz_out and int_out. 
    */
    template <typename IterT>
    bool filter(
        IterT mz_in_start,
        IterT mz_in_end,
        IterT int_in_start,
        IterT& mz_out,
        IterT& int_out)
    {
      bool use_ppm_tolerance(param_.getValue("use_ppm_tolerance").toBool());
      DoubleReal ppm_tolerance((DoubleReal)param_.getValue("ppm_tolerance"));
      bool found_signal = false;

      IterT mz_it = mz_in_start;
      IterT int_it = int_in_start;
      for (; mz_it != mz_in_end; mz_it++, int_it++)
      {
        // if ppm tolerance is used, calculate a reasonable width value for this m/z
        if (use_ppm_tolerance)
        {
          param_.setValue("gaussian_width", (*mz_it) * ppm_tolerance * 10e-6);
          updateMembers_();
        }

        DoubleReal new_int = integrate_(mz_it, int_it, mz_in_start, mz_in_end);
        
        // store new intensity and m/z into output iterator
        *mz_out = *mz_it;
        *int_out = new_int;
        ++mz_out;
        ++int_out;

        if (fabs(new_int) > 0) found_signal = true;
      }
      return found_signal;
    }

    template <typename PeakType>
    void filter(MSChromatogram<PeakType> & chromatogram)
    {

      if (param_.getValue("use_ppm_tolerance").toBool())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, 
          "GaussFilter: Cannot use ppm tolerance on chromatograms");
      }

      MSSpectrum<PeakType> filter_spectra;
      for (typename MSChromatogram<PeakType>::const_iterator it = chromatogram.begin(); it != chromatogram.end(); ++it)
      {
        filter_spectra.push_back(*it);
      }
      filter(filter_spectra);
      chromatogram.clear(false);
      for (typename MSSpectrum<PeakType>::const_iterator it = filter_spectra.begin(); it != filter_spectra.end(); ++it)
      {
        chromatogram.push_back(*it);
      }

    }

    /**
      @brief Smoothes an MSExperiment containing profile data.

        @exception Exception::IllegalArgument is thrown, if the @em gaussian_width parameter is too small.
          */
    template <typename PeakType>
    void filterExperiment(MSExperiment<PeakType> & map)
    {
      Size progress = 0;
      startProgress(0, map.size() + map.getChromatograms().size(), "smoothing data");
      for (Size i = 0; i < map.size(); ++i)
      {
        filter(map[i]);
        setProgress(++progress);
      }
      for (Size i = 0; i < map.getChromatograms().size(); ++i)
      {
        filter(map.getChromatogram(i));
        setProgress(++progress);
      }
      endProgress();
    }

protected:

    ///Coefficients
    std::vector<DoubleReal> coeffs_;
    /// The standard derivation  \f$ \sigma \f$.
    DoubleReal sigma_;
    /// The spacing of the pre-tabulated kernel coefficients
    DoubleReal spacing_;

    // Docu in base class
    virtual void updateMembers_();

    /// Computes the convolution of the raw data at position x and the gaussian kernel
    template <typename InputPeakIterator>
    DoubleReal integrate_(InputPeakIterator x /* mz */, InputPeakIterator y /* int */, InputPeakIterator first, InputPeakIterator last)
    {
      DoubleReal v = 0.;
      // norm the gaussian kernel area to one
      DoubleReal norm = 0.;
      Size middle = coeffs_.size();

      DoubleReal start_pos = (( (*x) - (middle * spacing_)) > (*first)) ? ((*x) - (middle * spacing_)) : (*first);
      DoubleReal end_pos = (( (*x) + (middle * spacing_)) < (*(last - 1))) ? ((*x) + (middle * spacing_)) : (*(last - 1));

      InputPeakIterator help_x = x;
      InputPeakIterator help_y = y;
#ifdef DEBUG_FILTERING

      std::cout << "integrate from middle to start_pos " << *help_x << " until " << start_pos << std::endl;
#endif

      //integrate from middle to start_pos
      while ((help_x != first) && (*(help_x - 1) > start_pos))
      {
        // search for the corresponding datapoint of help in the gaussian (take the left most adjacent point)
        DoubleReal distance_in_gaussian = fabs(*x - *help_x);
        Size left_position = (Size)floor(distance_in_gaussian / spacing_);

        // search for the true left adjacent data point (because of rounding errors)
        for (int j = 0; ((j < 3) &&  (distance(first, help_x - j) >= 0)); ++j)
        {
          if (((left_position - j) * spacing_ <= distance_in_gaussian) && ((left_position - j + 1) * spacing_ >= distance_in_gaussian))
          {
            left_position -= j;
            break;
          }

          if (((left_position + j) * spacing_ < distance_in_gaussian) && ((left_position + j + 1) * spacing_ < distance_in_gaussian))
          {
            left_position += j;
            break;
          }
        }

        // interpolate between the left and right data points in the gaussian to get the true value at position distance_in_gaussian
        Size right_position = left_position + 1;
        DoubleReal d = fabs((left_position * spacing_) - distance_in_gaussian) / spacing_;
        // check if the right data point in the gaussian exists
        DoubleReal coeffs_right = (right_position < middle) ? (1 - d) * coeffs_[left_position] + d * coeffs_[right_position]
                                  : coeffs_[left_position];
#ifdef DEBUG_FILTERING

        std::cout << "distance_in_gaussian " << distance_in_gaussian << std::endl;
        std::cout << " right_position " << right_position << std::endl;
        std::cout << " left_position " << left_position << std::endl;
        std::cout << "coeffs_ at left_position "  <<  coeffs_[left_position] << std::endl;
        std::cout << "coeffs_ at right_position "  <<  coeffs_[right_position] << std::endl;
        std::cout << "interpolated value left " << coeffs_right << std::endl;
#endif


        // search for the corresponding datapoint for (help-1) in the gaussian (take the left most adjacent point)
        distance_in_gaussian = fabs((*x) - (*(help_x - 1)));
        left_position = (Size)floor(distance_in_gaussian / spacing_);

        // search for the true left adjacent data point (because of rounding errors)
        for (UInt j = 0; ((j < 3) && (distance(first, help_x - j) >= 0)); ++j)
        {
          if (((left_position - j) * spacing_ <= distance_in_gaussian) && ((left_position - j + 1) * spacing_ >= distance_in_gaussian))
          {
            left_position -= j;
            break;
          }

          if (((left_position + j) * spacing_ < distance_in_gaussian) && ((left_position + j + 1) * spacing_ < distance_in_gaussian))
          {
            left_position += j;
            break;
          }
        }

        // start the interpolation for the true value in the gaussian
        right_position = left_position + 1;
        d = fabs((left_position * spacing_) - distance_in_gaussian) / spacing_;
        DoubleReal coeffs_left = (right_position < middle) ? (1 - d) * coeffs_[left_position] + d * coeffs_[right_position]
                                 : coeffs_[left_position];
#ifdef DEBUG_FILTERING

        std::cout << " help_x-1 " << *(help_x - 1) << " distance_in_gaussian " << distance_in_gaussian << std::endl;
        std::cout << " right_position " << right_position << std::endl;
        std::cout << " left_position " << left_position << std::endl;
        std::cout << "coeffs_ at left_position " <<  coeffs_[left_position] << std::endl;
        std::cout << "coeffs_ at right_position " <<   coeffs_[right_position] << std::endl;
        std::cout << "interpolated value right " << coeffs_left << std::endl;

        std::cout << " intensity " << fabs(*(help_x - 1) - (*help_x)) / 2. << " * " << *(help_y - 1) << " * " << coeffs_left << " + " << *help_y << "* " << coeffs_right
                  << std::endl;
#endif


        norm += fabs((*(help_x - 1)) - (*help_x)) / 2. * (coeffs_left + coeffs_right);

        v += fabs((*(help_x - 1)) - (*help_x)) / 2. * (*(help_y - 1) * coeffs_left + (*help_y) * coeffs_right);
        --help_x;
        --help_y;
      }


      //integrate from middle to end_pos
      help_x = x;
      help_y = y;
#ifdef DEBUG_FILTERING

      std::cout << "integrate from middle to endpos " << *help_x << " until " << end_pos << std::endl;
#endif

      while ((help_x != (last - 1)) && (*(help_x + 1) < end_pos))
      {
        // search for the corresponding datapoint for help in the gaussian (take the left most adjacent point)
        DoubleReal distance_in_gaussian = fabs((*x) - (*help_x));
        int left_position = (UInt)floor(distance_in_gaussian / spacing_);

        // search for the true left adjacent data point (because of rounding errors)
        for (int j = 0; ((j < 3) && (distance(help_x + j, last - 1) >= 0)); ++j)
        {
          if (((left_position - j) * spacing_ <= distance_in_gaussian) && ((left_position - j + 1) * spacing_ >= distance_in_gaussian))
          {
            left_position -= j;
            break;
          }

          if (((left_position + j) * spacing_ < distance_in_gaussian) && ((left_position + j + 1) * spacing_ < distance_in_gaussian))
          {
            left_position += j;
            break;
          }
        }
        // start the interpolation for the true value in the gaussian
        Size right_position = left_position + 1;
        DoubleReal d = fabs((left_position * spacing_) - distance_in_gaussian) / spacing_;
        DoubleReal coeffs_left = (right_position < middle) ? (1 - d) * coeffs_[left_position] + d * coeffs_[right_position]
                                 : coeffs_[left_position];

#ifdef DEBUG_FILTERING

        std::cout << " help " << *help_x << " distance_in_gaussian " << distance_in_gaussian << std::endl;
        std::cout << " left_position " << left_position << std::endl;
        std::cout << "coeffs_ at right_position " <<  coeffs_[left_position] << std::endl;
        std::cout << "coeffs_ at left_position " <<  coeffs_[right_position] << std::endl;
        std::cout << "interpolated value left " << coeffs_left << std::endl;
#endif

        // search for the corresponding datapoint for (help+1) in the gaussian (take the left most adjacent point)
        distance_in_gaussian = fabs((*x) - (*(help_x + 1)));
        left_position = (UInt)floor(distance_in_gaussian / spacing_);

        // search for the true left adjacent data point (because of rounding errors)
        for (int j = 0; ((j < 3) && (distance(help_x + j, last - 1) >= 0)); ++j)
        {
          if (((left_position - j) * spacing_ <= distance_in_gaussian) && ((left_position - j + 1) * spacing_ >= distance_in_gaussian))
          {
            left_position -= j;
            break;
          }

          if (((left_position + j) * spacing_ < distance_in_gaussian) && ((left_position + j + 1) * spacing_ < distance_in_gaussian))
          {
            left_position += j;
            break;
          }
        }

        // start the interpolation for the true value in the gaussian
        right_position = left_position + 1;
        d = fabs((left_position * spacing_) - distance_in_gaussian) / spacing_;
        DoubleReal coeffs_right = (right_position < middle) ? (1 - d) * coeffs_[left_position] + d * coeffs_[right_position]
                                  : coeffs_[left_position];
#ifdef DEBUG_FILTERING

        std::cout << " (help + 1) " << *(help_x + 1) << " distance_in_gaussian " << distance_in_gaussian << std::endl;
        std::cout << " left_position " << left_position << std::endl;
        std::cout << "coeffs_ at right_position " <<   coeffs_[left_position] << std::endl;
        std::cout << "coeffs_ at left_position " <<  coeffs_[right_position] << std::endl;
        std::cout << "interpolated value right " << coeffs_right << std::endl;

        std::cout << " intensity " <<  fabs(*help_x - *(help_x + 1)) / 2.
                  << " * " << *help_y << " * " << coeffs_left << " + " << *(help_y + 1)
                  << "* " << coeffs_right
                  << std::endl;
#endif
        norm += fabs((*help_x) - (*(help_x + 1)) ) / 2. * (coeffs_left + coeffs_right);

        v += fabs((*help_x) - (*(help_x + 1)) ) / 2. * ((*help_y) * coeffs_left + (*(help_y + 1)) * coeffs_right);
        ++help_x;
        ++help_y;
      }

      if (v > 0)
      {
        return v / norm;
      }
      else
      {
        return 0;
      }
    }

  };

} // namespace OpenMS
#endif
