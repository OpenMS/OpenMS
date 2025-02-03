// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Eva Lange $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/INTERFACES/DataStructures.h>
#include <OpenMS/MATH/MathFunctions.h>

#include <cmath>
#include <vector>

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

    @ingroup SignalProcessing
  */

// #define DEBUG_FILTERING

  class OPENMS_DLLAPI GaussFilterAlgorithm 
  {
public:
    /// Constructor
    GaussFilterAlgorithm();

    /// Destructor
    virtual ~GaussFilterAlgorithm();

    /**
      @brief Smoothes an Spectrum containing profile data.
    */
    bool filter(OpenMS::Interfaces::SpectrumPtr spectrum)
    {
      // create new arrays for mz / intensity data and set their size
      OpenMS::Interfaces::BinaryDataArrayPtr intensity_array(new OpenMS::Interfaces::BinaryDataArray);
      OpenMS::Interfaces::BinaryDataArrayPtr mz_array(new OpenMS::Interfaces::BinaryDataArray);
      mz_array->data.resize(spectrum->getMZArray()->data.size());
      intensity_array->data.resize(spectrum->getMZArray()->data.size());

      // apply the filter
      bool ret_val = filter(
          spectrum->getMZArray()->data.begin(), 
          spectrum->getMZArray()->data.end(),
          spectrum->getIntensityArray()->data.begin(),
          mz_array->data.begin(), intensity_array->data.begin()
          );
      // set the data of the spectrum to the new mz / int arrays
      spectrum->setMZArray(mz_array);
      spectrum->setIntensityArray(intensity_array);
      return ret_val;
    }

    /**
      @brief Smoothes an Chromatogram containing profile data.
    */
    bool filter(OpenMS::Interfaces::ChromatogramPtr chromatogram)
    {
      // create new arrays for rt / intensity data and set their size
      OpenMS::Interfaces::BinaryDataArrayPtr intensity_array(new OpenMS::Interfaces::BinaryDataArray);
      OpenMS::Interfaces::BinaryDataArrayPtr rt_array(new OpenMS::Interfaces::BinaryDataArray);
      rt_array->data.resize(chromatogram->getTimeArray()->data.size());
      intensity_array->data.resize(chromatogram->getTimeArray()->data.size());

      // apply the filter
      bool ret_val = filter(
          chromatogram->getTimeArray()->data.begin(), 
          chromatogram->getTimeArray()->data.end(),
          chromatogram->getIntensityArray()->data.begin(),
          rt_array->data.begin(), intensity_array->data.begin()
          );
      // set the data of the chromatogram to the new rt / int arrays
      chromatogram->setTimeArray(rt_array);
      chromatogram->setIntensityArray(intensity_array);
      return ret_val;
    }

    /**
      @brief Smoothes two data arrays.

      Convolutes the filter and the profile data and writes the results into the output iterators mz_out and int_out. 
    */
    template <typename ConstIterT, typename IterT>
    bool filter(
        ConstIterT mz_in_start,
        ConstIterT mz_in_end,
        ConstIterT int_in_start,
        IterT mz_out,
        IterT int_out)
    {
      bool found_signal = false;

      ConstIterT mz_it = mz_in_start;
      ConstIterT int_it = int_in_start;
      for (; mz_it != mz_in_end; mz_it++, int_it++)
      {
        // if ppm tolerance is used, calculate a reasonable width value for this m/z
        if (use_ppm_tolerance_)
        {
          initialize(Math::ppmToMass(ppm_tolerance_, *mz_it), spacing_, ppm_tolerance_, use_ppm_tolerance_);
        }

        double new_int = integrate_(mz_it, int_it, mz_in_start, mz_in_end);
        
        // store new intensity and m/z into output iterator
        *mz_out = *mz_it;
        *int_out = new_int;
        ++mz_out;
        ++int_out;

        if (fabs(new_int) > 0) found_signal = true;
      }
      return found_signal;
    }

    void initialize(double gaussian_width, double spacing, double ppm_tolerance, bool use_ppm_tolerance);

protected:

    ///Coefficients
    std::vector<double> coeffs_;
    /// The standard derivation  \f$ \sigma \f$.
    double sigma_;
    /// The spacing of the pre-tabulated kernel coefficients
    double spacing_;

    // tolerance in ppm
    bool use_ppm_tolerance_;
    double ppm_tolerance_;

    /// Computes the convolution of the raw data at position x and the gaussian kernel
    template <typename InputPeakIterator>
    double integrate_(InputPeakIterator x /* mz */, InputPeakIterator y /* int */, InputPeakIterator first, InputPeakIterator last)
    {
      double v = 0.;
      // norm the gaussian kernel area to one
      double norm = 0.;
      Size middle = coeffs_.size();

      double start_pos = (( (*x) - (middle * spacing_)) > (*first)) ? ((*x) - (middle * spacing_)) : (*first);
      double end_pos = (( (*x) + (middle * spacing_)) < (*(last - 1))) ? ((*x) + (middle * spacing_)) : (*(last - 1));

      InputPeakIterator help_x = x;
      InputPeakIterator help_y = y;
#ifdef DEBUG_FILTERING

      std::cout << "integrate from middle to start_pos " << *help_x << " until " << start_pos << std::endl;
#endif

      //integrate from middle to start_pos
      while ((help_x != first) && (*(help_x - 1) > start_pos))
      {
        // search for the corresponding datapoint of help in the gaussian (take the left most adjacent point)
        double distance_in_gaussian = fabs(*x - *help_x);
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
        double d = fabs((left_position * spacing_) - distance_in_gaussian) / spacing_;
        // check if the right data point in the gaussian exists
        double coeffs_right = (right_position < middle) ? (1 - d) * coeffs_[left_position] + d * coeffs_[right_position]
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
        double coeffs_left = (right_position < middle) ? (1 - d) * coeffs_[left_position] + d * coeffs_[right_position]
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
        double distance_in_gaussian = fabs((*x) - (*help_x));
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
        double d = fabs((left_position * spacing_) - distance_in_gaussian) / spacing_;
        double coeffs_left = (right_position < middle) ? (1 - d) * coeffs_[left_position] + d * coeffs_[right_position]
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
        double coeffs_right = (right_position < middle) ? (1 - d) * coeffs_[left_position] + d * coeffs_[right_position]
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
