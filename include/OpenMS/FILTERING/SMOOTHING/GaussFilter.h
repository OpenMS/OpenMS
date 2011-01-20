// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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

  class OPENMS_DLLAPI GaussFilter
  	: public ProgressLogger,
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
      void filter(MSSpectrum<PeakType>& spectrum)
      {
        // make sure the right data type is set
        spectrum.setType(SpectrumSettings::RAWDATA);
        //create container for output peaks
        std::vector<DoubleReal> output(spectrum.size());

				bool use_ppm_tolerance(param_.getValue("use_ppm_tolerance").toBool());
				DoubleReal ppm_tolerance((DoubleReal)param_.getValue("ppm_tolerance"));

        bool found_signal = false;
        for (Size p=0; p<spectrum.size(); ++p)
        {
        	// if ppm tolerance is used, calculate a reasonable width value for this m/z
					if (use_ppm_tolerance)
          {
						param_.setValue("gaussian_width", spectrum[p].getMZ() * ppm_tolerance * 10e-6);
            updateMembers_();
          }

          DoubleReal new_int = integrate_(spectrum.begin()+p,spectrum.begin(),spectrum.end());
          output[p] = std::max(new_int, 0.0);
          if (fabs(new_int) > 0) found_signal = true;
        }

        // If all intensities are zero in the scan and the scan has a reasonable size, throw an exception.
        // This is the case if the gaussian filter is smaller than the spacing of raw data
        if (!found_signal && spectrum.size()>=3)
        {
        	String error_message = "Found no signal. The gaussian width is probably smaller than the spacing in your profile data. Try to use a bigger width.";
        	if (spectrum.getRT()>0.0)
        	{
        		error_message += String(" The error occured in the spectrum with retention time ") + spectrum.getRT() + ".\n";
        	}
          std::cerr << error_message;
        }
				else
				{
					// copy the new data into the spectrum
					for (Size p=0; p<spectrum.size(); ++p)
					{
        		spectrum[p].setIntensity(output[p]);
					}
				}
      }


      /**
      	@brief Smoothes an MSExperiment containing profile data.
 
	      @exception Exception::IllegalArgument is thrown, if the @em gaussian_width parameter is too small.
			*/
      template <typename PeakType>
      void filterExperiment(MSExperiment<PeakType>& map)
      {
        startProgress(0,map.size(),"smoothing data");
        for (Size i = 0; i < map.size(); ++i)
        {
          filter(map[i]);
          setProgress(i);
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
      template < typename InputPeakIterator >
      DoubleReal integrate_(InputPeakIterator x, InputPeakIterator first, InputPeakIterator last)
      {
        DoubleReal v = 0.;
        // norm the gaussian kernel area to one
        DoubleReal norm = 0.;
        Size middle = coeffs_.size();

        DoubleReal start_pos = ((x->getMZ()-(middle*spacing_)) > first->getMZ()) ? (x->getMZ()-(middle*spacing_))
                           : first->getMZ();
        DoubleReal end_pos = ((x->getMZ()+(middle*spacing_)) < (last-1)->getMZ()) ? (x->getMZ()+(middle*spacing_))
                         : (last-1)->getMZ();


        InputPeakIterator help = x;
#ifdef DEBUG_FILTERING

        std::cout << "integrate from middle to start_pos "<< help->getMZ() << " until " << start_pos << std::endl;
#endif

        //integrate from middle to start_pos
        while ((help != first) && ((help-1)->getMZ() > start_pos))
        {
          // search for the corresponding datapoint of help in the gaussian (take the left most adjacent point)
          DoubleReal distance_in_gaussian = fabs(x->getMZ() - help->getMZ());
          Size left_position = (Size)floor(distance_in_gaussian / spacing_);

          // search for the true left adjacent data point (because of rounding errors)
          for (int j=0; ((j<3) &&  (distance(first,help-j) >= 0)); ++j)
          {
            if (((left_position-j)*spacing_ <= distance_in_gaussian) && ((left_position-j+1)*spacing_ >= distance_in_gaussian))
            {
              left_position -= j;
              break;
            }

            if (((left_position+j)*spacing_ < distance_in_gaussian) && ((left_position+j+1)*spacing_ < distance_in_gaussian))
            {
              left_position +=j;
              break;
            }
          }

          // interpolate between the left and right data points in the gaussian to get the true value at position distance_in_gaussian
          Size right_position = left_position+1;
          DoubleReal d = fabs((left_position*spacing_)-distance_in_gaussian) / spacing_;
          // check if the right data point in the gaussian exists
          DoubleReal coeffs_right = (right_position < middle) ? (1-d)*coeffs_[left_position]+d*coeffs_[right_position]
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
          distance_in_gaussian = fabs(x->getMZ() - (help-1)->getMZ());
          left_position = (Size)floor(distance_in_gaussian / spacing_);

          // search for the true left adjacent data point (because of rounding errors)
          for (UInt j=0; ((j<3) && (distance(first,help-j) >= 0)); ++j)
          {
            if (((left_position-j)*spacing_ <= distance_in_gaussian) && ((left_position-j+1)*spacing_ >= distance_in_gaussian))
            {
              left_position -= j;
              break;
            }

            if (((left_position+j)*spacing_ < distance_in_gaussian) && ((left_position+j+1)*spacing_ < distance_in_gaussian))
            {
              left_position +=j;
              break;
            }
          }

          // start the interpolation for the true value in the gaussian
          right_position = left_position+1;
          d = fabs((left_position*spacing_)-distance_in_gaussian) / spacing_;
          DoubleReal coeffs_left= (right_position < middle) ? (1-d)*coeffs_[left_position]+d*coeffs_[right_position]
                              : coeffs_[left_position];
#ifdef DEBUG_FILTERING

          std::cout << " help-1 " << (help-1)->getMZ() << " distance_in_gaussian " << distance_in_gaussian << std::endl;
          std::cout << " right_position " << right_position << std::endl;
          std::cout << " left_position " << left_position << std::endl;
          std::cout << "coeffs_ at left_position " <<  coeffs_[left_position]<<std::endl;
          std::cout << "coeffs_ at right_position " <<   coeffs_[right_position]<<std::endl;
          std::cout << "interpolated value right " << coeffs_left << std::endl;

          std::cout << " intensity " << fabs((help-1)->getMZ()-help->getMZ()) / 2. << " * " << (help-1)->getIntensity() << " * " << coeffs_left <<" + " << (help)->getIntensity()<< "* " << coeffs_right
          << std::endl;
#endif


          norm += fabs((help-1)->getMZ()-help->getMZ()) / 2. * (coeffs_left + coeffs_right);

          v+= fabs((help-1)->getMZ()-help->getMZ()) / 2. * ((help-1)->getIntensity()*coeffs_left + help->getIntensity()*coeffs_right);
          --help;
        }


        //integrate from middle to end_pos
        help = x;
#ifdef DEBUG_FILTERING

        std::cout << "integrate from middle to endpos "<< (help)->getMZ() << " until " << end_pos << std::endl;
#endif

        while ((help != (last-1)) && ((help+1)->getMZ() < end_pos))
        {
          // search for the corresponding datapoint for help in the gaussian (take the left most adjacent point)
          DoubleReal distance_in_gaussian = fabs(x->getMZ() - help->getMZ());
          int left_position = (UInt)floor(distance_in_gaussian / spacing_);

          // search for the true left adjacent data point (because of rounding errors)
          for (int j=0; ((j<3) && (distance(help+j,last-1) >= 0)); ++j)
          {
            if (((left_position-j)*spacing_ <= distance_in_gaussian) && ((left_position-j+1)*spacing_ >= distance_in_gaussian))
            {
              left_position -= j;
              break;
            }

            if (((left_position+j)*spacing_ < distance_in_gaussian) && ((left_position+j+1)*spacing_ < distance_in_gaussian))
            {
              left_position +=j;
              break;
            }
          }
          // start the interpolation for the true value in the gaussian
          Size right_position = left_position+1;
          DoubleReal d = fabs((left_position*spacing_)-distance_in_gaussian) / spacing_;
          DoubleReal coeffs_left= (right_position < middle) ? (1-d)*coeffs_[left_position]+d*coeffs_[right_position]
                              : coeffs_[left_position];

#ifdef DEBUG_FILTERING

          std::cout << " help " << (help)->getMZ() << " distance_in_gaussian " << distance_in_gaussian << std::endl;
          std::cout << " left_position " << left_position << std::endl;
          std::cout << "coeffs_ at right_position " <<  coeffs_[left_position]<<std::endl;
          std::cout << "coeffs_ at left_position " <<  coeffs_[right_position]<<std::endl;
          std::cout << "interpolated value left " << coeffs_left << std::endl;
#endif

          // search for the corresponding datapoint for (help+1) in the gaussian (take the left most adjacent point)
          distance_in_gaussian = fabs(x->getMZ() - (help+1)->getMZ());
          left_position = (UInt)floor(distance_in_gaussian / spacing_);

          // search for the true left adjacent data point (because of rounding errors)
          for (int j=0; ((j<3) && (distance(help+j,last-1) >= 0)); ++j)
          {
            if (((left_position-j)*spacing_ <= distance_in_gaussian) && ((left_position-j+1)*spacing_ >= distance_in_gaussian))
            {
              left_position -= j;
              break;
            }

            if (((left_position+j)*spacing_ < distance_in_gaussian) && ((left_position+j+1)*spacing_ < distance_in_gaussian))
            {
              left_position +=j;
              break;
            }
          }

          // start the interpolation for the true value in the gaussian
          right_position = left_position+1;
          d = fabs((left_position*spacing_)-distance_in_gaussian) / spacing_;
          DoubleReal coeffs_right = (right_position < middle) ? (1-d)*coeffs_[left_position]+d*coeffs_[right_position]
                                : coeffs_[left_position];
#ifdef DEBUG_FILTERING

          std::cout << " (help + 1) " << (help+1)->getMZ() << " distance_in_gaussian " << distance_in_gaussian << std::endl;
          std::cout << " left_position " << left_position << std::endl;
          std::cout << "coeffs_ at right_position " <<   coeffs_[left_position]<<std::endl;
          std::cout << "coeffs_ at left_position " <<  coeffs_[right_position]<<std::endl;
          std::cout << "interpolated value right " << coeffs_right << std::endl;

          std::cout << " intensity " <<  fabs(help->getMZ() - (help+1)->getMZ()) / 2.
          << " * " << help->getIntensity() << " * " << coeffs_left <<" + " << (help+1)->getIntensity()
          << "* " << coeffs_right
          << std::endl;
#endif
          norm += fabs(help->getMZ() - (help+1)->getMZ()) / 2. * (coeffs_left + coeffs_right);

          v+= fabs(help->getMZ() - (help+1)->getMZ()) / 2. * (help->getIntensity()*coeffs_left + (help+1)->getIntensity()*coeffs_right);
          ++help;
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
