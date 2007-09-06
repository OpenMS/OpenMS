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
// $Maintainer: Eva Lange  $
// --------------------------------------------------------------------------

#ifndef OPENMS_FILTERING_SMOOTHING_GAUSSFILTER_H
#define OPENMS_FILTERING_SMOOTHING_GAUSSFILTER_H

#include <OpenMS/FILTERING/SMOOTHING/SmoothFilter.h>
#include <OpenMS/KERNEL/MSExperimentExtern.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>


#include <math.h>

namespace OpenMS
{
  /**
    @brief This class represents a Gaussian lowpass-filter which works on uniform as well as on non-uniform raw data.
   
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
		 
		@ref GaussFilter_Parameters are explained on a separate page.
          
    @ingroup Filtering
  */
//#define DEBUG_FILTERING

  class GaussFilter : public SmoothFilter, public DefaultParamHandler
  {
    public:
      using SmoothFilter::coeffs_;

      /// Constructor
      inline GaussFilter()
      : SmoothFilter(),
        DefaultParamHandler("GaussFilter")
      {
      	//Parameter settings
      	defaults_.setValue("gaussian_width",0.8,"Use a gaussian filter kernel which has approximately the same width as your mass peaks");
        
         //members
        sigma_ = .1;
        spacing_ = 0.01;

        //compute the filter kernel coefficients
        init(sigma_,spacing_);
        
        defaultsToParam_();
      }

      /// Destructor
      virtual ~GaussFilter()
      {
      }

      /// Non-mutable access to the sigma
      inline DoubleReal getSigma() const
      {
        return sigma_;
      }
      /// Mutable access to the sigma
      inline void setSigma(DoubleReal sigma)
      {
        sigma_ = sigma;
        spacing_ = 4*sigma_ / 50;
        init(sigma_,spacing_);
        
        param_.setValue("gaussian_width",8*sigma_);
      }
      
      /// Non-mutable access to the kernel width
      inline DoubleReal getKernelWidth() const
      {
        return (sigma_ * 8.);
      }
      /// Mutable access to the kernel width
      inline void setKernelWidth(DoubleReal kernel_width) throw (Exception::InvalidValue)
      {
        if (kernel_width <= 0)
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"The kernel width should be greater than zero!",String(kernel_width));
        }
        
        sigma_ = kernel_width / 8.;
        init(sigma_,spacing_);
        param_.setValue("gaussian_width",kernel_width);
      }
      
      /// Non-mutable access to the spacing
      inline DoubleReal getSpacing() const
      {
        return spacing_;
      }
      /// Mutable access to the spacing
      inline void setSpacing(DoubleReal spacing)
      {
        spacing_=spacing;
        OPENMS_PRECONDITION((4*sigma_ > spacing), "You have to choose a smaller spacing for the kernel coefficients!" );
        init(sigma_,spacing_);
      }
      
      /** @brief Build a gaussian distribution for the current spacing and standard deviation.
          
          We store the coefficiens of gaussian in the vector<DoubleReal> coeffs_;

          We only need a finite amount of points since the gaussian distribution
          decays fast. We take 4*sigma (99.993666% of the area is within four standard deviations), since at that point the function
          has dropped to ~ -10^-4
      */
      void init(DoubleReal sigma, DoubleReal spacing);


      /** @brief Applies the convolution with the filter coefficients to an given iterator range.

      Convolutes the filter and the raw data in the iterator intervall [first,last) and writes the
      resulting data to the smoothed_data_container.

      @note This method assumes that the InputPeakIterator (e.g. of type MSSpectrum<RawDataPoint1D >::const_iterator)
            points to a data point of type RawDataPoint1D or any other class derived from RawDataPoint1D.

            The resulting peaks in the smoothed_data_container (e.g. of type MSSpectrum<RawDataPoint1D >)
            can be of type RawDataPoint1D or any other class derived from DRawDataPoint. 
       
            If you use MSSpectrum iterators you have to set the SpectrumSettings by your own.
       */
      template <typename InputPeakIterator, typename OutputPeakContainer  >
      void filter(InputPeakIterator first, InputPeakIterator last, OutputPeakContainer& smoothed_data_container)
      {
        smoothed_data_container.resize(distance(first,last));

#ifdef DEBUG_FILTERING
        std::cout << "KernelWidth: " << 8*sigma_ << std::endl;
        std::cout << "Spacing: " << spacing_ << std::endl;
#endif

        InputPeakIterator help = first;
        typename OutputPeakContainer::iterator out_it = smoothed_data_container.begin();
        while (help != last)
        {
          out_it->setPosition(help->getMZ());
          out_it->setIntensity(integrate_(help,first,last));
          ++out_it;
          ++help;
        }
      }


      /** @brief Convolutes the filter coefficients and the input raw data.

         Convolutes the filter and the raw data in the input_peak_container and writes the
          resulting data to the smoothed_data_container.

      @note This method assumes that the elements of the InputPeakContainer (e.g. of type MSSpectrum<RawDataPoint1D >)
            are of type RawDataPoint1D or any other class derived from RawDataPoint1D.

            The resulting peaks in the smoothed_data_container (e.g. of type MSSpectrum<RawDataPoint1D >)
            can be of type RawDataPoint1D or any other class derived from DRawDataPoint. 
       
            If you use MSSpectrum iterators you have to set the SpectrumSettings by your own.
         */
      template <typename InputPeakContainer, typename OutputPeakContainer >
      void filter(const InputPeakContainer& input_peak_container, OutputPeakContainer& smoothed_data_container)
      {
        filter(input_peak_container.begin(), input_peak_container.end(), smoothed_data_container);
      }


       /** @brief Filters every MSSpectrum in a given iterator range.
      		
      	Filters the data successive in every scan in the intervall [first,last).
      	The filtered data are stored in a MSExperiment.
      					
      	@note The InputSpectrumIterator should point to a MSSpectrum. Elements of the input spectren should be of type RawDataPoint1D 
                or any other derived class of DRawDataPoint.

          @note You have to copy the ExperimentalSettings of the raw data by your own. 	
      */
      template <typename InputSpectrumIterator, typename OutputPeakType >
      void filterExperiment(InputSpectrumIterator first,
                           	InputSpectrumIterator last,
                           	MSExperiment<OutputPeakType>& ms_exp_filtered)
      {
        UInt n = distance(first,last);
        startProgress(0,n,"smoothing data");
        
        // pick peaks on each scan
        for (UInt i = 0; i < n; ++i)
        {
          MSSpectrum< OutputPeakType > spectrum;
          InputSpectrumIterator input_it = first+i;
          
          // filter scan i
          filter(*input_it,spectrum);
          setProgress(i);
          
          // copy the spectrum settings
          static_cast<SpectrumSettings&>(spectrum) = *input_it;
          spectrum.setType(SpectrumSettings::RAWDATA);

          // copy the spectrum information
          spectrum.getPrecursorPeak() = input_it->getPrecursorPeak();
          spectrum.setRT(input_it->getRT());
          spectrum.setMSLevel(input_it->getMSLevel());
          spectrum.getName() = input_it->getName();

          ms_exp_filtered.push_back(spectrum);
        }
        endProgress();
      }
	  
	     /** @brief Filters every MSSpectrum in a given iterator range.
      		
      	Filters the data successive in every scan in the intervall [first,last).
      	The filtered data are stored in a MSExperiment.
      					
      	@note The InputSpectrumIterator should point to a MSSpectrum. Elements of the input spectren should be of type RawDataPoint1D 
                or any other derived class of DRawDataPoint.

          @note You have to copy the ExperimentalSettings of the raw data by your own. 	
      */
      template <typename InputSpectrumIterator, typename OutputPeakType >
      void filterExperiment(InputSpectrumIterator first,
                            		InputSpectrumIterator last,
                            		MSExperimentExtern<OutputPeakType>& ms_exp_filtered)
      {
        UInt n = distance(first,last);
        startProgress(0,n,"smoothing data");
        // pick peaks on each scan
        for (UInt i = 0; i < n; ++i)
        {
          MSSpectrum< OutputPeakType > spectrum;
          InputSpectrumIterator input_it = first+i;

          // pick the peaks in scan i
          filter(*input_it,spectrum);

          // if any peaks are found copy the spectrum settings
          if (spectrum.size() > 0)
          {
            // copy the spectrum settings
            spectrum.setType(SpectrumSettings::RAWDATA);

            // copy the spectrum information
            spectrum.getPrecursorPeak() = input_it->getPrecursorPeak();
            spectrum.setRT(input_it->getRT());
            spectrum.setMSLevel(input_it->getMSLevel());
            spectrum.getName() = input_it->getName();

            ms_exp_filtered.push_back(spectrum);
          }
        }
        endProgress();
      }



      /** @brief Filters a MSExperiment.
      	
      Filters the data every scan in the MSExperiment.
      The filtered data are stored in a MSExperiment.
      				
      @note The InputPeakType as well as the OutputPeakType should be of type RawDataPoint1D 
               or any other derived class of DRawDataPoint.
      */
      template <typename InputPeakType, typename OutputPeakType >
      void filterExperiment(const MSExperiment< InputPeakType >& ms_exp_raw,
                            MSExperiment<OutputPeakType>& ms_exp_filtered)
      {
        // copy the experimental settings
        static_cast<ExperimentalSettings&>(ms_exp_filtered) = ms_exp_raw;

        filterExperiment(ms_exp_raw.begin(), ms_exp_raw.end(), ms_exp_filtered);
      }
	  
	  /** @brief Smoothes an instance of MSExperimentExtern
      	
      Filters the data every scan in the MSExperimentExtern.
      The filtered data are stored in a MSExperimentExtern.
      				
      @note The InputPeakType as well as the OutputPeakType should be of type RawDataPoint1D 
               or any other derived class of DRawDataPoint.
      */
      template <typename InputPeakType, typename OutputPeakType >
      void filterExperiment(const MSExperimentExtern< InputPeakType >& ms_exp_raw,
                            MSExperimentExtern<OutputPeakType>& ms_exp_filtered)
      {
      	filterExperiment(ms_exp_raw.begin(), ms_exp_raw.end(), ms_exp_filtered);
      }

    protected:
      /// The standard derivation  \f$ \sigma \f$.
      DoubleReal sigma_;
      /// The spacing of the pre-tabulated kernel coefficients
      DoubleReal spacing_;
     
      virtual void updateMembers_() 
      {
        DoubleReal kernel_width = (DoubleReal)param_.getValue("gaussian_width"); 
        
        sigma_ = kernel_width / 8.;
        init(sigma_,spacing_);
      }


      /// Computes the value of the gaussian distribution (mean=0 and standard deviation=sigma) at position x
      inline DoubleReal gauss_(DoubleReal x)
      {
        return (1.0/(sigma_ * sqrt(2.0 * M_PI)) * exp(-(x*x) / (2 * sigma_ * sigma_)));
      }

      /// Computes the convolution of the raw data at position x and the gaussian kernel
      template < typename InputPeakIterator >
      DoubleReal integrate_(InputPeakIterator x, InputPeakIterator first, InputPeakIterator last)
      {
        DoubleReal v = 0.;
        // norm the gaussian kernel area to one
        DoubleReal norm = 0.;
        int middle = coeffs_.size();

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
          UInt left_position = (UInt)floor(distance_in_gaussian / spacing_);

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
          int right_position = left_position+1;
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
          left_position = (UInt)floor(distance_in_gaussian / spacing_);

          // search for the true left adjacent data point (because of rounding errors)
          for (int j=0; ((j<3) && (distance(first,help-j) >= 0)); ++j)
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
          int right_position = left_position+1;
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
