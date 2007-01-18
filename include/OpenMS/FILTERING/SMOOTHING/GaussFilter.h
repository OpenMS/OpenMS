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
#include <OpenMS/FORMAT/Param.h>

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
          
    @ingroup Smoothers
  */
//#define DEBUG_FILTERING

  class GaussFilter : public SmoothFilter
  {
    public:
      using SmoothFilter::coeffs_;

      /// Constructor
      inline GaussFilter()
      {
      	//Parameter settings
      	defaults_.setValue("gaussian_width",0.8);
        
        //members
        sigma_ = .1;
        spacing_ = 0.01;

        //compute the filter kernel coefficients
        init(sigma_,spacing_);
      }

      /** @brief Constructor given a param object.
          
          @note Please note that the frame_size must be odd.
      */
      GaussFilter(const Param& parameters) throw (Exception::InvalidValue);

      /// Copy constructor
      inline GaussFilter(const GaussFilter& g)
          : SmoothFilter(g),
          sigma_(g.sigma_),
          spacing_(g.spacing_),
          defaults_(g.defaults_)
      {
      }

      /// Destructor
      virtual ~GaussFilter()
      {
      }

      /// Assignment operator
      inline GaussFilter& operator=(const GaussFilter& s)
      {
        // take care of self assignments
        if (this == &s) return *this;

        defaults_ = s.defaults_;

        spacing_=s.spacing_;
        coeffs_=s.coeffs_;
        sigma_=s.sigma_;

        return *this;
      }

      /// Non-mutable access to the sigma
      inline const double& getSigma() const
      {
        return sigma_;
      }
      /// Mutable access to the sigma
      inline void setSigma(const double& sigma)
      {
        sigma_ = sigma;
        spacing_ = 4*sigma_ / 50;
        init(sigma_,spacing_);
      }
      
      /// Non-mutable access to the kernel width
      inline double getKernelWidth() const
      {
        return (sigma_ * 8.);
      }
      /// Mutable access to the kernel width
      inline void setKernelWidth(const double& kernel_width)
      {
        sigma_ = kernel_width / 8.;
        init(sigma_,spacing_);
      }
      
      /// Non-mutable access to the spacing
      inline const double& getSpacing() const
      {
        return spacing_;
      }
      /// Mutable access to the spacing
      inline void setSpacing(const double& spacing)
      {
        spacing_=spacing;
        OPENMS_PRECONDITION((4*sigma_ > spacing), "You have to choose a smaller spacing for the kernel coefficients!" );
        init(sigma_,spacing_);
      }
      
      /// Sets the parameters through a Param
      inline void setParam(Param param) throw (Exception::InvalidValue)
      {
      	param.setDefaults(defaults_);
      	param.checkDefaults("GaussFilter", defaults_);
      	
        DataValue dv = param.getValue("gaussian_width");
        if (dv.toString() != "")
        {
          double kernel_width = (float)dv;

          if (kernel_width <= 0)
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"The kernel width should be greater than zero!",String(kernel_width));
          }

          sigma_ = kernel_width / 8.;
          init(sigma_,spacing_);
        }
      }


      /** @brief Build a gaussian distribution for the current spacing and standard deviation.
          
          We store the coefficiens of gaussian in the vector<double> coeffs_;

          We only need a finite amount of points since the gaussian distribution
          decays fast. We take 4*sigma (99.993666% of the area is within four standard deviations), since at that point the function
          has dropped to ~ -10^-4
      */
      void init(float sigma, float spacing);


      /** @brief Applies the convolution with the filter coefficients to an given iterator range.

      Convolutes the filter and the raw data in the iterator intervall [first,last) and writes the
      resulting data to the smoothed_data_container.

      @note This method assumes that the InputPeakIterator (e.g. of type MSSpectrum<DRawDataPoint<1> >::const_iterator)
            points to a data point of type DRawDataPoint<1> or any other class derived from DRawDataPoint<1>.

            The resulting peaks in the smoothed_data_container (e.g. of type MSSpectrum<DRawDataPoint<1> >)
            can be of type DRawDataPoint<1> or any other class derived from DRawDataPoint. 
       
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
          out_it->setPosition(help->getPos());
          out_it->setIntensity(integrate_(help,first,last));
          ++out_it;
          ++help;
        }
      }


      /** @brief Convolutes the filter coefficients and the input raw data.

         Convolutes the filter and the raw data in the input_peak_container and writes the
          resulting data to the smoothed_data_container.

      @note This method assumes that the elements of the InputPeakContainer (e.g. of type MSSpectrum<DRawDataPoint<1> >)
            are of type DRawDataPoint<1> or any other class derived from DRawDataPoint<1>.

            The resulting peaks in the smoothed_data_container (e.g. of type MSSpectrum<DRawDataPoint<1> >)
            can be of type DRawDataPoint<1> or any other class derived from DRawDataPoint. 
       
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
      					
      	@note The InputSpectrumIterator should point to a MSSpectrum. Elements of the input spectren should be of type DRawDataPoint<1> 
                or any other derived class of DRawDataPoint.

          @note You have to copy the ExperimentalSettings of the raw data by your own. 	
      */
      template <typename InputSpectrumIterator, typename OutputPeakType >
      void filterExperiment(InputSpectrumIterator first,
                            		InputSpectrumIterator last,
                            		MSExperiment<OutputPeakType>& ms_exp_filtered)
      {
        unsigned int n = distance(first,last);
        // pick peaks on each scan
        for (unsigned int i = 0; i < n; ++i)
        {
          MSSpectrum< OutputPeakType > spectrum;
          InputSpectrumIterator input_it = first+i;

          // pick the peaks in scan i
          filter(*input_it,spectrum);

          // if any peaks are found copy the spectrum settings
          if (spectrum.size() > 0)
          {
            // copy the spectrum settings
            static_cast<SpectrumSettings&>(spectrum) = *input_it;
            spectrum.setType(SpectrumSettings::RAWDATA);

            // copy the spectrum information
            spectrum.getPrecursorPeak() = input_it->getPrecursorPeak();
            spectrum.setRetentionTime(input_it->getRetentionTime());
            spectrum.setMSLevel(input_it->getMSLevel());
            spectrum.getName() = input_it->getName();

            ms_exp_filtered.push_back(spectrum);
          }
        }
      }
	  
	     /** @brief Filters every MSSpectrum in a given iterator range.
      		
      	Filters the data successive in every scan in the intervall [first,last).
      	The filtered data are stored in a MSExperiment.
      					
      	@note The InputSpectrumIterator should point to a MSSpectrum. Elements of the input spectren should be of type DRawDataPoint<1> 
                or any other derived class of DRawDataPoint.

          @note You have to copy the ExperimentalSettings of the raw data by your own. 	
      */
      template <typename InputSpectrumIterator, typename OutputPeakType >
      void filterExperiment(InputSpectrumIterator first,
                            		InputSpectrumIterator last,
                            		MSExperimentExtern<OutputPeakType>& ms_exp_filtered)
      {
        unsigned int n = distance(first,last);
        // pick peaks on each scan
        for (unsigned int i = 0; i < n; ++i)
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
            spectrum.setRetentionTime(input_it->getRetentionTime());
            spectrum.setMSLevel(input_it->getMSLevel());
            spectrum.getName() = input_it->getName();

            ms_exp_filtered.push_back(spectrum);
          }
        }
      }



      /** @brief Filters a MSExperiment.
      	
      Filters the data every scan in the MSExperiment.
      The filtered data are stored in a MSExperiment.
      				
      @note The InputPeakType as well as the OutputPeakType should be of type DRawDataPoint<1> 
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
      				
      @note The InputPeakType as well as the OutputPeakType should be of type DRawDataPoint<1> 
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
      double sigma_;
      /// The spacing of the pre-tabulated kernel coefficients
      double spacing_;
      /// Parameter object
      Param defaults_;


      /// Computes the value of the gaussian distribution (mean=0 and standard deviation=sigma) at position x
      inline double gauss_(double x)
      {
        return (1.0/(sigma_ * sqrt(2.0 * M_PI)) * exp(-(x*x) / (2 * sigma_ * sigma_)));
      }

      /// Computes the convolution of the raw data at position x and the gaussian kernel
      template < typename InputPeakIterator >
      double integrate_(InputPeakIterator x, InputPeakIterator first, InputPeakIterator last)
      {
        double v = 0.;
        // norm the gaussian kernel area to one
        double norm = 0.;
        int middle = coeffs_.size();

        double start_pos = ((x->getPos()-(middle*spacing_)) > first->getPos()) ? (x->getPos()-(middle*spacing_))
                           : first->getPos();
        double end_pos = ((x->getPos()+(middle*spacing_)) < (last-1)->getPos()) ? (x->getPos()+(middle*spacing_))
                         : (last-1)->getPos();


        InputPeakIterator help = x;
#ifdef DEBUG_FILTERING

        std::cout << "integrate from middle to start_pos "<< help->getPos() << " until " << start_pos << std::endl;
#endif

        //integrate from middle to start_pos
        while ((help != first) && ((help-1)->getPos() > start_pos))
        {
          // search for the corresponding datapoint of help in the gaussian (take the left most adjacent point)
          double distance_in_gaussian = fabs(x->getPos() - help->getPos());
          unsigned int left_position = (unsigned int)floor(distance_in_gaussian / spacing_);

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
          double d = fabs((left_position*spacing_)-distance_in_gaussian) / spacing_;
          // check if the right data point in the gaussian exists
          double coeffs_right = (right_position < middle) ? (1-d)*coeffs_[left_position]+d*coeffs_[right_position]
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
          distance_in_gaussian = fabs(x->getPos() - (help-1)->getPos());
          left_position = (unsigned int)floor(distance_in_gaussian / spacing_);

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
          double coeffs_left= (right_position < middle) ? (1-d)*coeffs_[left_position]+d*coeffs_[right_position]
                              : coeffs_[left_position];
#ifdef DEBUG_FILTERING

          std::cout << " help-1 " << (help-1)->getPos() << " distance_in_gaussian " << distance_in_gaussian << std::endl;
          std::cout << " right_position " << right_position << std::endl;
          std::cout << " left_position " << left_position << std::endl;
          std::cout << "coeffs_ at left_position " <<  coeffs_[left_position]<<std::endl;
          std::cout << "coeffs_ at right_position " <<   coeffs_[right_position]<<std::endl;
          std::cout << "interpolated value right " << coeffs_left << std::endl;

          std::cout << " intensity " << fabs((help-1)->getPos()-help->getPos()) / 2. << " * " << (help-1)->getIntensity() << " * " << coeffs_left <<" + " << (help)->getIntensity()<< "* " << coeffs_right
          << std::endl;
#endif


          norm += fabs((help-1)->getPos()-help->getPos()) / 2. * (coeffs_left + coeffs_right);

          v+= fabs((help-1)->getPos()-help->getPos()) / 2. * ((help-1)->getIntensity()*coeffs_left + help->getIntensity()*coeffs_right);
          --help;
        }


        //integrate from middle to end_pos
        help = x;
#ifdef DEBUG_FILTERING

        std::cout << "integrate from middle to endpos "<< (help)->getPos() << " until " << end_pos << std::endl;
#endif

        while ((help != (last-1)) && ((help+1)->getPos() < end_pos))
        {
          // search for the corresponding datapoint for help in the gaussian (take the left most adjacent point)
          double distance_in_gaussian = fabs(x->getPos() - help->getPos());
          int left_position = (unsigned int)floor(distance_in_gaussian / spacing_);

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
          double d = fabs((left_position*spacing_)-distance_in_gaussian) / spacing_;
          double coeffs_left= (right_position < middle) ? (1-d)*coeffs_[left_position]+d*coeffs_[right_position]
                              : coeffs_[left_position];

#ifdef DEBUG_FILTERING

          std::cout << " help " << (help)->getPos() << " distance_in_gaussian " << distance_in_gaussian << std::endl;
          std::cout << " left_position " << left_position << std::endl;
          std::cout << "coeffs_ at right_position " <<  coeffs_[left_position]<<std::endl;
          std::cout << "coeffs_ at left_position " <<  coeffs_[right_position]<<std::endl;
          std::cout << "interpolated value left " << coeffs_left << std::endl;
#endif

          // search for the corresponding datapoint for (help+1) in the gaussian (take the left most adjacent point)
          distance_in_gaussian = fabs(x->getPos() - (help+1)->getPos());
          left_position = (unsigned int)floor(distance_in_gaussian / spacing_);

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
          double coeffs_right = (right_position < middle) ? (1-d)*coeffs_[left_position]+d*coeffs_[right_position]
                                : coeffs_[left_position];
#ifdef DEBUG_FILTERING

          std::cout << " (help + 1) " << (help+1)->getPos() << " distance_in_gaussian " << distance_in_gaussian << std::endl;
          std::cout << " left_position " << left_position << std::endl;
          std::cout << "coeffs_ at right_position " <<   coeffs_[left_position]<<std::endl;
          std::cout << "coeffs_ at left_position " <<  coeffs_[right_position]<<std::endl;
          std::cout << "interpolated value right " << coeffs_right << std::endl;

          std::cout << " intensity " <<  fabs(help->getPos() - (help+1)->getPos()) / 2.
          << " * " << help->getIntensity() << " * " << coeffs_left <<" + " << (help+1)->getIntensity()
          << "* " << coeffs_right
          << std::endl;
#endif
          norm += fabs(help->getPos() - (help+1)->getPos()) / 2. * (coeffs_left + coeffs_right);

          v+= fabs(help->getPos() - (help+1)->getPos()) / 2. * (help->getIntensity()*coeffs_left + (help+1)->getIntensity()*coeffs_right);
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
