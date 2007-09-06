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
//

#ifndef OPENMS_FILTERING_SMOOTHING_SMOOTHFILTER_H
#define OPENMS_FILTERING_SMOOTHING_SMOOTHFILTER_H

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>


namespace OpenMS
{
	/**
		@brief Base class for all noise filter implementations
	*/
  class SmoothFilter : public ProgressLogger
  {
    public:

      /// Constructor
      inline SmoothFilter()
      : coeffs_(0)
      {}

      /// Destructor
      virtual ~SmoothFilter()
      {}

      /// Non-mutable access to the coefficients of the filter
      inline const std::vector<DoubleReal>& getCoeffs() const
      {
        return coeffs_;
      }
      /// Mutable access access to the coefficients of the filter
      inline std::vector<DoubleReal>& getCoeffs()
      {
        return coeffs_;
      }
      /// Mutable access to the coefficients of the filter
      inline void setCoeffs(std::vector<DoubleReal>& coeffs)
      {
        coeffs_ = coeffs;
      }


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

        // needed for multiply the signal with the filter coefficients
        InputPeakIterator it_back;
        typename OutputPeakContainer::iterator out_it = smoothed_data_container.begin();

        int m,i,j;
        DoubleReal help;

        int frame_size = coeffs_.size();
        // compute the transient on
        for (i=0; i<frame_size;++i)
        {
          it_back=first;
          help=0;
          m=0;

          for (j=i; j>=0; --j)
          {
            help+=it_back->getIntensity()*coeffs_[m];
            --it_back;
            ++m;
          }

          out_it->setPosition(first->getPosition());
          out_it->setIntensity(help);
          ++out_it;
          ++first;
        }

        // compute the steady state output
        while (first!=last)
        {
          it_back=first;
          help=0;

          for (j=0; j<frame_size; ++j)
          {
            help+=it_back->getIntensity()*coeffs_[j];
            --it_back;
          }

          out_it->setPosition(first->getPosition());
          out_it->setIntensity(help);
          ++out_it;
          ++first;
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
          InputSpectrumIterator input_it(first+i);

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
            spectrum.setRT(input_it->getRT());
            spectrum.setMSLevel(input_it->getMSLevel());
            spectrum.getName() = input_it->getName();

            ms_exp_filtered.push_back(spectrum);
          }
        }
      }



      /** @brief Filters an MSExperiment.
      	
      Filters the data every scan in the MSExperiment.
      The filtered data are stored in an MSExperiment.
      				
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

    protected:
      /// The coefficient matrix.
      std::vector<DoubleReal> coeffs_;
  };

}// namespace OpenMS


#endif
