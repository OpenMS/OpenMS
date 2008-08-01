// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_FILTERING_BASELINE_TOPHATFILTER_H
#define OPENMS_FILTERING_BASELINE_TOPHATFILTER_H

#include <OpenMS/FILTERING/BASELINE/MorphFilter.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>

#include <algorithm>

namespace OpenMS
{

  /**
  	@brief This class represents a Top Hat baseline filter.
   
    This filter can be used by supposing that the required lineaments are brighter than the environment.
    The main advantage of this filter is to be able to detect an over brightness even if the environment is not uniform.
    Moreover it is possible to regulate the size or the width of the over brightnesses very easily.
    The principle is based on the subtraction of an signal \f$ s \f$  from its opening  \f$ \gamma \f$.
    The opening consists of an erosion followed by a dilation,
    the size (the frameSize) of the structuring element (here a flat line) being conditioned by the width of the lineament
    to be detected.
    
    @note This filter works only for uniform raw data!
    
    @ref TopHatFilter_Parameters are explained on a separate page.
    
		@ingroup SignalProcessing
  */
  class TopHatFilter 
  	: public MorphFilter
  {
    public:
      typedef MorphFilter BaseClass;
      using BaseClass::struc_size_;

      /// Constructor
      inline TopHatFilter() 
      	: MorphFilter()
      {
      }

      /// Destructor
      virtual ~TopHatFilter()
      {
      }

      /** 
      	@brief Applies the baseline removal algorithm to an given iterator range.

        Removes the baseline in the given iterator intervall [first,last) and writes the
        resulting data to the baseline_filtered_container.
        
        @note This method assumes that the InputPeakIterator (e.g. of type MSSpectrum<Peak1D >::const_iterator)
              points to a data point of type Peak1D or any other class derived from Peak1D.
        
        @note The resulting peaks in the baseline_filtered_container (e.g. of type MSSpectrum<Peak1D >)
              can be of type Peak1D or any other class derived from Peak1D.
      */
      template <typename InputPeakIterator, typename OutputPeakContainer  >
      void filter(InputPeakIterator first, InputPeakIterator last, OutputPeakContainer& baseline_filtered_container)
      {
      	baseline_filtered_container.clear();
      	
				//filter the baseline only if the scan contains enough raw data points
        if ( distance(first,last)==0 || struc_size_ > fabs(first->getMZ()-(last-1)->getMZ()))
        {
					return;
        }
      
        // compute the number of data points of the structuring element given the spacing of the raw data
        // and the size (in Th) of the structuring element
        DoubleReal spacing= ((last-1)->getMZ() - first->getMZ()) / (distance(first,last)-1);
        int struc_elem_number_of_points = (int) ceil(struc_size_ / spacing );

        // the number has to be odd
        if (!Math::isOdd(struc_elem_number_of_points))
        {
          struc_elem_number_of_points += 1;
        }

        std::vector<typename InputPeakIterator::value_type> erosion_result;
        // compute the erosion of raw data
        this->template erosion(first, last, erosion_result, struc_elem_number_of_points);
        // compute the dilation of erosion_result
        this->template dilatation(erosion_result.begin(),erosion_result.end(), baseline_filtered_container, struc_elem_number_of_points);
        // subtract the result from the original data
        this->template minusIntensities_(first,last,baseline_filtered_container);
      }


      /** 
      	@brief Applies the baseline removal algorithm to to a raw data point container.

        Removes the baseline in the the input container (e.g. of type MSSpectrum<Peak1D >) and writes the 
        resulting data to the baseline_filtered_container.
        
        @note This method assumes that the InputPeakIterator (e.g. of type MSSpectrum<Peak1D >::const_iterator)
              points to a data point of type Peak1D or any other class derived from Peak1D.
        
        @note The resulting peaks in the baseline_filtered_container (e.g. of type MSSpectrum<Peak1D >)
              can be of type Peak1D or any other class derived from DPeak. 
      */
      template <typename InputPeakContainer, typename OutputPeakContainer >
      void filter(const InputPeakContainer& input_peak_container, OutputPeakContainer& baseline_filtered_container)
      {
        // copy the experimental settings
        static_cast<SpectrumSettings&>(baseline_filtered_container) = input_peak_container;
        baseline_filtered_container.setType(SpectrumSettings::RAWDATA);
        filter(input_peak_container.begin(), input_peak_container.end(), baseline_filtered_container);
      }


      /** 
      	@brief Convenience method that removes the baseline from an MSExperiment containing raw data.
      */
      template <typename PeakType>
      void filterExperiment(MSExperiment<PeakType>& map)
      {
        startProgress(0,map.size(),"filtering baseline");
        for (UInt i = 0; i < map.size(); ++i)
        {
          typename MSExperiment<PeakType>::SpectrumType spectrum;
          filter(map[i],spectrum);
          map[i].getContainer() = spectrum.getContainer();
          setProgress(i);
        }
        endProgress();
      }
  };

}// namespace OpenMS
#endif
