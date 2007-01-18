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

#ifndef OPENMS_FILTERING_SMOOTHING_SAVITZKYGOLAYSVDFILTER_H
#define OPENMS_FILTERING_SMOOTHING_SAVITZKYGOLAYSVDFILTER_H

#include <OpenMS/FILTERING/SMOOTHING/SmoothFilter.h>

#include <OpenMS/KERNEL/MSExperimentExtern.h>

#include <OpenMS/FORMAT/Param.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_pow_int.h>

namespace OpenMS
{
  /**
    @brief Computes the Savitzky-Golay filter coefficients using QR decomposition.
   
    This class represents a Savitzky-Golay lowpass-filter. The idea of the Savitzky-Golay filter
    is to find filtercoefficients that preserve higher moments, which means to approximate the underlying
    function within the moving window by a polynomial of higher order (typically quadratic or quartic).
    Therefore we least-squares fit for each data point a polynomial to all points \f$ f_i \f$ in the window
    and set \f$g_i\f$ to be the value of that polynomial at position \f$ i \f$. This method is superior
    to adjacent averaging because it tends to preserve features of the data such as peak height and width, which
    are usually 'washed out' by adjacent averaging.
   
    Because of the linearity of the problem, we can reduce the work computing by fitting in advance, for fictious data
    consisiting of all zeros except for a singe 1 and then do the fits on the real data just by taking linear
    combinations. There are a particular sets of filter coefficients \f$ c_n \f$ which accomplish the process of
    polynomial least-squares fit inside a moving window. To get the symmetric coefficient-matrix
    \f$C \in R^{frameSize \times frameSize}\f$ with
    \f[ C= \left[ \begin{array}{cccc} c_{0,0} & c_{0,1} & \cdots & c_{0,frameSize-1} \\
                               \vdots  &         &         & \vdots            \\
             c_{frameSize-1,0} & c_{frameSize-1,2} & \ldots  & c_{frameSize-1,frameSize-1} \end{array} \right]\f]
    The first (last) \f$ \frac{frameSize}{2} \f$ rows of \f$ C \f$ we need to smooth the first (last)
    \f$ frameSize \f$ data points of the signal. So we use for the smoothing of the first data point the data
    point itself and the next \f$ frameSize-1 \f$ future points. For the second point we take the first datapoint,
    the data point itself and \f$ frameSize-2 \f$ of rightward data points... .
    We compute the Matrix \f$ C \f$ by solving the underlying least-squares problems with the singular value decomposition.
    Here we demonstrate the computation of the first row of a coefficient-matrix \f$ C \f$ for a Savitzky-Golay Filter
    of order=3 and frameSize=5:
    The design-matrix for the least-squares fit of a linear combination of 3 basis functions to 5 data points is:
    \f[ A=\left[ \begin{array}{ccc} x_0^0 & x_0^1 & x_0^2 \\ \\
                             x_1^0 & x_1^1 & x_1^2 \\ \\
           x_2^0 & x_2^1 & x_2^2 \\ \\
           x_3^0 & x_3^1 & x_3^2 \\ \\
           x_4^0 & x_4^1 & x_4^2 \end{array} \right]. \f]
    To smooth the first data point we have to create a design-matrix with \f$ x=[0,\ldots, frameSize-1] \f$.
    Now we have to solve the over-determined set of \f$ frameSize \f$ linear equations
    \f[ Ac=b \f]
    where \f$ b=[1,0,\ldots,0] \f$ represents the fictious data.
    Therefore we solve the normal equations of the least-squares problem
    \f[ A^TAc=A^Tb. \f]
    Now, it is possible to get
    \f[ c_n=\sum_{m=0}^8 \{(A^TA)^{-1}\}_{0,m} n^m, \f]
    with \f$ 0\le n \le 8\f$. Because we only need one row of the inverse matrix, it is possible to use LU decomposition with
    only a single backsubstitution.
    The vector \f$c=[c_0,\ldots,c_8] \f$ represents the wanted coefficients.
    Note that the solution of a least-squares problem directly from the normal equations is faster than the singular value
    decomposition but rather susceptible to roundoff error!
   
    @note This filter works only for uniform raw data!
          A polynom order of 4 is recommended.
          The bigger the frame size the smoother the signal (the more detail information get lost!). The frame size corresponds to the number
          of filter coefficients, so the width of the smoothing intervall is given by frame_size*spacing of the raw data.
    
    @ingroup Smoothers
   
  */
  class SavitzkyGolaySVDFilter : public SmoothFilter
  {
    public:
      using SmoothFilter::coeffs_;

      /// Constructor
      SavitzkyGolaySVDFilter();

      /// Copy constructor
      inline SavitzkyGolaySVDFilter(const SavitzkyGolaySVDFilter& s)
          : SmoothFilter(s),
          defaults_(s.defaults_),
          frame_size_(s.frame_size_),
          order_(s.order_)
      {
      }

      /// Destructor
      virtual ~SavitzkyGolaySVDFilter()
      {
      }

      /// Assignment operator
      inline SavitzkyGolaySVDFilter& operator=(const SavitzkyGolaySVDFilter& s)
      {
        // take care of self assignments
        if (this == &s)
        {
          return *this;
        }

        defaults_ = s.defaults_;
        frame_size_=s.frame_size_;
        coeffs_=s.coeffs_;
        order_=s.order_;

        return *this;
      }

      
      /// Non-mutable access to the order
      inline const unsigned int& getOrder() const
      {
        return order_;
      }
      
      /// Mutable access to the order
      void setOrder(const unsigned int& order);

      /// Non-mutable access to length of the smoothing window
      inline const unsigned int& getWindowSize() const
      {
        return frame_size_;
      }
      /// Mutable access to the length of the window
      void setWindowSize(const unsigned int& frame_size);

      /// Sets the parameters through a Param
      void setParam(Param param) throw (Exception::InvalidValue);


      /** @brief Applies the convolution with the filter coefficients to an given iterator range.

        Convolutes the filter and the raw data in the iterator intervall [first,last) and writes the
        resulting data to the smoothed_data_container.

        @note This method assumes that the InputPeakIterator (e.g. of type MSSpectrum<DRawDataPoint<1> >::const_iterator)
              points to a data point of type DRawDataPoint<1> or any other class derived from DRawDataPoint<1>.

              The resulting peaks in the smoothed_data_container (e.g. of type MSSpectrum<DRawDataPoint<1> >)
              can be of type DRawDataPoint<1> or any other class derived from DRawDataPoint. 
         
              If you use MSSpectrum iterators you have to set the SpectrumSettings by your own.
         */
      template < typename InputPeakIterator, typename OutputPeakContainer  >
      void filter(InputPeakIterator first, InputPeakIterator last, OutputPeakContainer& smoothed_data_container) throw (Exception::InvalidSize)
      {
        if (distance(first,last) <= (int)frame_size_)
        {
          throw Exception::InvalidSize(__FILE__,__LINE__,__PRETTY_FUNCTION__,distance(first,last));
        }

        smoothed_data_container.resize(distance(first,last));

        int i;
        unsigned int j;
        int mid=(frame_size_/2);
        double help;

        InputPeakIterator it_forward;
        InputPeakIterator it_help;
        typename OutputPeakContainer::iterator out_it = smoothed_data_container.begin();

        // compute the transient on
        for (i=0; i <= mid; ++i)
        {
          it_forward=(first-i);
          help=0;

          for (j=0; j < frame_size_; ++j)
          {
            help+=it_forward->getIntensity()*coeffs_[(i+1)*frame_size_-1-j];
            ++it_forward;
          }


          out_it->setPosition(first->getPos());
          out_it->setIntensity(help);
          ++out_it;
          ++first;
        }

        // compute the steady state output
        it_help=(last-mid);
        while (first!=it_help)
        {
          it_forward=(first-mid);
          help=0;

          for (j=0; j < frame_size_; ++j)
          {
            help+=it_forward->getIntensity()*coeffs_[mid*frame_size_+j];
            ++it_forward;
          }


          out_it->setPosition(first->getPos());
          out_it->setIntensity(help);
          ++out_it;
          ++first;
        }

        // compute the transient off
        for (i=(mid-1); i >= 0; --i)
        {
          it_forward=(first-(frame_size_-i-1));
          help=0;

          for (j=0; j < frame_size_; ++j)
          {
            help+=it_forward->getIntensity()*coeffs_[i*frame_size_+j];
            ++it_forward;
          }

          out_it->setPosition(first->getPos());
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
      void filter(const InputPeakContainer& input_peak_container, OutputPeakContainer& baseline_filtered_container)
      {
        filter(input_peak_container.begin(), input_peak_container.end(), baseline_filtered_container);
      }

      /** @brief Filters every MSSpectrum in a given iterator range.
          		
          	Filters the data successive in every scan in the intervall [first,last).
          	The filtered data are stored in a MSExperiment.
          					
          	@note The InputSpectrumIterator should point to a MSSpectrum. Elements of the input spectren should be of type DRawDataPoint<1> 
                    or any other derived class of DRawDataPoint.

              @note You have to copy the ExperimentalSettings of the raw data on your own. 	
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

          // smooth the peaks in scan i
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

	  /** @brief Filters a MSExperimentExtern.
           	
           Filters the data every scan in MSExperimentExtern.
             
       */
      template <typename InputPeakType, typename OutputPeakType >
      void filterExperiment(const MSExperimentExtern< InputPeakType >& ms_exp_raw,
                            MSExperimentExtern<OutputPeakType>& ms_exp_filtered)
      {
        // copy the experimental settings
        //static_cast<ExperimentalSettings&>(ms_exp_filtered) = ms_exp_raw;

        filterExperiment(ms_exp_raw.begin(), ms_exp_raw.end(), ms_exp_filtered);
      }
	  
	  /** @brief Filters every MSSpectrum in a given iterator range.
          		
          	Filters the data successive in every scan in the intervall [first,last).
          	The filtered data are stored in a MSExperiment.
          					
          	@note The InputSpectrumIterator should point to a MSSpectrum. Elements of the input spectren should be of type DRawDataPoint<1> 
                    or any other derived class of DRawDataPoint.

              @note You have to copy the ExperimentalSettings of the raw data on your own. 	
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
          InputSpectrumIterator input_it(first+i);

          // smooth the peaks in scan i
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

    protected:
      /// parameter defaults
      Param defaults_;
      /// Size of the filter kernel (number of pre-tabulated coefficients)
      unsigned int frame_size_;
      /// The order of the smoothing polynomial.
      unsigned int order_;

      /// Compute the coefficient-matrix \f$ C \f$ of the filter.
      void computeCoeffs_() throw (Exception::InvalidValue);
  };

} // namespace OpenMS
#endif

