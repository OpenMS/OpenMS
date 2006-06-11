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
// $Id: DTopHatFilter.h,v 1.21 2006/05/30 15:46:38 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_FILTERING_BASELINE_DTOPHATFILTER_H
#define OPENMS_FILTERING_BASELINE_DTOPHATFILTER_H

#include <OpenMS/FILTERING/BASELINE/DMorphFilter.h>

#include <algorithm>

namespace OpenMS
{

  /**
  @brief This class represents a Top Hat baseline filter.

   @note This filter works only for uniform raw data!

      This filter can be used by supposing that the required lineaments are brighter than the environment.
      The main advantage of this filter is to be able to detect an over brightness even if the environment is not uniform.
      Moreover it is possible to regulate the size or the width of the over brightnesses very easily.
      The principle is based on the subtraction of an signal \f$ s \f$  from its opening  \f$ \gamma \f$.
      The opening consists of an erosion followed by a dilation,
      the size (the frameSize) of the structuring element (here a flat line) being conditioned by the width of the lineament
      to be detected.

      <i>The DTop Hat-Filter filters the data between two given iterators.
      In the 2-dimensional case that should be the begin and end iterators of the raw data container or a begin iterator of a scan
      and an end iterator of a consecutive scan. Otherwise the correctness of the filtering method is not warrented.
      <b> Only complete consecutive scans can be filtered! A choice of certain m/z chromatograms is at present not possible! <b> </i>
  */
  template <Size D, typename MapType = MSExperiment<DRawDataPoint<1> > >
  class DTopHatFilter : public DMorphFilter<D, MapType>
  {
  public:

    /** @name Type definitions
     */
    //@{
    ///
    typedef DMorphFilter<D, MapType> BaseClass;
    typedef typename BaseClass::DimensionDescription DimensionDescription;
    typedef typename BaseClass::RawDataConstIterator RawDataConstIterator;
    typedef typename BaseClass::RawDataIterator RawDataIterator;
    typedef typename BaseClass::RawData RawData;
    typedef typename BaseClass::MSExperimentFilteredData MSExperimentFilteredData;
    typedef typename BaseClass::MSExperimentRawData MSExperimentRawData;
    typedef typename MapType::const_iterator SpectrumConstIterator;
	typedef typename MapType::SpectrumType SpectrumType;
    ///
    using BaseClass::struc_size_;
    using BaseClass::raw_filtered_;
    using BaseClass::mz_dim_;
    using BaseClass::rt_dim_;
    using BaseClass::param_;
    using BaseClass::ms_exp_filtered_;
    ///
    //@}

    /** @name Constructors and Destructor
     */
    //@{
    ///
    inline DTopHatFilter() : DMorphFilter<D, MapType>()
    {}

    ///
    inline DTopHatFilter(const Param& parameters) : DMorphFilter<D, MapType>(parameters)
    {}
    ///
    inline DTopHatFilter(const DTopHatFilter& t) : DMorphFilter<D, MapType>(t)
    {}
    ///
    inline ~DTopHatFilter()
    {}
    //@}

    /** @name Assignment
     */
    //@{
    ///
    inline DTopHatFilter& operator=(const DTopHatFilter& t)
    {
      param_=t.param_;
      struc_size_=t.struc_size_;
      raw_filtered_=t.raw_filtered_;
      mz_dim_=t.mz_dim_;
      rt_dim_=t.rt_dim_;

      return *this;
    }
    //@}



    /** @name Implements the Filtering method.
     */
    //@{
    /// Compute tophat_filtered_data=signal-opening(signal)
    void tophatMSExperiment(typename SpectrumType::const_iterator first,
                            			  typename SpectrumType::const_iterator last,
                            			  typename SpectrumType::iterator it_new_data,
                            			  DTopHatFilter<1, MapType> const*)
    {
      int number_points = distance(first,last);
      DPeakArrayNonPolymorphic<1,DRawDataPoint<1> > dpa(number_points);

      SpectrumType spectrum;
      spectrum.setContainer(dpa);

      float spacing= ((last-1)->getPosition()[mz_dim_] - first->getPosition()[mz_dim_])/(distance(first,last)+1);
      int struc_elem_number_of_points = (int) ceil(struc_size_ / spacing + 1);
      if (!isOdd(struc_elem_number_of_points))
      {
        struc_elem_number_of_points += 1;
      }

      this->template erosion<SpectrumType >(first,last,spectrum.begin(),struc_elem_number_of_points);
      this->template dilatation<SpectrumType >(spectrum.begin(),spectrum.end(),it_new_data, struc_elem_number_of_points);
      this->template minusIntensities_<SpectrumType >(first,last,it_new_data);
    }

    /// Compute tophat_filtered_data=signal-opening(signal)
    void tophatMSExperiment(typename SpectrumType::const_iterator,
                            			  typename SpectrumType::const_iterator,
                            			  typename SpectrumType::iterator,
                            			  DTopHatFilter<2, MapType> const*) throw (Exception::InvalidValue)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Use the one dimensional tophat filter for MSExperiments","1");
    }

    void tophat(RawDataConstIterator first, RawDataConstIterator last, RawDataIterator new_first)
    {
      RawData erosion_result(distance(first,last));

      float spacing= ((last-1)->getPosition()[mz_dim_] - first->getPosition()[mz_dim_])/(distance(first,last)+1);
      int struc_elem_number_of_points = (int) ceil(struc_size_ / spacing + 1);
      if (!isOdd(struc_elem_number_of_points))
      {
        struc_elem_number_of_points += 1;
      }

      this->template erosion<RawData>(first,last,erosion_result.begin(),struc_elem_number_of_points);
      this->template dilatation<RawData>(erosion_result.begin(),erosion_result.end(),new_first, struc_elem_number_of_points);
      this->template minusIntensities_<RawData>(first,last,new_first);
    }

    ///
    virtual void filter(RawDataConstIterator first, RawDataConstIterator last, RawDataIterator new_first)
    {
      double precision = 1e-5;

      if (D==1)
      {
        tophat(first,last,new_first);
      }
      else
      {
        RawDataConstIterator scan_first=first;
        RawDataConstIterator scan_last=first+1;

        while (scan_first != last)
        {
          // new scan
          if (fabs(scan_first->getPosition()[rt_dim_]-scan_last->getPosition()[rt_dim_]) > precision)
          {
            int scan_length=distance(scan_first,scan_last);
            tophat(scan_first,scan_last,new_first);
            new_first+=scan_length;
            scan_first=scan_last;
          }
          ++scan_last;
        }
      }
    }


    ///
    virtual void filter(const MapType& ms_exp_raw)
    {
      SpectrumConstIterator first_scan=ms_exp_raw.begin();
      SpectrumConstIterator last_scan=ms_exp_raw.end();

      // copy the experimental settings
      //*static_cast<ExperimentalSettings*>(ms_exp_filtered_) = ms_exp_raw;

      while (first_scan != last_scan)
      {
        typename SpectrumType::const_iterator first_data_point = first_scan->begin();
        typename SpectrumType::const_iterator last_data_point = first_scan->end();

        // create the filtered spectrum
        DPeakArrayNonPolymorphic<1,DRawDataPoint<1> > filtered_data(distance(first_data_point,last_data_point));
        tophatMSExperiment(first_data_point,last_data_point,filtered_data.begin(),this);


       	SpectrumType spectrum;
        spectrum.setContainer(filtered_data);
        
        spectrum.setRetentionTime(first_scan->getRetentionTime(), first_scan->getRetentionTimeStart(), first_scan->getRetentionTimeStop());
        spectrum.setMSLevel(first_scan->getMSLevel());
        spectrum.setName(first_scan->getName());

        ms_exp_filtered_->push_back(spectrum);
        
        ++first_scan;
      }
    }
    //@}


  };

}// namespace OpenMS
#endif
