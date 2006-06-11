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
// $Id: DSmoothFilter.h,v 1.21 2006/05/30 15:46:38 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_FILTERING_SMOOTHING_DSMOOTHFILTER_H
#define OPENMS_FILTERING_SMOOTHING_DSMOOTHFILTER_H

#include <OpenMS/CONCEPT/Macros.h>

#include <OpenMS/KERNEL/DPeakArrayNonPolymorphic.h>

#include <OpenMS/MATH/MISC/MathFunctions.h>

#include <OpenMS/KERNEL/DimensionDescription.h>

#include <OpenMS/KERNEL/MSExperiment.h>

#include <iostream>

namespace OpenMS
{

  template <Size D, typename MapType = MSExperiment<DRawDataPoint<1> > >
  class DSmoothFilter
  {
  public:

    /** @name Type definitions
     */
    //@{
    ///
    typedef MapType MSExperimentFilteredData;
    ///
    typedef MapType MSExperimentRawData;
	///
	typedef typename MapType::SpectrumType SpectrumType;
    ///
    typedef typename MapType::const_iterator SpectrumConstIterator;
    ///
    typedef OpenMS::DimensionDescription < DimensionDescriptionTagLCMS > DimensionDescription;
    ///
    typedef DPeakArrayNonPolymorphic< D,DRawDataPoint<D> > RawData;
    ///
    typedef typename RawData::Iterator RawDataIterator;
    ///
    typedef typename RawData::ConstIterator RawDataConstIterator;
    ///
    //@}


    /** @name Constructors and Destructor
     */
    //@{
    ///
    inline DSmoothFilter()
        : coeffs_(0),
        raw_filtered_(0)
    {
      if (D == 1)
      {
        rt_dim_ = -1;
        mz_dim_ = 0;
      }
      else
        if (D == 2)
        {
          rt_dim_ = DimensionDescription::RT;
          mz_dim_ = DimensionDescription::MZ;
        }
    }

    ///
    inline DSmoothFilter(const DSmoothFilter& fir)
        : mz_dim_(fir.mz_dim_),
        rt_dim_(fir.rt_dim_)
    {
      coeffs_ = fir.coeffs_;

      // Note: you have to set the member raw_filtered_ by your own
      raw_filtered_ = 0;
    }
    ///
    virtual ~DSmoothFilter(){}

    //@}

    /** @name Assignment
     */
    //@{
    ///
    inline DSmoothFilter& operator=(const DSmoothFilter& fir)
    {
      coeffs_=fir.coeffs_;

      // Note: you have to set the member raw_filtered_ by your own
      raw_filtered_ = 0;
      mz_dim_=fir.mz_dim_;
      rt_dim_=fir.rt_dim_;

      return *this;
    }
    //@}

    /** Accessors
     */
    //@{
    /// Non-mutable access to the filtered raw data
    inline const RawData* getFilteredDataPointer() const { return raw_filtered_; }
    /// Mutable access to the filtered raw data
    inline RawData* getFilteredDataPointer() { return raw_filtered_; }
    /// Mutable access to the filtered raw data
    inline void setFilteredDataPointer(RawData& raw_filtered) { raw_filtered_ = &raw_filtered; }

    /// Non-mutable access to the coefficients of the filter
    inline const std::vector<double>& getCoeffs() const { return coeffs_; }
    /// Mutable access access to the coefficients of the filter
    inline std::vector<double>& getCoeffs() { return coeffs_; }
    /// Mutable access to the coefficients of the filter
    inline void setCoeffs(std::vector<double>& coeffs)
    {
      coeffs_ = coeffs;
    }

    /// Non-mutable access to the the mz dimension
    inline const int& getMZdim() const { return mz_dim_; }
    /// Mutable access to the the mz dimension
    inline int& getMZdim() { return mz_dim_; }
    /// Mutable access to the the mz dimension
    inline void setMZdim(int mz_dim) { mz_dim_=mz_dim; }
    /// Non-mutable access to he rt dimension

    inline const int& getRTdim() const { return rt_dim_; }
    /// Mutable access to he rt dimension
    inline int& getRTdim() { return rt_dim_; }
    /// Mutable access to the the rt dimension
    inline void setRTdim(int rt_dim) { rt_dim_=rt_dim; }
    /// Mutable access to the filtered data

    inline const MSExperimentFilteredData& getMSExperimentFilteredData() const { return *ms_exp_filtered_; }
    /// Mutable access to the filtered data
    inline MSExperimentFilteredData& getMSExperimentFilteredData() { return *ms_exp_filtered_; }
    /// Mutable access to the filtered data
    inline void setMSExperimentFilteredData(const MSExperimentFilteredData& ms_exp_filtered) { ms_exp_filtered_ = &ms_exp_filtered; }
    //@}


    /** @name The Filtering.
     */
    //@{
    /// This method convolutes the signal with the filter-coefficients.
    /// NOTE: Resize the resulting raw data to the length of [first,last)
    void filter(RawDataConstIterator first, RawDataConstIterator last, RawDataIterator new_first)
    {
      double precision = 1e-5;

      if (D==1)
      {
        convolute_(first,last,new_first);
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
            convolute_(scan_first,scan_last,new_first);
            new_first+=scan_length;
            scan_first=scan_last;
          }
          ++scan_last;
        }
      }
    }

    ///
    void filter(const MapType& ms_exp_raw)
    {
      SpectrumConstIterator first_scan=ms_exp_raw.begin();
      SpectrumConstIterator last_scan=ms_exp_raw.end();

      //static_cast<ExperimentalSettings&>(*ms_exp_filtered_) = ms_exp_raw;

      while (first_scan != last_scan)
      {
        typename SpectrumType::const_iterator first_data_point = first_scan->begin();
        typename SpectrumType::const_iterator last_data_point = first_scan->end();


        SpectrumType spectrum;
        // create the filtered spectrum
        // the spectrum should contain at least 2 data points
        if (first_data_point != (last_data_point-1))
        {
          DPeakArrayNonPolymorphic<1,DRawDataPoint<1> > filtered_data(distance(first_data_point,last_data_point));

          startConvolution_(first_data_point,last_data_point,filtered_data.begin(),this);

          spectrum.setContainer(filtered_data);

	        spectrum.setRetentionTime(first_scan->getRetentionTime(), first_scan->getRetentionTimeStart(), first_scan->getRetentionTimeStop());
	        spectrum.setMSLevel(first_scan->getMSLevel());
	        spectrum.setName(first_scan->getName());

        }
        else
        {
          std::cout << "only one data point " << std::endl;
          spectrum = *first_scan;
        }

        ms_exp_filtered_->push_back(spectrum);
        ++first_scan;
      }
    }
    //@}

    ///
    DSmoothFilter& operator()(RawData& raw)
    {
      raw_filtered_ = &raw;
      return *this;
    }

    DSmoothFilter& operator()(MSExperimentFilteredData& ms_exp_filtered)
    {
      ms_exp_filtered_ = &ms_exp_filtered;
      return *this;
    }

  protected:
    /// The coefficient matrix.
    std::vector<double> coeffs_;
    /// MZ dimension
    int mz_dim_;
    /// RT dimension
    int rt_dim_;
    // Pointer to the filtered raw data
    RawData* raw_filtered_;
    // ms_exp_raw_ points to MSExperiment containing the filtered data
    MSExperimentFilteredData* ms_exp_filtered_;

    void startConvolution_(typename SpectrumType::const_iterator first,
                           			typename SpectrumType::const_iterator last,
                           			typename SpectrumType::iterator it_new_data,
                           			DSmoothFilter<1, MapType> const*)
    {
      convolute_(first,last,it_new_data);
    }

    /// Compute tophat_filtered_data=signal-opening(signal)
    void startConvolution_(typename SpectrumType::const_iterator first,
                           			typename SpectrumType::const_iterator last,
                           			typename SpectrumType::iterator it_new_data,
                           			DSmoothFilter<2, MapType> const*) throw (Exception::InvalidValue)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Use the one dimensional smoothing filter for MSExperiments","1");
    }

    /// The discrete convolution
    virtual void convolute_(RawDataConstIterator it_begin, RawDataConstIterator it_end, RawDataIterator new_raw_first)
    {
      // needed for multiply the signal with the filter coefficients
      RawDataConstIterator it_back;

      int m,i,j;
      float help;

      int frame_size = coeffs_.size();
      // compute the transient on
      for (i=0; i<frame_size;++i)
      {
        it_back=it_begin;
        help=0;
        m=0;

        for (j=i; j>=0; --j)
        {
          help+=it_back->getIntensity()*coeffs_[m];
          --it_back;
          ++m;
        }

        new_raw_first->setPosition(it_begin->getPosition());
        new_raw_first->setIntensity(help);
        ++new_raw_first;
        ++it_begin;
      }

      // compute the steady state output
      while (it_begin!=it_end)
      {
        it_back=it_begin;
        help=0;

        for (j=0; j<frame_size; ++j)
        {
          help+=it_back->getIntensity()*coeffs_[j];
          --it_back;
        }

        new_raw_first->setPosition(it_begin->getPosition());
        new_raw_first->setIntensity(help);
        ++new_raw_first;
        ++it_begin;
      }
    }
  };

  //Usage: "ms_exp_raw >> filter(ms_experiment_raw)"
  template <typename MapType>
  const MapType&
  operator>>(const MapType& ms_exp_raw, DSmoothFilter<1, MapType>& m)
  {	
// 		(m.getMSExperimentFilteredData()).setSample(ms_exp_raw.getSample());
// 		(m.getMSExperimentFilteredData()).setSourceFile(ms_exp_raw.getSourceFile());
// 		(m.getMSExperimentFilteredData()).setContacts(ms_exp_raw.getContacts());
// 		(m.getMSExperimentFilteredData()).setInstrument(ms_exp_raw.getInstrument());
// 		(m.getMSExperimentFilteredData()).setSoftware(ms_exp_raw.getSoftware());
// 		(m.getMSExperimentFilteredData()).setProcessingMethod(ms_exp_raw.getProcessingMethod());
// 		(m.getMSExperimentFilteredData()).setHPLC(ms_exp_raw.getHPLC());
// 		(m.getMSExperimentFilteredData()).setType(ms_exp_raw.getType());
// 		(m.getMSExperimentFilteredData()).setDate(ms_exp_raw.getDate());

		m.filter(ms_exp_raw);
		
		return m.getMSExperimentFilteredData();
  }

  /// Usage: "input_raw_data >> filter(output_raw_data)"
  template <Size D, typename MapType>
  const typename DSmoothFilter<D, MapType>::RawData&
  operator>>(const typename DSmoothFilter<D,MapType>::RawData& raw, DSmoothFilter<D, MapType>& m)
  {
    // Resize the resulting raw data
    int array_length=distance(raw.begin(),raw.end());
    (m.getFilteredDataPointer())->resize(array_length);

    m.filter(raw.begin(),raw.end(),(m.getFilteredDataPointer())->begin());

    return *(m.getFilteredDataPointer());
  }


}// namespace OpenMS


#endif
