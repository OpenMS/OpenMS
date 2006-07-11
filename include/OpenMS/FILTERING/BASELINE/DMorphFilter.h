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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#ifndef OPENMS_FILTERING_BASELINE_DMORPHFILTER_H
#define OPENMS_FILTERING_BASELINE_DMORPHFILTER_H

#include <OpenMS/CONCEPT/Macros.h>

#include <OpenMS/KERNEL/DPeakArrayNonPolymorphic.h>

#include <OpenMS/MATH/MISC/MathFunctions.h>

#include <OpenMS/KERNEL/DimensionDescription.h>

#include <OpenMS/FORMAT/Param.h>

#include <OpenMS/KERNEL/MSExperiment.h>


#include <iostream>

namespace OpenMS
{
  using namespace Math;

  /** This class represents a morphological Filter.
       The basic idea of a morphological filter is to inhibit selected signal structures.
       Such structures could be noise or some irrelevant signal structures like the baseline.
       A morphological filter is an increasing operator and has the feature of idempotenz,
       which means that structures which should be received will not be modified by more
       applications of the filter.
       This class provides the basic morphological operations Erosion and Dilatation with a
       structuring element (a flat line) of length frameSize_.

       Erosion and dilatation are implemented using van Herk's method.
   */

  template <Size D , typename MapType = MSExperiment<DRawDataPoint<1> > >
  class DMorphFilter
  {
  public:

    /** @name Type definitions
      */
    //@{
    ///
    typedef DimensionDescription < DimensionDescriptionTagLCMS > DimensionDescription;
    ///
    typedef DPeakArrayNonPolymorphic< D,DRawDataPoint<D> > RawData;
    ///
    typedef typename RawData::iterator RawDataIterator;
    ///
    typedef typename RawData::const_iterator RawDataConstIterator;
    ///
    typedef MapType MSExperimentFilteredData;
    ///
    typedef MapType MSExperimentRawData;
    ///
    //@}


    /** @name Constructors and Destructor
     */
    //@{
    inline DMorphFilter()
        : struc_size_(0),
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
    inline DMorphFilter(const Param& parameters)
    {
      param_ = parameters;

      if (D == 1) // one dimensional picking
      {
        mz_dim_ = 0;
        rt_dim_ = -1;
      }
      else // due to our precondition, we know that D == 2
      {
        rt_dim_= DimensionDescription::RT;
        mz_dim_ = DimensionDescription::MZ;
      }

      // if a parameter is missed in the param object the value should be substituted by a dv value
      DataValue dv = param_.getValue("StrucElementLength");
      if (dv.isEmpty() || dv.toString() == "") struc_size_ = 3;
      else struc_size_ = (float)dv;

      raw_filtered_=0;
    }
    ///
    inline DMorphFilter(const DMorphFilter& m)
        : struc_size_(m.struc_size_),
        mz_dim_(m.mz_dim_),
        rt_dim_(m.rt_dim_)
    {

      // NOTE: The raw_filtered_ pointer points to the same raw data (use setFilteredData())
      raw_filtered_ = m.raw_filtered_;
      param_ = m.param_;

    }
    ///
    virtual ~DMorphFilter()
    { }
    //@}

    /** @name Assignment
      */
    //@{
    ///
    inline DMorphFilter& operator=(const DMorphFilter& m)
    {
      struc_size_=m.struc_size_;
      raw_filtered_ = m.raw_filtered_;
      mz_dim_=m.mz_dim_;
      rt_dim_=m.rt_dim_;
      param_=m.param_;

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

    /// Non-mutable access to length of the structuring element
    inline const float& getStrucElemSize() const { return struc_size_; }
    /// Mutable access to length of the structuring element
    inline float& getStrucElemSize() { return struc_size_; }
    /// Mutable access to the length of the structuring element
    inline void setStrucElemSize(const float& struc_size) { struc_size_=struc_size; }

    /// Non-mutable access to the the mz dimension
    inline const int& getMZdim() const { return mz_dim_; }
    /// Mutable access to the the mz dimension
    inline int& getMZdim() { return mz_dim_; }
    /// Mutable access to the the mz dimension
    inline void setMZdim(const int mz_dim) { mz_dim_=mz_dim; }

    /// Non-mutable access to he rt dimension
    inline const int& getRTdim() const { return rt_dim_; }
    /// Mutable access to he rt dimension
    inline int& getRTdim() { return rt_dim_; }
    /// Mutable access to the the rt dimension
    inline void setRTdim(const int rt_dim) { rt_dim_=rt_dim; }

    /// Non-mutable access to the parameter object
    inline const Param& getParam() const { return param_; }
    /// Mutable access to the parameter object
    inline void setParam(const Param& param)
    {
      param_ = param;

      // read the new parameter
      DataValue dv = param_.getValue("StrucElementLength");
      if (!(dv.isEmpty() || dv.toString() == "")) struc_size_ = (float)dv;
    }

    /// Non-mutable access to the filtered data
  inline const MSExperimentFilteredData& getMSExperimentFilteredData() const { return *ms_exp_filtered_; }
    /// Mutable access to the filtered data
    inline MSExperimentFilteredData& getMSExperimentFilteredData() { return *ms_exp_filtered_; }
    /// Mutable access to the filtered data
    inline void setMSExperimentFilteredData(const MSExperimentFilteredData& ms_exp_filtered) { ms_exp_filtered_ = &ms_exp_filtered; }
    //@}


    ///
    DMorphFilter& operator()(RawData& raw)
    {
      raw_filtered_ = &raw;
      return *this;
    }

    DMorphFilter& operator()(MSExperimentFilteredData& ms_exp_filtered)
    {
      ms_exp_filtered_ = &ms_exp_filtered;
      return *this;
    }

    virtual void filter(RawDataConstIterator first, RawDataConstIterator last, RawDataIterator new_first)=0;
    virtual void filter(const MapType& ms_exp_raw)=0;

    /** @name The Filtering.
     */
    //@{
    /** Van Herk's method of the Dilatation. The algorithm requires 3 min/max comparisons for every data point
        independent from the length of the structuring element.
        Basic idea of the dilatation is "Does the structuring element touch a given set?". The value of a data point
        \f$ x \f$ in the signal \f$s \f$ after a dilatation is the maximal data point in a window which is represented by
        the structuring element \f$ B\f$, when the \f$ B\f$'s point of reference is at \f$ x \f$:
        \f[ [\delta_B(s)](x)=max_{b \in B} s(x+b). \f]
        \image html Dilatation.jpg "Erosion with a structuring element of length 3"
        \image latex Dilatation.eps "Erosion with a structuring element of length 3"
    */
    template<typename ContainerType>
    void dilatation(typename ContainerType::const_iterator itBegin,
                    typename ContainerType::const_iterator itEnd,
                    typename ContainerType::iterator itNew, int l);
    //@}

    //@{
    /**Van Herk's method of the Erosion. The algorithm requires 3 min/max comparisons for every data point
       independent from the length of the structuring element.
       Basic idea of the erosion is "Does the structuring element fit completely in a given set?". The value of a data point
       \f$ x \f$ in the signal \f$s \f$ after an erosion is the minimal data point in a window which is represented by the
       structuring element \f$ B\f$, when the \f$ B\f$'s point of reference is at \f$ x \f$:
       \f[ [\epsilon_B(s)](x)=min_{b \in B} s(x+b). \f]
       \image html Erosion.jpg "Erosion with a structuring element of length 3"
       \image latex Erosion.eps "Erosion with a structuring element of length 3"
    */
    template<typename ContainerType>
    void erosion(typename ContainerType::const_iterator itBegin,
                 typename ContainerType::const_iterator itEnd,
                 typename ContainerType::iterator itNew, int l);
    //@}


  protected:
    ///The length of the structuring element.
    float struc_size_;
    /// MZ dimension
    int mz_dim_;
    /// RT dimension
    int rt_dim_;
    // Pointer to the filtered raw data
    RawData* raw_filtered_;
    // ms_exp_raw_ points to MSExperiment containing the filtered data
    MSExperimentFilteredData* ms_exp_filtered_;
    /// Parameter object
    Param param_;

    ///
    template<typename ContainerType>
    inline void minusIntensities_(typename ContainerType::const_iterator itBegin,
                                  typename ContainerType::const_iterator itEnd,
                                  typename ContainerType::iterator itNew)
    {
      while (itBegin!=itEnd)
      {
        itNew->setIntensity(itBegin->getIntensity()-itNew->getIntensity());
        ++itBegin;
        ++itNew;
      }
    }

    // compute the auxiliary fields g and h
    template<typename ContainerType>
    void calcGErosion_(typename ContainerType::const_iterator itBegin,
                       typename ContainerType::const_iterator itEnd,
                       int l, double* g, bool b);
    //
    template<typename ContainerType>
    void calcHErosion_(typename ContainerType::const_iterator itBegin,
                       int l, double* h, bool b);
    // compute the auxiliary fields g and h
    template<typename ContainerType>
    void calcGDilatation_(typename ContainerType::const_iterator itBegin,
                          typename ContainerType::const_iterator itEnd,
                          int l, double* g, bool b);
    //
    template<typename ContainerType>
    void calcHDilatation_(typename ContainerType::const_iterator itBegin,
                          typename ContainerType::const_iterator itEnd,
                          int l, double* h, bool b);
  };


  template <Size D, typename MapType>
  template <typename ContainerType>
  void  DMorphFilter<D, MapType>::dilatation(typename ContainerType::const_iterator itBegin,
                                    							 typename ContainerType::const_iterator itEnd,
                                    							 typename ContainerType::iterator itNew, int l)
  {
    //--------------van Herk's method of the dilatation --------------------
    RawDataIterator itForward, itBack;

    int middle=l/2;
    int i,k,m,n;
    int length=distance(itBegin,itEnd);

    double *g = new double[l];
    double *h = new double[l];
    k=length-(length%l)-1;


    calcGDilatation_<ContainerType>(itBegin,itEnd,l,g,true);
    calcHDilatation_<ContainerType>(itBegin,itBegin+l-1,l,h,true);

    for (i=0; i<middle; ++i)
    {
      itNew->setIntensity(g[i+middle]);
      itNew->setPosition(itBegin->getPosition());
      ++itNew;
      ++itBegin;
    }


    m=l-1;
    n=0;
    for (i=middle; i<length; ++i)
    {
      if ((i%l)==(middle+1))
      {
        if (i==k)
        {
          calcGDilatation_<ContainerType>((itBegin+middle),itEnd,l,g,false);
        }
        else
        {
          calcGDilatation_<ContainerType>((itBegin+middle),itEnd,l,g,true);
        }
        m=0;
      }
      if ((i%l)==middle && (i>middle))
      {
        if (i>k)
        {
          calcHDilatation_<ContainerType>(itBegin,itEnd,l,h,false);
        }
        else
        {
          calcHDilatation_<ContainerType>((itBegin-middle),(itBegin+middle),l,h,true);
        }
        n=0;
      }

      itNew->setIntensity(std::max(g[m],h[n]));
      itNew->setPosition(itBegin->getPosition());
      ++itBegin;
      ++itNew;
      ++m;
      ++n;
    }
    delete [] g;
    delete [] h;
  }

  template <Size D, typename MapType>
  template <typename ContainerType>
  void  DMorphFilter<D, MapType>::erosion(typename ContainerType::const_iterator itBegin,
                                 							   typename ContainerType::const_iterator itEnd,
                                	 						   typename ContainerType::iterator itNew, int l)
  {
    //-------------- van Herk's method of the erosion --------------------
    RawDataConstIterator itForward=itBegin;
    int middle=l/2;
    int i,k,m,n;
    int length=distance(itBegin,itEnd);

    double *g = new double[l];
    double *h = new double[l];
    k=length-(length%l)-1;

    calcGErosion_<ContainerType>(itBegin,itEnd,l,g,true);
    calcHErosion_<ContainerType>(itBegin+l-1,l,h,true);

    for (i=0; i<middle; ++i)
    {
      itNew->setIntensity(0);
      itNew->setPosition(itBegin->getPosition());
      ++itBegin;
      ++itNew;
    }

    m=l-1;
    n=0;
    for (i=middle; i<length; ++i)
    {
      if ((i%l)==(middle+1))
      {
        if (i==k)
        {
          calcGErosion_<ContainerType>((itBegin+middle),itEnd,l,g,false);
        }
        else
        {
          calcGErosion_<ContainerType>((itBegin+middle),itEnd,l,g,true);
        }
        m=0;
      }
      if ((i%l)==middle && (i>middle) )
      {
        if (i>k)
        {
          calcHErosion_<ContainerType>((itBegin+middle),l,h,false);
        }
        else
        {
          calcHErosion_<ContainerType>((itBegin+middle),l,h,true);
        }
        n=0;
      }

      itNew->setIntensity(std::min(g[m],h[n]));
      itNew->setPosition(itBegin->getPosition());
      ++itNew;
      ++itBegin;
      ++m;
      ++n;
    }

    delete [] g;
    delete [] h;
  }


  template <Size D, typename MapType>
  template<typename ContainerType>
  void  DMorphFilter<D, MapType>::calcGErosion_(typename ContainerType::const_iterator itBegin,
                                  										  typename ContainerType::const_iterator itEnd,
                                  										  int l, double* g, bool b)
  {
    int i,j;

    if (b)
    {
      for (j=0; j<l; ++j)
      {
        if (itBegin < itEnd)
        {
          if (j==0)
          {
            g[j]=itBegin->getIntensity();
          }
          else
          {
            g[j]=std::min(itBegin->getIntensity(),g[j-1]);
          }
          ++itBegin;
        }
        else
        {
          break;
        }
      }
    }
    else
    {
      j=0;
      while (itBegin!=itEnd)
      {
        if (j==0)
        {
          g[j]=itBegin->getIntensity();
        }
        else
        {
          g[j]=std::min(itBegin->getIntensity(),g[j-1]);
        }
        ++itBegin;
        ++j;
      }

      for (i=j; i<l; ++i)
      {
        g[i]=0;
      }
    }
  }

  template <Size D, typename MapType>
  template <typename ContainerType>
  void  DMorphFilter<D, MapType>::calcHErosion_(typename ContainerType::const_iterator itBegin,
                                                                          int l, double* h, bool b)
  {
    int j;
    if (b)
    {
      for (j=l-1; j>=0; --j)
      {
        if (j==(l-1))
        {
          h[j]=itBegin->getIntensity();
        }
        else
        {
          h[j]=std::min(itBegin->getIntensity(),h[j+1]);
        }
        --itBegin;
      }
    }
    else
    {
      for (j=0;j<l;++j)
      {
        h[j]=0;
      }
    }
  }


  /// compute the auxiliary fields g and h
  template <Size D, typename MapType>
  template <typename ContainerType>
  void  DMorphFilter<D, MapType >::calcGDilatation_(typename ContainerType::const_iterator itBegin,
                                          									 typename ContainerType::const_iterator itEnd,
                                          									 int l, double* g, bool b)
  {
    int i,j;

    if (b)
    {
      for (j=0; j<l; ++j)
      {
        if (itBegin < itEnd)
        {
          if (j==0)
          {
            g[j]=itBegin->getIntensity();
          }
          else
          {
            g[j]=std::max(itBegin->getIntensity(),g[j-1]);
          }
          ++itBegin;
        }
        else
        {
          break;
        }
      }
    }
    else
    {
      j=0;
      while (itBegin!=itEnd)
      {
        if (j==0)
        {
          g[j]=itBegin->getIntensity();
        }
        else
        {
          g[j]=std::max(itBegin->getIntensity(),g[j-1]);
        }
        ++itBegin;
        ++j;
      }
      for (i=j; i<l; ++i)
      {
        g[i]=g[j-1];
      }
    }
  }

  template <Size D, typename MapType>
  template <typename ContainerType>
  void  DMorphFilter<D, MapType>::calcHDilatation_(typename ContainerType::const_iterator itBegin,
                                          									typename ContainerType::const_iterator itEnd,
                                         									int l, double* h, bool b)
  {
    int j;

    if (b)
    {
      for (j=l-1; j>=0; --j)
      {
        if (j==(l-1))
        {
          h[j]=itEnd->getIntensity();
        }
        else
        {
          h[j]=std::max(itEnd->getIntensity(),h[j+1]);
        }
        --itEnd;
      }
    }
    else
    {
      j=(itEnd-itBegin)-1;
      h[j]=(--itEnd)->getIntensity();
      while (itEnd!=itBegin)
      {
        --j;
        h[j]=std::max(itBegin->getIntensity(),h[j+1]);;
        --itEnd;
      }
    }
  }

  /**
  @brief Start baseline filtering of a whole MS experiment

  Usage: "ms_raw_data >> filter(ms_filtered_data)"
  */
  ///
  template <typename MapType>
  const MapType &
  operator>>(const MapType& ms_exp_raw, DMorphFilter<1, MapType>& m)
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

  /**
  @brief Start baseline filtering of a DPeakArrayNonPolymorphic of DRawDataPoints

  Usage: "ms_raw_data >> filter(ms_filtered_data)"
  */
  template <Size D>
  const typename DMorphFilter<D>::RawData&
  operator>>(const typename DMorphFilter<D>::RawData& raw, DMorphFilter<D>& m)
  {
    // Resize the resulting raw data
    int array_length=distance(raw.begin(),raw.end());
    (m.getFilteredDataPointer())->resize(array_length);

    m.filter(raw.begin(),raw.end(),(m.getFilteredDataPointer())->begin());

    return *(m.getFilteredDataPointer());
  }





}

#endif
