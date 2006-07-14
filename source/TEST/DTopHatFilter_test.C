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
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FILTERING/BASELINE/DTopHatFilter.h>

#include <OpenMS/KERNEL/MSExperiment.h>

///////////////////////////

START_TEST(DTopHatFilter<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

DTopHatFilter<2>* tophat_ptr = 0;

CHECK(DTopHatFilter())
  tophat_ptr = new DTopHatFilter<2>;
  TEST_NOT_EQUAL(tophat_ptr, 0) // ???
RESULT

CHECK(~DTopHatFilter())
  delete tophat_ptr;
RESULT

CHECK(DTopHatFilter(const Param& p))
  Param p;
  p.setValue("StrucElementLength",3);
  DTopHatFilter<2> tophat(p);
  TEST_EQUAL(tophat.getStrucElemSize(), 3)
RESULT

CHECK(DTopHatFilter(const DTopHatFilter& t))
    DPeakArrayNonPolymorphic<2,DRawDataPoint<2> > tophat_data;
    DTopHatFilter<2> tophat;
    tophat.setStrucElemSize(3);
    tophat(tophat_data);

    DTopHatFilter<2> tophat_copy(tophat);
    TEST_EQUAL(tophat_copy.getFilteredDataPointer(), tophat.getFilteredDataPointer())
    TEST_EQUAL(tophat_copy.getStrucElemSize(), 3)
    TEST_EQUAL(tophat_copy.getRTdim(), 0)
    TEST_EQUAL(tophat_copy.getMZdim(), 1)
RESULT

CHECK(DTopHatFilter& operator=(const DTopHatFilter& t))
    DPeakArrayNonPolymorphic<2,DRawDataPoint<2> > tophat_data;
    DTopHatFilter<2> tophat;
    tophat.setStrucElemSize(3);
    tophat(tophat_data);

    DTopHatFilter<2> tophat_copy;
    tophat_copy = tophat;
    TEST_EQUAL(tophat_copy.getFilteredDataPointer(), tophat.getFilteredDataPointer())
    TEST_EQUAL(tophat_copy.getStrucElemSize(), 3)
    TEST_EQUAL(tophat_copy.getRTdim(), 0)
    TEST_EQUAL(tophat_copy.getMZdim(), 1)
RESULT

CHECK(const RawData&  operator>>(const RawData& raw, DMorphFilter& m))
RESULT


CHECK((void filter(RawDataConstIterator first, RawDataConstIterator last, RawDataIterator new_first)))
    DPeakArrayNonPolymorphic<1,DRawDataPoint<1> > raw_data;
    int i;
    for (i=0; i < 24; ++i)
    {
      DRawDataPoint<1> p;
      DPosition<1> pos = i;
      if ((1<i) && (i<5))
        {
          p.setIntensity(1);
        }
      else
        {
          p.setIntensity(0);
        }
      p.setPosition(pos);
      raw_data.push_back(p);
    }

    DPeakArrayNonPolymorphic<1,DRawDataPoint<1> > tophat_data(24);
    DTopHatFilter<1> tophat;
    tophat.setStrucElemSize(3);
    tophat.filter(raw_data.begin(),raw_data.end(),tophat_data.begin());

    DPeakArrayNonPolymorphic<1,DRawDataPoint<1> >::ConstIterator it=tophat_data.begin();
    for (int i=0; i<24; ++i)
    {
      TEST_REAL_EQUAL(it->getIntensity(), 0)
    }
RESULT


CHECK((void tophat(RawDataConstIterator scan_beg, RawDataConstIterator scan_end, RawDataArray& it_ero, RawDataIterator it_new_data)))
    DPeakArrayNonPolymorphic<1,DRawDataPoint<1> > raw_data;
    int i;
    for (i=0; i<8; ++i)
    {
      DRawDataPoint<1> p;
      DPosition<1> pos;
      pos=i;

      if ( (1<i) && (i<5))
        {
          p.setIntensity(1);
        }
      else
        {
          p.setIntensity(0);
        }
      p.setPosition(pos);
      raw_data.push_back(p);
     }

     DPeakArrayNonPolymorphic<1,DRawDataPoint<1> > tophat_data(8);

     DTopHatFilter<1> tophat;
     tophat.setStrucElemSize(3);
     tophat.tophat(raw_data.begin(),raw_data.end(),tophat_data.begin());

     DPeakArrayNonPolymorphic<1,DRawDataPoint<1> >::ConstIterator it=tophat_data.begin();
     for (int i=0; i < 8; ++i)
       {
         TEST_REAL_EQUAL(it->getIntensity(), 0)
       }
RESULT

CHECK(void filter(const MSExperiment<DRawDataPoint<1> >& ms_exp_raw))
  MSExperiment<DRawDataPoint<1> > ms_exp_raw;
  MSExperiment<DRawDataPoint<1> > ms_exp_filtered;

  DPeakArrayNonPolymorphic<2,DRawDataPoint<2> > raw_data;
  int i;
  for (i=0; i < 8; ++i)
  {
    DRawDataPoint<2> p;
    DPosition<2> pos;
    pos[0]=10.;
    pos[1]=i;
    p.setPosition(pos);

    if ( (1<i) && (i<5))
    {
      p.setIntensity(1);
    }
    else
    {
      p.setIntensity(0);
    }
    raw_data.push_back(p);
  }

  ms_exp_raw.set2DData(raw_data);

  DTopHatFilter<1> tophat;
  tophat.setStrucElemSize(3);

  tophat(ms_exp_filtered);
  tophat.filter(ms_exp_raw);

  DPeakArrayNonPolymorphic<2> dpeak_arra_filtered;
  ms_exp_filtered.get2DData(dpeak_arra_filtered);
  DPeakArrayNonPolymorphic<2>::iterator it = dpeak_arra_filtered.begin();
  for (int i=0; i<8; ++i)
  {
    TEST_REAL_EQUAL(it->getIntensity(), 0)
   }
RESULT

CHECK((void tophatMSExperiment(typename MSSpectrum<DRawDataPoint<1> >::const_it
erator first, typename MSSpectrum<DRawDataPoint<1> >::const_iterator last, type
name MSSpectrum<DRawDataPoint<1> >::iterator it_new_data, DTopHatFilter<1> cons
t*)))
  MSExperiment<DRawDataPoint<1> > ms_exp_raw;
  MSExperiment<DRawDataPoint<1> > ms_exp_filtered;

  DPeakArrayNonPolymorphic<2,DRawDataPoint<2> > raw_data;
    int i;
    for (i=0; i<8; ++i)
    {
      DRawDataPoint<2> p;
      DPosition<2> pos;
      pos[0]=10.;
      pos[1]=i;
      if ( (1<i) && (i<5))
        {
          p.setIntensity(1);
        }
      else
        {
          p.setIntensity(0);
        }
      p.setPosition(pos);
      raw_data.push_back(p);
    }

    ms_exp_raw.set2DData(raw_data);

    DPeakArrayNonPolymorphic<2> dpeak_arra_filtered(8);
    ms_exp_filtered.set2DData(dpeak_arra_filtered);

    DTopHatFilter<1> tophat;
    tophat.setStrucElemSize(3);

    tophat(ms_exp_filtered);
    tophat.tophatMSExperiment(ms_exp_raw[0].begin(),ms_exp_raw[0].end(), ms_exp_filtered[0].begin(),&tophat);

    ms_exp_filtered.get2DData(dpeak_arra_filtered);
    DPeakArrayNonPolymorphic<2>::iterator it = dpeak_arra_filtered.begin();
    for (int i=0; i<8; ++i)
       {
         TEST_REAL_EQUAL(it->getIntensity(), 0)
       }
RESULT

CHECK((void tophatMSExperiment(typename MSSpectrum<DRawDataPoint<1> >::const_it
erator, typename MSSpectrum<DRawDataPoint<1> >::const_iterator, typename MSSpec
trum<DRawDataPoint<1> >::iterator, DTopHatFilter<2> const*)))
  MSExperiment<DRawDataPoint<1> > ms_exp_raw;
  MSExperiment<DRawDataPoint<1> > ms_exp_filtered;

  DPeakArrayNonPolymorphic<2,DRawDataPoint<2> > raw_data;
  int i;
  for (i=0; i<8; ++i)
    {
      DRawDataPoint<2> p;
      DPosition<2> pos;
      p.setIntensity(0);
      if ( (1<i) && (i<5))
        {
          pos[0]=10.;
          pos[1]=1;
        }
      else
        {
          pos[0]=10.;
          pos[1]=0;
        }
      p.setPosition(pos);
      raw_data.push_back(p);
    }

    ms_exp_raw.set2DData(raw_data);

    DPeakArrayNonPolymorphic<2> dpeak_arra_filtered(8);
    ms_exp_filtered.set2DData(dpeak_arra_filtered);

    DTopHatFilter<2> tophat;
    tophat(ms_exp_filtered);

    TEST_EXCEPTION(Exception::InvalidValue,tophat.tophatMSExperiment(ms_exp_raw[0].begin(),ms_exp_raw[0].end(), ms_exp_filtered[0].begin(),&tophat))
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

