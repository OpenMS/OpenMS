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

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FILTERING/BASELINE/TopHatFilter.h>

#include <OpenMS/KERNEL/MSExperiment.h>

///////////////////////////

START_TEST(TopHatFilter<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

TopHatFilter* tophat_ptr = 0;
CHECK((TopHatFilter()))
  tophat_ptr = new TopHatFilter;
  TEST_NOT_EQUAL(tophat_ptr, 0) 
RESULT

CHECK((~TopHatFilter()))
  delete tophat_ptr;
RESULT

CHECK((TopHatFilter(const Param& parameters)))
  Param p;
  p.setValue("struc_elem_length",3);
  TopHatFilter tophat(p);
  
  TEST_EQUAL(tophat.getStrucElemSize(), 3)
RESULT

CHECK((TopHatFilter(const TopHatFilter& t)))
  TopHatFilter tophat;
  tophat.setStrucElemSize(3);

  TopHatFilter tophat_copy(tophat);
  TEST_EQUAL(tophat_copy.getStrucElemSize(), 3)
RESULT

CHECK((TopHatFilter& operator=(const TopHatFilter& t)))
    TopHatFilter tophat;
    tophat.setStrucElemSize(3);

    TopHatFilter tophat_copy;
    tophat_copy = tophat;
    TEST_EQUAL(tophat_copy.getStrucElemSize(), 3)
RESULT

CHECK((TopHatFilter& operator=(const TopHatFilter& t)))
RESULT


CHECK((template<typename InputPeakIterator, typename OutputPeakContainer  > void filter(InputPeakIterator first, InputPeakIterator last, OutputPeakContainer& baseline_filtered_container)))
    DPeakArray<1,DRawDataPoint<1> > raw_data;
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

    DPeakArray<1,DRawDataPoint<1> > tophat_data;
    TopHatFilter tophat;
    tophat.setStrucElemSize(3);
    tophat.filter(raw_data,tophat_data);

    DPeakArray<1,DRawDataPoint<1> >::ConstIterator it=tophat_data.begin();
    for (int i=0; i<24; ++i)
    {
      TEST_REAL_EQUAL(it->getIntensity(), 0)
    }
RESULT


CHECK((template<typename InputPeakIterator, typename OutputPeakContainer  > void filter(InputPeakIterator first, InputPeakIterator last, OutputPeakContainer& baseline_filtered_container)))
    DPeakArray<1,DRawDataPoint<1> > raw_data;
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

     DPeakArray<1,DRawDataPoint<1> > tophat_data;

     TopHatFilter tophat;
     tophat.setStrucElemSize(3);
     tophat.filter(raw_data.begin(),raw_data.end(),tophat_data);

     DPeakArray<1,DRawDataPoint<1> >::ConstIterator it=tophat_data.begin();
     for (int i=0; i < 8; ++i)
       {
         TEST_REAL_EQUAL(it->getIntensity(), 0)
       }
RESULT

CHECK((template<typename InputPeakType, typename OutputPeakType > void filterExperiment(const MSExperiment< InputPeakType >& ms_exp_raw, MSExperiment<OutputPeakType>& ms_exp_filtered)))
  MSExperiment< DRawDataPoint<1> > ms_exp_raw;
  MSExperiment< DRawDataPoint<1> > ms_exp_filtered;

  DPeakArray<2,DRawDataPoint<2> > raw_data;
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

  TopHatFilter tophat;
  tophat.setStrucElemSize(3);
  
  tophat.filterExperiment(ms_exp_raw, ms_exp_filtered);

  DPeakArray<2> dpeak_arra_filtered;
  ms_exp_filtered.get2DData(dpeak_arra_filtered);
  DPeakArray<2>::iterator it = dpeak_arra_filtered.begin();
  for (int i=0; i < 8; ++i)
  {
    TEST_REAL_EQUAL(it->getIntensity(), 0)
   }
RESULT

CHECK((template<typename InputSpectrumIterator, typename OutputPeakType > void filterExperiment(InputSpectrumIterator first, InputSpectrumIterator last, MSExperiment<OutputPeakType>& ms_exp_filtered)))
  MSExperiment<DRawDataPoint<1> > ms_exp_raw;
  MSExperiment<DRawDataPoint<1> > ms_exp_filtered;

  DPeakArray<2,DRawDataPoint<2> > raw_data;
  DPeakArray<2,DRawDataPoint<2> > filtered_data;

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

    TopHatFilter tophat;
    tophat.setStrucElemSize(3);
    tophat.filterExperiment(ms_exp_raw.begin(),ms_exp_raw.end(), ms_exp_filtered);

    ms_exp_filtered.get2DData(filtered_data);
    DPeakArray<2,DRawDataPoint<2> >::iterator it = filtered_data.begin();
    for (int i=0; i<8; ++i)
       {
         TEST_REAL_EQUAL(it->getIntensity(), 0)
       }
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

