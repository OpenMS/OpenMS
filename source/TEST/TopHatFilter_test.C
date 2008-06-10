// -*- Mode: C++; tab-width: 2; -*-
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

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FILTERING/BASELINE/TopHatFilter.h>

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak2D.h>

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

CHECK((virtual ~TopHatFilter()))
  delete tophat_ptr;
RESULT

CHECK((template <typename InputPeakContainer, typename OutputPeakContainer> void filter(const InputPeakContainer &input_peak_container, OutputPeakContainer &baseline_filtered_container)))
    MSSpectrum<Peak1D > raw_data;
    int i;
    for (i=0; i < 24; ++i)
    {
      Peak1D p;
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

    MSSpectrum<Peak1D > tophat_data;
    TopHatFilter tophat;
    tophat.filter(raw_data,tophat_data);

    DPeakArray<Peak1D >::ConstIterator it=tophat_data.begin();
    for (int i=0; i<24; ++i)
    {
      TEST_REAL_EQUAL(it->getIntensity(), 0)
    }
RESULT


CHECK((template<typename InputPeakIterator, typename OutputPeakContainer  > void filter(InputPeakIterator first, InputPeakIterator last, OutputPeakContainer& baseline_filtered_container)))
    DPeakArray<Peak1D > raw_data;
    int i;
    for (i=0; i<8; ++i)
    {
      Peak1D p;
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

     DPeakArray<Peak1D > tophat_data;

     TopHatFilter tophat;
     tophat.filter(raw_data.begin(),raw_data.end(),tophat_data);

     DPeakArray<Peak1D >::ConstIterator it=tophat_data.begin();
     for (int i=0; i < 8; ++i)
       {
         TEST_REAL_EQUAL(it->getIntensity(), 0)
       }
RESULT

CHECK((template <typename InputSpectrumIterator, typename OutputPeakType, typename OutputAllocType> void filterExperiment(InputSpectrumIterator first, InputSpectrumIterator last, MSExperiment< OutputPeakType, OutputAllocType > &ms_exp_filtered)))
  MSExperiment<Peak1D > ms_exp_raw;
  MSExperiment<Peak1D > ms_exp_filtered;

  DPeakArray<Peak2D > raw_data;
  DPeakArray<Peak2D > filtered_data;

    int i;
    for (i=0; i<8; ++i)
    {
      Peak2D p;
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
    tophat.filterExperiment(ms_exp_raw.begin(),ms_exp_raw.end(), ms_exp_filtered);

    ms_exp_filtered.get2DData(filtered_data);
    DPeakArray<Peak2D >::iterator it = filtered_data.begin();
    for (int i=0; i<8; ++i)
       {
         TEST_REAL_EQUAL(it->getIntensity(), 0)
       }
RESULT



CHECK((template <typename InputPeakType, typename InputAllocType, typename OutputPeakType, typename OutputAllocType> void filterExperiment(const MSExperiment< InputPeakType, InputAllocType > &ms_exp_raw, MSExperiment< OutputPeakType, OutputAllocType > &ms_exp_filtered)))
  MSExperiment< Peak1D > ms_exp_raw;
  MSExperiment< Peak1D > ms_exp_filtered;

  DPeakArray<Peak2D > raw_data;
  int i;
  for (i=0; i < 8; ++i)
  {
    Peak2D p;
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
  tophat.filterExperiment(ms_exp_raw, ms_exp_filtered);

  DPeakArray<Peak2D> dpeak_arra_filtered;
  ms_exp_filtered.get2DData(dpeak_arra_filtered);
  DPeakArray<Peak2D>::iterator it = dpeak_arra_filtered.begin();
  for (int i=0; i < 8; ++i)
  {
    TEST_REAL_EQUAL(it->getIntensity(), 0)
   }
RESULT

CHECK((template <typename InputSpectrumIterator, typename OutputPeakType, typename OutputAllocType> void filterExperiment(InputSpectrumIterator first, InputSpectrumIterator last, MSExperiment< OutputPeakType, OutputAllocType > &ms_exp_filtered)))
  MSExperiment<Peak1D > ms_exp_raw;
  MSExperiment<Peak1D > ms_exp_filtered;

  DPeakArray<Peak2D > raw_data;
  DPeakArray<Peak2D > filtered_data;

    int i;
    for (i=0; i<8; ++i)
    {
      Peak2D p;
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
    tophat.filterExperiment(ms_exp_raw.begin(),ms_exp_raw.end(), ms_exp_filtered);

    ms_exp_filtered.get2DData(filtered_data);
    DPeakArray<Peak2D >::iterator it = filtered_data.begin();
    for (int i=0; i<8; ++i)
       {
         TEST_REAL_EQUAL(it->getIntensity(), 0)
       }
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

