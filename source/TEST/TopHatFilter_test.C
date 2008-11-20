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
START_SECTION((TopHatFilter()))
  tophat_ptr = new TopHatFilter;
  TEST_NOT_EQUAL(tophat_ptr, 0) 
END_SECTION

START_SECTION((virtual ~TopHatFilter()))
  delete tophat_ptr;
END_SECTION


Param param;
param.setValue("struc_elem_length",3.0);

START_SECTION((template <typename InputPeakContainer, typename OutputPeakContainer> void filter(const InputPeakContainer &input_peak_container, OutputPeakContainer &baseline_filtered_container)))
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
    tophat.setParameters(param);
    tophat.filter(raw_data,tophat_data);

    std::vector<Peak1D >::const_iterator it=tophat_data.begin();
    for (int i=0; i<24; ++i)
    {
      TEST_REAL_SIMILAR(it->getIntensity(), 0)
    }
END_SECTION


START_SECTION((template<typename InputPeakIterator, typename OutputPeakContainer  > void filter(InputPeakIterator first, InputPeakIterator last, OutputPeakContainer& baseline_filtered_container)))
    std::vector<Peak1D > raw_data;
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

     std::vector<Peak1D > tophat_data;

     TopHatFilter tophat;
     tophat.setParameters(param);
     tophat.filter(raw_data.begin(),raw_data.end(),tophat_data);

     std::vector<Peak1D >::const_iterator it=tophat_data.begin();
     for (int i=0; i < 8; ++i)
       {
         TEST_REAL_SIMILAR(it->getIntensity(), 0)
       }
END_SECTION

START_SECTION((template <typename PeakType> void filterExperiment(MSExperiment<PeakType>& map)))
  MSExperiment<Peak1D> exp;
  exp.resize(4);
	
	Peak1D p;
	for (int i=0; i<8; ++i)
	{
		p.setMZ(i);
		p.setIntensity(0);
		if (i>1 && i<5)
		{
			p.setIntensity(1);
		}
		exp[0].push_back(p);
		exp[1].push_back(p);
	}
	exp[2].push_back(p);
	
	TopHatFilter tophat;
  tophat.setParameters(param);
	tophat.filterExperiment(exp);
	
	TEST_EQUAL(exp.size(),4);
	TEST_EQUAL(exp[0].size(),8);
	TEST_EQUAL(exp[1].size(),8);
	TEST_EQUAL(exp[2].size(),0);
	TEST_EQUAL(exp[3].size(),0);

	for (UInt i=0; i<exp[0].size(); ++i)
	{
		TEST_REAL_SIMILAR(exp[0][i].getIntensity(), 0.0)
		TEST_REAL_SIMILAR(exp[1][i].getIntensity(), 0.0)
	}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

