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
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  021-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Eva Lange  $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FILTERING/SMOOTHING/SavitzkyGolayFilter.h>
#include <OpenMS/KERNEL/Peak2D.h>

///////////////////////////

START_TEST(SavitzkyGolayFilter<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

SavitzkyGolayFilter* dsg_ptr = 0;
CHECK((SavitzkyGolayFilter()))
  dsg_ptr = new SavitzkyGolayFilter;
  TEST_NOT_EQUAL(dsg_ptr, 0)
RESULT

CHECK((virtual ~SavitzkyGolayFilter()))
  delete dsg_ptr;
RESULT

Param param;
param.setValue("polynomial_order",2);
param.setValue("frame_length",3);

CHECK((template <typename InputPeakIterator, typename OutputPeakContainer> void filter(InputPeakIterator first, InputPeakIterator last, OutputPeakContainer &smoothed_data_container)))
  MSSpectrum<Peak1D> raw;
  raw.resize(5);
  MSSpectrum<Peak1D> filtered;

  MSSpectrum<Peak1D>::Iterator it = raw.begin();
  for (int i=0; i<5; ++i, ++it)
  {
    if (i==2)
    {
      it->setIntensity(1);
    }
    else
    {
      it->setIntensity(0);
    }
  }

  SavitzkyGolayFilter sgolay;
	sgolay.setParameters(param);
  sgolay.filter(raw.begin(),raw.end(),filtered);
  it=filtered.begin();
  TEST_REAL_EQUAL(it->getIntensity(),0.)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),0.)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),1)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),0)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),5.55112e-17)
RESULT 

CHECK((template<typename InputSpectrumIterator, typename OutputPeakType > void filterExperiment(InputSpectrumIterator first, InputSpectrumIterator last, MSExperiment<OutputPeakType>& ms_exp_filtered)))
	MSExperiment< Peak1D > raw_exp;
	MSExperiment< Peak1D > filtered_exp;
	MSSpectrum< Peak1D > raw_spectrum;
	raw_spectrum.resize(5);
	
  
  MSSpectrum< Peak1D >::iterator it=raw_spectrum.begin();
  for (int i=0; i<5; ++i, ++it)
  {
    if (i==2)
    {
      it->setIntensity(1);
    }
    else
    {
      it->setIntensity(0);
    }
  }

  SavitzkyGolayFilter sgolay;
	sgolay.setParameters(param);
  raw_exp.resize(1);
  raw_exp[0] = raw_spectrum;
  sgolay.filterExperiment(raw_exp.begin(),raw_exp.end(),filtered_exp);

  MSExperiment< Peak1D >::SpectrumType::iterator it2 = filtered_exp[0].begin();
  TEST_REAL_EQUAL(it2->getIntensity(),0.)
  ++it2;
  TEST_REAL_EQUAL(it2->getIntensity(),0.)
  ++it2;
  TEST_REAL_EQUAL(it2->getIntensity(),1)
  ++it2;
  TEST_REAL_EQUAL(it2->getIntensity(),0)
  ++it2;
  TEST_REAL_EQUAL(it2->getIntensity(),5.55112e-17)
RESULT

CHECK((template<typename InputPeakType, typename OutputPeakType > void filterExperiment(const MSExperiment< InputPeakType >& ms_exp_raw, MSExperiment<OutputPeakType>& ms_exp_filtered)))
	MSExperiment< Peak1D > raw_exp;
	MSExperiment< Peak1D > filtered_exp;
	MSSpectrum< Peak1D > raw_spectrum;
	raw_spectrum.resize(5);
	
  
  MSSpectrum< Peak1D >::iterator it=raw_spectrum.begin();
  for (int i=0; i<5; ++i, ++it)
  {
    if (i==2)
    {
      it->setIntensity(1);
    }
    else
    {
      it->setIntensity(0);
    }
  }

  SavitzkyGolayFilter sgolay;
	sgolay.setParameters(param);
  raw_exp.resize(1);
  raw_exp[0] = raw_spectrum;
  sgolay.filterExperiment(raw_exp,filtered_exp);

  MSExperiment< Peak1D >::SpectrumType::iterator it2 = filtered_exp[0].begin();
  TEST_REAL_EQUAL(it2->getIntensity(),0.)
  ++it2;
  TEST_REAL_EQUAL(it2->getIntensity(),0.)
  ++it2;
  TEST_REAL_EQUAL(it2->getIntensity(),1)
  ++it2;
  TEST_REAL_EQUAL(it2->getIntensity(),0)
  ++it2;
  TEST_REAL_EQUAL(it2->getIntensity(),5.55112e-17)
RESULT

CHECK((template <typename InputPeakContainer, typename OutputPeakContainer> void filter(const InputPeakContainer &input_peak_container, OutputPeakContainer &baseline_filtered_container)))
  MSSpectrum<Peak1D> raw;
  raw.resize(5);
  MSSpectrum<Peak1D> filtered;

  MSSpectrum<Peak1D>::Iterator it = raw.begin();
  for (int i=0; i<5; ++i, ++it)
  {
    if (i==2)
    {
      it->setIntensity(1);
    }
    else
    {
      it->setIntensity(0);
    }
  }

  SavitzkyGolayFilter sgolay;
	sgolay.setParameters(param);
  sgolay.filter(raw,filtered);
  it=filtered.begin();
  TEST_REAL_EQUAL(it->getIntensity(),0.)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),0.)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),1)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),0)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),5.55112e-17)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
