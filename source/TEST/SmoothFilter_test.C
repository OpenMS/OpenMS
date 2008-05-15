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
// $Maintainer: Eva Lange  $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FILTERING/SMOOTHING/SmoothFilter.h>

///////////////////////////

START_TEST(SmoothFilter<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

class SmoothFilterDummy
  : public SmoothFilter
{
	public:
	  const std::vector<DoubleReal>& getCoeffs() const
	  {
	    return coeffs_;
	  }
	  std::vector<DoubleReal>& getCoeffs()
	  {
	    return coeffs_;
	  }
	  void setCoeffs(std::vector<DoubleReal>& coeffs)
	  {
	    coeffs_ = coeffs;
	  }
};

SmoothFilter* dsmooth_ptr = 0;
CHECK((SmoothFilter()))
  dsmooth_ptr = new SmoothFilter;
  TEST_NOT_EQUAL(dsmooth_ptr, 0)
RESULT

CHECK((~SmoothFilter()))
  delete dsmooth_ptr;
RESULT

CHECK((template<typename InputPeakIterator, typename OutputPeakContainer  > void filter(InputPeakIterator first, InputPeakIterator last, OutputPeakContainer& smoothed_data_container)))
  MSSpectrum<RawDataPoint1D> raw;
  raw.resize(5);
  MSSpectrum<RawDataPoint1D> filtered;

  MSSpectrum<RawDataPoint1D>::Iterator it=raw.begin();
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

  std::vector<double> coeffs(4);
  coeffs[0]=0;
  coeffs[1]=1;
  coeffs[2]=1;
  coeffs[3]=0;

  SmoothFilterDummy smooth;
  smooth.setCoeffs(coeffs);
  smooth.filter(raw.begin(),raw.end(),filtered);

  it=filtered.begin();
  TEST_REAL_EQUAL(it->getIntensity(),0)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),0)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),0)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),1)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),1)
RESULT

CHECK((template<typename InputSpectrumIterator, typename OutputPeakType > void filterExperiment(InputSpectrumIterator first, InputSpectrumIterator last, MSExperiment<OutputPeakType>& ms_exp_filtered)))
	MSExperiment< RawDataPoint1D > raw_exp;
	MSExperiment< RawDataPoint1D > filtered_exp;
	MSSpectrum< RawDataPoint1D > raw_spectrum;
	raw_spectrum.resize(5);
	
  
  MSSpectrum< RawDataPoint1D >::iterator it=raw_spectrum.begin();
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

  std::vector<double> coeffs(4);
  coeffs[0]=0;
  coeffs[1]=1;
  coeffs[2]=1;
  coeffs[3]=0;

  SmoothFilterDummy smooth;
  smooth.setCoeffs(coeffs);
  raw_exp.resize(1);
  raw_exp[0] = raw_spectrum;
  smooth.filterExperiment(raw_exp.begin(),raw_exp.end(),filtered_exp);


	MSExperiment< RawDataPoint1D >::SpectrumType::iterator it2 = filtered_exp[0].begin();
  TEST_REAL_EQUAL(it2->getIntensity(),0)
  ++it2;
  TEST_REAL_EQUAL(it2->getIntensity(),0)
  ++it2;
  TEST_REAL_EQUAL(it2->getIntensity(),0)
  ++it2;
  TEST_REAL_EQUAL(it2->getIntensity(),1)
  ++it2;
  TEST_REAL_EQUAL(it2->getIntensity(),1)
RESULT

CHECK((template<typename InputPeakType, typename OutputPeakType > void filterExperiment(const MSExperiment< InputPeakType >& ms_exp_raw, MSExperiment<OutputPeakType>& ms_exp_filtered)))
	MSExperiment< RawDataPoint1D > raw_exp;
	MSExperiment< RawDataPoint1D > filtered_exp;
	MSSpectrum< RawDataPoint1D > raw_spectrum;
	raw_spectrum.resize(5);
	
  
  MSSpectrum< RawDataPoint1D >::iterator it=raw_spectrum.begin();
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

  std::vector<double> coeffs(4);
  coeffs[0]=0;
  coeffs[1]=1;
  coeffs[2]=1;
  coeffs[3]=0;

  SmoothFilterDummy smooth;
  smooth.setCoeffs(coeffs);
  raw_exp.resize(1);
  raw_exp[0] = raw_spectrum;
  smooth.filterExperiment(raw_exp ,filtered_exp);


	MSExperiment< RawDataPoint1D >::SpectrumType::iterator it2 = filtered_exp[0].begin();
  TEST_REAL_EQUAL(it2->getIntensity(),0)
  ++it2;
  TEST_REAL_EQUAL(it2->getIntensity(),0)
  ++it2;
  TEST_REAL_EQUAL(it2->getIntensity(),0)
  ++it2;
  TEST_REAL_EQUAL(it2->getIntensity(),1)
  ++it2;
  TEST_REAL_EQUAL(it2->getIntensity(),1)
RESULT

CHECK((template<typename InputPeakContainer, typename OutputPeakContainer > void filter(const InputPeakContainer& input_peak_container, OutputPeakContainer& smoothed_data_container)))
 MSSpectrum<RawDataPoint1D> raw;
  raw.resize(5);
  MSSpectrum<RawDataPoint1D> filtered;
  filtered.resize(5);

  MSSpectrum<RawDataPoint1D>::Iterator it=raw.begin();
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

  std::vector<double> coeffs(4);
  coeffs[0]=0;
  coeffs[1]=1;
  coeffs[2]=1;
  coeffs[3]=0;

  SmoothFilterDummy smooth;
  smooth.setCoeffs(coeffs);
  smooth.filter(raw,filtered);

  it=filtered.begin();
  TEST_REAL_EQUAL(it->getIntensity(),0)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),0)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),0)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),1)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),1)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

