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
// $Maintainer: Eva Lange  $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>
#include <OpenMS/KERNEL/RawDataPoint2D.h>

///////////////////////////

START_TEST(GaussFilter<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;


typedef DPeakArray<RawDataPoint1D > RawDataArray1D;
typedef RawDataArray1D::Iterator RawDataIterator1D;
typedef RawDataArray1D::ConstIterator RawDataConstIterator1D;


GaussFilter* dgauss_ptr = 0;
CHECK((GaussFilter()))
  dgauss_ptr = new GaussFilter;
  TEST_NOT_EQUAL(dgauss_ptr, 0)
RESULT

CHECK((virtual ~GaussFilter()))
    delete dgauss_ptr;
RESULT

CHECK((DoubleReal getSigma() const))
  const GaussFilter gaussian;

  TEST_REAL_EQUAL(gaussian.getSigma(),.1);
RESULT

CHECK((DoubleReal getSpacing() const))
  const GaussFilter gaussian;

  TEST_REAL_EQUAL(gaussian.getSpacing(),0.01);
RESULT

CHECK((DoubleReal getKernelWidth() const))
  const GaussFilter gaussian;

  TEST_REAL_EQUAL(gaussian.getKernelWidth(),.8);
RESULT

CHECK((void init(DoubleReal sigma, DoubleReal spacing)))
  GaussFilter gaussian;
  gaussian.init(0.2,0.001);

  TEST_REAL_EQUAL(gaussian.getSpacing(),0.001);
  TEST_REAL_EQUAL(gaussian.getSigma(),0.2);
  TEST_REAL_EQUAL(gaussian.getKernelWidth(),1.6);
RESULT

CHECK((void setKernelWidth(DoubleReal kernel_width) throw (Exception::InvalidValue)))
  GaussFilter gaussian;
  gaussian.setKernelWidth(1.6);

  TEST_REAL_EQUAL(gaussian.getKernelWidth(),1.6);
RESULT

CHECK((void setSigma(DoubleReal sigma)))
  GaussFilter gauss;
  gauss.setSigma(2.434);
  
  TEST_REAL_EQUAL(gauss.getSigma(), 2.434)
RESULT

CHECK((void setSpacing(DoubleReal spacing)))
  GaussFilter gaussian;
  gaussian.setSpacing(0.0001);

  TEST_REAL_EQUAL(gaussian.getSpacing(),0.0001);
RESULT

CHECK((template<typename InputPeakIterator, typename OutputPeakContainer  > void filter(InputPeakIterator first, InputPeakIterator last, OutputPeakContainer& smoothed_data_container)))
  RawDataArray1D raw(5);
  RawDataArray1D filtered;

  RawDataIterator1D it = raw.begin();
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

  GaussFilter gauss;
  gauss.filter(raw.begin(),raw.end(),filtered);
  it=filtered.begin();
  TEST_REAL_EQUAL(it->getIntensity(),0.)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),0.)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),0)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),0)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),0)
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

  GaussFilter gauss;
  raw_exp.resize(1);
  raw_exp[0] = raw_spectrum;
  gauss.filterExperiment(raw_exp.begin(),raw_exp.end(),filtered_exp);

  it = filtered_exp[0].begin();
  TEST_REAL_EQUAL(it->getIntensity(),0.)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),0.)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),0)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),0)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),0)
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

  GaussFilter gauss;
  raw_exp.resize(1);
  raw_exp[0] = raw_spectrum;
  gauss.filterExperiment(raw_exp,filtered_exp);

  it = filtered_exp[0].begin();
  TEST_REAL_EQUAL(it->getIntensity(),0.)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),0.)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),0)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),0)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),0)
RESULT

CHECK((template<typename InputPeakContainer, typename OutputPeakContainer > void filter(const InputPeakContainer& input_peak_container, OutputPeakContainer& smoothed_data_container)))
  RawDataArray1D raw(5);
  RawDataArray1D filtered;

  RawDataIterator1D it = raw.begin();
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

  GaussFilter gauss;
  gauss.filter(raw,filtered);
  it=filtered.begin();
  TEST_REAL_EQUAL(it->getIntensity(),0.)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),0.)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),0)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),0)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),0)
RESULT

CHECK((template <typename InputSpectrumIterator, typename OutputPeakType> void filterExperiment(InputSpectrumIterator first, InputSpectrumIterator last, MSExperimentExtern< OutputPeakType > &ms_exp_filtered)))
	MSExperiment< RawDataPoint1D > raw_exp;
	MSExperiment< RawDataPoint1D > filtered_exp;
	
	DPeakArray< RawDataPoint2D > raw_data(5);
	DPeakArray< RawDataPoint2D >::iterator it_2 = raw_data.begin();
	for (int i=0; i<5; ++i, ++it_2)
  {
    if (i==2)
    {
      it_2->setIntensity(1);
    }
    else
    {
      it_2->setIntensity(0);
    }
  }
  
  raw_exp.set2DData(raw_data);
  GaussFilter gauss;
  gauss.filterExperiment(raw_exp.begin(),raw_exp.end(),filtered_exp);

  MSSpectrum< RawDataPoint1D >::iterator it = filtered_exp[0].begin();;
  it = filtered_exp[0].begin();
  TEST_REAL_EQUAL(it->getIntensity(),0.)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),0.)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),0)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),0)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),0)
   
//  MSExperiment< RawDataPoint1D > raw_exp;
//	MSExperiment< RawDataPoint1D > filtered_exp;  
//	MSSpectrum< RawDataPoint1D > raw_spectrum;
//	raw_spectrum.resize(5);
//	
//  
//  MSSpectrum< RawDataPoint1D >::iterator it=raw_spectrum.begin();
//  for (int i=0; i<5; ++i, ++it)
//  {
//    if (i==2)
//    {
//      it->setIntensity(1);
//    }
//    else
//    {
//      it->setIntensity(0);
//    }
//  }
//
//  GaussFilter gauss;
//  raw_exp.resize(1);
//  raw_exp[0] = raw_spectrum;
//  gauss.filterExperiment(raw_exp.begin(),raw_exp.end(),filtered_exp);
//
//  it = filtered_exp[0].begin();
//  TEST_REAL_EQUAL(it->getIntensity(),0.)
//  ++it;
//  TEST_REAL_EQUAL(it->getIntensity(),0.)
//  ++it;
//  TEST_REAL_EQUAL(it->getIntensity(),0)
//  ++it;
//  TEST_REAL_EQUAL(it->getIntensity(),0)
//  ++it;
//  TEST_REAL_EQUAL(it->getIntensity(),0)
RESULT

CHECK((template<typename InputPeakType, typename OutputPeakType > void filterExperiment(const MSExperimentExtern< InputPeakType >& ms_exp_raw, MSExperimentExtern<OutputPeakType>& ms_exp_filtered)))
	MSExperiment< RawDataPoint1D > raw_exp;
	MSExperiment< RawDataPoint1D > filtered_exp;
	
	DPeakArray< RawDataPoint2D > raw_data(5);
	DPeakArray< RawDataPoint2D >::iterator it_2 = raw_data.begin();
	for (int i=0; i<5; ++i, ++it_2)
  {
    if (i==2)
    {
      it_2->setIntensity(1);
    }
    else
    {
      it_2->setIntensity(0);
    }
  }
  
  raw_exp.set2DData(raw_data);
  GaussFilter gauss;
  gauss.filterExperiment(raw_exp,filtered_exp);

  MSSpectrum< RawDataPoint1D >::iterator it = filtered_exp[0].begin();
  it = filtered_exp[0].begin();
  TEST_REAL_EQUAL(it->getIntensity(),0.)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),0.)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),0)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),0)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),0)
//	MSExperimentExtern< RawDataPoint1D > raw_exp;
//	MSExperimentExtern< RawDataPoint1D > filtered_exp;
//	MSSpectrum< RawDataPoint1D > raw_spectrum;
//	raw_spectrum.resize(5);
//	
//  
//  MSSpectrum< RawDataPoint1D >::iterator it=raw_spectrum.begin();
//  for (int i=0; i<5; ++i, ++it)
//  {
//    if (i==2)
//    {
//      it->setIntensity(1);
//    }
//    else
//    {
//      it->setIntensity(0);
//    }
//  }
//
//  GaussFilter gauss;
//  raw_exp.resize(1);
//  raw_exp[0] = raw_spectrum;
//  gauss.filterExperiment(raw_exp,filtered_exp);
//
//  it = filtered_exp[0].begin();
//  TEST_REAL_EQUAL(it->getIntensity(),0.)
//  ++it;
//  TEST_REAL_EQUAL(it->getIntensity(),0.)
//  ++it;
//  TEST_REAL_EQUAL(it->getIntensity(),0)
//  ++it;
//  TEST_REAL_EQUAL(it->getIntensity(),0)
//  ++it;
//  TEST_REAL_EQUAL(it->getIntensity(),0)
RESULT
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
