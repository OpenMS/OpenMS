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

#include <OpenMS/FILTERING/SMOOTHING/SmoothFilter.h>

///////////////////////////

START_TEST(SmoothFilter<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;


typedef DPeakArray<1,DRawDataPoint<1> > RawDataArray1D;
typedef RawDataArray1D::Iterator RawDataIterator1D;
typedef RawDataArray1D::ConstIterator RawDataConstIterator1D;


SmoothFilter* dsmooth_ptr = 0;
CHECK((SmoothFilter()))
  dsmooth_ptr = new SmoothFilter;
  TEST_NOT_EQUAL(dsmooth_ptr, 0)
RESULT

CHECK((~SmoothFilter()))
  delete dsmooth_ptr;
RESULT


CHECK((SmoothFilter& operator=(const SmoothFilter& fir)))
  std::vector<double> coeffs(1);
  coeffs[0]=1.23423;
  SmoothFilter smooth;
  smooth.setCoeffs(coeffs);
  SmoothFilter smooth_copy;

  smooth_copy=smooth;

  TEST_EQUAL(smooth_copy.getCoeffs()[0], 1.23423)
RESULT

CHECK((SmoothFilter(const SmoothFilter& fir)))
  std::vector<double> coeffs(1);
  coeffs[0]=1.23423;
  SmoothFilter smooth;
  smooth.setCoeffs(coeffs);
  
  SmoothFilter smooth_copy(smooth);
  TEST_EQUAL(smooth_copy.getCoeffs()[0], 1.23423)
RESULT

CHECK((const std::vector<double>& getCoeffs() const))
  const SmoothFilter smooth;
  TEST_EQUAL(smooth.getCoeffs().size(),0)
RESULT

CHECK((std::vector<double>& getCoeffs()))
  std::vector<double> coeffs(1);
  coeffs[0]=1;
  SmoothFilter smooth;
  smooth.setCoeffs(coeffs);
  TEST_EQUAL(smooth.getCoeffs()[0],1);
RESULT


CHECK((template<typename InputPeakIterator, typename OutputPeakContainer  > void filter(InputPeakIterator first, InputPeakIterator last, OutputPeakContainer& smoothed_data_container)))
  RawDataArray1D raw(5);
  RawDataArray1D filtered;

  RawDataIterator1D it=raw.begin();
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

  SmoothFilter smooth;
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
	MSExperiment< DRawDataPoint<1> > raw_exp;
	MSExperiment< DRawDataPoint<1> > filtered_exp;
	MSSpectrum< DRawDataPoint<1> > raw_spectrum;
	raw_spectrum.resize(5);
	
  
  MSSpectrum< DRawDataPoint<1> >::iterator it=raw_spectrum.begin();
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

  SmoothFilter smooth;
  smooth.setCoeffs(coeffs);
  raw_exp.resize(1);
  raw_exp[0] = raw_spectrum;
  smooth.filterExperiment(raw_exp.begin(),raw_exp.end(),filtered_exp);


	it = filtered_exp[0].begin();
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


CHECK((void setCoeffs(std::vector<double>& coeffs)))
  std::vector<double> coeffs(3);
  coeffs[0]=1.2;
  coeffs[1]=2.4;
  coeffs[2]= 1.4;
  SmoothFilter smooth;
  smooth.setCoeffs(coeffs);
  TEST_EQUAL(smooth.getCoeffs()[0],1.2)
  TEST_EQUAL(smooth.getCoeffs()[1],2.4)
  TEST_EQUAL(smooth.getCoeffs()[2],1.4)
RESULT

CHECK((template<typename InputPeakType, typename OutputPeakType > void filterExperiment(const MSExperiment< InputPeakType >& ms_exp_raw, MSExperiment<OutputPeakType>& ms_exp_filtered)))
	MSExperiment< DRawDataPoint<1> > raw_exp;
	MSExperiment< DRawDataPoint<1> > filtered_exp;
	MSSpectrum< DRawDataPoint<1> > raw_spectrum;
	raw_spectrum.resize(5);
	
  
  MSSpectrum< DRawDataPoint<1> >::iterator it=raw_spectrum.begin();
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

  SmoothFilter smooth;
  smooth.setCoeffs(coeffs);
  raw_exp.resize(1);
  raw_exp[0] = raw_spectrum;
  smooth.filterExperiment(raw_exp ,filtered_exp);


	it = filtered_exp[0].begin();
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


CHECK((void setCoeffs(std::vector<double>& coeffs)))
  std::vector<double> coeffs(3);
  coeffs[0]=1.2;
  coeffs[1]=2.4;
  coeffs[2]= 1.4;
  SmoothFilter smooth;
  smooth.setCoeffs(coeffs);
  TEST_EQUAL(smooth.getCoeffs()[0],1.2)
  TEST_EQUAL(smooth.getCoeffs()[1],2.4)
  TEST_EQUAL(smooth.getCoeffs()[2],1.4)
RESULT

CHECK((template<typename InputPeakContainer, typename OutputPeakContainer > void filter(const InputPeakContainer& input_peak_container, OutputPeakContainer& smoothed_data_container)))
 RawDataArray1D raw(5);
  RawDataArray1D filtered(5);

  RawDataIterator1D it=raw.begin();
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

  SmoothFilter smooth;
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

CHECK((void setCoeffs(std::vector<double>& coeffs)))
  std::vector<double> coeffs(3);
  coeffs[0]=1.2;
  coeffs[1]=2.4;
  coeffs[2]= 1.4;
  SmoothFilter smooth;
  smooth.setCoeffs(coeffs);
  TEST_EQUAL(smooth.getCoeffs()[0],1.2)
  TEST_EQUAL(smooth.getCoeffs()[1],2.4)
  TEST_EQUAL(smooth.getCoeffs()[2],1.4)
RESULT



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

