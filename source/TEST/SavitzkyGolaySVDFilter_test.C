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

#include <OpenMS/FILTERING/SMOOTHING/SavitzkyGolaySVDFilter.h>

///////////////////////////

START_TEST(SavitzkyGolaySVDFilter<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

typedef DPeakArray<1,DRawDataPoint<1> > RawDataArray1D;
typedef RawDataArray1D::Iterator RawDataIterator1D;
typedef RawDataArray1D::ConstIterator RawDataConstIterator1D;


SavitzkyGolaySVDFilter* dsg_ptr = 0;
CHECK((SavitzkyGolaySVDFilter()))
  dsg_ptr = new SavitzkyGolaySVDFilter;
  TEST_NOT_EQUAL(dsg_ptr, 0)
RESULT

CHECK((~SavitzkyGolaySVDFilter()))
  delete dsg_ptr;
RESULT

CHECK((SavitzkyGolaySVDFilter& operator=(const SavitzkyGolaySVDFilter& s)))
  SavitzkyGolaySVDFilter sgolay;
  sgolay.setOrder(4);
  sgolay.setWindowSize(5);

  SavitzkyGolaySVDFilter sgolay_copy;
  sgolay_copy = sgolay;
  TEST_REAL_EQUAL(sgolay_copy.getOrder(),sgolay.getOrder())
  TEST_EQUAL(sgolay_copy.getWindowSize(),sgolay.getWindowSize())
RESULT

CHECK((SavitzkyGolaySVDFilter(const SavitzkyGolaySVDFilter& s)))
   SavitzkyGolaySVDFilter sgolay;
  sgolay.setOrder(4);
  sgolay.setWindowSize(5);

  SavitzkyGolaySVDFilter sgolay_copy(sgolay);
  TEST_REAL_EQUAL(sgolay_copy.getOrder(),sgolay.getOrder())
  TEST_EQUAL(sgolay_copy.getWindowSize(),sgolay.getWindowSize())
RESULT

CHECK((const unsigned int& getOrder() const))
  SavitzkyGolaySVDFilter sgolay;

  TEST_EQUAL(sgolay.getOrder(),4);
RESULT

CHECK((const unsigned int& getWindowSize() const))
  SavitzkyGolaySVDFilter sgolay;

  TEST_EQUAL(sgolay.getWindowSize(),17);
RESULT

CHECK((const unsigned int& getOrder() const))
  SavitzkyGolaySVDFilter sgolay;

  TEST_EQUAL(sgolay.getOrder(),4);
RESULT

CHECK((void setOrder(const unsigned int& order)))
  SavitzkyGolaySVDFilter sgolay;
  sgolay.setOrder(3);

  TEST_EQUAL(sgolay.getOrder(),3);
RESULT

CHECK((void setParam(Param param) throw(Exception::InvalidValue)))
  Param p;
  p.setValue("polynomial_order",2);
  p.setValue("frame_length",3);
  SavitzkyGolaySVDFilter sgolay;
  sgolay.setParam(p);

  TEST_REAL_EQUAL(sgolay.getOrder(),2);
  TEST_REAL_EQUAL(sgolay.getWindowSize(),3);
RESULT

CHECK((void setWindowSize(const unsigned int& frame_size)))
  SavitzkyGolaySVDFilter sgolay;
  sgolay.setWindowSize(7);

  TEST_EQUAL(sgolay.getWindowSize(),7);
RESULT

CHECK((template< typename InputPeakIterator, typename OutputPeakContainer  > void filter(InputPeakIterator first, InputPeakIterator last, OutputPeakContainer& smoothed_data_container) throw(Exception::InvalidSize)))
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

  SavitzkyGolaySVDFilter sgolay;
  sgolay.setOrder(2);
  sgolay.setWindowSize(3);
  sgolay.filter(raw.begin(),raw.end(),filtered);
  it=filtered.begin();
  TEST_REAL_EQUAL(it->getIntensity(),0.)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),0.)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),sgolay.getCoeffs()[2])
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),sgolay.getCoeffs()[1])
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),sgolay.getCoeffs()[0])
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

  SavitzkyGolaySVDFilter sgolay;
  sgolay.setOrder(2);
  sgolay.setWindowSize(3);
  raw_exp.resize(1);
  raw_exp[0] = raw_spectrum;
  sgolay.filterExperiment(raw_exp.begin(),raw_exp.end(),filtered_exp);

  it = filtered_exp[0].begin();
  TEST_REAL_EQUAL(it->getIntensity(),0.)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),0.)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),sgolay.getCoeffs()[2])
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),sgolay.getCoeffs()[1])
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),sgolay.getCoeffs()[0])
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

  SavitzkyGolaySVDFilter sgolay;
  sgolay.setOrder(2);
  sgolay.setWindowSize(3);
  raw_exp.resize(1);
  raw_exp[0] = raw_spectrum;
  sgolay.filterExperiment(raw_exp,filtered_exp);

  it = filtered_exp[0].begin();
  TEST_REAL_EQUAL(it->getIntensity(),0.)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),0.)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),sgolay.getCoeffs()[2])
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),sgolay.getCoeffs()[1])
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),sgolay.getCoeffs()[0])
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

  SavitzkyGolaySVDFilter sgolay;
  sgolay.setOrder(2);
  sgolay.setWindowSize(3);
  sgolay.filter(raw,filtered);
  it=filtered.begin();
  TEST_REAL_EQUAL(it->getIntensity(),0.)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),0.)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),sgolay.getCoeffs()[2])
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),sgolay.getCoeffs()[1])
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),sgolay.getCoeffs()[0])
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
