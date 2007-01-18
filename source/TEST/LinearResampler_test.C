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

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/KERNEL/DimensionDescription.h>

#include <OpenMS/FILTERING/TRANSFORMERS/LinearResampler.h>

#include <OpenMS/KERNEL/MSExperiment.h>

///////////////////////////

START_TEST(LinearResampler, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

LinearResampler* lr_ptr = 0;
CHECK((LinearResampler()))
  lr_ptr = new LinearResampler;
  TEST_NOT_EQUAL(lr_ptr,0);
RESULT

CHECK((~LinearResampler()))
  delete lr_ptr;
RESULT

CHECK((LinearResampler(const Param& parameters)))
  Param p;
  p.setValue("resampling_width",0.5);
  LinearResampler lr(p);

  TEST_REAL_EQUAL(lr.getSpacing(),0.5);
RESULT

CHECK((LinearResampler( LinearResampler const & lr )))
  Param p;
  p.setValue("resampling_width",0.5);
  LinearResampler tmp(p);

  LinearResampler tmp2(tmp);
  TEST_REAL_EQUAL(tmp2.getSpacing(),0.5);
  TEST_EQUAL(tmp2.getParam(),p);
RESULT

CHECK((LinearResampler& operator= (const LinearResampler& source)))
  Param p;
  p.setValue("resampling_width",0.5);
  LinearResampler tmp(p);

  LinearResampler tmp2;
  tmp2 = tmp;
  TEST_REAL_EQUAL(tmp2.getSpacing(),0.5);
  TEST_EQUAL(tmp2.getParam(),p);
RESULT

CHECK((template<typename InputSpectrumIterator, typename OutputPeakType > void rasterExperiment(InputSpectrumIterator first, InputSpectrumIterator last, MSExperiment<OutputPeakType>& ms_exp_filtered)))
  MSExperiment< DRawDataPoint <1> > raw;
  raw.resize(1);
  MSExperiment< DRawDataPoint <1> > resampled;
  MSSpectrum< DRawDataPoint<1> > spec;
  spec.getContainer().resize(5);
  spec.getContainer()[0].getPos() = 0;
  spec.getContainer()[0].getIntensity() = 3;
  spec.getContainer()[1].getPos() = 0.5;
  spec.getContainer()[1].getIntensity() = 6;
  spec.getContainer()[2].getPos() = 1.;
  spec.getContainer()[2].getIntensity() = 8;
  spec.getContainer()[3].getPos() = 1.6;
  spec.getContainer()[3].getIntensity() = 2;
  spec.getContainer()[4].getPos() = 1.8;
  spec.getContainer()[4].getIntensity() = 1;
  raw[0] = spec;

  LinearResampler lr;
  lr.setSpacing(0.5);
  lr.rasterExperiment(raw.begin(),raw.end(),resampled);

  double sum = 0.;
  MSSpectrum< DRawDataPoint<1> >::const_iterator it = resampled[0].begin();
  while(it != resampled[0].end())
  {
    sum += it->getIntensity();
    ++it;
  }

  TEST_REAL_EQUAL(sum, 20);
RESULT

CHECK((template<typename InputPeakType, typename OutputPeakType > void rasterExperiment(const MSExperiment< InputPeakType >& ms_exp_raw, MSExperiment<OutputPeakType>& ms_exp_filtered)))
  MSExperiment< DRawDataPoint <1> > raw;
  raw.resize(1);
  MSExperiment< DRawDataPoint <1> > resampled;
  MSSpectrum< DRawDataPoint<1> > spec;
  spec.getContainer().resize(5);
  spec.getContainer()[0].getPos() = 0;
  spec.getContainer()[0].getIntensity() = 3;
  spec.getContainer()[1].getPos() = 0.5;
  spec.getContainer()[1].getIntensity() = 6;
  spec.getContainer()[2].getPos() = 1.;
  spec.getContainer()[2].getIntensity() = 8;
  spec.getContainer()[3].getPos() = 1.6;
  spec.getContainer()[3].getIntensity() = 2;
  spec.getContainer()[4].getPos() = 1.8;
  spec.getContainer()[4].getIntensity() = 1;
  raw[0] = spec;

  LinearResampler lr;
  lr.setSpacing(0.5);
  lr.rasterExperiment(raw,resampled);

  double sum = 0.;
  MSSpectrum< DRawDataPoint<1> >::const_iterator it = resampled[0].begin();
  while(it != resampled[0].end())
  {
    sum += it->getIntensity();
    ++it;
  }

  TEST_REAL_EQUAL(sum, 20);
RESULT

CHECK((const Param& getParam() const))
  Param p;
  p.setValue("resampling_width",0.5);
  const LinearResampler tmp(p);

  TEST_EQUAL(tmp.getParam(),p);
RESULT

CHECK((const double& getSpacing() const))
  const LinearResampler tmp;

  TEST_EQUAL(tmp.getSpacing(),0.05);
RESULT

CHECK((double& getSpacing()))
  LinearResampler tmp;
  tmp.getSpacing() = 0.1;

  TEST_EQUAL(tmp.getSpacing(),0.1);
RESULT

CHECK((void setParam(const Param& param)))
  Param p;
  p.setValue("resampling_width",0.5);
  LinearResampler tmp;
  tmp.setParam(p);

  TEST_EQUAL(tmp.getParam(),p);
RESULT

CHECK((template< typename InputPeakIterator, typename OutputPeakContainer > void raster(InputPeakIterator first, InputPeakIterator last, OutputPeakContainer& resampled_peak_container)))
  MSSpectrum< DRawDataPoint<1> > spec;
  spec.getContainer().resize(5);
  spec.getContainer()[0].getPos() = 0;
  spec.getContainer()[0].getIntensity() = 3;
  spec.getContainer()[1].getPos() = 0.5;
  spec.getContainer()[1].getIntensity() = 6;
  spec.getContainer()[2].getPos() = 1.;
  spec.getContainer()[2].getIntensity() = 8;
  spec.getContainer()[3].getPos() = 1.6;
  spec.getContainer()[3].getIntensity() = 2;
  spec.getContainer()[4].getPos() = 1.8;
  spec.getContainer()[4].getIntensity() = 1;

  LinearResampler lr;
  lr.setSpacing(0.5);

  MSSpectrum< DRawDataPoint<1> > spec_resampled;
  lr.raster(spec.begin(),spec.end(),spec_resampled);

  double sum = 0.;
  MSSpectrum< DRawDataPoint<1> >::const_iterator it = spec_resampled.begin();
  while(it != spec_resampled.end())
  {
    sum += it->getIntensity();
    ++it;
  }

  TEST_REAL_EQUAL(sum, 20);
RESULT

CHECK((void setSpacing(const double& spacing)))
  LinearResampler tmp;
  tmp.setSpacing(0.1);

  TEST_EQUAL(tmp.getSpacing(),0.1);
RESULT

CHECK((template<typename InputPeakContainer, typename OutputPeakContainer > void raster(const InputPeakContainer& input_peak_container, OutputPeakContainer& baseline_filtered_container)))
  MSSpectrum< DRawDataPoint<1> > spec;
  spec.getContainer().resize(5);
  spec.getContainer()[0].getPos() = 0;
  spec.getContainer()[0].getIntensity() = 3;
  spec.getContainer()[1].getPos() = 0.5;
  spec.getContainer()[1].getIntensity() = 6;
  spec.getContainer()[2].getPos() = 1.;
  spec.getContainer()[2].getIntensity() = 8;
  spec.getContainer()[3].getPos() = 1.6;
  spec.getContainer()[3].getIntensity() = 2;
  spec.getContainer()[4].getPos() = 1.8;
  spec.getContainer()[4].getIntensity() = 1;

  LinearResampler lr;
  lr.setSpacing(0.5);

  MSSpectrum< DRawDataPoint<1> > spec_resampled;
  lr.raster(spec,spec_resampled);

  double sum = 0.;
  MSSpectrum< DRawDataPoint<1> >::const_iterator it = spec_resampled.begin();
  while(it != spec_resampled.end())
  {
    sum += it->getIntensity();
    ++it;
  }

  TEST_REAL_EQUAL(sum, 20);
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
