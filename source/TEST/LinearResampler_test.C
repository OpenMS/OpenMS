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

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/KERNEL/DimensionDescription.h>

#include <OpenMS/FILTERING/TRANSFORMERS/LinearResampler.h>

#include <OpenMS/KERNEL/MSExperiment.h>

///////////////////////////

START_TEST(LinearResampler< PeakType >, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

LinearResampler< DRawDataPoint<1> >* lr_ptr = 0;
CHECK(LinearResampler())
  lr_ptr = new LinearResampler< DRawDataPoint<1> >;
  TEST_NOT_EQUAL(lr_ptr,0);
RESULT

CHECK(~LinearResampler())
  delete lr_ptr;
RESULT

CHECK(LinearResampler( const Param& parameters))
  Param p;
  p.setValue("ResamplingWidth",0.5);
  LinearResampler< DRawDataPoint<1> > lr(p);

  TEST_REAL_EQUAL(lr.getSpacing(),0.5);
RESULT

CHECK(LinearResampler( LinearResampler const & lr ))
  Param p;
  p.setValue("ResamplingWidth",0.5);
  LinearResampler< DRawDataPoint<1> > tmp(p);

  LinearResampler< DRawDataPoint<1> > tmp2(tmp);
  TEST_REAL_EQUAL(tmp2.getSpacing(),0.5);
  TEST_EQUAL(tmp2.getParam(),p);
RESULT

CHECK(LinearResampler& operator= (const LinearResampler& source))
  Param p;
  p.setValue("ResamplingWidth",0.5);
  LinearResampler< DRawDataPoint<1> > tmp(p);

  LinearResampler< DRawDataPoint<1> > tmp2;
  tmp2 = tmp;
  TEST_REAL_EQUAL(tmp2.getSpacing(),0.5);
  TEST_EQUAL(tmp2.getParam(),p);
RESULT

CHECK(LinearResampler& operator()(MSExperiment< PeakType >& ms_exp))
  LinearResampler< DRawDataPoint<1> > lr;
  MSExperiment< DRawDataPoint<1> > ms_exp;
  lr(ms_exp);

  TEST_EQUAL(lr.getResampledDataPointer(), &ms_exp);
RESULT

CHECK(const MSExperiment<PeakType>* getResampledDataPointer() const)
  const LinearResampler< DRawDataPoint<1> > tmp;

  TEST_EQUAL(tmp.getResampledDataPointer(), 0);
RESULT

CHECK(const Param& getParam() const)
  Param p;
  p.setValue("ResamplingWidth",0.5);
  const LinearResampler< DRawDataPoint<1> > tmp(p);

  TEST_EQUAL(tmp.getParam(),p);
RESULT

CHECK(const double getSpacing() const)
  const LinearResampler< DRawDataPoint<1> > tmp;

  TEST_EQUAL(tmp.getSpacing(),0.05);
RESULT

CHECK(double getSpacing())
  LinearResampler< DRawDataPoint<1> > tmp;
  tmp.getSpacing() = 0.1;

  TEST_EQUAL(tmp.getSpacing(),0.1);
RESULT

CHECK(void setParam(const Param& param))
  Param p;
  p.setValue("ResamplingWidth",0.5);
  LinearResampler< DRawDataPoint<1> > tmp;
  tmp.setParam(p);

  TEST_EQUAL(tmp.getParam(),p);
RESULT

CHECK(void setResampledDataPointer(const MSExperiment<PeakType>& ms_exp))
  LinearResampler< DRawDataPoint<1> > lr;
  MSExperiment< DRawDataPoint<1> > ms_exp;
  lr.setResampledDataPointer(ms_exp);

  TEST_EQUAL(lr.getResampledDataPointer(), &ms_exp);
RESULT

CHECK(void setSpacing(const double spacing))
  LinearResampler< DRawDataPoint<1> > tmp;
  tmp.setSpacing(0.1);

  TEST_EQUAL(tmp.getSpacing(),0.1);
RESULT

CHECK((void start(ConstPeakIterator first, ConstPeakIterator last, PeakIterator resampled_first)))
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

  LinearResampler< DRawDataPoint<1> > lr;
  lr.setSpacing(0.5);

  int number_resampled_points = (int)ceil(((spec.end()-1)->getPos() -spec.begin()->getPos()) / lr.getSpacing() + 1);
  MSSpectrum< DRawDataPoint<1> > spec_resampled;
  spec_resampled.getContainer().resize(number_resampled_points);
  lr.start(spec.begin(),spec.end(),spec_resampled.begin());

  double sum = 0.;
  for (int i=0; i < 5; ++i)
  {
    sum += spec_resampled.getContainer()[i].getIntensity();
  }

  TEST_REAL_EQUAL(sum, 20);
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
