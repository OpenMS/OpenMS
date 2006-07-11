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
// $Maintainer: Eva Lange  $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FILTERING/SMOOTHING/DGaussFilter.h>

///////////////////////////

START_TEST(DGaussFilter<D>, "$Id: DGaussFilter_test.C,v 1.9 2006/04/13 12:57:40 elange Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;


typedef DPeakArrayNonPolymorphic<2,DRawDataPoint<2> > RawDataArray2D;
typedef RawDataArray2D::Iterator RawDataIterator2D;
typedef RawDataArray2D::ConstIterator RawDataConstIterator2D;
typedef DPeakArrayNonPolymorphic<1,DRawDataPoint<1> > RawDataArray1D;
typedef RawDataArray1D::Iterator RawDataIterator1D;
typedef RawDataArray1D::ConstIterator RawDataConstIterator1D;


DGaussFilter<1>* dgauss_ptr = 0;
CHECK(DGaussFilter())
  dgauss_ptr = new DGaussFilter<1>;
  TEST_NOT_EQUAL(dgauss_ptr, 0)
RESULT

CHECK(~DGaussFilter())
    delete dgauss_ptr;
RESULT

CHECK(DGaussFilter& operator=(const DGaussFilter& s))
  RawDataArray2D gauss_data;
  DGaussFilter<2> gauss;
  gauss(gauss_data);
  DGaussFilter<2> gauss_copy;
  gauss_copy = gauss;

  TEST_EQUAL(gauss_copy.getFilteredDataPointer(),gauss.getFilteredDataPointer())
  TEST_REAL_EQUAL(gauss_copy.getSigma(),gauss.getSigma())
  TEST_EQUAL(gauss_copy.getRTdim(), gauss.getRTdim())
  TEST_EQUAL(gauss_copy.getMZdim(),gauss.getMZdim())
RESULT

CHECK(DGaussFilter(const DGaussFilter& g))
  RawDataArray2D gauss_data;
  DGaussFilter<2> gauss;
  gauss(gauss_data);
  DGaussFilter<2> gauss_copy(gauss);

  TEST_EQUAL(gauss_copy.getFilteredDataPointer(),0)
  TEST_REAL_EQUAL(gauss_copy.getSigma(),gauss.getSigma())
  TEST_EQUAL(gauss_copy.getRTdim(), gauss.getRTdim())
  TEST_EQUAL(gauss_copy.getMZdim(),gauss.getMZdim())
RESULT

CHECK(DGaussFilter(const Param& p))
  Param p;
  p.setValue("GaussianWidth",1.6);
  DGaussFilter<1> gauss(p);
  TEST_REAL_EQUAL(gauss.getSigma(),.2);
  TEST_REAL_EQUAL(gauss.getKernelWidth(),1.6);
RESULT

CHECK(const Param& getParam() const)
  Param p;
  p.setValue("GaussianWidth",1.6);
  DGaussFilter<1> gauss(p);

  TEST_REAL_EQUAL(gauss.getParam().getValue("GaussianWidth"),1.6);
RESULT

CHECK(const double& getSigma() const)
  const DGaussFilter<1> gaussian;

  TEST_REAL_EQUAL(gaussian.getSigma(),.1);
RESULT

CHECK(const double& getSpacing() const)
  const DGaussFilter<1> gaussian;

  TEST_REAL_EQUAL(gaussian.getSpacing(),0.01);
RESULT

CHECK(double getKernelWidth() const)
  const DGaussFilter<1> gaussian;

  TEST_REAL_EQUAL(gaussian.getKernelWidth(),.8);
RESULT

CHECK((void init(float sigma, float spacing)))
  DGaussFilter<1> gaussian;
  gaussian.init(0.2,0.001);

  TEST_REAL_EQUAL(gaussian.getSpacing(),0.001);
  TEST_REAL_EQUAL(gaussian.getSigma(),0.2);
  TEST_REAL_EQUAL(gaussian.getKernelWidth(),1.6);
RESULT

CHECK(void setKernelWidth(const double kernel_width))
  DGaussFilter<1> gaussian;
  gaussian.setKernelWidth(1.6);

  TEST_REAL_EQUAL(gaussian.getKernelWidth(),1.6);
RESULT

CHECK(void setParam(const Param& param))
  Param p;
  p.setValue("GaussianWidth",1.6);
  DGaussFilter<1> gaussian(p);

  TEST_REAL_EQUAL(gaussian.getSigma(),0.2);
  TEST_REAL_EQUAL(gaussian.getKernelWidth(),1.6);
RESULT

CHECK(void setSigma(float sigma))
  DGaussFilter<2> gauss;
  gauss.setSigma(2.434);
  TEST_REAL_EQUAL(gauss.getSigma(), 2.434)
RESULT

CHECK(void setSpacing(const double spacing))
  DGaussFilter<1> gaussian;
  gaussian.setSpacing(0.0001);

  TEST_REAL_EQUAL(gaussian.getSpacing(),0.0001);
RESULT

CHECK(void setKernelWidth(float kernel_width))
  DGaussFilter<2> gauss;
  gauss.setSigma(0.2);

  TEST_REAL_EQUAL(gauss.getSigma(), 0.2)
  TEST_REAL_EQUAL(gauss.getKernelWidth(), 1.6)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
