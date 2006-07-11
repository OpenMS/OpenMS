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

#include <OpenMS/FILTERING/SMOOTHING/DSmoothFilter.h>

///////////////////////////

START_TEST(DSmoothFilter<D>, "$Id: DSmoothFilter_test.C,v 1.5 2006/04/13 12:57:40 elange Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;


typedef DPeakArrayNonPolymorphic<2,DRawDataPoint<2> > RawDataArray2D;
typedef RawDataArray2D::Iterator RawDataIterator2D;
typedef RawDataArray2D::ConstIterator RawDataConstIterator2D;
typedef DPeakArrayNonPolymorphic<1,DRawDataPoint<1> > RawDataArray1D;
typedef RawDataArray1D::Iterator RawDataIterator1D;
typedef RawDataArray1D::ConstIterator RawDataConstIterator1D;


DSmoothFilter<1>* dsmooth_ptr = 0;
CHECK(DSmoothFilter())
  dsmooth_ptr = new DSmoothFilter<1>;
  TEST_NOT_EQUAL(dsmooth_ptr, 0)
RESULT

CHECK(~DSmoothFilter())
  delete dsmooth_ptr;
RESULT


CHECK(DSmoothFilter& operator=(const DSmoothFilter& fir))
  RawDataArray2D smooth_data;
  std::vector<double> coeffs(1);
  coeffs[0]=1.23423;
  DSmoothFilter<2> smooth;
  smooth.setCoeffs(coeffs);
  smooth(smooth_data);

  DSmoothFilter<2> smooth_copy;

  smooth_copy=smooth;

  TEST_EQUAL(smooth_copy.getFilteredDataPointer(), 0)
  TEST_EQUAL(smooth_copy.getCoeffs()[0], 1.23423)
  TEST_EQUAL(smooth_copy.getRTdim(), 0)
  TEST_EQUAL(smooth_copy.getMZdim(), 1)
RESULT

CHECK(DSmoothFilter(const DSmoothFilter& fir))
  RawDataArray2D smooth_data;
  std::vector<double> coeffs(1);
  coeffs[0]=1.23423;
  DSmoothFilter<2> smooth;
  smooth.setCoeffs(coeffs);
  smooth(smooth_data);

  DSmoothFilter<2> smooth_copy(smooth);
  TEST_EQUAL(smooth_copy.getFilteredDataPointer(), 0)
  TEST_EQUAL(smooth_copy.getCoeffs()[0], 1.23423)
  TEST_EQUAL(smooth_copy.getRTdim(), 0)
  TEST_EQUAL(smooth_copy.getMZdim(), 1)
RESULT

CHECK(const RawDataArray getFilteredDataPointer() const)
  const DSmoothFilter<1> smooth;

  TEST_EQUAL(smooth.getFilteredDataPointer(),0)
RESULT

CHECK(const std::vector<double>& getCoeffs() const)
  const DSmoothFilter<1> smooth;
  TEST_EQUAL(smooth.getCoeffs().size(),0)
RESULT

CHECK(const int& getMZdim() const)
  const DSmoothFilter<1> smooth;
  TEST_EQUAL(smooth.getMZdim(),0)
RESULT

CHECK(const int& getRTdim() const)
  const DSmoothFilter<2> smooth;
  TEST_EQUAL(smooth.getRTdim(),0)
RESULT

CHECK(std::vector<double>& getCoeffs())
  std::vector<double> coeffs(1);
  coeffs[0]=1;
  DSmoothFilter<1> smooth;
  smooth.setCoeffs(coeffs);
  TEST_EQUAL(smooth.getCoeffs()[0],1);
RESULT

CHECK(int& getMZdim())
  DSmoothFilter<1> smooth;
  TEST_EQUAL(smooth.getMZdim(),0)
RESULT

CHECK(int& getRTdim())
  DSmoothFilter<2> smooth;
  TEST_EQUAL(smooth.getRTdim(),0)
RESULT

CHECK((void filter(RawDataConstIterator first, RawDataConstIterator last, RawDataIterator new_first)))
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

  DSmoothFilter<1> smooth;
  smooth.setCoeffs(coeffs);
  smooth.filter(raw.begin(),raw.end(),filtered.begin());

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

CHECK(void setCoeffs(std::vector<double>& coeffs))
  std::vector<double> coeffs(3);
  coeffs[0]=1.2;
  coeffs[1]=2.4;
  coeffs[2]= 1.4;
  DSmoothFilter<1> smooth;
  smooth.getCoeffs()=coeffs;
  TEST_EQUAL(smooth.getCoeffs()[0],1.2)
  TEST_EQUAL(smooth.getCoeffs()[1],2.4)
  TEST_EQUAL(smooth.getCoeffs()[2],1.4)
RESULT

CHECK(void setFilteredDataPointer(RawDataArray& raw_filtered))
  RawDataArray2D raw_filtered;
  DSmoothFilter<2> smooth;
  smooth.setFilteredDataPointer(raw_filtered);

  TEST_EQUAL(*smooth.getFilteredDataPointer(),raw_filtered)
RESULT

CHECK(void setMZdim(int mz_dim))
  DSmoothFilter<2> smooth;
  smooth.setMZdim(0);
  TEST_EQUAL(smooth.getMZdim(),0)
RESULT

CHECK(void setRTdim(int rt_dim))
  DSmoothFilter<2> smooth;
  smooth.setRTdim(1);
  TEST_EQUAL(smooth.getRTdim(),1)
RESULT

CHECK(DSmoothFilter& operator()(RawDataArray& raw))
  RawDataArray2D raw_filtered;
  DSmoothFilter<2> smooth;
  smooth(raw_filtered);

  TEST_EQUAL(smooth.getFilteredDataPointer(),&raw_filtered)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

