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

#include <OpenMS/FILTERING/SMOOTHING/DSavitzkyGolayQRFilter.h>

///////////////////////////

START_TEST(DSavitzkyGolayQRFilter<D>, "$Id: DSavitzkyGolayQRFilter_test.C,v 1.10 2006/04/13 08:40:46 elange Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

typedef DPeakArrayNonPolymorphic<2,DRawDataPoint<2> > RawDataArray2D;
typedef RawDataArray2D::Iterator RawDataIterator2D;
typedef RawDataArray2D::ConstIterator RawDataConstIterator2D;
typedef DPeakArrayNonPolymorphic<1,DRawDataPoint<1> > RawDataArray1D;
typedef RawDataArray1D::Iterator RawDataIterator1D;
typedef RawDataArray1D::ConstIterator RawDataConstIterator1D;


DSavitzkyGolayQRFilter<1>* dsg_ptr = 0;
CHECK(DSavitzkyGolayQRFilter())
  dsg_ptr = new DSavitzkyGolayQRFilter<1>;
  TEST_NOT_EQUAL(dsg_ptr, 0)
RESULT

CHECK(~DSavitzkyGolayQRFilter())
  delete dsg_ptr;
RESULT

CHECK(DSavitzkyGolayQRFilter(const Param& p))
  Param p;
  p.setValue("PolynomOrder",2);
  p.setValue("FrameLength",3);
  DSavitzkyGolayQRFilter<1> sgolay(p);

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


  std::vector<double> coeffs = sgolay.getCoeffs();

  sgolay.filter(raw.begin(),raw.end(),filtered.begin());
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

  TEST_EQUAL(sgolay.getCoeffs().size(), 6)
  TEST_REAL_EQUAL(sgolay.getOrder(), 2)
  TEST_EQUAL(sgolay.getRTdim(), -1)
  TEST_EQUAL(sgolay.getMZdim(), 0)
RESULT

CHECK(DSavitzkyGolayQRFilter& operator=(const DSavitzkyGolayQRFilter& s))
  RawDataArray2D* sgolay_data=0;
  DSavitzkyGolayQRFilter<2> sgolay;
  sgolay.setOrder(4);
  sgolay.setWindowSize(5);
  sgolay(*sgolay_data);

  DSavitzkyGolayQRFilter<2> sgolay_copy=sgolay;
  TEST_EQUAL(sgolay_copy.getFilteredDataPointer(), 0)
  TEST_EQUAL(sgolay_copy.getFilteredDataPointer(), 0)
  TEST_REAL_EQUAL(sgolay_copy.getOrder(),sgolay.getOrder())
  TEST_EQUAL(sgolay_copy.getWindowSize(),sgolay.getWindowSize())
  TEST_EQUAL(sgolay_copy.getRTdim(), sgolay.getRTdim())
  TEST_EQUAL(sgolay_copy.getMZdim(), sgolay.getMZdim())
RESULT

CHECK(DSavitzkyGolayQRFilter(const DSavitzkyGolayQRFilter& s))
  RawDataArray2D* sgolay_data=0;
  DSavitzkyGolayQRFilter<2> sgolay;
  sgolay.setOrder(4);
  sgolay.setWindowSize(5);
  sgolay(*sgolay_data);

  DSavitzkyGolayQRFilter<2> sgolay_copy(sgolay);
  TEST_EQUAL(sgolay_copy.getFilteredDataPointer(), 0)
  TEST_REAL_EQUAL(sgolay_copy.getOrder(),sgolay.getOrder())
  TEST_EQUAL(sgolay_copy.getWindowSize(),sgolay.getWindowSize())
  TEST_EQUAL(sgolay_copy.getRTdim(), sgolay.getRTdim())
  TEST_EQUAL(sgolay_copy.getMZdim(), sgolay.getMZdim())
RESULT

CHECK(const Param& getParam() const)
  Param p;
  p.setValue("PolynomOrder",2);
  p.setValue("FrameLength",3);
  const DSavitzkyGolayQRFilter<1> sgolay(p);

  TEST_REAL_EQUAL(sgolay.getParam().getValue("PolynomOrder"),2);
  TEST_REAL_EQUAL(sgolay.getParam().getValue("FrameLength"),3);
RESULT

CHECK(const unsigned int& getOrder() const)
  DSavitzkyGolayQRFilter<1> sgolay;

  TEST_EQUAL(sgolay.getOrder(),4);
RESULT

CHECK(const unsigned int& getWindowSize() const)
  DSavitzkyGolayQRFilter<1> sgolay;

  TEST_EQUAL(sgolay.getWindowSize(),17);
RESULT

CHECK(unsigned int& getOrder())
  DSavitzkyGolayQRFilter<1> sgolay;
  sgolay.getOrder() = 3;

  TEST_EQUAL(sgolay.getOrder(),3);
RESULT

CHECK(void setOrder(const unsigned int order))
  DSavitzkyGolayQRFilter<1> sgolay;
  sgolay.setOrder(3);

  TEST_EQUAL(sgolay.getOrder(),3);
RESULT

CHECK(void setParam(const Param& param))
  Param p;
  p.setValue("PolynomOrder",2);
  p.setValue("FrameLength",3);
  DSavitzkyGolayQRFilter<1> sgolay;
  sgolay.setParam(p);

  TEST_REAL_EQUAL(sgolay.getParam().getValue("PolynomOrder"),2);
  TEST_REAL_EQUAL(sgolay.getParam().getValue("FrameLength"),3);
RESULT

CHECK(void setWindowSize(const unsigned int frame_size))
  DSavitzkyGolayQRFilter<1> sgolay;
  sgolay.setWindowSize(7);

  TEST_EQUAL(sgolay.getWindowSize(),7);
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
