// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Marcel Grunert $
// --------------------------------------------------------------------------


#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/LmaIsotopeModel.h>

///////////////////////////

START_TEST(LmaIsotopeModel, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using std::stringstream;

// default ctor
LmaIsotopeModel* ptr = 0;
CHECK((LmaIsotopeModel()))
	ptr = new LmaIsotopeModel();
  	TEST_EQUAL(ptr->getName(), "LmaIsotopeModel")
	TEST_NOT_EQUAL(ptr, 0)
RESULT

// destructor
CHECK((virtual ~LmaIsotopeModel()))
	delete ptr;
RESULT

CHECK((static const String getProductName()))
	TEST_EQUAL(LmaIsotopeModel::getProductName(),"LmaIsotopeModel")
	TEST_EQUAL(LmaIsotopeModel().getName(),"LmaIsotopeModel")
RESULT

CHECK((static BaseModel<1>* create()))
	BaseModel<1>* ptr = LmaIsotopeModel::create();
	TEST_EQUAL(ptr->getName(), "LmaIsotopeModel")
	TEST_NOT_EQUAL(ptr, 0)
RESULT

// assignment operator
CHECK((virtual LmaIsotopeModel& operator=(const LmaIsotopeModel &source)))
  LmaIsotopeModel lim1;
	
  Param tmp;
  tmp.setValue("charge", 3);
  tmp.setValue("isotope:stdev",0.8);
  tmp.setValue("statistics:mean", 670.5);
  lim1.setParameters(tmp);

  LmaIsotopeModel lim2;
  lim2 = lim1;

  LmaIsotopeModel lim3;
  lim3.setParameters(tmp);

  lim1 = LmaIsotopeModel();
  TEST_EQUAL(lim3.getParameters(), lim2.getParameters())
RESULT

// copy constructor
CHECK((LmaIsotopeModel(const LmaIsotopeModel& source)))
	LmaIsotopeModel lim1;
	
	Param tmp;
  tmp.setValue("charge", 3);
  tmp.setValue("isotope:stdev",0.8);
  tmp.setValue("statistics:mean", 670.5);
  lim1.setParameters(tmp);

  LmaIsotopeModel lim2;
  lim2 = lim1;

  LmaIsotopeModel lim3;
  lim3.setParameters(tmp);

  lim1 = LmaIsotopeModel();
  TEST_EQUAL(lim3.getParameters(), lim2.getParameters())
RESULT

      
CHECK([EXTRA] DefaultParamHandler::setParameters(...))
  PRECISION(0.001)
  LmaIsotopeModel im1;
  Param tmp;
  tmp.setValue("charge", 3);
  tmp.setValue("isotope:stdev",0.8);
  tmp.setValue("statistics:mean", 670.5);
  im1.setParameters(tmp);

  LmaIsotopeModel im2;
  im2.setParameters(im1.getParameters());

  DPeakArray<DPeak<1> > dpa1;
  DPeakArray<DPeak<1> > dpa2;
  im1.getSamples(dpa1);
  im2.getSamples(dpa2);

  PRECISION(0.00001)
  TEST_EQUAL(dpa1.size(),dpa2.size())
  ABORT_IF(dpa1.size()!=dpa2.size());
  for (UInt i=0; i<dpa1.size(); ++i)
  {
    TEST_REAL_EQUAL(dpa1[i].getPosition()[0],dpa2[i].getPosition()[0])
    TEST_REAL_EQUAL(dpa1[i].getIntensity(),dpa2[i].getIntensity())
  }
RESULT

CHECK(UInt getCharge() )
  // can only reliably be tested after fitting, only sanity check here
  LmaIsotopeModel im1;
  TEST_EQUAL(im1.getCharge() == 1, true)		// default charge is 1
RESULT    
    
CHECK((CoordinateType getCenter() const))
  // can only reliably be tested after fitting, only sanity check here
  LmaIsotopeModel im1;
	TEST_EQUAL(im1.getCenter() == 0, true)
RESULT

CHECK( void setOffset(CoordinateType offset) )
  LmaIsotopeModel im1;
  Param tmp;
  tmp.setValue("charge", 3);
  tmp.setValue("isotope:stdev",0.8);
  tmp.setValue("statistics:mean", 670.5);
  im1.setParameters(tmp);
  im1.setOffset( 673.5 );
	
  LmaIsotopeModel im2;
  im2.setParameters(tmp);
  im2.setOffset( 673.5 );
	
  DPeakArray<DPeak<1> > dpa1;
  DPeakArray<DPeak<1> > dpa2;
  im1.getSamples(dpa1);
  im2.getSamples(dpa2);

  TEST_EQUAL(dpa1.size(),dpa2.size())
  ABORT_IF(dpa1.size()!=dpa2.size());
  for (UInt i=0; i<dpa1.size(); ++i)
  {
    TEST_REAL_EQUAL(dpa1[i].getPosition()[0],dpa2[i].getPosition()[0])
    TEST_REAL_EQUAL(dpa1[i].getIntensity(),dpa2[i].getIntensity())
  }
RESULT

CHECK( CoordinateType getOffset() )
  LmaIsotopeModel im1;
  Param tmp;
  tmp.setValue("charge", 3);
  tmp.setValue("isotope:stdev",0.8);
  tmp.setValue("statistics:mean", 670.5);
  im1.setParameters(tmp);
  im1.setOffset( 673.5 );
	
  LmaIsotopeModel im2;
  im2.setParameters(tmp);
  im2.setOffset( 673.5 );
	
  DPeakArray<DPeak<1> > dpa1;
  DPeakArray<DPeak<1> > dpa2;
  im1.getSamples(dpa1);
  im2.getSamples(dpa2);

  TEST_EQUAL(dpa1.size(),dpa2.size())
  ABORT_IF(dpa1.size()!=dpa2.size());
  for (UInt i=0; i<dpa1.size(); ++i)
  {
    TEST_REAL_EQUAL(dpa1[i].getPosition()[0],dpa2[i].getPosition()[0])
    TEST_REAL_EQUAL(dpa1[i].getIntensity(),dpa2[i].getIntensity())
  }	
RESULT

CHECK((void setSamples()))
{
  // dummy subtest
	TEST_EQUAL(1,1)
}
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
