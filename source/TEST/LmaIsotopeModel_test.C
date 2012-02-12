// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// $Authors: $
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
LmaIsotopeModel* nullPointer = 0;
START_SECTION((LmaIsotopeModel()))
	ptr = new LmaIsotopeModel();
  	TEST_EQUAL(ptr->getName(), "LmaIsotopeModel")
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

// destructor
START_SECTION((virtual ~LmaIsotopeModel()))
	delete ptr;
END_SECTION

START_SECTION((static const String getProductName()))
	TEST_EQUAL(LmaIsotopeModel::getProductName(),"LmaIsotopeModel")
	TEST_EQUAL(LmaIsotopeModel().getName(),"LmaIsotopeModel")
END_SECTION

START_SECTION((static BaseModel<1>* create()))
	BaseModel<1>* ptr = LmaIsotopeModel::create();
	TEST_EQUAL(ptr->getName(), "LmaIsotopeModel")
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

// assignment operator
START_SECTION((virtual LmaIsotopeModel& operator=(const LmaIsotopeModel &source)))
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
END_SECTION

// copy constructor
START_SECTION((LmaIsotopeModel(const LmaIsotopeModel& source)))
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
END_SECTION

      
START_SECTION([EXTRA] DefaultParamHandler::setParameters(...))
  TOLERANCE_ABSOLUTE(0.001)
  LmaIsotopeModel im1;
  Param tmp;
  tmp.setValue("charge", 3);
  tmp.setValue("isotope:stdev",0.8);
  tmp.setValue("statistics:mean", 670.5);
  im1.setParameters(tmp);

  LmaIsotopeModel im2;
  im2.setParameters(im1.getParameters());

	std::vector<Peak1D> dpa1;
	std::vector<Peak1D> dpa2;
  im1.getSamples(dpa1);
  im2.getSamples(dpa2);

  TOLERANCE_ABSOLUTE(0.00001)
  TEST_EQUAL(dpa1.size(),dpa2.size())
  ABORT_IF(dpa1.size()!=dpa2.size());
  for (Size i=0; i<dpa1.size(); ++i)
  {
    TEST_REAL_SIMILAR(dpa1[i].getPosition()[0],dpa2[i].getPosition()[0])
    TEST_REAL_SIMILAR(dpa1[i].getIntensity(),dpa2[i].getIntensity())
  }
END_SECTION

START_SECTION(UInt getCharge() )
  // can only reliably be tested after fitting, only sanity check here
  LmaIsotopeModel im1;
  TEST_EQUAL(im1.getCharge() == 1, true)		// default charge is 1
END_SECTION    
    
START_SECTION((CoordinateType getCenter() const))
  // can only reliably be tested after fitting, only sanity check here
  LmaIsotopeModel im1;
	TEST_EQUAL(im1.getCenter() == 0, true)
END_SECTION

START_SECTION( void setOffset(CoordinateType offset) )
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
	
	std::vector<Peak1D> dpa1;
	std::vector<Peak1D> dpa2;
  im1.getSamples(dpa1);
  im2.getSamples(dpa2);

  TEST_EQUAL(dpa1.size(),dpa2.size())
  ABORT_IF(dpa1.size()!=dpa2.size());
  for (Size i=0; i<dpa1.size(); ++i)
  {
    TEST_REAL_SIMILAR(dpa1[i].getPosition()[0],dpa2[i].getPosition()[0])
    TEST_REAL_SIMILAR(dpa1[i].getIntensity(),dpa2[i].getIntensity())
  }
END_SECTION

START_SECTION( CoordinateType getOffset() )
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
	
	std::vector<Peak1D> dpa1;
	std::vector<Peak1D> dpa2;
  im1.getSamples(dpa1);
  im2.getSamples(dpa2);

  TEST_EQUAL(dpa1.size(),dpa2.size())
  ABORT_IF(dpa1.size()!=dpa2.size());
  for (Size i=0; i<dpa1.size(); ++i)
  {
    TEST_REAL_SIMILAR(dpa1[i].getPosition()[0],dpa2[i].getPosition()[0])
    TEST_REAL_SIMILAR(dpa1[i].getIntensity(),dpa2[i].getIntensity())
  }	
END_SECTION

START_SECTION((void setSamples()))
{
  // dummy subtest
	TEST_EQUAL(1,1)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
