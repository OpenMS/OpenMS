// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Stephan Aiche$
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/SIMULATION/EGHModel.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(EGHModel, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

EGHModel* ptr = 0;
EGHModel* nullPointer = 0;
START_SECTION(EGHModel())
{
	ptr = new EGHModel();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((EGHModel(const EGHModel &source)))
{
  EGHModel egh1;
  egh1.setInterpolationStep(0.2);

  Param tmp;
  tmp.setValue("statistics:mean", 680.1 );
  tmp.setValue("statistics:variance",  2.0);
  tmp.setValue("egh:retention", 680.1);
  tmp.setValue("egh:height",100000.0);
  tmp.setValue("egh:A",150.0);
  tmp.setValue("egh:B",100.0);
  tmp.setValue("egh:alpha",0.4);
  egh1.setParameters(tmp);

  EGHModel egh2(egh1);
  EGHModel egh3;
  egh3.setInterpolationStep(0.2);
  egh3.setParameters(tmp);

  egh1 = EGHModel();
  TEST_EQUAL(egh3.getParameters(), egh2.getParameters())
  TEST_EQUAL(egh3==egh2,true)
}
END_SECTION

START_SECTION((virtual ~EGHModel()))
{
  delete ptr;
}
END_SECTION

START_SECTION((virtual EGHModel& operator=(const EGHModel &source)))
{
  EGHModel egh1;
  egh1.setInterpolationStep(0.2);

  Param tmp;
  tmp.setValue("statistics:mean", 680.1 );
  tmp.setValue("statistics:variance",  2.0);
  tmp.setValue("egh:retention", 680.1);
  tmp.setValue("egh:height",100000.0);
  tmp.setValue("egh:A",150.0);
  tmp.setValue("egh:B",100.0);
  tmp.setValue("egh:alpha",0.4);
  egh1.setParameters(tmp);

  EGHModel egh2;
  egh2 = egh1;

  EGHModel egh3;
  egh3.setInterpolationStep(0.2);
  egh3.setParameters(tmp);

  egh1 = EGHModel();
  TEST_EQUAL(egh3.getParameters(), egh2.getParameters())
  TEST_EQUAL(egh3==egh2,true)
}
END_SECTION

START_SECTION((void setOffset(CoordinateType offset)))
{
  //
  EGHModel egh1;
  egh1.setInterpolationStep(0.2);

  Param tmp;
  tmp.setValue("statistics:mean", 680.1 );
  tmp.setValue("statistics:variance",  2.0);
  tmp.setValue("egh:retention", 680.1);
  tmp.setValue("egh:height",100000.0);
  tmp.setValue("egh:A",150.0);
  tmp.setValue("egh:B",100.0);
  tmp.setValue("egh:alpha",0.4);
  egh1.setParameters(tmp);

  //
  EGHModel::CoordinateType current_offset = egh1.getInterpolation().getOffset();
  EGHModel::CoordinateType current_mean = egh1.getCenter();
  EGHModel::CoordinateType new_offset = current_offset + 10.0;
  egh1.setOffset(new_offset);

  TEST_REAL_SIMILAR(egh1.getInterpolation().getOffset(), new_offset)
  TEST_REAL_SIMILAR(egh1.getCenter(), current_mean + 10.0)
}
END_SECTION

START_SECTION((void setSamples()))
{
  EGHModel egh1;

  Param tmp;
  tmp.setValue("statistics:mean", 1000.0 );
  tmp.setValue("statistics:variance",  2.0);
  tmp.setValue("egh:retention", 1000.0);
  tmp.setValue("egh:height",100.0);
  tmp.setValue("egh:A",10.0);
  tmp.setValue("egh:B",20.0);
  tmp.setValue("egh:alpha",0.5);
  egh1.setInterpolationStep(0.2);
  egh1.setParameters(tmp); // setSamples() is called here

  TEST_REAL_SIMILAR(egh1.getInterpolation().value(1000.0), 100.0)
  TEST_REAL_SIMILAR(egh1.getInterpolation().value(990.0), 50.0) // corresponds to A_
  TEST_REAL_SIMILAR(egh1.getInterpolation().value(1020.0), 50.0) // corresponds to B_
}
END_SECTION

START_SECTION((CoordinateType getCenter() const ))
{
  //
  EGHModel egh1;
  egh1.setInterpolationStep(0.2);

  Param tmp;
  tmp.setValue("statistics:mean", 680.1 );
  tmp.setValue("statistics:variance",  2.0);
  tmp.setValue("egh:retention", 680.1);
  tmp.setValue("egh:height",100000.0);
  tmp.setValue("egh:A",150.0);
  tmp.setValue("egh:B",100.0);
  tmp.setValue("egh:alpha",0.4);
  egh1.setParameters(tmp);

  //
  EGHModel::CoordinateType current_offset = egh1.getInterpolation().getOffset();
  EGHModel::CoordinateType current_mean = egh1.getCenter();
  EGHModel::CoordinateType new_offset = current_offset + 10.0;
  egh1.setOffset(new_offset);

  TEST_REAL_SIMILAR(egh1.getInterpolation().getOffset(), new_offset)
  TEST_REAL_SIMILAR(egh1.getCenter(), current_mean + 10.0)
}
END_SECTION

START_SECTION((static BaseModel<1>* create()))
{
  BaseModel<1>* ptr = EGHModel::create();
  TEST_EQUAL(ptr->getName(), "EGHModel")
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((static const String getProductName()))
{
  TEST_EQUAL(EGHModel::getProductName(),"EGHModel")
  TEST_EQUAL(EGHModel().getName(),"EGHModel")
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



