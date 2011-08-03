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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/VISUAL/MultiGradient.h>
#include <QtGui/QColor>
#include <OpenMS/CONCEPT/Types.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MultiGradient, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MultiGradient* d10_ptr = 0;
MultiGradient* d10_nullPointer = 0;
START_SECTION((MultiGradient()))
	d10_ptr = new MultiGradient();
  TEST_NOT_EQUAL(d10_ptr, d10_nullPointer)
END_SECTION

START_SECTION((~MultiGradient()))
	delete d10_ptr;
END_SECTION

START_SECTION((InterpolationMode getInterpolationMode() const))
	TEST_EQUAL(MultiGradient().getInterpolationMode(),MultiGradient::IM_LINEAR)
END_SECTION

START_SECTION((void setInterpolationMode(InterpolationMode mode)))
	MultiGradient mg;
	mg.setInterpolationMode(MultiGradient::IM_STAIRS);
	TEST_EQUAL(mg.getInterpolationMode(),MultiGradient::IM_STAIRS)
END_SECTION

START_SECTION((Size size() const))
	MultiGradient mg;
	TEST_EQUAL(mg.size(),2);
END_SECTION

START_SECTION((UInt position(UInt index)))
	MultiGradient mg;
	TEST_EQUAL(mg.position(0),0);
	TEST_EQUAL(mg.position(1),100);
END_SECTION

START_SECTION((QColor color(UInt index)))
	MultiGradient mg;
	TEST_EQUAL(mg.color(0)==Qt::white,true);
	TEST_EQUAL(mg.color(1)==Qt::black,true);
END_SECTION

START_SECTION((void insert(DoubleReal position, QColor color)))
	MultiGradient mg;
	mg.insert(50,Qt::red);
	TEST_EQUAL(mg.size(),3);
	TEST_EQUAL(mg.position(0),0);
	TEST_EQUAL(mg.position(1),50);	
	TEST_EQUAL(mg.position(2),100);
	TEST_EQUAL(mg.color(0)==Qt::white,true);
	TEST_EQUAL(mg.color(1)==Qt::red,true);
	TEST_EQUAL(mg.color(2)==Qt::black,true);
	mg.insert(50,Qt::red);
	TEST_EQUAL(mg.size(),3);
	TEST_EQUAL(mg.position(0),0);
	TEST_EQUAL(mg.position(1),50);	
	TEST_EQUAL(mg.position(2),100);
	TEST_EQUAL(mg.color(0)==Qt::white,true);
	TEST_EQUAL(mg.color(1)==Qt::red,true);
	TEST_EQUAL(mg.color(2)==Qt::black,true);
	mg.insert(25,Qt::green);
	mg.insert(75,Qt::blue);
	TEST_EQUAL(mg.size(),5);
	TEST_EQUAL(mg.position(0),0);
	TEST_EQUAL(mg.position(1),25);	
	TEST_EQUAL(mg.position(2),50);
	TEST_EQUAL(mg.position(3),75);	
	TEST_EQUAL(mg.position(4),100);
	TEST_EQUAL(mg.color(0)==Qt::white,true);
	TEST_EQUAL(mg.color(1)==Qt::green,true);
	TEST_EQUAL(mg.color(2)==Qt::red,true);
	TEST_EQUAL(mg.color(3)==Qt::blue,true);
	TEST_EQUAL(mg.color(4)==Qt::black,true);
	mg.insert(76,Qt::magenta);
	TEST_EQUAL(mg.size(),6);
	TEST_EQUAL(mg.position(0),0);
	TEST_EQUAL(mg.position(1),25);	
	TEST_EQUAL(mg.position(2),50);
	TEST_EQUAL(mg.position(3),75);
	TEST_EQUAL(mg.position(4),76);
	TEST_EQUAL(mg.position(5),100);
	TEST_EQUAL(mg.color(0)==Qt::white,true);
	TEST_EQUAL(mg.color(1)==Qt::green,true);
	TEST_EQUAL(mg.color(2)==Qt::red,true);
	TEST_EQUAL(mg.color(3)==Qt::blue,true);
	TEST_EQUAL(mg.color(4)==Qt::magenta,true);
	TEST_EQUAL(mg.color(5)==Qt::black,true);
END_SECTION

START_SECTION((bool remove(DoubleReal position)))
	MultiGradient mg;
	mg.insert(25,Qt::green);
	mg.insert(50,Qt::red);
	mg.insert(75,Qt::blue);
	mg.remove(50);
	TEST_EQUAL(mg.size(),4);
	TEST_EQUAL(mg.position(0),0);
	TEST_EQUAL(mg.position(1),25);
	TEST_EQUAL(mg.position(2),75);	
	TEST_EQUAL(mg.position(3),100);
	TEST_EQUAL(mg.color(0)==Qt::white,true);
	TEST_EQUAL(mg.color(1)==Qt::green,true);
	TEST_EQUAL(mg.color(2)==Qt::blue,true);
	TEST_EQUAL(mg.color(3)==Qt::black,true);
	mg.remove(25);
	TEST_EQUAL(mg.size(),3);
	TEST_EQUAL(mg.position(0),0);
	TEST_EQUAL(mg.position(1),75);	
	TEST_EQUAL(mg.position(2),100);
	TEST_EQUAL(mg.color(0)==Qt::white,true);
	TEST_EQUAL(mg.color(1)==Qt::blue,true);
	TEST_EQUAL(mg.color(2)==Qt::black,true);
	mg.remove(75);
	TEST_EQUAL(mg.size(),2);
	TEST_EQUAL(mg.position(0),0);
	TEST_EQUAL(mg.position(1),100);
	TEST_EQUAL(mg.color(0)==Qt::white,true);
	TEST_EQUAL(mg.color(1)==Qt::black,true);
END_SECTION

START_SECTION((bool exists(DoubleReal position)))
	MultiGradient mg;
	mg.insert(25,Qt::green);
	mg.insert(50,Qt::red);
	mg.insert(75,Qt::blue);
	TEST_EQUAL(mg.exists(0),true);
	TEST_EQUAL(mg.exists(1),false);
	TEST_EQUAL(mg.exists(25),true);
	TEST_EQUAL(mg.exists(49),false);
	TEST_EQUAL(mg.exists(50),true);
	TEST_EQUAL(mg.exists(51),false);
	TEST_EQUAL(mg.exists(75),true);
	TEST_EQUAL(mg.exists(99),false);
	TEST_EQUAL(mg.exists(100),true);
END_SECTION

START_SECTION((QColor interpolatedColorAt(DoubleReal position) const))
	MultiGradient mg;
	TEST_EQUAL(mg.interpolatedColorAt(0)==Qt::white,true);
	TEST_EQUAL(mg.interpolatedColorAt(25)==QColor(191,191,191),true);
	TEST_EQUAL(mg.interpolatedColorAt(50)==QColor(127,127,127),true);
	TEST_EQUAL(mg.interpolatedColorAt(75)==QColor(63,63,63),true);
	TEST_EQUAL(mg.interpolatedColorAt(100)==Qt::black,true);
	mg.insert(50,Qt::red);
	TEST_EQUAL(mg.interpolatedColorAt(0)==Qt::white,true);
	TEST_EQUAL(mg.interpolatedColorAt(25)==QColor(255,127,127),true);
	TEST_EQUAL(mg.interpolatedColorAt(50)==Qt::red,true);
	TEST_EQUAL(mg.interpolatedColorAt(75)==QColor(127,0,0),true);
	TEST_EQUAL(mg.interpolatedColorAt(100)==Qt::black,true);	
	mg.insert(50,Qt::green);
	TEST_EQUAL(mg.interpolatedColorAt(0)==Qt::white,true);
	TEST_EQUAL(mg.interpolatedColorAt(25)==QColor(127,255,127),true);
	TEST_EQUAL(mg.interpolatedColorAt(50)==Qt::green,true);
	TEST_EQUAL(mg.interpolatedColorAt(75)==QColor(0,127,0),true);
	TEST_EQUAL(mg.interpolatedColorAt(100)==Qt::black,true);	
	mg.insert(50,Qt::blue);
	TEST_EQUAL(mg.interpolatedColorAt(0)==Qt::white,true);
	TEST_EQUAL(mg.interpolatedColorAt(25)==QColor(127,127,255),true);
	TEST_EQUAL(mg.interpolatedColorAt(50)==Qt::blue,true);
	TEST_EQUAL(mg.interpolatedColorAt(75)==QColor(0,0,127),true);
	TEST_EQUAL(mg.interpolatedColorAt(100)==Qt::black,true);	

	MultiGradient mg2;
	mg2.setInterpolationMode(MultiGradient::IM_STAIRS);
	TEST_EQUAL(mg2.interpolatedColorAt(0)==Qt::white,true);
	TEST_EQUAL(mg2.interpolatedColorAt(25)==Qt::white,true);
	TEST_EQUAL(mg2.interpolatedColorAt(100)==Qt::black,true);
	mg2.insert(50,Qt::red);
	TEST_EQUAL(mg2.interpolatedColorAt(0)==Qt::white,true);
	TEST_EQUAL(mg2.interpolatedColorAt(49)==Qt::white,true);
	TEST_EQUAL(mg2.interpolatedColorAt(50)==Qt::red,true);
	TEST_EQUAL(mg2.interpolatedColorAt(51)==Qt::red,true);
	TEST_EQUAL(mg2.interpolatedColorAt(99)==Qt::red,true);
	TEST_EQUAL(mg2.interpolatedColorAt(100)==Qt::black,true);	
END_SECTION

START_SECTION((QColor interpolatedColorAt(DoubleReal position, DoubleReal min, DoubleReal max) const))
	MultiGradient mg;
	mg.insert(50,Qt::red);
	TEST_EQUAL(mg.interpolatedColorAt(0,0,100)==Qt::white,true);
	TEST_EQUAL(mg.interpolatedColorAt(25,0,100)==QColor(255,127,127),true);
	TEST_EQUAL(mg.interpolatedColorAt(50,0,100)==Qt::red,true);
	TEST_EQUAL(mg.interpolatedColorAt(75,0,100)==QColor(127,0,0),true);
	TEST_EQUAL(mg.interpolatedColorAt(100,0,100)==Qt::black,true);	

	MultiGradient mg2;
	mg2.setInterpolationMode(MultiGradient::IM_STAIRS);
	mg2.insert(50,Qt::red);
	TEST_EQUAL(mg2.interpolatedColorAt(0)==Qt::white,true);
	TEST_EQUAL(mg2.interpolatedColorAt(49)==Qt::white,true);
	TEST_EQUAL(mg2.interpolatedColorAt(50)==Qt::red,true);
	TEST_EQUAL(mg2.interpolatedColorAt(51)==Qt::red,true);
	TEST_EQUAL(mg2.interpolatedColorAt(99)==Qt::red,true);
	TEST_EQUAL(mg2.interpolatedColorAt(100)==Qt::black,true);	
END_SECTION

START_SECTION((void activatePrecalculationMode(DoubleReal min, DoubleReal max, UInt steps)))
NOT_TESTABLE
END_SECTION

START_SECTION((QColor precalculatedColorAt(DoubleReal position) const ))
	MultiGradient mg;
	mg.insert(0,Qt::white);
	mg.insert(100,Qt::blue);
	mg.activatePrecalculationMode(-50.0,50.0,100);
	
	//Test precalclulated Values
	TEST_EQUAL(mg.precalculatedColorAt(-50.0).red(),255);
	TEST_EQUAL(mg.precalculatedColorAt(-50.0).green(),255);
	TEST_EQUAL(mg.precalculatedColorAt(-50.0).blue(),255);
	
	TEST_EQUAL(mg.precalculatedColorAt(-25.0).red(),193);
	TEST_EQUAL(mg.precalculatedColorAt(-25.0).green(),193);
	TEST_EQUAL(mg.precalculatedColorAt(-25.0).blue(),255);
	
	TEST_EQUAL(mg.precalculatedColorAt(0.0).red(),128);
	TEST_EQUAL(mg.precalculatedColorAt(0.0).green(),128);
	TEST_EQUAL(mg.precalculatedColorAt(0.0).blue(),255);
	
	TEST_EQUAL(mg.precalculatedColorAt(25.0).red(),64);
	TEST_EQUAL(mg.precalculatedColorAt(25.0).green(),64);
	TEST_EQUAL(mg.precalculatedColorAt(25.0).blue(),255);

	TEST_EQUAL(mg.precalculatedColorAt(50.0).red(),2);
	TEST_EQUAL(mg.precalculatedColorAt(50.0).green(),2);
	TEST_EQUAL(mg.precalculatedColorAt(50.0).blue(),255);
END_SECTION

START_SECTION((void deactivatePrecalculationMode()))
	MultiGradient mg;
	mg.activatePrecalculationMode(-50,50,100);
	mg.deactivatePrecalculationMode();
	NOT_TESTABLE
END_SECTION

START_SECTION((std::string toString() const))
	MultiGradient mg;
	TEST_EQUAL(mg.toString(),"Linear|0,#ffffff;100,#000000")
	mg.setInterpolationMode(MultiGradient::IM_STAIRS);
	mg.insert(50,Qt::red);
	TEST_EQUAL(mg.toString(),"Stairs|0,#ffffff;50,#ff0000;100,#000000")
END_SECTION

START_SECTION((void fromString(const std::string& gradient)))
	MultiGradient mg;
	mg.fromString("Linear|0,#ff0000;100,#000000");
	TEST_EQUAL(mg.getInterpolationMode(),MultiGradient::IM_LINEAR)
	TEST_EQUAL(mg.size(),2)
	TEST_EQUAL(mg.color(0)==Qt::red, true);
	TEST_EQUAL(mg.color(1)==Qt::black, true);
	TEST_EQUAL(mg.position(0), 0);
	TEST_EQUAL(mg.position(1), 100);
	mg.fromString("Stairs|0,#ffffff;50,#ff0000;100,#000000");
	TEST_EQUAL(mg.getInterpolationMode(),MultiGradient::IM_STAIRS)
	TEST_EQUAL(mg.size(),3)
	TEST_EQUAL(mg.color(0)==Qt::white, true);
	TEST_EQUAL(mg.color(1)==Qt::red, true);
	TEST_EQUAL(mg.color(2)==Qt::black, true);
	TEST_EQUAL(mg.position(0), 0);
	TEST_EQUAL(mg.position(1), 50);
	TEST_EQUAL(mg.position(2), 100);
END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



