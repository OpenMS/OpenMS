// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Stephan Aiche $
// $Authors: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/SIMULATION/IsotopeModelGeneral.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(IsotopeModelGeneral, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IsotopeModelGeneral* ptr = 0;
START_SECTION(IsotopeModelGeneral())
{
	ptr = new IsotopeModelGeneral();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION((virtual ~IsotopeModelGeneral()))
{
	delete ptr;
}
END_SECTION

START_SECTION((IsotopeModelGeneral(const IsotopeModelGeneral &source)))
{
	IsotopeModelGeneral im1;
	
	Param tmp;
	tmp.setValue("charge", 3);
	tmp.setValue("isotope:stdev",0.8);
	tmp.setValue("statistics:mean", 670.5);
	im1.setParameters(tmp);

  IsotopeModelGeneral im2;
  im2 = im1;

  IsotopeModelGeneral im3;
	im3.setParameters(tmp);

  im1 = IsotopeModelGeneral();
	TEST_EQUAL(im3.getParameters(), im2.getParameters())
	TEST_EQUAL(im3==im2,true)
}
END_SECTION

START_SECTION((virtual IsotopeModelGeneral& operator=(const IsotopeModelGeneral &source)))
{
	IsotopeModelGeneral im1;
	
	Param tmp;
	tmp.setValue("charge", 3);
	tmp.setValue("isotope:stdev",0.8);
	tmp.setValue("statistics:mean", 670.5);
	im1.setParameters(tmp);

	IsotopeModelGeneral im2(im1);
  IsotopeModelGeneral im3;
	im3.setParameters(tmp);

  im1 = IsotopeModelGeneral();
	TEST_EQUAL(im3.getParameters(), im2.getParameters())
	TEST_EQUAL(im3==im2,true)
}
END_SECTION

START_SECTION((void setSamples(EmpiricalFormula formula)))
{
	TOLERANCE_ABSOLUTE(0.02)
	IsotopeModelGeneral img;
	EmpiricalFormula ef("C36H24");

	ef += "H1"; // adducts are modeled via ef, not internally by the class
	Param p;
	p.setValue("interpolation_step",0.0005);
  p.setValue("statistics:mean",ef.getAverageWeight());
  p.setValue("isotope:stdev",0.5);
  p.setValue("charge",1);

	// init model
	img.setSamples(ef);
	img.setParameters(p);

	TEST_REAL_SIMILAR(img.getIntensity( DPosition<1>(ef.getAverageWeight()-1.0) ),0.211)
	TEST_REAL_SIMILAR(img.getIntensity( DPosition<1>(ef.getAverageWeight())),0.59)
	TEST_REAL_SIMILAR(img.getIntensity( DPosition<1>(ef.getAverageWeight()+1.0) ),0.158)
}
END_SECTION

START_SECTION((static BaseModel<1>* create()))
{
	BaseModel<1>* ptr = IsotopeModelGeneral::create();
	TEST_EQUAL(ptr->getName(), "IsotopeModelGeneral")
	TEST_NOT_EQUAL(ptr, 0)
	delete ptr;
}
END_SECTION

START_SECTION((static const String getProductName()))
{
	TEST_EQUAL(IsotopeModelGeneral::getProductName(),"IsotopeModelGeneral")
	TEST_EQUAL(IsotopeModelGeneral().getName(),"IsotopeModelGeneral")
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



