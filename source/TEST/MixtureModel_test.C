// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/SIMULATION/MixtureModel.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MixtureModel, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MixtureModel* ptr = 0;
START_SECTION(MixtureModel())
{
	ptr = new MixtureModel();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION((virtual ~MixtureModel()))
{
	delete ptr;
}
END_SECTION

START_SECTION((MixtureModel(const MixtureModel &source)))
{
	MixtureModel mm1;
	
	Param p;
	p.setValue("mix", 0.5);
	p.setValue("statistics:variance1",0.8);
	p.setValue("statistics:variance2",0.8);
	p.setValue("statistics:mean1", 670.5);
	p.setValue("statistics:mean2", 672.5);
	mm1.setParameters(p);

  MixtureModel mm2(mm1);

  MixtureModel mm3;
	mm3.setParameters(p);

  mm1 = MixtureModel();
	TEST_EQUAL(mm3.getParameters(), mm2.getParameters())
	TEST_EQUAL(mm3==mm2,true)
}
END_SECTION

START_SECTION((virtual MixtureModel& operator=(const MixtureModel &source)))
{
	MixtureModel mm1;
	
	Param p;
	p.setValue("mix", 0.5);
	p.setValue("statistics:variance1",0.8);
	p.setValue("statistics:variance2",0.8);
	p.setValue("statistics:mean1", 670.5);
	p.setValue("statistics:mean2", 672.5);
	mm1.setParameters(p);

  MixtureModel mm2;
  mm2 = mm1;

  MixtureModel mm3;
	mm3.setParameters(p);

  mm1 = MixtureModel();
	TEST_EQUAL(mm3.getParameters(), mm2.getParameters())
	TEST_EQUAL(mm3==mm2,true)
}
END_SECTION

START_SECTION((void setOffset(double offset)))
{
  // Shamelessly copied from IsotopeModel_test !!
	MixtureModel mm1;
	Param p;
	p.setValue("mix", 0.5);
	p.setValue("statistics:variance1",0.8);
	p.setValue("statistics:variance2",0.8);
	p.setValue("statistics:mean1", 670.5);
	p.setValue("statistics:mean2", 672.5);
	mm1.setParameters(p);
	mm1.setSamples();
	mm1.setOffset( 673.5 );
	
	MixtureModel mm2;
	mm2.setParameters(mm1.getParameters());
	mm2.setSamples();
	mm2.setOffset( 673.5 );
	
	std::vector<Peak1D> v1;
	std::vector<Peak1D> v2;
	mm1.getSamples(v1);
	mm2.getSamples(v2);

	TEST_EQUAL(v1.size(),v2.size())
	ABORT_IF(v1.size()!=v2.size());
	for (Size i=0; i<v1.size(); ++i)
	{
		TEST_REAL_SIMILAR(v1[i].getPosition()[0],v2[i].getPosition()[0])
		TEST_REAL_SIMILAR(v1[i].getIntensity(),v2[i].getIntensity())
	}
}
END_SECTION

START_SECTION((void setSamples()))
{
	// already tested above, but well...
	MixtureModel mm1;
	Param p;
	p.setValue("mix", 0.5);
	p.setValue("statistics:variance1",0.8);
	p.setValue("statistics:variance2",0.8);
	p.setValue("statistics:mean1", 670.5);
	p.setValue("statistics:mean2", 672.5);
	mm1.setParameters(p);
	mm1.setSamples();
	mm1.setOffset( 673.5 );
	
	MixtureModel mm2;
	mm2.setParameters(mm1.getParameters());
	mm2.setSamples();
	mm2.setOffset( 673.5 );
	
	std::vector<Peak1D> v1;
	std::vector<Peak1D> v2;
	mm1.getSamples(v1);
	mm2.getSamples(v2);

	TEST_EQUAL(v1.size(),v2.size())
	ABORT_IF(v1.size()!=v2.size());
	for (Size i=0; i<v1.size(); ++i)
	{
		TEST_REAL_SIMILAR(v1[i].getPosition()[0],v2[i].getPosition()[0])
		TEST_REAL_SIMILAR(v1[i].getIntensity(),v2[i].getIntensity())
	}
}
END_SECTION

START_SECTION((CoordinateType getCenter() const ))
{
	TOLERANCE_ABSOLUTE(0.001)
	MixtureModel mm;
	
	Param p;
	p.setValue("mix", 0.5);
	p.setValue("statistics:variance1",0.8);
	p.setValue("statistics:variance2",0.8);
	p.setValue("statistics:mean1", 670.5);
	p.setValue("statistics:mean2", 672.5);
	mm.setParameters(p);
	mm.setOffset(680.0);
	TEST_REAL_SIMILAR(mm.getCenter(),1013.375)
}
END_SECTION

START_SECTION((static BaseModel<1>* create()))
{
	BaseModel<1>* ptr = MixtureModel::create();
	TEST_EQUAL(ptr->getName(), "MixtureModel")
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION((static const String getProductName()))
{
	TEST_EQUAL(MixtureModel::getProductName(),"MixtureModel")
	TEST_EQUAL(MixtureModel().getName(),"MixtureModel")
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



