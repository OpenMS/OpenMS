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
// $Maintainer: Clemens Groepl $
// $Authors: Clemens Groepl, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeModel.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <sstream>


///////////////////////////

START_TEST(IsotopeModel, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using std::stringstream;

// default ctor
IsotopeModel* ptr = 0;
START_SECTION((IsotopeModel()))
	ptr = new IsotopeModel();
  TEST_EQUAL(ptr->getName(), "IsotopeModel")
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

// destructor
START_SECTION((virtual ~IsotopeModel()))
	delete ptr;
END_SECTION

START_SECTION(static BaseModel<1>* create())
	BaseModel<1>* ptr = IsotopeModel::create();
	TEST_EQUAL(ptr->getName(), "IsotopeModel")
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((static const String getProductName()))
	TEST_EQUAL(IsotopeModel::getProductName(),"IsotopeModel")
	TEST_EQUAL(IsotopeModel().getName(),"IsotopeModel")
END_SECTION

// assignment operator
START_SECTION((virtual IsotopeModel& operator=(const IsotopeModel &source)))
	IsotopeModel im1;
	
	Param tmp;
	tmp.setValue("charge", 3);
	tmp.setValue("isotope:mode:GaussianSD",0.8);
	tmp.setValue("statistics:mean", 670.5);
	im1.setParameters(tmp);

  IsotopeModel im2;
  im2 = im1;

  IsotopeModel im3;
	im3.setParameters(tmp);

  im1 = IsotopeModel();
	TEST_EQUAL(im3.getParameters(), im2.getParameters())
END_SECTION

// copy ctor
START_SECTION((IsotopeModel(const IsotopeModel& source)))
	IsotopeModel im1;
	
	Param tmp;
	tmp.setValue("charge", 3);
	tmp.setValue("isotope:mode:GaussianSD",0.8);
	tmp.setValue("statistics:mean", 670.5);
	im1.setParameters(tmp);

	IsotopeModel im2(im1);
  IsotopeModel im3;
	im3.setParameters(tmp);

  im1 = IsotopeModel();
	TEST_EQUAL(im3.getParameters(), im2.getParameters())
END_SECTION

START_SECTION([EXTRA] DefaultParamHandler::setParameters(...))
	TOLERANCE_ABSOLUTE(0.001)
	IsotopeModel im1;
	Param tmp;
	tmp.setValue("charge", 3);
	tmp.setValue("isotope:mode:GaussianSD",0.8);
	tmp.setValue("statistics:mean", 670.5);
	im1.setParameters(tmp);

	IsotopeModel im2;
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
	IsotopeModel im1;
	TEST_EQUAL(im1.getCharge() == 1, true)		// default charge is 1
END_SECTION

START_SECTION( CoordinateType getCenter() const )
	// can only reliably be tested after fitting, only sanity check here
	IsotopeModel im1;
	TEST_EQUAL(im1.getCenter() == 0, true)
END_SECTION

START_SECTION( void setSamples(const EmpiricalFormula &formula) )
  IsotopeModel im1;
  Param tmp;
  EmpiricalFormula ef("C66H129O3");
  tmp.setValue("statistics:mean", ef.getAverageWeight() / 1);		
  tmp.setValue("interpolation_step", 0.03);
  tmp.setValue("charge", 1);
  tmp.setValue("isotope:maximum",100);
  tmp.setValue("isotope:mode:mode","Gaussian");
  tmp.setValue("isotope:mode:GaussianSD",0.15);

  im1.setParameters(tmp);
  im1.setSamples(ef);

  {
    DoubleReal data[] = {0.000429512, 0.00093697, 0.00196383, 0.00395466, 0.00765145, 0.0142235, 0.0254037, 0.043593, 0.0718726, 0.113852, 0.173278, 0.253381, 0.355987, 0.480533, 0.623218, 0.776577, 0.929731, 1.06945, 1.18192, 1.25501, 1.28036, 1.25501, 1.18192, 1.06945, 0.929731, 0.776577, 0.623218, 0.480533, 0.355987, 0.253381, 0.173278, 0.113852, 0.0718726, 0.0439064, 0.0260875, 0.0156567, 0.0105376, 0.00953883, 0.0123444, 0.019477, 0.0318149, 0.0524539, 0.0830909, 0.126461, 0.184922, 0.259806, 0.350701, 0.454835, 0.566759, 0.678534, 0.7805, 0.862586, 0.915925, 0.934428, 0.915925, 0.862586, 0.7805, 0.678534, 0.566759, 0.454835, 0.350701, 0.259806, 0.184922, 0.126461, 0.0830909, 0.0524539, 0.0318149, 0.0186555, 0.0106322, 0.00611169, 0.00394849, 0.00348857, 0.00450454, 0.00682396, 0.01171, 0.0193065, 0.0305829, 0.0465459, 0.0680634, 0.0956255, 0.129081, 0.167409, 0.208604, 0.249745, 0.287275, 0.317488, 0.33712, 0.34393, 0.33712, 0.317488, 0.287275, 0.249745, 0.208604, 0.167409, 0.129081, 0.0956255, 0.0680634, 0.0465459, 0.0305829, 0.0193065, 0.0117385, 0.00688626, 0.00395131, 0.0023183, 0.00157108, 0.00147331, 0.00194089, 0.00289868, 0.00477912, 0.00757047, 0.011522, 0.0168484, 0.0236711, 0.0319527, 0.0414404, 0.0516379, 0.0618218, 0.071112, 0.0785909, 0.0834507, 0.0851365, 0.0834507, 0.0785909, 0.071112, 0.0618218, 0.0516379, 0.0414404, 0.0319527, 0.0236711, 0.0168484, 0.011522, 0.00757047, 0.00477912, 0.00290403, 0.00170087, 0.000970228, 0.000558009, 0.000358215, 0.000307651, 0.000378553, 0.000542686, 0.000894739, 0.00141733, 0.00215713, 0.00315433, 0.00443167, 0.00598213, 0.0077584, 0.00966757, 0.0115742, 0.0133135, 0.0147137, 0.0156235, 0.0159391, 0.0156235, 0.0147137, 0.0133135, 0.0115742, 0.00966757, 0.0077584, 0.00598213, 0.00443167, 0.00315433, 0.00215713, 0.00141733, 0.000894739, 0.000542686, 0.00031625, 0.000177068, 9.52526e-005, 4.92314e-005, 2.44476e-005, 1.16643e-005};
    int size = sizeof( data ) / sizeof( data[0] );
    std::vector<DoubleReal> dpa2( data, &data[ size ] );

    std::vector<Peak1D> dpa1;
    im1.getSamples(dpa1);

    TEST_EQUAL(dpa1.size(),dpa2.size())
      ABORT_IF(dpa1.size()!=dpa2.size());
    for (Size i=0; i<dpa1.size(); ++i)
    {
      TEST_REAL_SIMILAR(dpa1[i].getIntensity(),dpa2[i])
    }
  }

  {
    // lorentzian
    tmp.setValue("isotope:mode:mode","Lorentzian");
    tmp.setValue("isotope:mode:LorentzFWHM",0.05);

    im1.setParameters(tmp);
    im1.setSamples(ef);

    DoubleReal data[] = {0.134448975324631, 0.171693190932274, 0.226403713226318, 0.311002969741821, 0.450473368167877, 0.699134647846222, 1.18097066879272, 2.13150811195374, 3.36122441291809, 3.01351141929626, 1.74783658981323, 0.981930673122406, 0.598574221134186, 0.39543816447258, 0.278317928314209, 0.205627843737602, 0.157746985554695, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0981232225894928, 0.125304698944092, 0.16523340344429, 0.226975426077843, 0.328763365745544, 0.510240733623505, 0.861893117427826, 1.55561196804047, 2.45308041572571, 2.19931364059448, 1.27560186386108, 0.716630280017853, 0.43684995174408, 0.288597702980042, 0.203121319413185, 0.150070801377296, 0.115126520395279, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0361157357692719, 0.0461202934384346, 0.060816653072834, 0.0835417360067368, 0.121006332337856, 0.187801823019981, 0.317232817411423, 0.572566568851471, 0.902893424034119, 0.809490621089935, 0.469504565000534, 0.26376661658287, 0.160789236426353, 0.106222756206989, 0.074761874973774, 0.0552358329296112, 0.0423740595579147, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00894008204340935, 0.0114166075363755, 0.0150545425713062, 0.0206799041479826, 0.0299538820981979, 0.0464884266257286, 0.0785277485847473, 0.141733005642891, 0.223502054810524, 0.200381144881248, 0.116221062839031, 0.0652927309274673, 0.0398017354309559, 0.0262943580746651, 0.0185065399855375, 0.0136730661615729, 0.010489265434444, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00167374848388135, 0.00213739974424243, 0.00281848828308284, 0.00387166021391749, 0.00560792023316026, 0.00870349165052176, 0.0147018451243639, 0.0265350360423326, 0.0418437123298645, 0.0375150516629219, 0.0217587295919657, 0.0122240055352449, 0.0074516199529171, 0.00492278952151537, 0.00346476584672928, 0.00255985069088638, 0.00196378421969712};
    int size = sizeof( data ) / sizeof( data[0] );
    std::vector<DoubleReal> dpa2( data, &data[ size ] );

    std::vector<Peak1D> dpa1;
    im1.getSamples(dpa1);

    TEST_EQUAL(dpa1.size(),dpa2.size())
    ABORT_IF(dpa1.size()!=dpa2.size());
    for (Size i=0; i<dpa1.size(); ++i)
    {
      TEST_REAL_SIMILAR(dpa1[i].getIntensity(),dpa2[i])
    }
  }

END_SECTION

START_SECTION( void setOffset(CoordinateType offset) )
	TOLERANCE_ABSOLUTE(0.1)
	IsotopeModel im1;
	Param tmp;
	tmp.setValue("charge", 3);
	tmp.setValue("isotope:mode:GaussianSD",0.8);
	tmp.setValue("statistics:mean", 670.5);
	im1.setParameters(tmp);
	im1.setOffset( 673.5 );
	
	IsotopeModel im2;
	im2.setParameters(im1.getParameters());
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
	TOLERANCE_ABSOLUTE(0.1)
	IsotopeModel im1;
	Param tmp;
	tmp.setValue("charge", 3);
	tmp.setValue("isotope:mode:GaussianSD",0.8);
	tmp.setValue("statistics:mean", 670.5);
	im1.setParameters(tmp);
	im1.setOffset( 673.5 );
	
	IsotopeModel im2;
	im2.setParameters(im1.getParameters());
	im2.setOffset( im1.getOffset() );
	
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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
