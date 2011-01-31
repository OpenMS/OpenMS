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

    DoubleReal data[] = {0.014159276150167, 0.0153580168262124, 0.0167154017835855, 0.0182606503367424, 0.0200300235301256, 0.0220689419656992, 0.0244349148124456, 0.0272016022354364, 0.0304645504802465, 0.0343494676053524, 0.0390243455767632, 0.0447176694869995, 0.0517463870346546, 0.0605601146817207, 0.071813203394413, 0.086486391723156, 0.106100678443909, 0.133111327886581, 0.171673625707626, 0.229226112365723, 0.319999635219574, 0.473372250795364, 0.754716157913208, 1.31147396564484, 2.35293865203857, 3.19999647140503, 2.35293865203857, 1.31147396564484, 0.754716157913208, 0.473372250795364, 0.319999635219574, 0.229226112365723, 0.171673625707626, 0.143445029854774, 0.117309227585793, 0.098685584962368, 0.0851401463150978, 0.0751783773303032, 0.067852683365345, 0.0625507012009621, 0.0588765554130077, 0.0565830320119858, 0.0555333979427814, 0.0556822568178177, 0.057070653885603, 0.059834361076355, 0.064227856695652, 0.0706711858510971, 0.0798346847295761, 0.0927921533584595, 0.0971469879150391, 0.125290423631668, 0.167293235659599, 0.23354135453701, 0.345475375652313, 0.55080509185791, 0.957136690616608, 1.71721589565277, 2.33541369438171, 1.71721589565277, 0.957136690616608, 0.55080509185791, 0.345475375652313, 0.23354135453701, 0.167293235659599, 0.125290423631668, 0.0971469879150391, 0.0812376067042351, 0.0672447606921196, 0.0569006353616714, 0.0491030178964138, 0.0431458912789822, 0.0385639071464539, 0.0350443683564663, 0.0323757492005825, 0.030416963621974, 0.0290791746228933, 0.0283157657831907, 0.0281183794140816, 0.0285183973610401, 0.0295946262776852, 0.031489685177803, 0.0344405584037304, 0.0285008065402508, 0.0357564203441143, 0.0461150407791138, 0.0615748092532158, 0.0859584361314774, 0.12715744972229, 0.202732160687447, 0.352288663387299, 0.632047295570374, 0.859584331512451, 0.632047295570374, 0.352288663387299, 0.202732160687447, 0.12715744972229, 0.0859584361314774, 0.0615748092532158, 0.0461150407791138, 0.0366979315876961, 0.0295220259577036, 0.0243434868752956, 0.0205047205090523, 0.0175995640456676, 0.015367591753602, 0.0136368591338396, 0.0122914854437113, 0.0112526854500175, 0.010467441752553, 0.00990179926156998, 0.00953718367964029, 0.00936900451779366, 0.00940737128257751, 0.00968034751713276, 0.0102409441024065, 0.0111805601045489, 0.00885113701224327, 0.0114153074100614, 0.0152422152459621, 0.0212781336158514, 0.0314765274524689, 0.0501842759549618, 0.0872054621577263, 0.156456857919693, 0.212781324982643, 0.156456857919693, 0.0872054621577263, 0.0501842759549618, 0.0314765274524689, 0.0212781336158514, 0.0152422152459621, 0.0114153074100614, 0.00902740471065044, 0.00724627496674657, 0.00595893571153283, 0.00500249024480581, 0.00427625142037869, 0.00371557171456516, 0.00327765638940036, 0.00293352571316063, 0.00266329338774085, 0.00245333183556795, 0.0022945620585233, 0.00218146899715066, 0.00211164564825594, 0.0020857909694314, 0.00210822699591517, 0.00218814262188971, 0.00234206160530448, 0.00165709631983191, 0.00213715643621981, 0.00285362428985536, 0.00398365966975689, 0.00589298736304045, 0.00939542334526777, 0.0163264740258455, 0.0292916130274534, 0.0398365966975689, 0.0292916130274534, 0.0163264740258455, 0.00939542334526777, 0.00589298736304045, 0.00398365966975689, 0.00285362428985536, 0.00213715643621981, 0.00165709631983191, 0.00132084195502102, 0.00107666477560997, 0.00089399900753051, 0.000753909815102816, 0.000644188141450286, 0.000556687999051064, 0.000485812139231712, 0.000427614810178056, 0.000379251665435731, 0.000338631361955777, 0.000304189015878364, 0.000274735124548897, 0.000249352742685005, 0.000227325930609368, 0.000208089186344296, 0.000191191182238981};
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
