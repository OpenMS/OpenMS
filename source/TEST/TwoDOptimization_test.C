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
// $Maintainer: Alexandra Zerck $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/TwoDOptimization.h>
#include <OpenMS/FORMAT/MzMLFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(TwoDOptimization, "$Id$")
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TwoDOptimization* ptr = 0;
TwoDOptimization* nullPointer = 0;
START_SECTION((TwoDOptimization   ( )))
	ptr = new TwoDOptimization();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~TwoDOptimization()))
	delete ptr;
END_SECTION

START_SECTION((TwoDOptimization& operator=(const TwoDOptimization& opt)))
  TwoDOptimization opt_2d;
  struct OptimizationFunctions::PenaltyFactorsIntensity penalties;
  penalties.pos = 2;
  penalties.height = 3;
  penalties.lWidth = 4;
  penalties.rWidth = 5;
  opt_2d.setPenalties(penalties);
  opt_2d.setMaxIterations(10);
  opt_2d.setMaxAbsError(0.01);
  opt_2d.setMaxRelError(0.001);
  
  TwoDOptimization opt_2d_copy;
  opt_2d_copy = opt_2d;
  struct OptimizationFunctions::PenaltyFactorsIntensity penalties_copy = opt_2d_copy.getPenalties();
  unsigned int number = opt_2d_copy.getMaxIterations();
  DoubleReal abs_err = opt_2d_copy.getMaxAbsError();
  DoubleReal rel_err = opt_2d_copy.getMaxRelError();
  TEST_REAL_SIMILAR(penalties.pos,penalties_copy.pos)
  TEST_REAL_SIMILAR(penalties.lWidth,penalties_copy.lWidth)
  TEST_REAL_SIMILAR(penalties.rWidth,penalties_copy.rWidth)
  TEST_REAL_SIMILAR(penalties.height,penalties_copy.height)
    
  TEST_EQUAL(number == 10, true)
  TEST_REAL_SIMILAR(abs_err, 0.01)
  TEST_REAL_SIMILAR(rel_err, 0.001)
END_SECTION

START_SECTION((TwoDOptimization(const TwoDOptimization& opt)))
  TOLERANCE_ABSOLUTE(0.001)
  TwoDOptimization opt_2d;
  struct OptimizationFunctions::PenaltyFactorsIntensity penalties;
  opt_2d.setPenalties(penalties);
  opt_2d.setMaxIterations(10);
  opt_2d.setMaxAbsError(0.01);
  opt_2d.setMaxRelError(0.001);
  
  TwoDOptimization opt_2d_copy(opt_2d);
  struct OptimizationFunctions::PenaltyFactorsIntensity penalties_copy = opt_2d_copy.getPenalties();
  unsigned int number = opt_2d_copy.getMaxIterations();
  DoubleReal abs_err = opt_2d_copy.getMaxAbsError();
  DoubleReal rel_err = opt_2d_copy.getMaxRelError();
  TEST_REAL_SIMILAR(penalties.pos,penalties_copy.pos)
  TEST_REAL_SIMILAR(penalties.lWidth,penalties_copy.lWidth)
  TEST_REAL_SIMILAR(penalties.rWidth,penalties_copy.rWidth)
  TEST_REAL_SIMILAR(penalties.height,penalties_copy.height)
    
  TEST_EQUAL(number == 10, true)
  TEST_REAL_SIMILAR(abs_err, 0.01)
  TEST_REAL_SIMILAR(rel_err, 0.001)
END_SECTION


START_SECTION(( template <typename InputSpectrumIterator,typename OutputPeakType>  void optimize(InputSpectrumIterator first,InputSpectrumIterator last,MSExperiment< OutputPeakType >& ms_exp,bool real2D=true)))
  //******************************************************************
  //test exception with unequal number of scans
  {
  	MSExperiment<Peak1D> exp_in;
	  exp_in.resize(1);
	  MSExperiment<Peak1D>::const_iterator first1 = exp_in.begin();
	  MSExperiment<Peak1D>::const_iterator last1 = exp_in.end();
	  MSExperiment<> exp_out;
		TwoDOptimization opt1;
		TEST_EXCEPTION(Exception::IllegalArgument, opt1.optimize(first1,last1,exp_out));
  }

  //******************************************************************
  //test exception when meta data is missing
  {
  	MSExperiment<Peak1D> exp_in;
	  exp_in.resize(1);
	  MSExperiment<Peak1D>::const_iterator first1 = exp_in.begin();
	  MSExperiment<Peak1D>::const_iterator last1 = exp_in.end();
	  MSExperiment<> exp_out;
	  exp_out.resize(1);
		TwoDOptimization opt1;
		TEST_EXCEPTION(Exception::IllegalArgument, opt1.optimize(first1,last1,exp_out));
  }
  
  //******************************************************************
// test for 2D optimization
  TOLERANCE_ABSOLUTE(0.04)
  MSSpectrum<> peaks;
  peaks.getFloatDataArrays().resize(6);
  peaks.getFloatDataArrays()[1].setName("maximumIntensity");
  peaks.getFloatDataArrays()[1].push_back(700.); //intensity
  peaks.getFloatDataArrays()[3].setName("leftWidth");
	peaks.getFloatDataArrays()[3].push_back(12.5f); //left width
  peaks.getFloatDataArrays()[4].setName("rightWidth");
	peaks.getFloatDataArrays()[4].push_back(12.5f); //right width
  peaks.getFloatDataArrays()[5].setName("peakShape");
	peaks.getFloatDataArrays()[5].push_back(0); //shape
  peaks.getFloatDataArrays()[1].push_back(700.); //intensity
	peaks.getFloatDataArrays()[3].push_back(12.5f); //left width
	peaks.getFloatDataArrays()[4].push_back(12.5f); //right width
	peaks.getFloatDataArrays()[5].push_back(0); //shape
  MSSpectrum<> peaks2;
  peaks2.getFloatDataArrays().resize(6);
  peaks2.getFloatDataArrays()[1].setName("maximumIntensity");
  peaks2.getFloatDataArrays()[1].push_back(700.); //intensity
  peaks2.getFloatDataArrays()[3].setName("leftWidth");
	peaks2.getFloatDataArrays()[3].push_back(12.5f); //left width
  peaks2.getFloatDataArrays()[4].setName("rightWidth");
	peaks2.getFloatDataArrays()[4].push_back(12.5f); //right width
  peaks2.getFloatDataArrays()[5].setName("peakShape");
	peaks2.getFloatDataArrays()[5].push_back(0); //shape
  peaks2.getFloatDataArrays()[1].push_back(700.); //intensity
	peaks2.getFloatDataArrays()[3].push_back(12.5f); //left width
	peaks2.getFloatDataArrays()[4].push_back(12.5f); //right width
	peaks2.getFloatDataArrays()[5].push_back(0); //shape
	
	Peak1D peak;
  PeakShape peak_shape,peak_shape2;
  peak.setMZ(500.);
  peak.setIntensity(171.69f);
  peak_shape.mz_position = 500;
  peak_shape.left_width = 12.5;
  peak_shape.right_width = 12.5;
  peak_shape.area = 171.69;
  peak_shape.height = 700;
  peak_shape.type = PeakShape::LORENTZ_PEAK;  
	peaks.push_back(peak);
  peak.setMZ(501.);
  peak.setIntensity(171.69f);
  peak_shape2.mz_position = 501;
  peak_shape2.left_width = 12.5;
  peak_shape2.right_width = 12.5;
  peak_shape2.area = 171.69;
  peak_shape2.height = 700;
  peak_shape2.type = PeakShape::LORENTZ_PEAK;  
	peaks.push_back(peak);

  MSExperiment<> ms_exp;
	ms_exp.push_back(peaks);
	ms_exp.begin()->setRT(100);
			
  float origin = 499;
  float spacing = 0.1f;

	MSSpectrum<Peak1D >	 raw_spec;
  for (Size i = 0; i < 20 ;++i)
  {
		Peak1D data_point;
		data_point.setMZ(origin +i*spacing);
		data_point.setIntensity(peak_shape(origin +i*spacing)+peak_shape2(origin +i*spacing));
    raw_spec.push_back(data_point);
  }
  peak.setMZ(500.02);
  peak.setIntensity(171.69f);
  peak_shape.mz_position = 500;
  peak_shape.left_width = 12.5;
  peak_shape.right_width = 12.5;
  peak_shape.area = 171.69;
  peak_shape.height = 700;
  peak_shape.type = PeakShape::LORENTZ_PEAK;  
	peaks2.push_back(peak);
  peak.setMZ(501);
  peak.setIntensity(171.69f);
  peak_shape2.mz_position = 501;
  peak_shape2.left_width = 12.5;
  peak_shape2.right_width = 12.5;
  peak_shape2.area = 171.69;
  peak_shape2.height = 700;
  peak_shape2.type = PeakShape::LORENTZ_PEAK;  
	peaks2.push_back(peak);


	ms_exp.push_back(peaks2);
 (ms_exp.begin()+1)->setRT(101);

  MSSpectrum<Peak1D >	 raw_spec2;
  for (Size i = 0; i < 20 ;++i)
  {
		Peak1D data_point;
		data_point.setMZ(origin +i*spacing);
		data_point.setIntensity(peak_shape(origin +i*spacing)+peak_shape2(origin +i*spacing));
    raw_spec2.push_back(data_point);
  }

  MSExperiment<Peak1D > raw_exp;
  raw_exp.push_back(raw_spec);
  raw_exp.push_back(raw_spec2);
	raw_exp.begin()->setRT(100);
  (raw_exp.begin()+1)->setRT(101);
  String file = OPENMS_GET_TEST_DATA_PATH("TwoDOptimization.xml");	
  Param param;
	param.load(file);
  MSExperiment<Peak1D >::const_iterator first,last;
  first = raw_exp.begin();
  last = raw_exp.end();
  TwoDOptimization opt_2d;
 	opt_2d.setParameters(param);
  opt_2d.optimize(first,last,ms_exp,true);
  TEST_REAL_SIMILAR(ms_exp[0][0].getMZ(),500)
 	TEST_REAL_SIMILAR(ms_exp[0].getFloatDataArrays()[3][0],12.5)
 	TEST_REAL_SIMILAR(ms_exp[0].getFloatDataArrays()[4][0],12.5)
 	TEST_REAL_SIMILAR(ms_exp[0][0].getIntensity(),171.69)
	TEST_REAL_SIMILAR(ms_exp[0][1].getMZ(),501)
 	TEST_REAL_SIMILAR(ms_exp[0].getFloatDataArrays()[3][1],12.5)
 	TEST_REAL_SIMILAR(ms_exp[0].getFloatDataArrays()[4][1],12.5)
 	TEST_REAL_SIMILAR(ms_exp[0][1].getIntensity(),171.69)
	TEST_REAL_SIMILAR(ms_exp[1][0].getMZ(),500)
 	TEST_REAL_SIMILAR(ms_exp[1].getFloatDataArrays()[3][0],12.5)
 	TEST_REAL_SIMILAR(ms_exp[1].getFloatDataArrays()[4][0],12.5)
 	TEST_REAL_SIMILAR(ms_exp[1][0].getIntensity(),171.69)
	TEST_REAL_SIMILAR(ms_exp[1][1].getMZ(),501)
 	TEST_REAL_SIMILAR(ms_exp[1].getFloatDataArrays()[3][1],12.5)
 	TEST_REAL_SIMILAR(ms_exp[1].getFloatDataArrays()[4][1],12.5)
 	TEST_REAL_SIMILAR(ms_exp[1][1].getIntensity(),171.69)

  
  peaks2.clear(false);
  peaks2.getFloatDataArrays()[3][0]=12.5; //left width
  peaks2.getFloatDataArrays()[4][0]=12.5; //right width
  peaks2.getFloatDataArrays()[3][1]=12.5; //left width
	peaks2.getFloatDataArrays()[4][1]=12.5; //right width
  peaks2.getFloatDataArrays()[1][0]=700.; //intensity
  peaks2.getFloatDataArrays()[1][1]=800.; //intensity


  peak.setMZ(900);
  peak.setIntensity(105.79f);
  peak_shape.mz_position = 900;
  peak_shape.left_width = 12.5;
  peak_shape.right_width = 12.5;
  peak_shape.area = 171.69;
  peak_shape.height = 700;
  peak_shape.type = PeakShape::LORENTZ_PEAK;  
	peaks2.push_back(peak);
  peak.setMZ(901);
  peak.setIntensity(105.79f);
  peak_shape2.mz_position = 901;
  peak_shape2.left_width = 12.5;
  peak_shape2.right_width = 12.5;
  peak_shape2.area = 171.69;
  peak_shape2.height = 700;
  peak_shape2.type = PeakShape::LORENTZ_PEAK;  
	peaks2.push_back(peak);
	ms_exp.push_back(peaks2);
  (ms_exp.begin()+2)->setRT(102);

  raw_spec2.clear(true);
  origin = 899;
  for (Size i = 0; i < 20 ;++i)
  {
		Peak1D data_point;
		data_point.setMZ(origin +i*spacing);
		data_point.setIntensity(peak_shape(origin +i*spacing)+peak_shape2(origin +i*spacing));
    raw_spec2.push_back(data_point);
  }


  raw_exp.push_back(raw_spec2);
  (raw_exp.begin()+2)->setRT(102);
  first = raw_exp.begin();
  last = raw_exp.end();
  TwoDOptimization opt_1d;
 	opt_1d.setParameters(param);
  opt_1d.optimize(first,last,ms_exp,false); // test 1D optimization
  TEST_REAL_SIMILAR(ms_exp[0][0].getMZ(),500)
 	TEST_REAL_SIMILAR(ms_exp[0].getFloatDataArrays()[3][0],12.5)
 	TEST_REAL_SIMILAR(ms_exp[0].getFloatDataArrays()[4][0],12.5)
 	TEST_REAL_SIMILAR(ms_exp[0][0].getIntensity(),171.69)
	TEST_REAL_SIMILAR(ms_exp[0][1].getMZ(),501)
 	TEST_REAL_SIMILAR(ms_exp[0].getFloatDataArrays()[3][1],12.5)
 	TEST_REAL_SIMILAR(ms_exp[0].getFloatDataArrays()[4][1],12.5)
 	TEST_REAL_SIMILAR(ms_exp[0][1].getIntensity(),171.69)
	TEST_REAL_SIMILAR(ms_exp[1][0].getMZ(),500)
 	TEST_REAL_SIMILAR(ms_exp[1].getFloatDataArrays()[3][0],12.5)
 	TEST_REAL_SIMILAR(ms_exp[1].getFloatDataArrays()[4][0],12.5)
 	TEST_REAL_SIMILAR(ms_exp[1][0].getIntensity(),171.69)
	TEST_REAL_SIMILAR(ms_exp[1][1].getMZ(),501)
 	TEST_REAL_SIMILAR(ms_exp[1].getFloatDataArrays()[3][1],12.5)
 	TEST_REAL_SIMILAR(ms_exp[1].getFloatDataArrays()[4][1],12.5)
 	TEST_REAL_SIMILAR(ms_exp[1][1].getIntensity(),171.69)

  TEST_REAL_SIMILAR(ms_exp[2][0].getMZ(),900)
 	TEST_REAL_SIMILAR(ms_exp[2].getFloatDataArrays()[3][0],12.5)
 	TEST_REAL_SIMILAR(ms_exp[2].getFloatDataArrays()[4][0],12.5)
 	TEST_REAL_SIMILAR(ms_exp[2][0].getIntensity(),171.69)
	TEST_REAL_SIMILAR(ms_exp[2][1].getMZ(),901)
 	TEST_REAL_SIMILAR(ms_exp[2].getFloatDataArrays()[3][1],12.5)
 	TEST_REAL_SIMILAR(ms_exp[2].getFloatDataArrays()[4][1],12.5)
 	TEST_REAL_SIMILAR(ms_exp[2][1].getIntensity(),171.69)

END_SECTION

START_SECTION((void setMaxAbsError(DoubleReal eps_abs)))
  TOLERANCE_ABSOLUTE(0.0001)
  DoubleReal abs_err = 0.01;
   
  TwoDOptimization opt_2d;
  opt_2d.setMaxAbsError(abs_err);
    
 	TEST_REAL_SIMILAR(abs_err, opt_2d.getMaxAbsError())
END_SECTION

START_SECTION((void setMaxRelError(DoubleReal eps_rel)))
  TOLERANCE_ABSOLUTE(0.0001)
  DoubleReal rel_err = 0.01;
   
  TwoDOptimization opt_2d;
  opt_2d.setMaxRelError(rel_err);
    
 	TEST_REAL_SIMILAR(rel_err, opt_2d.getMaxRelError())
END_SECTION
	
START_SECTION((DoubleReal getMaxAbsError() const))
  TOLERANCE_ABSOLUTE(0.0001)
  DoubleReal abs_err = 0.01;
   
  TwoDOptimization opt_2d;
  opt_2d.setMaxAbsError(abs_err);
    
 	TEST_REAL_SIMILAR(abs_err, opt_2d.getMaxAbsError())
END_SECTION

START_SECTION((DoubleReal getMaxRelError() const))
  TOLERANCE_ABSOLUTE(0.0001)
  DoubleReal rel_err = 0.01;
   
  TwoDOptimization opt_2d;
  opt_2d.setMaxRelError(rel_err);
    
 	TEST_REAL_SIMILAR(rel_err, opt_2d.getMaxRelError())
END_SECTION
	
START_SECTION((void setMaxPeakDistance(DoubleReal max_peak_distance)))
  TOLERANCE_ABSOLUTE(0.0001)
  DoubleReal max_peak_distance = 0.01;
   
  TwoDOptimization opt_2d;
  opt_2d.setMaxPeakDistance(max_peak_distance);
    
 	TEST_REAL_SIMILAR(max_peak_distance, opt_2d.getMaxPeakDistance())
END_SECTION

START_SECTION((DoubleReal getMaxPeakDistance() const))
  TOLERANCE_ABSOLUTE(0.0001)
  DoubleReal max_peak_distance = 0.01;
   
  TwoDOptimization opt_2d;
  opt_2d.setMaxPeakDistance(max_peak_distance);
    
 	TEST_REAL_SIMILAR(max_peak_distance, opt_2d.getMaxPeakDistance())
END_SECTION



START_SECTION((void setMZTolerance(DoubleReal tolerance_mz)))
  TOLERANCE_ABSOLUTE(0.0001)
  DoubleReal mz_tolerance = 0.01;
   
  TwoDOptimization opt_2d;
  opt_2d.setMZTolerance(mz_tolerance);
    
 	TEST_REAL_SIMILAR(mz_tolerance, opt_2d.getMZTolerance())
END_SECTION

START_SECTION((DoubleReal getMZTolerance() const))
  TOLERANCE_ABSOLUTE(0.0001)
  DoubleReal mz_tolerance = 0.01;
   
  TwoDOptimization opt_2d;
  opt_2d.setMZTolerance(mz_tolerance);
    
 	TEST_REAL_SIMILAR(mz_tolerance, opt_2d.getMZTolerance())
END_SECTION

START_SECTION((void setMaxIterations(UInt max_iteration)))
	UInt number = 20;
   
  TwoDOptimization opt_2d;
  opt_2d.setMaxIterations(number);
    
 	TEST_EQUAL(number == opt_2d.getMaxIterations(), true)
END_SECTION

START_SECTION((UInt getMaxIterations() const))
  UInt number = 20;
   
  TwoDOptimization opt_2d;
  opt_2d.setMaxIterations(number);
    
 	TEST_EQUAL(number == opt_2d.getMaxIterations(), true)
END_SECTION

	
START_SECTION((void setPenalties(OptimizationFunctions::PenaltyFactorsIntensity& penalties)))
  TOLERANCE_ABSOLUTE(0.0001)
  struct OptimizationFunctions::PenaltyFactorsIntensity penalties;
  penalties.pos = 0;
  penalties.lWidth = 1;
  penalties.rWidth = 2;
  penalties.height = 3;

  TwoDOptimization opt_2d;
  opt_2d.setPenalties(penalties);
  TEST_REAL_SIMILAR(penalties.pos,opt_2d.getPenalties().pos)
  TEST_REAL_SIMILAR(penalties.lWidth,opt_2d.getPenalties().lWidth)
  TEST_REAL_SIMILAR(penalties.rWidth,opt_2d.getPenalties().rWidth)
	TEST_REAL_SIMILAR(penalties.height,opt_2d.getPenalties().height)
END_SECTION

START_SECTION((const OptimizationFunctions::PenaltyFactorsIntensity& getPenalties() const))
  TOLERANCE_ABSOLUTE(0.0001)
  struct OptimizationFunctions::PenaltyFactorsIntensity penalties;
  penalties.pos = 0;
  penalties.lWidth = 1;
  penalties.rWidth = 2;
  penalties.height = 3;

  TwoDOptimization opt_2d;
  opt_2d.setPenalties(penalties);
  TEST_REAL_SIMILAR(penalties.pos,opt_2d.getPenalties().pos)
  TEST_REAL_SIMILAR(penalties.lWidth,opt_2d.getPenalties().lWidth)
  TEST_REAL_SIMILAR(penalties.rWidth,opt_2d.getPenalties().rWidth)
	TEST_REAL_SIMILAR(penalties.height,opt_2d.getPenalties().height)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
