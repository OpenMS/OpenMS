// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Alexandra Zerck $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/TwoDOptimization.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(TwoDOptimization, "$Id$")
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TwoDOptimization* ptr = 0;
CHECK((TwoDOptimization   ( )))
	ptr = new TwoDOptimization();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~TwoDOptimization()))
	delete ptr;
RESULT

CHECK((TwoDOptimization& operator=(const TwoDOptimization& opt)))
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
  double abs_err = opt_2d_copy.getMaxAbsError();
  double rel_err = opt_2d_copy.getMaxRelError();
  TEST_REAL_EQUAL(penalties.pos,penalties_copy.pos)
  TEST_REAL_EQUAL(penalties.lWidth,penalties_copy.lWidth)
  TEST_REAL_EQUAL(penalties.rWidth,penalties_copy.rWidth)
  TEST_REAL_EQUAL(penalties.height,penalties_copy.height)
    
  TEST_EQUAL(number == 10, true)
  TEST_REAL_EQUAL(abs_err, 0.01)
  TEST_REAL_EQUAL(rel_err, 0.001)
RESULT

CHECK((TwoDOptimization(const TwoDOptimization& opt)))
  PRECISION(0.0001)
  TwoDOptimization opt_2d;
  struct OptimizationFunctions::PenaltyFactorsIntensity penalties;
  opt_2d.setPenalties(penalties);
  opt_2d.setMaxIterations(10);
  opt_2d.setMaxAbsError(0.01);
  opt_2d.setMaxRelError(0.001);
  
  TwoDOptimization opt_2d_copy(opt_2d);
  struct OptimizationFunctions::PenaltyFactorsIntensity penalties_copy = opt_2d_copy.getPenalties();
  unsigned int number = opt_2d_copy.getMaxIterations();
  double abs_err = opt_2d_copy.getMaxAbsError();
  double rel_err = opt_2d_copy.getMaxRelError();
  TEST_REAL_EQUAL(penalties.pos,penalties_copy.pos)
  TEST_REAL_EQUAL(penalties.lWidth,penalties_copy.lWidth)
  TEST_REAL_EQUAL(penalties.rWidth,penalties_copy.rWidth)
  TEST_REAL_EQUAL(penalties.height,penalties_copy.height)
    
  TEST_EQUAL(number == 10, true)
  TEST_REAL_EQUAL(abs_err, 0.01)
  TEST_REAL_EQUAL(rel_err, 0.001)
RESULT


CHECK(( template <typename InputSpectrumIterator,typename OutputPeakType>  void optimize(InputSpectrumIterator& first,InputSpectrumIterator& last,MSExperiment< OutputPeakType >& ms_exp,bool real2D=true)))
  //******************************************************************
  //test exception with unequal number of scans
  {
  	MSExperiment<RawDataPoint1D> exp_in;
	  exp_in.resize(1);
	  MSExperiment<RawDataPoint1D>::const_iterator first1 = exp_in.begin();
	  MSExperiment<RawDataPoint1D>::const_iterator last1 = exp_in.end();
	  MSExperiment<> exp_out;
		TwoDOptimization opt1;
		TEST_EXCEPTION(Exception::IllegalArgument, opt1.optimize(first1,last1,exp_out));
  }

  //******************************************************************
  //test exception when meta data is missing
  {
  	MSExperiment<RawDataPoint1D> exp_in;
	  exp_in.resize(1);
	  MSExperiment<RawDataPoint1D>::const_iterator first1 = exp_in.begin();
	  MSExperiment<RawDataPoint1D>::const_iterator last1 = exp_in.end();
	  MSExperiment<> exp_out;
	  exp_out.resize(1);
		TwoDOptimization opt1;
		TEST_EXCEPTION(Exception::IllegalArgument, opt1.optimize(first1,last1,exp_out));
  }
  
  //******************************************************************
  //actual test
  MSSpectrum<> peaks;
  peaks.getMetaDataArrays().resize(6);
  peaks.getMetaDataArrays()[1].setName("area");
  peaks.getMetaDataArrays()[1].push_back(100.); //area
  peaks.getMetaDataArrays()[3].setName("leftWidth");
	peaks.getMetaDataArrays()[3].push_back(2.5); //left width
  peaks.getMetaDataArrays()[4].setName("rightWidth");
	peaks.getMetaDataArrays()[4].push_back(2.6); //right width
  peaks.getMetaDataArrays()[5].setName("peakShape");
	peaks.getMetaDataArrays()[5].push_back(0); //shape
	
	Peak1D peak;
	PeakShape peak_shape;
  peak.setMZ(500);
  peak.setIntensity(400);
  peak_shape.mz_position = 500;
  peak_shape.left_width = 2.5;
  peak_shape.right_width = 2.5;
  peak_shape.area = 100;
  peak_shape.height = 400;
  peak_shape.type = PeakShape::LORENTZ_PEAK;  
	peaks.push_back(peak);
	MSExperiment<> ms_exp;
	ms_exp.push_back(peaks);
	ms_exp.begin()->setRT(100);
			
  float origin = 499;
  float spacing = 0.1;

	MSSpectrum<RawDataPoint1D >	 raw_spec;
  for (unsigned int i = 0; i < 20 ;++i)
  {
		RawDataPoint1D data_point;
		data_point.setMZ(origin +i*spacing);
		data_point.setIntensity(peak_shape(origin +i*spacing));
    raw_spec.push_back(data_point);
  }
  MSExperiment<RawDataPoint1D > raw_exp;
  raw_exp.push_back(raw_spec);
	raw_exp.begin()->setRT(100);
  String file = "data/TwoDOptimization.xml";	
  Param param;
	param.load(file);

 	TwoDOptimization opt_2d;
 	opt_2d.setParameters(param);
  MSExperiment<RawDataPoint1D >::const_iterator first,last;
  first = raw_exp.begin();
  last = raw_exp.end();
 	opt_2d.optimize(first,last,ms_exp);
 	TEST_REAL_EQUAL(peak_shape.mz_position,500)
 	TEST_REAL_EQUAL(peak_shape.left_width,2.5)
 	TEST_REAL_EQUAL(peak_shape.right_width,2.5)
 	TEST_REAL_EQUAL(peak_shape.area,100)
 	TEST_REAL_EQUAL(peak_shape.height,400)
RESULT

CHECK((void setMaxAbsError(double eps_abs)))
  PRECISION(0.0001)
  double abs_err = 0.01;
   
  TwoDOptimization opt_2d;
  opt_2d.setMaxAbsError(abs_err);
    
 	TEST_REAL_EQUAL(abs_err, opt_2d.getMaxAbsError())
RESULT

CHECK((void setMaxRelError(double eps_rel)))
  PRECISION(0.0001)
  double rel_err = 0.01;
   
  TwoDOptimization opt_2d;
  opt_2d.setMaxRelError(rel_err);
    
 	TEST_REAL_EQUAL(rel_err, opt_2d.getMaxRelError())
RESULT
	
CHECK((DoubleReal getMaxAbsError() const))
  PRECISION(0.0001)
  double abs_err = 0.01;
   
  TwoDOptimization opt_2d;
  opt_2d.setMaxAbsError(abs_err);
    
 	TEST_REAL_EQUAL(abs_err, opt_2d.getMaxAbsError())
RESULT

CHECK((DoubleReal getMaxRelError() const))
  PRECISION(0.0001)
  double rel_err = 0.01;
   
  TwoDOptimization opt_2d;
  opt_2d.setMaxRelError(rel_err);
    
 	TEST_REAL_EQUAL(rel_err, opt_2d.getMaxRelError())
RESULT
	
CHECK((void setMaxPeakDistance(double max_peak_distance)))
  PRECISION(0.0001)
  double max_peak_distance = 0.01;
   
  TwoDOptimization opt_2d;
  opt_2d.setMaxPeakDistance(max_peak_distance);
    
 	TEST_REAL_EQUAL(max_peak_distance, opt_2d.getMaxPeakDistance())
RESULT

CHECK((DoubleReal getMaxPeakDistance() const))
  PRECISION(0.0001)
  double max_peak_distance = 0.01;
   
  TwoDOptimization opt_2d;
  opt_2d.setMaxPeakDistance(max_peak_distance);
    
 	TEST_REAL_EQUAL(max_peak_distance, opt_2d.getMaxPeakDistance())
RESULT



CHECK((void setMZTolerance(double tolerance_mz)))
  PRECISION(0.0001)
  double mz_tolerance = 0.01;
   
  TwoDOptimization opt_2d;
  opt_2d.setMZTolerance(mz_tolerance);
    
 	TEST_REAL_EQUAL(mz_tolerance, opt_2d.getMZTolerance())
RESULT

CHECK((DoubleReal getMZTolerance() const))
  PRECISION(0.0001)
  double mz_tolerance = 0.01;
   
  TwoDOptimization opt_2d;
  opt_2d.setMZTolerance(mz_tolerance);
    
 	TEST_REAL_EQUAL(mz_tolerance, opt_2d.getMZTolerance())
RESULT

CHECK((void setMaxIterations(int max_iteration)))
  int number = 20;
   
  TwoDOptimization opt_2d;
  opt_2d.setMaxIterations(number);
    
 	TEST_EQUAL(number == opt_2d.getMaxIterations(), true)
RESULT

CHECK((Int getMaxIterations() const))
  int number = 20;
   
  TwoDOptimization opt_2d;
  opt_2d.setMaxIterations(number);
    
 	TEST_EQUAL(number == opt_2d.getMaxIterations(), true)
RESULT

	
CHECK((void setPenalties(OptimizationFunctions::PenaltyFactorsIntensity& penalties)))
  PRECISION(0.0001)
  struct OptimizationFunctions::PenaltyFactorsIntensity penalties;
  penalties.pos = 0;
  penalties.lWidth = 1;
  penalties.rWidth = 2;
  penalties.height = 3;

  TwoDOptimization opt_2d;
  opt_2d.setPenalties(penalties);
  TEST_REAL_EQUAL(penalties.pos,opt_2d.getPenalties().pos)
  TEST_REAL_EQUAL(penalties.lWidth,opt_2d.getPenalties().lWidth)
  TEST_REAL_EQUAL(penalties.rWidth,opt_2d.getPenalties().rWidth)
	TEST_REAL_EQUAL(penalties.height,opt_2d.getPenalties().height)
RESULT

CHECK((const OptimizationFunctions::PenaltyFactorsIntensity& getPenalties() const))
  PRECISION(0.0001)
  struct OptimizationFunctions::PenaltyFactorsIntensity penalties;
  penalties.pos = 0;
  penalties.lWidth = 1;
  penalties.rWidth = 2;
  penalties.height = 3;

  TwoDOptimization opt_2d;
  opt_2d.setPenalties(penalties);
  TEST_REAL_EQUAL(penalties.pos,opt_2d.getPenalties().pos)
  TEST_REAL_EQUAL(penalties.lWidth,opt_2d.getPenalties().lWidth)
  TEST_REAL_EQUAL(penalties.rWidth,opt_2d.getPenalties().rWidth)
	TEST_REAL_EQUAL(penalties.height,opt_2d.getPenalties().height)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
