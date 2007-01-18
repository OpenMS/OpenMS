// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
  PRECISION(0.0001)
  TwoDOptimization opt_2d;
  struct OptimizationFunctions::PenaltyFactorsInt penalties;
  opt_2d.getPenalties() = penalties;
  opt_2d.getMaxIterations() = 10;
  opt_2d.getMaxAbsError() = 0.01;
  opt_2d.getMaxRelError() = 0.001;
  
  TwoDOptimization opt_2d_copy;
  opt_2d_copy = opt_2d;
  struct OptimizationFunctions::PenaltyFactorsInt penalties_copy = opt_2d_copy.getPenalties();
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
  struct OptimizationFunctions::PenaltyFactorsInt penalties;
  opt_2d.getPenalties() = penalties;
  opt_2d.getMaxIterations() = 10;
  opt_2d.getMaxAbsError() = 0.01;
  opt_2d.getMaxRelError() = 0.001;
  
  TwoDOptimization opt_2d_copy(opt_2d);
  struct OptimizationFunctions::PenaltyFactorsInt penalties_copy = opt_2d_copy.getPenalties();
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

CHECK((TwoDOptimization(const Param& param)))
  PRECISION(0.0001)
  String file = "data/TwoDOptimization.xml";	
  Param param;
	param.load(file);


  TwoDOptimization opt_2d(param);
  TEST_REAL_EQUAL(10,opt_2d.getPenalties().pos)
  TEST_REAL_EQUAL(opt_2d.getPenalties().lWidth,0)
	TEST_REAL_EQUAL(opt_2d.getPenalties().rWidth,0)
  TEST_REAL_EQUAL(opt_2d.getPenalties().height,1)		
 	TEST_EQUAL(20 == opt_2d.getMaxIterations(), true)
 	TEST_REAL_EQUAL( opt_2d.getMaxAbsError(),0.001)
 	TEST_REAL_EQUAL(opt_2d.getMaxRelError(),0.002)
	TEST_REAL_EQUAL(opt_2d.getMZTolerance(),0.2)
	TEST_REAL_EQUAL(opt_2d.getMaxPeakDistance(),1)			
RESULT

CHECK(( template <typename InputSpectrumIterator,typename OutputPeakType>
			 void twoDOptimize(InputSpectrumIterator& first,
												 InputSpectrumIterator& last,
												 MSExperiment< OutputPeakType >& ms_exp)   ))
  MSSpectrum<DPickedPeak<1> > peaks;
	
	DPickedPeak<1> peak;
	PeakShape peak_shape;
  peak.getPos() = 500;
  peak.getLeftWidthParameter() = 2.5;
  peak.getRightWidthParameter() = 2.6;
  peak.getArea() = 100;
  peak.getIntensity() = 400;
  peak_shape.mz_position = 500;
  peak_shape.left_width = 2.5;
  peak_shape.right_width = 2.5;
  peak_shape.area = 100;
  peak_shape.height = 400;
  peak_shape.type = PeakShapeType::LORENTZ_PEAK;  
	peaks.push_back(peak);
	MSExperiment<DPickedPeak<1> > ms_exp;
	ms_exp.push_back(peaks);
	ms_exp.begin()->setRetentionTime(100);
			
  float origin = 499;
  float spacing = 0.1;

	MSSpectrum<DRawDataPoint<1> >	 raw_spec;
  for (unsigned int i = 0; i < 20 ;++i)
  {
		DRawDataPoint<1> data_point;
		data_point.setPos(origin +i*spacing);
		data_point.setIntensity(peak_shape(origin +i*spacing));
    raw_spec.push_back(data_point);
  }
  MSExperiment<DRawDataPoint<1> > raw_exp;
  raw_exp.push_back(raw_spec);
	raw_exp.begin()->setRetentionTime(100);
  String file = "data/TwoDOptimization.xml";	
  Param param;
	param.load(file);

 	TwoDOptimization opt_2d(param);
  MSExperiment<DRawDataPoint<1> >::const_iterator first,last;
  first = raw_exp.begin();
  last = raw_exp.end();
 	opt_2d.twoDOptimize(first,last,ms_exp);
 	TEST_REAL_EQUAL(peak_shape.mz_position,500)
 	TEST_REAL_EQUAL(peak_shape.left_width,2.5)
 	TEST_REAL_EQUAL(peak_shape.right_width,2.5)
 	TEST_REAL_EQUAL(peak_shape.area,100)
 	TEST_REAL_EQUAL(peak_shape.height,400)
RESULT

CHECK((void setMaxAbsError(const double eps_abs)))
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

CHECK((void setMaxPeakDistance(double max_peak_distance)))
  PRECISION(0.0001)
  double max_peak_distance = 0.01;
   
  TwoDOptimization opt_2d;
  opt_2d.setMaxPeakDistance(max_peak_distance);
    
 	TEST_REAL_EQUAL(max_peak_distance, opt_2d.getMaxPeakDistance())
RESULT

CHECK((void setMZTolerance(double mz_tolerance)))
  PRECISION(0.0001)
  double mz_tolerance = 0.01;
   
  TwoDOptimization opt_2d;
  opt_2d.setMZTolerance(mz_tolerance);
    
 	TEST_REAL_EQUAL(mz_tolerance, opt_2d.getMZTolerance())
RESULT
			
CHECK((void setMaxIterations(const int max_iteration)))
  int number = 20;
   
  TwoDOptimization opt_2d;
  opt_2d.setMaxIterations(number);
    
 	TEST_EQUAL(number == opt_2d.getMaxIterations(), true)
RESULT


	
CHECK((void setPenalties(const struct OptimizationFunctions::PenaltyFactorsInt& penalties)))
  PRECISION(0.0001)
  struct OptimizationFunctions::PenaltyFactorsInt penalties;
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
