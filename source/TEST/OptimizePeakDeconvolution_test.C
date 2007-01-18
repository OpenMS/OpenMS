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
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePeakDeconvolution.h>
#include <OpenMS/FORMAT/Param.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(OptimizePeakDeconvolution, "$Id$")
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

OptimizePeakDeconvolution* ptr = 0;
CHECK((OptimizePeakDeconvolution   ( )))
	ptr = new OptimizePeakDeconvolution();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~OptimizePeakDeconvolution()))
	delete ptr;
RESULT

CHECK((OptimizePeakDeconvolution& operator=(const OptimizePeakDeconvolution& opt)))
  PRECISION(0.0001)
  OptimizePeakDeconvolution opt_deconv;
  struct OptimizationFunctions::PenaltyFactorsInt penalties;
  opt_deconv.getPenalties() = penalties;
  opt_deconv.getNumberIterations() = 10;
  opt_deconv.getMaxAbsError() = 0.01;
  opt_deconv.getMaxRelError() = 0.001;
  opt_deconv.getCharge() = 2;
  
  OptimizePeakDeconvolution opt_deconv_copy;
  opt_deconv_copy = opt_deconv;
  struct OptimizationFunctions::PenaltyFactorsInt penalties_copy = opt_deconv_copy.getPenalties();
  unsigned int number = opt_deconv_copy.getNumberIterations();
  double abs_err = opt_deconv_copy.getMaxAbsError();
  double rel_err = opt_deconv_copy.getMaxRelError();
  double charge = opt_deconv_copy.getCharge();
  TEST_REAL_EQUAL(penalties.pos,penalties_copy.pos)
  TEST_REAL_EQUAL(penalties.lWidth,penalties_copy.lWidth)
  TEST_REAL_EQUAL(penalties.rWidth,penalties_copy.rWidth)
  TEST_REAL_EQUAL(penalties.height,penalties_copy.height)
    
 	TEST_EQUAL(number == 10, true)
	TEST_EQUAL(charge == 2, true)
 	TEST_REAL_EQUAL(abs_err, 0.01)
 	TEST_REAL_EQUAL(rel_err, 0.001)
RESULT

CHECK((OptimizePeakDeconvolution( )))
	OptimizePeakDeconvolution opt_deconv;
  int max_iteration = opt_deconv.getNumberIterations();
  double eps_abs = opt_deconv.getMaxAbsError();
  double eps_rel = opt_deconv.getMaxRelError();
  int charge = opt_deconv.getCharge();

  TEST_EQUAL(max_iteration == 0, true)
	TEST_EQUAL(charge == 1, true)
 	TEST_REAL_EQUAL(eps_abs, 0.0)
 	TEST_REAL_EQUAL(eps_rel, 0.0)		
RESULT

CHECK((OptimizePeakDeconvolution(const OptimizePeakDeconvolution& opt)))
  PRECISION(0.0001)
  OptimizePeakDeconvolution opt_deconv;
  struct OptimizationFunctions::PenaltyFactorsInt penalties;
  opt_deconv.getPenalties() = penalties;
  opt_deconv.getNumberIterations() = 10;
  opt_deconv.getMaxAbsError() = 0.01;
  opt_deconv.getMaxRelError() = 0.001;
  opt_deconv.getCharge() = 2;
  
  OptimizePeakDeconvolution opt_deconv_copy(opt_deconv);
  struct OptimizationFunctions::PenaltyFactorsInt penalties_copy = opt_deconv_copy.getPenalties();
  unsigned int number = opt_deconv_copy.getNumberIterations();
  double abs_err = opt_deconv_copy.getMaxAbsError();
  double rel_err = opt_deconv_copy.getMaxRelError();
  double charge = opt_deconv_copy.getCharge();
  TEST_REAL_EQUAL(penalties.pos,penalties_copy.pos)
  TEST_REAL_EQUAL(penalties.lWidth,penalties_copy.lWidth)
  TEST_REAL_EQUAL(penalties.rWidth,penalties_copy.rWidth)
  TEST_REAL_EQUAL(penalties.height,penalties_copy.height)
    
 	TEST_EQUAL(number == 10, true)
	TEST_EQUAL(charge == 2, true)
 	TEST_REAL_EQUAL(abs_err, 0.01)
 	TEST_REAL_EQUAL(rel_err, 0.001)
RESULT

CHECK((OptimizePeakDeconvolution(const struct OptimizationFunctions::PenaltyFactorsInt& penalties, const int max_iteration, const double eps_abs, const double eps_rel, const int charge)))
  PRECISION(0.0001)
  struct OptimizationFunctions::PenaltyFactorsInt penalties;
  penalties.pos = 0;
  penalties.lWidth = 1;
  penalties.rWidth = 2;
  penalties.height = 3;
  int number = 10;
  double abs_err = 0.01;
  double rel_err = 0.001;
  int charge = 5;

  OptimizePeakDeconvolution opt_deconv(penalties,number,abs_err,rel_err,charge);
  TEST_REAL_EQUAL(penalties.pos,opt_deconv.getPenalties().pos)
  TEST_REAL_EQUAL(penalties.lWidth,opt_deconv.getPenalties().lWidth)
  TEST_REAL_EQUAL(penalties.rWidth,opt_deconv.getPenalties().rWidth)
  TEST_REAL_EQUAL(penalties.height,opt_deconv.getPenalties().height)		
 	TEST_EQUAL(number == opt_deconv.getNumberIterations(), true)
 	TEST_REAL_EQUAL(abs_err, opt_deconv.getMaxAbsError())
 	TEST_REAL_EQUAL(rel_err, opt_deconv.getMaxRelError())
  TEST_EQUAL(charge == opt_deconv.getCharge(), true)		
RESULT

CHECK((bool optimize(std::vector<PeakShape>& peaks,Param& param,int failure)))
	std::vector<PeakShape> peak_shapes(1);
	PeakShape peak_shape;
  peak_shape.mz_position = 500;
  peak_shape.left_width = 2.5;
  peak_shape.right_width = 2.5;
  peak_shape.area = 100;
  peak_shape.height = 400;
  peak_shape.type = PeakShapeType::LORENTZ_PEAK;
  peak_shapes[0] = peak_shape;
//  peak_shapes[1] = peak_shape;

std::cout << "peaks ready"<<std::endl;
  float origin = 499;
  float spacing = 0.1;
 
	OptimizationFunctions::positions_DC_.resize(20);
  OptimizationFunctions::signal_DC_.resize(20);
  for (unsigned int i = 0; i < 20 ;++i)
  {
  	OptimizationFunctions::positions_DC_[i] = origin +i*spacing;
    OptimizationFunctions::signal_DC_[i] = peak_shape(origin +i*spacing);
   }
std::cout << "positions ready"<<std::endl;
  String file = "data/OptimizePeakDeconvolution.xml";
  Param param;
  param.load(file);

std::cout << "param loaded"<<std::endl;

 	OptimizePeakDeconvolution opt_deconv;
 	opt_deconv.optimize(peak_shapes,param,1);
std::cout << "opt finished"<<std::endl;
 	TEST_REAL_EQUAL(peak_shape.mz_position,500)
 	TEST_REAL_EQUAL(peak_shape.left_width,2.5)
 	TEST_REAL_EQUAL(peak_shape.right_width,2.5)
 	TEST_REAL_EQUAL(peak_shape.area,100)
 	TEST_REAL_EQUAL(peak_shape.height,400)
RESULT

CHECK((void setMaxAbsError(const double eps_abs)))
  PRECISION(0.0001)
  double abs_err = 0.01;
   
  OptimizePeakDeconvolution opt_deconv;
  opt_deconv.setMaxAbsError(abs_err);
    
 	TEST_REAL_EQUAL(abs_err, opt_deconv.getMaxAbsError())
RESULT

CHECK((void setMaxRelError(const double eps_rel)))
  PRECISION(0.0001)
  double rel_err = 0.01;
   
  OptimizePeakDeconvolution opt_deconv;
  opt_deconv.setMaxRelError(rel_err);
    
 	TEST_REAL_EQUAL(rel_err, opt_deconv.getMaxRelError())
RESULT

CHECK((void setNumberIterations(const int max_iteration)))
  int number = 20;
   
  OptimizePeakDeconvolution opt_deconv;
  opt_deconv.setNumberIterations(number);
    
 	TEST_EQUAL(number == opt_deconv.getNumberIterations(), true)
RESULT

CHECK((void setCharge(const int charge)))
  int charge = 2;
   
  OptimizePeakDeconvolution opt_deconv;
  opt_deconv.setCharge(charge);
    
 	TEST_EQUAL(charge == opt_deconv.getCharge(), true)
RESULT

	
CHECK((void setPenalties(const struct OptimizationFunctions::PenaltyFactorsInt& penalties)))
  PRECISION(0.0001)
  struct OptimizationFunctions::PenaltyFactorsInt penalties;
  penalties.pos = 0;
  penalties.lWidth = 1;
  penalties.rWidth = 2;
  penalties.height = 3;

  OptimizePeakDeconvolution opt_deconv;
  opt_deconv.setPenalties(penalties);
  TEST_REAL_EQUAL(penalties.pos,opt_deconv.getPenalties().pos)
  TEST_REAL_EQUAL(penalties.lWidth,opt_deconv.getPenalties().lWidth)
  TEST_REAL_EQUAL(penalties.rWidth,opt_deconv.getPenalties().rWidth)
	TEST_REAL_EQUAL(penalties.height,opt_deconv.getPenalties().height)
RESULT

	

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



