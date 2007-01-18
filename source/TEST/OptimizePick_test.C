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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePick.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(OptimizePick, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

OptimizePick* ptr = 0;
CHECK((OptimizePick( )))
	ptr = new OptimizePick();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~OptimizePick()))
	delete ptr;
RESULT

CHECK((OptimizePick& operator=(const OptimizePick& opt)))
	PRECISION(0.0001)
  OptimizePick opt_pick;
  struct OptimizationFunctions::PenaltyFactors penalties;
  opt_pick.getPenalties() = penalties;
  opt_pick.getNumberIterations() = 10;
  opt_pick.getMaxAbsError() = 0.01;
  opt_pick.getMaxRelError() = 0.001;
  
  OptimizePick opt_pick_copy;
  opt_pick_copy = opt_pick;
  struct OptimizationFunctions::PenaltyFactors penalties_copy = opt_pick_copy.getPenalties();
  unsigned int number = opt_pick_copy.getNumberIterations();
  double abs_err = opt_pick_copy.getMaxAbsError();
  double rel_err = opt_pick_copy.getMaxRelError();
  TEST_REAL_EQUAL(penalties.pos,penalties_copy.pos)
  TEST_REAL_EQUAL(penalties.lWidth,penalties_copy.lWidth)
  TEST_REAL_EQUAL(penalties.rWidth,penalties_copy.rWidth)
 	TEST_EQUAL(number == 10, true)
 	TEST_REAL_EQUAL(abs_err, 0.01)
 	TEST_REAL_EQUAL(rel_err, 0.001)
RESULT

CHECK((OptimizePick( )))
  // ???
RESULT

CHECK((OptimizePick(const OptimizePick& opt)))
  PRECISION(0.0001)
  OptimizePick opt_pick;
  struct OptimizationFunctions::PenaltyFactors penalties;
  opt_pick.getPenalties() = penalties;
  opt_pick.getNumberIterations() = 10;
  opt_pick.getMaxAbsError() = 0.01;
  opt_pick.getMaxRelError() = 0.001;
  
  OptimizePick opt_pick_copy(opt_pick);
  struct OptimizationFunctions::PenaltyFactors penalties_copy = opt_pick_copy.getPenalties();
  unsigned int number = opt_pick_copy.getNumberIterations();
  double abs_err = opt_pick_copy.getMaxAbsError();
  double rel_err = opt_pick_copy.getMaxRelError();
  TEST_REAL_EQUAL(penalties.pos,penalties_copy.pos)
  TEST_REAL_EQUAL(penalties.lWidth,penalties_copy.lWidth)
  TEST_REAL_EQUAL(penalties.rWidth,penalties_copy.rWidth)
 	TEST_EQUAL(number == 10, true)
 	TEST_REAL_EQUAL(abs_err, 0.01)
 	TEST_REAL_EQUAL(rel_err, 0.001)
RESULT

CHECK((OptimizePick(const struct OptimizationFunctions::PenaltyFactors& penalties_, const int max_iteration_, const double eps_abs_, const double eps_rel_ )))
  PRECISION(0.0001)
  struct OptimizationFunctions::PenaltyFactors penalties;
  penalties.pos = 0;
  penalties.lWidth = 1;
  penalties.rWidth = 2;
  unsigned int number = 10;
  double abs_err = 0.01;
  double rel_err = 0.001;
  
  OptimizePick opt_pick(penalties,number,abs_err,rel_err);
  TEST_REAL_EQUAL(penalties.pos,opt_pick.getPenalties().pos)
  TEST_REAL_EQUAL(penalties.lWidth,opt_pick.getPenalties().lWidth)
  TEST_REAL_EQUAL(penalties.rWidth,opt_pick.getPenalties().rWidth)
 	TEST_EQUAL(number == opt_pick.getNumberIterations(), true)
 	TEST_REAL_EQUAL(abs_err, opt_pick.getMaxAbsError())
 	TEST_REAL_EQUAL(rel_err, opt_pick.getMaxRelError())
RESULT

CHECK((void optimize(std::vector<PeakShape>& peaks)))
	std::vector<PeakShape> peak_shapes(1);
	PeakShape peak_shape;
  peak_shape.mz_position = 500;
  peak_shape.left_width = 0.1;
  peak_shape.right_width = 0.1;
  peak_shape.area = 100;
  peak_shape.height = 400;
  peak_shape.type = PeakShapeType::LORENTZ_PEAK;
  peak_shapes[0] = peak_shape;
	
  float origin = 499;
  float spacing = 0.1;
 
	OptimizationFunctions::positions_.resize(20);
  OptimizationFunctions::signal_.resize(20);
  for (unsigned int i = 0; i < 20 ;++i)
  {
  	OptimizationFunctions::positions_[i] = origin +i*spacing;
    OptimizationFunctions::signal_[i] = peak_shape(origin +i*spacing);
   }
 	OptimizePick opt_pick;
 	opt_pick.optimize(peak_shapes);
 	TEST_REAL_EQUAL(peak_shape.mz_position,500)
 	TEST_REAL_EQUAL(peak_shape.left_width,0.1)
 	TEST_REAL_EQUAL(peak_shape.right_width,0.1)
 	TEST_REAL_EQUAL(peak_shape.area,100)
 	TEST_REAL_EQUAL(peak_shape.height,400)
RESULT

CHECK((void setMaxAbsError(const double eps_abs)))
  PRECISION(0.0001)
  double abs_err = 0.01;
   
  OptimizePick opt_pick;
  opt_pick.setMaxAbsError(abs_err);
    
 	TEST_REAL_EQUAL(abs_err, opt_pick.getMaxAbsError())
RESULT

CHECK((void setMaxRelError(const double eps_rel)))
  PRECISION(0.0001)
  double rel_err = 0.01;
   
  OptimizePick opt_pick;
  opt_pick.setMaxRelError(rel_err);
    
 	TEST_REAL_EQUAL(rel_err, opt_pick.getMaxRelError())
RESULT

CHECK((void setNumberIterations(const int max_iteration)))
  unsigned int number = 20;
   
  OptimizePick opt_pick;
  opt_pick.setNumberIterations(number);
    
 	TEST_EQUAL(number == opt_pick.getNumberIterations(), true)
RESULT

CHECK((void setPenalties(const struct OptimizationFunctions::PenaltyFactors& penalties)))
  PRECISION(0.0001)
  struct OptimizationFunctions::PenaltyFactors penalties;
  penalties.pos = 0;
  penalties.lWidth = 1;
  penalties.rWidth = 2;
    
  OptimizePick opt_pick;
  opt_pick.setPenalties(penalties);
  TEST_REAL_EQUAL(penalties.pos,opt_pick.getPenalties().pos)
  TEST_REAL_EQUAL(penalties.lWidth,opt_pick.getPenalties().lWidth)
  TEST_REAL_EQUAL(penalties.rWidth,opt_pick.getPenalties().rWidth)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



