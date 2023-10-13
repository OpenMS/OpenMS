// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePick.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(OptimizePick, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

OptimizePick* ptr = nullptr;
OptimizePick* nullPointer = nullptr;
START_SECTION((OptimizePick( )))
	ptr = new OptimizePick();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~OptimizePick()))
	delete ptr;
END_SECTION

START_SECTION((OptimizePick(const struct OptimizationFunctions::PenaltyFactors& penalties_, const int max_iteration_, const double eps_abs_, const double eps_rel_ )))
  TOLERANCE_ABSOLUTE(0.0001)
	struct OptimizationFunctions::PenaltyFactors penalties;
  penalties.pos = 0;
  penalties.lWidth = 1;
  penalties.rWidth = 2;
  unsigned int number = 10;
  OptimizePick opt_pick(penalties,number);
  TEST_REAL_SIMILAR(penalties.pos,opt_pick.getPenalties().pos)
  TEST_REAL_SIMILAR(penalties.lWidth,opt_pick.getPenalties().lWidth)
  TEST_REAL_SIMILAR(penalties.rWidth,opt_pick.getPenalties().rWidth)
 	TEST_EQUAL(number == opt_pick.getNumberIterations(), true)
END_SECTION

START_SECTION((void optimize(std::vector< PeakShape > &peaks, Data &data)))
	std::vector<PeakShape> peak_shapes(1);
	PeakShape peak_shape;
  peak_shape.mz_position = 500;
  peak_shape.left_width = 0.1;
  peak_shape.right_width = 0.1;
  peak_shape.area = 100;
  peak_shape.height = 400;
  peak_shape.type = PeakShape::LORENTZ_PEAK;
  peak_shapes[0] = peak_shape;

  float origin = 499.f;
  float spacing = 0.1f;
  OptimizePick::Data data;
	data.positions.resize(20);
  data.signal.resize(20);
  for (Size i = 0; i < 20 ;++i)
  {
  	data.positions[i] = origin +i*spacing;
    data.signal[i] = peak_shape(origin +i*spacing);
   }
 	OptimizePick opt_pick;
  opt_pick.optimize(peak_shapes,data);
 	TEST_REAL_SIMILAR(peak_shape.mz_position,500)
 	TEST_REAL_SIMILAR(peak_shape.left_width,0.1)
 	TEST_REAL_SIMILAR(peak_shape.right_width,0.1)
 	TEST_REAL_SIMILAR(peak_shape.area,100)
 	TEST_REAL_SIMILAR(peak_shape.height,400)
END_SECTION

START_SECTION((void setNumberIterations(const int max_iteration)))
  unsigned int number = 20;

  OptimizePick opt_pick;
  opt_pick.setNumberIterations(number);

 	TEST_EQUAL(number == opt_pick.getNumberIterations(), true)
END_SECTION

START_SECTION((unsigned int& getNumberIterations()))
  unsigned int number = 20;

  OptimizePick opt_pick;
  opt_pick.setNumberIterations(number);

 	TEST_EQUAL(number == opt_pick.getNumberIterations(), true)
END_SECTION


START_SECTION((UInt getNumberIterations() const))
  unsigned int number = 20;

  OptimizePick opt_pick;
  opt_pick.setNumberIterations(number);

 	TEST_EQUAL(number == opt_pick.getNumberIterations(), true)
END_SECTION

START_SECTION((void setPenalties(const struct OptimizationFunctions::PenaltyFactors& penalties)))
  TOLERANCE_ABSOLUTE(0.0001)
	struct OptimizationFunctions::PenaltyFactors penalties;
	penalties.pos = 0;
  penalties.lWidth = 1;
  penalties.rWidth = 2;

  OptimizePick opt_pick;
  opt_pick.setPenalties(penalties);
  TEST_REAL_SIMILAR(penalties.pos,opt_pick.getPenalties().pos)
  TEST_REAL_SIMILAR(penalties.lWidth,opt_pick.getPenalties().lWidth)
  TEST_REAL_SIMILAR(penalties.rWidth,opt_pick.getPenalties().rWidth)
END_SECTION

START_SECTION((struct OptimizationFunctions::PenaltyFactors& getPenalties() const ))
  TOLERANCE_ABSOLUTE(0.0001)
	struct OptimizationFunctions::PenaltyFactors penalties;
	penalties.pos = 0;
  penalties.lWidth = 1;
  penalties.rWidth = 2;

  OptimizePick opt_pick;
  opt_pick.setPenalties(penalties);
  TEST_REAL_SIMILAR(penalties.pos,opt_pick.getPenalties().pos)
  TEST_REAL_SIMILAR(penalties.lWidth,opt_pick.getPenalties().lWidth)
  TEST_REAL_SIMILAR(penalties.rWidth,opt_pick.getPenalties().rWidth)
END_SECTION

START_SECTION((struct OptimizationFunctions::PenaltyFactors& getPenalties()))
  TOLERANCE_ABSOLUTE(0.0001)
	struct OptimizationFunctions::PenaltyFactors penalties;
	penalties.pos = 0;
  penalties.lWidth = 1;
  penalties.rWidth = 2;

  OptimizePick opt_pick;
  opt_pick.setPenalties(penalties);
  TEST_REAL_SIMILAR(penalties.pos,opt_pick.getPenalties().pos)
  TEST_REAL_SIMILAR(penalties.lWidth,opt_pick.getPenalties().lWidth)
  TEST_REAL_SIMILAR(penalties.rWidth,opt_pick.getPenalties().rWidth)
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



