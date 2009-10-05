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
// $Maintainer: Eva Lange $
// $Authors: $
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
START_SECTION((OptimizePick( )))
	ptr = new OptimizePick();
	TEST_NOT_EQUAL(ptr, 0)
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
  double abs_err = 0.01;
  double rel_err = 0.001;
  OptimizePick opt_pick(penalties,number,abs_err,rel_err);
  TEST_REAL_SIMILAR(penalties.pos,opt_pick.getPenalties().pos)
  TEST_REAL_SIMILAR(penalties.lWidth,opt_pick.getPenalties().lWidth)
  TEST_REAL_SIMILAR(penalties.rWidth,opt_pick.getPenalties().rWidth)
 	TEST_EQUAL(number == opt_pick.getNumberIterations(), true)
 	TEST_REAL_SIMILAR(abs_err, opt_pick.getMaxAbsError())
 	TEST_REAL_SIMILAR(rel_err, opt_pick.getMaxRelError())
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

START_SECTION((void setMaxAbsError(double eps_abs)))
  TOLERANCE_ABSOLUTE(0.0001)
  double abs_err = 0.01;

  OptimizePick opt_pick;
  opt_pick.setMaxAbsError(abs_err);

 	TEST_REAL_SIMILAR(abs_err, opt_pick.getMaxAbsError())
END_SECTION


START_SECTION((DoubleReal getMaxAbsError() const))
  TOLERANCE_ABSOLUTE(0.0001)
  double abs_err = 0.01;

  OptimizePick opt_pick;
  opt_pick.setMaxAbsError(abs_err);

 	TEST_REAL_SIMILAR(abs_err, opt_pick.getMaxAbsError())
END_SECTION

START_SECTION((double& getMaxAbsError()))
  TOLERANCE_ABSOLUTE(0.0001)
  double abs_err = 0.01;

  OptimizePick opt_pick;
  opt_pick.setMaxAbsError(abs_err);

 	TEST_REAL_SIMILAR(abs_err, opt_pick.getMaxAbsError())
END_SECTION

START_SECTION((void setMaxRelError(double eps_rel)))
  TOLERANCE_ABSOLUTE(0.0001)
  double rel_err = 0.01;

  OptimizePick opt_pick;
  opt_pick.setMaxRelError(rel_err);

 	TEST_REAL_SIMILAR(rel_err, opt_pick.getMaxRelError())
END_SECTION

START_SECTION((DoubleReal getMaxRelError() const))
  TOLERANCE_ABSOLUTE(0.0001)
  double rel_err = 0.01;

  OptimizePick opt_pick;
  opt_pick.setMaxRelError(rel_err);

 	TEST_REAL_SIMILAR(rel_err, opt_pick.getMaxRelError())
END_SECTION

START_SECTION((double& getMaxRelError()))
  TOLERANCE_ABSOLUTE(0.0001)
  double rel_err = 0.01;

  OptimizePick opt_pick;
  opt_pick.setMaxRelError(rel_err);

 	TEST_REAL_SIMILAR(rel_err, opt_pick.getMaxRelError())
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



