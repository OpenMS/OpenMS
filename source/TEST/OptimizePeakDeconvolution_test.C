// -*- mode: C++; tab-width: 2; -*-
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
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePeakDeconvolution.h>
#include <OpenMS/DATASTRUCTURES/Param.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(OptimizePeakDeconvolution, "$Id$")
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

OptimizePeakDeconvolution* ptr = 0;
START_SECTION((OptimizePeakDeconvolution   ( )))
	ptr = new OptimizePeakDeconvolution();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((~OptimizePeakDeconvolution()))
	delete ptr;
END_SECTION

START_SECTION((OptimizePeakDeconvolution& operator=(const OptimizePeakDeconvolution& opt)))
  TOLERANCE_ABSOLUTE(0.0001)
  OptimizePeakDeconvolution opt_deconv;
  struct OptimizationFunctions::PenaltyFactorsIntensity penalties;
  opt_deconv.setPenalties(penalties);

  opt_deconv.setCharge(2);
  
  OptimizePeakDeconvolution opt_deconv_copy;
  opt_deconv_copy = opt_deconv;
  struct OptimizationFunctions::PenaltyFactorsIntensity penalties_copy = opt_deconv_copy.getPenalties();
 
  double charge = opt_deconv_copy.getCharge();
  TEST_REAL_SIMILAR(penalties.pos,penalties_copy.pos)
  TEST_REAL_SIMILAR(penalties.lWidth,penalties_copy.lWidth)
  TEST_REAL_SIMILAR(penalties.rWidth,penalties_copy.rWidth)
  TEST_REAL_SIMILAR(penalties.height,penalties_copy.height)
    

	TEST_EQUAL(charge == 2, true)
 
END_SECTION


START_SECTION((OptimizePeakDeconvolution(const OptimizePeakDeconvolution& opt)))
  TOLERANCE_ABSOLUTE(0.0001)
  OptimizePeakDeconvolution opt_deconv;
  struct OptimizationFunctions::PenaltyFactorsIntensity penalties;
  opt_deconv.setPenalties(penalties);
  opt_deconv.setCharge(2);
  
  OptimizePeakDeconvolution opt_deconv_copy(opt_deconv);
  struct OptimizationFunctions::PenaltyFactorsIntensity penalties_copy = opt_deconv_copy.getPenalties();
  double charge = opt_deconv_copy.getCharge();
  TEST_REAL_SIMILAR(penalties.pos,penalties_copy.pos)
  TEST_REAL_SIMILAR(penalties.lWidth,penalties_copy.lWidth)
  TEST_REAL_SIMILAR(penalties.rWidth,penalties_copy.rWidth)
  TEST_REAL_SIMILAR(penalties.height,penalties_copy.height)
    
	TEST_EQUAL(charge == 2, true)
 
END_SECTION


START_SECTION((bool optimize(std::vector<PeakShape>& peaks,int failure)))
	std::vector<PeakShape> peak_shapes(1);
	PeakShape peak_shape;
  peak_shape.mz_position = 500;
  peak_shape.left_width = 2.5;
  peak_shape.right_width = 2.5;
  peak_shape.area = 100;
  peak_shape.height = 400;
  peak_shape.type = PeakShape::LORENTZ_PEAK;
  peak_shapes[0] = peak_shape;
//  peak_shapes[1] = peak_shape;
  float origin = 499;
  float spacing = 0.1;
 
	OptimizationFunctions::positions_DC_.resize(20);
  OptimizationFunctions::signal_DC_.resize(20);
  for (unsigned int i = 0; i < 20 ;++i)
  {
  	OptimizationFunctions::positions_DC_[i] = origin +i*spacing;
    OptimizationFunctions::signal_DC_[i] = peak_shape(origin +i*spacing);
   }
  String file = "data/OptimizePeakDeconvolution.ini";
  Param param;
  param.load(file);


 	OptimizePeakDeconvolution opt_deconv;
  opt_deconv.setParameters(param.copy("deconvolution:fitting:",true));
 	opt_deconv.optimize(peak_shapes,1);
 	TEST_REAL_SIMILAR(peak_shape.mz_position,500)
 	TEST_REAL_SIMILAR(peak_shape.left_width,2.5)
 	TEST_REAL_SIMILAR(peak_shape.right_width,2.5)
 	TEST_REAL_SIMILAR(peak_shape.area,100)
 	TEST_REAL_SIMILAR(peak_shape.height,400)
END_SECTION


START_SECTION((void setCharge(const int charge)))
  int charge = 2;
   
  OptimizePeakDeconvolution opt_deconv;
  opt_deconv.setCharge(charge);
    
 	TEST_EQUAL(charge == opt_deconv.getCharge(), true)
END_SECTION

START_SECTION((int getCharge() const))
  int charge = 2;
   
  OptimizePeakDeconvolution opt_deconv;
  opt_deconv.setCharge(charge);
    
 	TEST_EQUAL(charge == opt_deconv.getCharge(), true)
END_SECTION
	

START_SECTION((void setPenalties(const OptimizationFunctions::PenaltyFactorsIntensity& penalties)))
  TOLERANCE_ABSOLUTE(0.0001)
  struct OptimizationFunctions::PenaltyFactorsIntensity penalties;
  penalties.pos = 0;
  penalties.lWidth = 1;
  penalties.rWidth = 2;
  penalties.height = 3;

  OptimizePeakDeconvolution opt_deconv;
  opt_deconv.setPenalties(penalties);
  TEST_REAL_SIMILAR(penalties.pos,opt_deconv.getPenalties().pos)
  TEST_REAL_SIMILAR(penalties.lWidth,opt_deconv.getPenalties().lWidth)
  TEST_REAL_SIMILAR(penalties.rWidth,opt_deconv.getPenalties().rWidth)
	TEST_REAL_SIMILAR(penalties.height,opt_deconv.getPenalties().height)
END_SECTION

START_SECTION(( const OptimizationFunctions::PenaltyFactorsIntensity& getPenalties() const))
  TOLERANCE_ABSOLUTE(0.0001)
  struct OptimizationFunctions::PenaltyFactorsIntensity penalties;
  penalties.pos = 0;
  penalties.lWidth = 1;
  penalties.rWidth = 2;
  penalties.height = 3;

  OptimizePeakDeconvolution opt_deconv;
  opt_deconv.setPenalties(penalties);
  TEST_REAL_SIMILAR(penalties.pos,opt_deconv.getPenalties().pos)
  TEST_REAL_SIMILAR(penalties.lWidth,opt_deconv.getPenalties().lWidth)
  TEST_REAL_SIMILAR(penalties.rWidth,opt_deconv.getPenalties().rWidth)
	TEST_REAL_SIMILAR(penalties.height,opt_deconv.getPenalties().height)
END_SECTION
	

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



