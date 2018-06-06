// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePeakDeconvolution.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(OptimizePeakDeconvolution, "$Id$")
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

OptimizePeakDeconvolution* ptr = nullptr;
OptimizePeakDeconvolution* nullPointer = nullptr;
START_SECTION((OptimizePeakDeconvolution   ( )))
	ptr = new OptimizePeakDeconvolution();
	TEST_NOT_EQUAL(ptr, nullPointer)
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


START_SECTION((bool optimize(std::vector<PeakShape>& peaks,Data& data)))
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
  float spacing = 0.1f;
  OptimizePeakDeconvolution::Data data;
	data.positions.resize(20);
  data.signal.resize(20);
  int scale = 1;
  for (Size i = 0; i < 20 ;++i)
  {
  	data.positions[i] = origin +i*spacing;
    data.signal[i] = peak_shape(origin +i*spacing) + scale * 0.1;
    scale *= -1;
   }
  String file = OPENMS_GET_TEST_DATA_PATH("OptimizePeakDeconvolution.ini");
  Param param;
  ParamXMLFile paramFile;
  paramFile.load(file, param);


 	OptimizePeakDeconvolution opt_deconv;
  opt_deconv.setParameters(param.copy("deconvolution:fitting:",true));
  opt_deconv.optimize(peak_shapes,data);
 	TEST_REAL_SIMILAR(peak_shape.mz_position,500)
 	TEST_REAL_SIMILAR(peak_shape.left_width,2.5)
 	TEST_REAL_SIMILAR(peak_shape.right_width,2.5)
 	TEST_REAL_SIMILAR(peak_shape.area,100)
 	TEST_REAL_SIMILAR(peak_shape.height,400)
END_SECTION


START_SECTION((void setCharge(const Int charge)))
  Int charge = 2;
   
  OptimizePeakDeconvolution opt_deconv;
  opt_deconv.setCharge(charge);
    
 	TEST_EQUAL(charge == opt_deconv.getCharge(), true)
END_SECTION

START_SECTION((Int getCharge() const))
  Int charge = 2;
   
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



