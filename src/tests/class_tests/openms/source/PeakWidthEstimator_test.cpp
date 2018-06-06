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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/PeakWidthEstimator.h>

#include <OpenMS/FORMAT/MzMLFile.h>

using namespace OpenMS;

START_TEST(PeakWidthEstimator, "$Id$")

PeakMap exp;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_orbitrap.mzML"), exp);

PeakPickerHiRes picker;
Param param = picker.getParameters();
param.setValue("ms_levels", ListUtils::create<Int>("1"));
param.setValue("signal_to_noise", 0.0);
picker.setParameters(param);

std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > boundaries_exp_s;
std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > boundaries_exp_c;

PeakMap exp_picked;
picker.pickExperiment(exp, exp_picked, boundaries_exp_s, boundaries_exp_c);

PeakWidthEstimator* nullPointer = nullptr;
PeakWidthEstimator* ptr;

START_SECTION(PeakWidthEstimator(const PeakMap & exp_picked, const std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > & boundaries))
{
  PeakWidthEstimator estimator(exp_picked, boundaries_exp_s);
  TEST_REAL_SIMILAR(estimator.getPeakWidth(365.3),0.00886469661896705);
  ptr = new PeakWidthEstimator(exp_picked, boundaries_exp_s);
  TEST_NOT_EQUAL(ptr, nullPointer);
  delete ptr;
}
END_SECTION

PeakWidthEstimator estimator2(exp_picked, boundaries_exp_s);

START_SECTION(double getPeakWidth(double mz))
{
  TEST_REAL_SIMILAR(estimator2.getPeakWidth(365.3),0.00886469661896705);
  TEST_REAL_SIMILAR(estimator2.getPeakWidth(305.1),0.00886699447290451);    // outside m/z range
  TEST_REAL_SIMILAR(estimator2.getPeakWidth(405.1),0.01184458329884600);    // outside m/z range
}
END_SECTION

END_TEST

