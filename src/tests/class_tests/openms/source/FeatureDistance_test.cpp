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
// $Maintainer: Hendrik Weisser
// $Authors: Clemens Groepl, Hendrik Weisser, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureDistance.h>

///////////////////////////

START_TEST(FeatureDistance, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

FeatureDistance* d_ptr = nullptr;
FeatureDistance* d_nullPointer = nullptr;
START_SECTION((FeatureDistance(double max_intensity=1.0, bool force_constraints=false)))
{
	d_ptr = new FeatureDistance();
  TEST_NOT_EQUAL(d_ptr, d_nullPointer);
}
END_SECTION

START_SECTION((~FeatureDistance()))
{
	delete d_ptr;
}
END_SECTION

START_SECTION((std::pair<bool, double> operator()(const BaseFeature& left, const BaseFeature& right)))
{
	FeatureDistance dist(1000.0, false);
	Param param = dist.getDefaults();
	param.setValue("distance_RT:max_difference", 100.0);
	param.setValue("distance_MZ:max_difference", 1.0);
	param.setValue("distance_MZ:exponent", 1.0);
	param.setValue("distance_intensity:weight", 1.0);
	dist.setParameters(param);
	BaseFeature left, right;
	left.setRT(100.0);
	left.setMZ(100.0);
	left.setIntensity(100.0);
	// all distance components vary by 10% of the maximum:
	right.setRT(110.0);
	right.setMZ(100.1);
	right.setIntensity(200.0);
	pair<bool, double> result = dist(left, right);
	TEST_EQUAL(result.first, true);
	TEST_REAL_SIMILAR(result.second, 0.1);
	// no differences:
	result = dist(left, left);
	TEST_EQUAL(result.first, true);
	TEST_REAL_SIMILAR(result.second, 0.0);
	// differences at maximum:
	right.setRT(200.0);
	right.setMZ(101.0);
	right.setIntensity(1000.0);
	left.setIntensity(0.0);
	result = dist(left, right);
	TEST_EQUAL(result.first, true);
	TEST_REAL_SIMILAR(result.second, 1.0);
	// differences beyond maximum:
	right.setRT(300.0);
	result = dist(left, right);
	TEST_EQUAL(result.first, false);
	TEST_REAL_SIMILAR(result.second, 1.33333333);
	FeatureDistance dist2(1000.0, true);
	result = dist2(left, right);
	TEST_EQUAL(result.first, false);
	TEST_EQUAL(result.second, FeatureDistance::infinity);

  // ppm for m/z
	param.setValue("distance_intensity:weight", 0.0);
	param.setValue("distance_RT:weight", 0.0);
	param.setValue("distance_MZ:weight", 1.0);
  param.setValue("distance_MZ:max_difference", 10.0);
  param.setValue("distance_MZ:unit", "ppm");
  dist.setParameters(param);
  left.setRT(100.0);
  left.setIntensity(100.0);
  left.setMZ(100.0);
  right.setRT(110.0);
  right.setIntensity(200.0);
  right.setMZ(100.0 + 100.0/1e6 * 5); // 5ppm off
  result = dist(left, right);
  TEST_EQUAL(result.first, true);
  TEST_REAL_SIMILAR(result.second, 0.5);

  right.setMZ(100.0 + 100.0/1e6 * 20); // 20ppm off
  result = dist(left, right);
  TEST_EQUAL(result.first, false);
  TEST_REAL_SIMILAR(result.second, 2);

  // charge
  param.setValue("ignore_charge", "false");
  dist.setParameters(param);
  right.setMZ(100.0 + 100.0/1e6 * 5); // 5ppm off --> valid in m/z
  // charges differ
  right.setCharge(1);
  left.setCharge(2);
  result = dist(left, right);
  TEST_EQUAL(result.first, false); // --> invalid
  TEST_REAL_SIMILAR(result.second, FeatureDistance::infinity);
  // one charge 0 -- pass filter
  right.setCharge(1);
  left.setCharge(0);
  result = dist(left, right);
  TEST_EQUAL(result.first, true); // --> valid
  TEST_REAL_SIMILAR(result.second, 0.5);
  // ignore charge
  param.setValue("ignore_charge", "true");
  dist.setParameters(param);
  // charges differ, but we don't care
  right.setCharge(1);
  left.setCharge(2);
  result = dist(left, right);
  TEST_EQUAL(result.first, true); // --> valid
  TEST_REAL_SIMILAR(result.second, 0.5);

}
END_SECTION

START_SECTION((FeatureDistance& operator=(const FeatureDistance& other)))
{
	FeatureDistance dist(1000.0, true);
	Param param = dist.getDefaults();
	param.setValue("distance_RT:max_difference", 100.0);
	param.setValue("distance_MZ:max_difference", 1.0);
	param.setValue("distance_MZ:exponent", 1.0);
	param.setValue("distance_intensity:weight", 1.0);
	dist.setParameters(param);
	FeatureDistance dist2;
	dist2 = dist;
	TEST_EQUAL(dist.getParameters(), dist2.getParameters());
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
