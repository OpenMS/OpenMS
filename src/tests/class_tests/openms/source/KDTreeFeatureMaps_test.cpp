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
// $Maintainer: Johannes Veit $
// $Authors: Johannes Veit $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureMaps.h>
#include <OpenMS/KERNEL/FeatureMap.h>

using namespace OpenMS;
using namespace std;

START_TEST(KDTreeFeatureMaps, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Feature f1;
f1.setCharge(2);
f1.setIntensity(100);
f1.setMZ(400);
f1.setRT(1000);

Feature f2;
f2.setCharge(3);
f2.setIntensity(1000);
f2.setMZ(500);
f2.setRT(2000);

FeatureMap fmap;
fmap.push_back(f1);
fmap.push_back(f2);

vector<FeatureMap> fmaps;
fmaps.push_back(fmap);

Param p;
p.setValue("rt_tol", 100);
p.setValue("mz_tol", 10);
p.setValue("mz_unit", "ppm");

KDTreeFeatureMaps* ptr = nullptr;
KDTreeFeatureMaps* nullPointer = nullptr;

START_SECTION((KDTreeFeatureMaps()))
  ptr = new KDTreeFeatureMaps();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~KDTreeFeatureMaps()))
  delete ptr;
END_SECTION

START_SECTION((KDTreeFeatureMaps(const std::vector<MapType>& maps, const Param& param)))
  ptr = new KDTreeFeatureMaps(fmaps, p);
  TEST_NOT_EQUAL(ptr, nullPointer);
  delete ptr;
END_SECTION

KDTreeFeatureMaps kd_data_1(fmaps, p);

START_SECTION((KDTreeFeatureMaps(const KDTreeFeatureMaps& rhs)))
  ptr = new KDTreeFeatureMaps(kd_data_1);
  TEST_NOT_EQUAL(ptr, nullPointer)
  TEST_EQUAL(ptr->size(), kd_data_1.size())
  TEST_EQUAL(ptr->size(), 2)
  TEST_EQUAL(ptr->mz(0), kd_data_1.mz(0))
  TEST_EQUAL(ptr->mz(1), kd_data_1.mz(1))
END_SECTION

START_SECTION((KDTreeFeatureMaps& operator=(const KDTreeFeatureMaps& rhs)))
  KDTreeFeatureMaps kd_data_2 = kd_data_1;
  TEST_EQUAL(kd_data_2.size(), kd_data_1.size())
  TEST_EQUAL(kd_data_2.size(), 2)
  TEST_EQUAL(kd_data_2.mz(0), kd_data_1.mz(0))
  TEST_EQUAL(kd_data_2.mz(1), kd_data_1.mz(1))
END_SECTION

KDTreeFeatureMaps kd_data_3;

START_SECTION((void addMaps(const std::vector<MapType>& maps)))
  kd_data_3.addMaps(fmaps);
  TEST_EQUAL(kd_data_3.size(), 2);
END_SECTION

START_SECTION((void addFeature(Size mt_map_index, const BaseFeature* feature)))
  Feature f3;
  f3.setMZ(300);
  f3.setRT(500);
  kd_data_3.addFeature(2, &f3);
  TEST_EQUAL(kd_data_3.size(), 3);
END_SECTION

START_SECTION((const BaseFeature* feature(Size i) const))
  TEST_EQUAL(kd_data_1.feature(0), &(fmaps[0][0]))
  TEST_EQUAL(kd_data_1.feature(1), &(fmaps[0][1]))
END_SECTION

START_SECTION((double rt(Size i) const))
  TEST_REAL_SIMILAR(kd_data_1.rt(0), 1000)
END_SECTION

START_SECTION((double mz(Size i) const))
  TEST_REAL_SIMILAR(kd_data_1.mz(0), 400)
END_SECTION

START_SECTION((float intensity(Size i) const))
  TEST_REAL_SIMILAR(kd_data_1.intensity(0), 100)
END_SECTION

START_SECTION((Int charge(Size i) const))
  TEST_EQUAL(kd_data_1.charge(0), 2)
END_SECTION

START_SECTION((Size mapIndex(Size i) const))
  TEST_EQUAL(kd_data_1.mapIndex(0), 0)
END_SECTION

START_SECTION((Size size() const))
  TEST_EQUAL(kd_data_1.size(), 2)
  TEST_EQUAL(kd_data_3.size(), 3)
END_SECTION

START_SECTION((Size treeSize() const))
  TEST_EQUAL(kd_data_1.treeSize(), 2)
  TEST_EQUAL(kd_data_3.treeSize(), 3)
END_SECTION

START_SECTION((Size numMaps() const))
  TEST_EQUAL(kd_data_1.numMaps(), 1)
END_SECTION

START_SECTION((void clear()))
  kd_data_3.clear();
  TEST_EQUAL(kd_data_3.size(), 0)
  TEST_EQUAL(kd_data_3.treeSize(), 0)
END_SECTION

START_SECTION((void optimizeTree()))
  NOT_TESTABLE;
END_SECTION

START_SECTION((void getNeighborhood(Size index, std::vector<Size>& result_indices, bool include_features_from_same_map = false) const))
  NOT_TESTABLE;
END_SECTION

START_SECTION((void queryRegion(double rt_low, double rt_high, double mz_low, double mz_high, std::vector<Size>& result_indices, Size ignored_map_index = std::numeric_limits<Size>::max()) const))
  NOT_TESTABLE;
END_SECTION

START_SECTION((void applyTransformations(const std::vector<TransformationModelLowess*>& trafos)))
  NOT_TESTABLE;
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
