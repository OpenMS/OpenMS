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

#include <OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureNode.h>
#include <OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureMaps.h>
#include <OpenMS/KERNEL/FeatureMap.h>

using namespace OpenMS;
using namespace std;

START_TEST(KDTreeFeatureNode, "$Id$")

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

KDTreeFeatureMaps* kd_data_ptr = new KDTreeFeatureMaps(fmaps, p);

KDTreeFeatureNode* ptr = nullptr;
KDTreeFeatureNode* nullPointer = nullptr;

START_SECTION((KDTreeFeatureNode(KDTreeFeatureMaps* data, Size idx)))
  ptr = new KDTreeFeatureNode(kd_data_ptr, 0);
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~KDTreeFeatureNode()))
  delete ptr;
END_SECTION

KDTreeFeatureNode node_1(kd_data_ptr, 1);

START_SECTION((KDTreeFeatureNode(const KDTreeFeatureNode& rhs)))
  ptr = new KDTreeFeatureNode(node_1);
  TEST_NOT_EQUAL(ptr, nullPointer)
  TEST_EQUAL(ptr->getIndex(), node_1.getIndex())
  TEST_REAL_SIMILAR((*ptr)[0], node_1[0])
  TEST_REAL_SIMILAR((*ptr)[1], node_1[1])
  delete ptr;
END_SECTION

START_SECTION((KDTreeFeatureNode& operator=(KDTreeFeatureNode const& rhs)))
  KDTreeFeatureNode node_2 = node_1;
  TEST_EQUAL(node_2.getIndex(), node_1.getIndex())
  TEST_REAL_SIMILAR(node_2[0], node_1[0])
  TEST_REAL_SIMILAR(node_2[1], node_1[1])
END_SECTION

START_SECTION((Size getIndex() const))
  TEST_EQUAL(node_1.getIndex(), 1)
END_SECTION

START_SECTION((value_type operator[](Size i) const))
  TEST_REAL_SIMILAR(node_1[0], 2000)
  TEST_REAL_SIMILAR(node_1[1], 500)
END_SECTION

delete kd_data_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
