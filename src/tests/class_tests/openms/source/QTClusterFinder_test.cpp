// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/QTClusterFinder.h>
///////////////////////////

#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>

using namespace OpenMS;
using namespace std;

START_TEST(QTClusterFinder, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

QTClusterFinder* ptr = 0;
QTClusterFinder* nullPointer = 0;
BaseGroupFinder* base_nullPointer = 0;

START_SECTION((QTClusterFinder()))
{
	ptr = new QTClusterFinder();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((virtual ~QTClusterFinder()))
	delete ptr;
END_SECTION

START_SECTION((static BaseGroupFinder* create()))
{
	BaseGroupFinder* base_ptr = 0;
	base_ptr = QTClusterFinder::create();
  TEST_NOT_EQUAL(base_ptr, base_nullPointer);
}
END_SECTION

START_SECTION((static const String getProductName()))
{
	QTClusterFinder finder;
  TEST_EQUAL(finder.getName() == "qt", true);
}
END_SECTION

START_SECTION((void run(const std::vector<FeatureMap >& input_maps, ConsensusMap& result_map)))
{
  vector<FeatureMap > input(2);
  Feature feat1;
  Feature feat2;
  DPosition<2> pos1(0,0);
  DPosition<2> pos2(100,200);
  feat1.setPosition(pos1);
  feat1.setUniqueId(0);
  feat2.setPosition(pos2);
  feat2.setUniqueId(1);
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("AAA"));
  feat1.getPeptideIdentifications().resize(1);
  feat1.getPeptideIdentifications()[0].insertHit(hit);
  hit.setSequence(AASequence::fromString("CCC"));
  feat2.getPeptideIdentifications().resize(1);
  feat2.getPeptideIdentifications()[0].insertHit(hit);
  input[0].push_back(feat1);
  input[0].push_back(feat2);

  Feature feat3;
	Feature feat4;
  Feature feat5;
  DPosition<2> pos3(4,0.04);
  DPosition<2> pos4(5,0.05);
  DPosition<2> pos5(104,200.04);
  feat3.setPosition(pos3);
  feat3.setUniqueId(0);
  feat4.setPosition(pos4);
  feat4.setUniqueId(1);
  feat5.setPosition(pos5);
  feat5.setUniqueId(2);
  hit.setSequence(AASequence::fromString("DDD"));
  feat3.getPeptideIdentifications().resize(1);
  feat3.getPeptideIdentifications()[0].insertHit(hit);
  hit.setSequence(AASequence::fromString("AAA"));
  feat4.getPeptideIdentifications().resize(1);
  feat4.getPeptideIdentifications()[0].insertHit(hit);
  // no peptide ID for "feat5"
  input[1].push_back(feat3);
  input[1].push_back(feat4);
  input[1].push_back(feat5);

  input[0].updateRanges();
  input[1].updateRanges();

  QTClusterFinder finder;
	Param param = finder.getDefaults();
	param.setValue("distance_RT:max_difference", 5.1);
	param.setValue("distance_MZ:max_difference", 0.1);
	finder.setParameters(param);
	ConsensusMap result;
	finder.run(input, result);
	TEST_EQUAL(result.size(), 3);
	ABORT_IF(result.size() != 3);

  ConsensusFeature::HandleSetType group1 = result[0].getFeatures();
  ConsensusFeature::HandleSetType group2 = result[1].getFeatures();
  ConsensusFeature::HandleSetType group3 = result[2].getFeatures();

  FeatureHandle ind1(0, feat1);
  FeatureHandle ind2(0, feat2);
  FeatureHandle ind3(1, feat3);
  FeatureHandle ind4(1, feat4);
  FeatureHandle ind5(1, feat5);

  ConsensusFeature::HandleSetType::const_iterator it;
	// don't know why the order is this way, but it shouldn't matter...
	it = group1.begin();
  STATUS(*it);
	STATUS(ind2);
	TEST_EQUAL(*(it) == ind2, true);
	++it;
  STATUS(*it);
	STATUS(ind5);
  TEST_EQUAL(*(it) == ind5, true);

	it = group2.begin();
  STATUS(*it);
	STATUS(ind1);
	TEST_EQUAL(*(it) == ind1, true);
	++it;
  STATUS(*it);
	STATUS(ind3);
  TEST_EQUAL(*(it) == ind3, true);

  it = group3.begin();
  STATUS(*it);
	STATUS(ind4);
  TEST_EQUAL(*(it) == ind4, true);


	// test annotation-specific matching (simple case):

	param.setValue("use_identifications", "true");
	finder.setParameters(param);
	finder.run(input, result);
	TEST_EQUAL(result.size(), 3);
	ABORT_IF(result.size() != 3);

  group1 = result[0].getFeatures();
  group2 = result[1].getFeatures();
  group3 = result[2].getFeatures();

	it = group1.begin();
  STATUS(*it);
	STATUS(ind2);
	TEST_EQUAL(*(it) == ind2, true);
	++it;
  STATUS(*it);
	STATUS(ind5);
  TEST_EQUAL(*(it) == ind5, true);

	it = group2.begin();
  STATUS(*it);
	STATUS(ind1);
	TEST_EQUAL(*(it) == ind1, true);
	++it;
  STATUS(*it);
	STATUS(ind3);
  TEST_EQUAL(*(it) == ind4, true);

  it = group3.begin();
  STATUS(*it);
	STATUS(ind4);
  TEST_EQUAL(*(it) == ind3, true);


	// test annotation-specific matching (complex case):

	input.resize(3);
  Feature feat6;
  Feature feat7;
  DPosition<2> pos6(104,200.04);
  DPosition<2> pos7(108,200.08);
  feat6.setPosition(pos6);
  feat6.setUniqueId(0);
  feat7.setPosition(pos7);
  feat7.setUniqueId(1);
  hit.setSequence(AASequence::fromString("EEE"));
  feat6.getPeptideIdentifications().resize(1);
  feat6.getPeptideIdentifications()[0].insertHit(hit);
  hit.setSequence(AASequence::fromString("CCC"));
  feat7.getPeptideIdentifications().resize(1);
  feat7.getPeptideIdentifications()[0].insertHit(hit);
  input[2].push_back(feat6);
  input[2].push_back(feat7);

	finder.run(input, result);
	TEST_EQUAL(result.size(), 4);
	ABORT_IF(result.size() != 4);

  FeatureHandle ind7(2, feat7);

  group1 = result[0].getFeatures();
	it = group1.begin();
  STATUS(*it);
	STATUS(ind2);
	TEST_EQUAL(*(it) == ind2, true);
	++it;
  STATUS(*it);
	STATUS(ind5);
  TEST_EQUAL(*(it) == ind5, true);
	++it;
  STATUS(*it);
	// "ind6" is closer, but its annotation doesn't match
	STATUS(ind7);
  TEST_EQUAL(*(it) == ind7, true);
}
END_SECTION

START_SECTION((void run(const std::vector<ConsensusMap>& input_maps, ConsensusMap& result_map)))
{
	NOT_TESTABLE; // same as "run" for feature maps (tested above)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
