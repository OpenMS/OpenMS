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
#include <OpenMS/KERNEL/StandardTypes.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/SimplePairFinder.h>
#include <OpenMS/KERNEL/FeatureMap.h>
///////////////////////////

#include <OpenMS/KERNEL/ConversionHelper.h>

using namespace OpenMS;
using namespace std;

typedef DPosition <2> PositionType;


START_TEST(SimplePairFinder, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SimplePairFinder* ptr = nullptr;
SimplePairFinder* nullPointer = nullptr;
BaseGroupFinder* base_nullPointer = nullptr;

START_SECTION((SimplePairFinder()))
	ptr = new SimplePairFinder();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~SimplePairFinder()))
	delete ptr;
END_SECTION

START_SECTION((static BaseGroupFinder* create()))
	BaseGroupFinder* base_ptr = nullptr;
	base_ptr = SimplePairFinder::create();
  TEST_NOT_EQUAL(base_ptr, base_nullPointer)
END_SECTION

START_SECTION((static const String getProductName()))
  SimplePairFinder spf;
  
  TEST_EQUAL(spf.getName() == "simple",true)
END_SECTION

START_SECTION((virtual void run(const std::vector< ConsensusMap > &input_maps, ConsensusMap &result_map)))
  FeatureMap scene;
  Feature feat1;
  Feature feat2;
  Feature feat3;
  PositionType pos1(0,0);
  PositionType pos2(200,300);
  PositionType pos3(400,500);
  feat1.setPosition(pos1);
  feat1.setIntensity(100.0f);
  feat1.setUniqueId(111111);
  feat2.setPosition(pos2);
  feat2.setIntensity(300.0f);
  feat2.setUniqueId(222222);
  feat3.setPosition(pos3);
  feat3.setIntensity(400.0f);
  feat3.setUniqueId(333333);
  scene.push_back(feat1);
  scene.push_back(feat2);
  scene.push_back(feat3);
  
  FeatureMap model;
  Feature feat4;
  Feature feat5;
  Feature feat6;
  PositionType pos4(4,4);
  PositionType pos5(204,304);
  PositionType pos6(404,504);
  feat4.setPosition(pos4);
  feat4.setIntensity(100.0f);
  feat4.setUniqueId(444444);
  feat5.setPosition(pos5);
  feat5.setIntensity(300.0f);
  feat5.setUniqueId(555555);
  feat6.setPosition(pos6);
  feat6.setIntensity(400.0f);
  feat6.setUniqueId(666666);
  model.push_back(feat4);
  model.push_back(feat5);
  model.push_back(feat6);
  
  SimplePairFinder spf;
	std::vector<ConsensusMap> input(2);
	MapConversion::convert(0,model,input[0]);
	MapConversion::convert(1,scene,input[1]);
	ConsensusMap result;
  spf.run(input,result);
	TEST_EQUAL(result.size(),3);
	ABORT_IF(result.size()!=3);

  ConsensusFeature::HandleSetType group1 = result[0].getFeatures();
  ConsensusFeature::HandleSetType group2 = result[1].getFeatures();
  ConsensusFeature::HandleSetType group3 = result[2].getFeatures();
  
  ConsensusFeature::HandleSetType::const_iterator it;
	it = group1.begin();
  STATUS(*it);
	TEST_EQUAL(it->getUniqueId(),444444);
	++it;
  STATUS(*it);
  TEST_EQUAL(it->getUniqueId(),111111);
	it = group2.begin();
  STATUS(*it);
  TEST_EQUAL(it->getUniqueId(),555555);
	++it;
  STATUS(*it);
  TEST_EQUAL(it->getUniqueId(),222222);
  it = group3.begin();
  STATUS(*it);
  TEST_EQUAL(it->getUniqueId(),666666);
	++it;
  STATUS(*it);
  TEST_EQUAL(it->getUniqueId(),333333);
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



