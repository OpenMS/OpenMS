// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/METADATA/MetaInfoInterfaceUtils.h>
#include <OpenMS/METADATA/PeptideHit.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MetaInfoInterfaceUtils, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((template<typename T_In, T_Out> static T_Out findCommonMetaKeys(const T::const_iterator& start, const T::const_iterator& end, const float min_frequency = 100.0)))
{
  vector<PeptideHit> hits; // some class derived from MetaInfoInterface
  for (Size i = 0; i < 10; ++i)
  {
    PeptideHit h;
    h.setMetaValue("commonMeta1", i);
    h.setMetaValue("commonMeta2", i);
    if (i % 2 == 0) 
    {
      h.setMetaValue("meta50pc", i);
    }
    hits.push_back(h);
  }
  hits.back().setMetaValue("metaSingle", "single");
  
  // common keys for ALL entries (i.e. 100% min_frequency)
  {
    std::vector<String> common = MetaInfoInterfaceUtils::findCommonMetaKeys<std::vector<PeptideHit>, std::vector<String> >(hits.begin(), hits.end(), 100.0);
    TEST_EQUAL(common.size(), 2);
    ABORT_IF(common.size() != 2);
    TEST_EQUAL(common[0], "commonMeta1");
    TEST_EQUAL(common[1], "commonMeta2");
    
    // exceeds 100% --> should be corrected to 100% internally
    std::vector<String> common2 = MetaInfoInterfaceUtils::findCommonMetaKeys<std::vector<PeptideHit>, std::vector<String> >(hits.begin(), hits.end(), 1110.0);
    TEST_TRUE(common == common2);
  }
  
  // occurrence of at least 50 (i.e. 50% min_frequency)
  {
    std::set<String> set50 = MetaInfoInterfaceUtils::findCommonMetaKeys<std::vector<PeptideHit>, std::set<String> >(hits.begin(), hits.end(), 50.0);
    TEST_EQUAL(set50.size(), 3);
    ABORT_IF(set50.size() != 3);
    std::set<String> set50_expected;
    set50_expected.insert("commonMeta1");
    set50_expected.insert("commonMeta2");
    set50_expected.insert("meta50pc");
    TEST_TRUE(set50 == set50_expected);
  }

  // ALL keys (i.e. 0% min_frequency)
  {
    std::set<String> set0 = MetaInfoInterfaceUtils::findCommonMetaKeys<std::vector<PeptideHit>, std::set<String> >(hits.begin(), hits.end(), 0.0);
    TEST_EQUAL(set0.size(), 4);
    ABORT_IF(set0.size() != 4);
    std::set<String> set0_expected;
    set0_expected.insert("commonMeta1");
    set0_expected.insert("commonMeta2");
    set0_expected.insert("meta50pc");
    set0_expected.insert("metaSingle");
    TEST_TRUE(set0 == set0_expected);

    // exceeds 0% --> should be corrected to 0% internally
    std::set<String> set0_2 = MetaInfoInterfaceUtils::findCommonMetaKeys<std::vector<PeptideHit>, std::set<String> >(hits.begin(), hits.end(), -10.0);
    TEST_TRUE(set0 == set0_2);
  }
  

}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



