// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/KERNEL/StandardTypes.h>

///////////////////////////
#include <OpenMS/KERNEL/ConversionHelper.h>
///////////////////////////

#include <algorithm>
#include <vector>

using namespace OpenMS;
using namespace std;

START_TEST(ConsensusMap, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((template < typename FeatureT > static void convert(UInt64 const input_map_index, FeatureMap< FeatureT > const &input_map, ConsensusMap &output_map, Size n=-1)))
{

  FeatureMap fm;
  Feature f;
  for ( UInt i = 0; i < 3; ++i )
  {
    f.setRT(i*77.7);
    f.setMZ(i+100.35);
    f.setUniqueId(i*33+17);
    fm.push_back(f);
  }
  ConsensusMap cm;
  MapConversion::convert(33,fm,cm);

  TEST_EQUAL(cm.size(),3);
  TEST_EQUAL(cm.getColumnHeaders()[33].size,3);
  for ( UInt i = 0; i < 3; ++i )
  {
    TEST_EQUAL(cm[i].size(),1);
    TEST_EQUAL(cm[i].begin()->getMapIndex(),33);
    TEST_EQUAL(cm[i].begin()->getUniqueId(),i*33+17);
    TEST_REAL_SIMILAR(cm[i].begin()->getRT(),i*77.7);
    TEST_REAL_SIMILAR(cm[i].begin()->getMZ(),i+100.35);
  }

cm.clear();
MapConversion::convert(33,fm,cm,2);
TEST_EQUAL(cm.size(),2);
TEST_EQUAL(cm.getColumnHeaders()[33].size,3);

}
END_SECTION

/////

ConsensusMap cm;
cm.setUniqueId();
for (UInt i = 0; i < 3; ++i)
{
  ConsensusFeature cf;
  cf.setRT(i * 5.0);
  cf.setMZ(i + 100.35);
  cf.setIntensity(900 + 5 * i);
  cf.setUniqueId();
  cm.push_back(cf);
}
cm.updateRanges();

START_SECTION((template < typename FeatureT > static void convert(ConsensusMap const &input_map, const bool keep_uids, FeatureMap< FeatureT > &output_map)))
{
    FeatureMap out_fm;
    MapConversion::convert(cm, true, out_fm);

    TEST_EQUAL(cm.getUniqueId(), out_fm.getUniqueId());
    TEST_EQUAL(cm.getProteinIdentifications().size(), out_fm.getProteinIdentifications().size());
    TEST_EQUAL(cm.getUnassignedPeptideIdentifications().size(), out_fm.getUnassignedPeptideIdentifications().size());
    TEST_EQUAL(cm.size(), out_fm.size());

    for (Size i = 0; i < cm.size(); ++i)
    {
        TEST_EQUAL(cm[i], out_fm[i]);
    }

    out_fm.clear();
    MapConversion::convert(cm, false, out_fm);
    TEST_NOT_EQUAL(cm.getUniqueId(), out_fm.getUniqueId());

    for (Size i = 0; i < cm.size(); ++i)
    {
        TEST_REAL_SIMILAR(cm[i].getRT(), out_fm[i].getRT());
        TEST_REAL_SIMILAR(cm[i].getMZ(), out_fm[i].getMZ());
        TEST_REAL_SIMILAR(cm[i].getIntensity(), out_fm[i].getIntensity());

        TEST_NOT_EQUAL(cm[i].getUniqueId(), out_fm[i].getUniqueId());
    }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



