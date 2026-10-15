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
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/NativeIdentificationTest.h>

///////////////////////////
#include <OpenMS/KERNEL/ConversionHelper.h>
///////////////////////////

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
    TEST_TRUE(cm.getIdentificationData() == out_fm.getIdentificationData());
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

START_SECTION(([EXTRA] conversions keep the identifications of the converted features))
{
  using namespace OpenMS::Internal::ClassTest;
  FeatureMap features;
  FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_1.featureXML"), features);
  auto input = identificationsOf(features);
  const auto with_subordinates = std::find_if(features.begin(), features.end(), [](const Feature& f) { return ! f.getSubordinates().empty(); });
  ABORT_IF(with_subordinates == features.end() || input.features.empty() || input.unassigned.empty())
  // the conversion to consensus features drops the identifications of subordinates
  input.features[{Size(with_subordinates - features.begin()), 0}] = input.features.begin()->second;
  setIdentifications(features, input);
  for (Size n : {Size(-1), Size(1), Size(0)})
  {
    // the identifications of the first n features, with the map index
    MapPeptideIdentifications expected;
    expected.proteins = input.proteins;
    expected.unassigned = input.unassigned;
    for (Size i = 0; i < std::min(n, features.size()); ++i)
    {
      PeptideIdentificationList peptides = input.at({i});
      if (peptides.empty()) continue;
      for (auto& peptide : peptides) peptide.setMetaValue("map_index", UInt64(7));
      expected.features[{i}] = peptides;
    }
    ConsensusMap out;
    MapConversion::convert(7, features, out, n);
    TEST_EQUAL(identificationDifference(expected, identificationsOf(out)), "")
  }

  // selected features, in the given order: the identifications of the first two identified features, swapped
  std::vector<Size> identified;
  for (Size i = 0; i < features.size(); ++i)
  {
    if (! input.at({i}).empty()) identified.push_back(i);
  }
  ABORT_IF(identified.size() < 2)
  const std::vector<Size> selected = {identified[1], identified[0]};
  MapPeptideIdentifications expected;
  expected.proteins = input.proteins;
  expected.unassigned = input.unassigned;
  for (Size i = 0; i < selected.size(); ++i)
  {
    PeptideIdentificationList peptides = input.at({selected[i]});
    for (auto& peptide : peptides) peptide.setMetaValue("map_index", UInt64(7));
    expected.features[{i}] = peptides;
  }
  ConsensusMap out;
  MapConversion::convert(7, features, out, selected);
  ABORT_IF(out.size() != 2)
  TEST_REAL_SIMILAR(out[0].getRT(), features[selected[0]].getRT())
  TEST_REAL_SIMILAR(out[1].getRT(), features[selected[1]].getRT())
  TEST_EQUAL(out.getColumnHeaders()[7].size, features.size())
  TEST_EQUAL(identificationDifference(expected, identificationsOf(out)), "")
  TEST_EXCEPTION(Exception::IndexOverflow, MapConversion::convert(7, features, out, std::vector<Size>{features.size()}))

  ConsensusMap consensus;
  ConsensusXMLFile().load(OPENMS_GET_TEST_DATA_PATH("ExperimentalDesign_ProteomicsLFQ_1_subset_out.consensusXML"), consensus);
  FeatureMap result;
  MapConversion::convert(consensus, true, result);
  TEST_EQUAL(identificationDifference(identificationsOf(consensus), identificationsOf(result)), "")
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



