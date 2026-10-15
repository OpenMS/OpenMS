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
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>

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

START_SECTION(([EXTRA] conversions treat the identification data like peptide identifications))
{
  FeatureMap features;
  FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_1.featureXML"), features);
  const auto with_subordinates = std::find_if(features.begin(), features.end(), [](const Feature& f) { return ! f.getSubordinates().empty(); });
  const auto with_ids = std::find_if(features.begin(), features.end(), [](const Feature& f) { return ! f.getPeptideIdentifications().empty(); });
  ABORT_IF(with_subordinates == features.end() || with_ids == features.end() || features.getUnassignedPeptideIdentifications().empty())
  // the conversion to consensus features drops the identifications of subordinates
  with_subordinates->getSubordinates()[0].setPeptideIdentifications(with_ids->getPeptideIdentifications());
  FeatureMap native = features;
  IdentificationDataConverter::importFeatureIDs(native);
  for (Size n : {Size(-1), Size(1), Size(0)})
  {
    // identifications linked by converted features get the map index; unlinked ones stay, others are removed
    ConsensusMap out;
    MapConversion::convert(7, native, out, n);
    Size marked = 0, unmarked = 0;
    for (const auto& run : out.getIdentificationData().getRuns())
      for (const auto& source : run.getSources())
        for (const auto& query : source.identifications)
        {
          if (! query.metaValueExists("map_index")) ++unmarked;
          else if (query.getMetaValue("map_index") == DataValue(7)) ++marked;
        }
    Size assigned = 0;
    for (Size i = 0; i < std::min(n, features.size()); ++i)
      assigned += features[i].getPeptideIdentifications().size();
    TEST_EQUAL(marked, assigned)
    TEST_EQUAL(unmarked, features.getUnassignedPeptideIdentifications().size())
    // the links of the converted features resolve
    for (const auto& feature : out)
    {
      TEST_EQUAL(feature.getLinkedIdentifications(out.getIdentificationData()).size(), feature.getIDQueries().size())
    }
  }

  // consensus features to features: identification data and links are copied
  ConsensusMap consensus;
  ConsensusXMLFile().load(OPENMS_GET_TEST_DATA_PATH("ExperimentalDesign_ProteomicsLFQ_1_subset_out.consensusXML"), consensus);
  IdentificationDataConverter::importConsensusIDs(consensus);
  FeatureMap result;
  MapConversion::convert(consensus, true, result);
  TEST_TRUE(result.getIdentificationData() == consensus.getIdentificationData())
  ABORT_IF(result.size() != consensus.size())
  for (Size i = 0; i < result.size(); ++i)
  {
    TEST_TRUE(result[i].getIDQueries() == consensus[i].getIDQueries() && result[i].getIDMatches() == consensus[i].getIDMatches())
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



