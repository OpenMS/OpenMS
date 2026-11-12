// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/TestFileValidation.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricNormalizer.h>
///////////////////////////

#include <OpenMS/ANALYSIS/QUANTITATION/ItraqFourPlexQuantitationMethod.h>

#include <OpenMS/FORMAT/ConsensusXMLFile.h>

using namespace OpenMS;
using namespace std;

START_TEST(IsobaricNormalizer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IsobaricNormalizer* ptr = nullptr;
IsobaricNormalizer* null_ptr = nullptr;

// 
ItraqFourPlexQuantitationMethod quant_meth;

START_SECTION((IsobaricNormalizer(const IsobaricQuantitationMethod *const quant_method)))
{
	ptr = new IsobaricNormalizer(&quant_meth);
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~IsobaricNormalizer())
{
	delete ptr;
}
END_SECTION

START_SECTION((IsobaricNormalizer(const IsobaricNormalizer &other)))
{
  IsobaricNormalizer normalizer(&quant_meth);
  IsobaricNormalizer * normalizer2 = new IsobaricNormalizer(normalizer);
  
	TEST_NOT_EQUAL(normalizer2, null_ptr)  
  delete normalizer2;
	
  // equality cannot be checked
  NOT_TESTABLE
}
END_SECTION

START_SECTION((IsobaricNormalizer& operator=(const IsobaricNormalizer &rhs)))
{
  IsobaricNormalizer normalizer(&quant_meth);
  IsobaricNormalizer normalizer2(&quant_meth);
  
  normalizer2 = normalizer;
  
  // equality cannot be checked
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void normalize(ConsensusMap &consensus_map)))
{
  IsobaricNormalizer normalizer(&quant_meth);
    
	ConsensusXMLFile cm_file;
	ConsensusMap cm_in;
	cm_file.load(OPENMS_GET_TEST_DATA_PATH("IsobaricNormalizer.consensusXML"),cm_in);

  normalizer.normalize(cm_in);

	std::string cm_file_out;
	NEW_TMP_FILE(cm_file_out);
	cm_file.store(cm_file_out,cm_in);

	WHITELIST("<?xml-stylesheet");
	TEST_FILE_SIMILAR(cm_file_out,OPENMS_GET_TEST_DATA_PATH("IsobaricNormalizer_out.consensusXML"));
}
END_SECTION

START_SECTION([EXTRA] void normalize(ConsensusMap &consensus_map) with map indices that do not start at 0)
{
  // The channels of a later file of a multi-file experiment (as in IsobaricWorkflow) do not start at
  // map index 0. The reference channel used to be looked up by its map index in the per-channel
  // vectors, which read out of bounds for such a map.
  ConsensusMap cm_expected;
  ConsensusXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IsobaricNormalizer.consensusXML"), cm_expected);
  IsobaricNormalizer(&quant_meth).normalize(cm_expected);

  ConsensusMap cm_in;
  ConsensusXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IsobaricNormalizer.consensusXML"), cm_in);
  const Size offset = 4;
  ConsensusMap cm_shifted = cm_in;
  cm_shifted.getColumnHeaders().clear();
  for (const auto& [index, header] : cm_in.getColumnHeaders())
  {
    cm_shifted.getColumnHeaders()[index + offset] = header;
  }
  for (ConsensusFeature& cf : cm_shifted)
  {
    ConsensusFeature::HandleSetType handles = cf.getFeatures();
    cf.clear();
    for (FeatureHandle handle : handles)
    {
      handle.setMapIndex(handle.getMapIndex() + offset);
      cf.insert(handle);
    }
  }

  IsobaricNormalizer(&quant_meth).normalize(cm_shifted);

  TEST_EQUAL(cm_shifted.size(), cm_expected.size())
  for (Size i = 0; i < cm_expected.size(); ++i)
  {
    TEST_EQUAL(cm_shifted[i].size(), cm_expected[i].size())
    auto it_shifted = cm_shifted[i].begin();
    for (auto it = cm_expected[i].begin(); it != cm_expected[i].end(); ++it, ++it_shifted)
    {
      TEST_EQUAL(it_shifted->getMapIndex(), it->getMapIndex() + offset)
      TEST_REAL_SIMILAR(it_shifted->getIntensity(), it->getIntensity())
    }
  }
}
END_SECTION

START_SECTION([EXTRA] void normalize(ConsensusMap &consensus_map) with a channel without signal)
{
  // An unused channel has a median ratio of 0 to the reference channel. Dividing by that turned its
  // intensities into NaN; the channel is now left as it is. The other channels are normalized as usual.
  ConsensusMap cm_expected;
  ConsensusXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IsobaricNormalizer.consensusXML"), cm_expected);
  IsobaricNormalizer(&quant_meth).normalize(cm_expected);

  const Size empty_map = 3;
  ConsensusMap cm_empty;
  ConsensusXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IsobaricNormalizer.consensusXML"), cm_empty);
  for (ConsensusFeature& cf : cm_empty)
  {
    ConsensusFeature::HandleSetType handles = cf.getFeatures();
    cf.clear();
    for (FeatureHandle handle : handles)
    {
      if (handle.getMapIndex() == empty_map)
      {
        handle.setIntensity(0);
      }
      cf.insert(handle);
    }
  }

  IsobaricNormalizer(&quant_meth).normalize(cm_empty);

  TEST_EQUAL(cm_empty.size(), cm_expected.size())
  for (Size i = 0; i < cm_expected.size(); ++i)
  {
    auto it_empty = cm_empty[i].begin();
    for (auto it = cm_expected[i].begin(); it != cm_expected[i].end(); ++it, ++it_empty)
    {
      if (it->getMapIndex() == empty_map)
      {
        TEST_EQUAL(it_empty->getIntensity(), 0)
      }
      else
      {
        TEST_REAL_SIMILAR(it_empty->getIntensity(), it->getIntensity())
      }
    }
  }
}
END_SECTION

START_SECTION([EXTRA] void normalize(ConsensusMap &consensus_map) with a reference channel without signal)
{
  // Without signal in the reference channel, the ratio of every other channel is the stand-in for x/0
  // (the largest float). Dividing by that turned all intensities into (almost) 0, without a warning.
  // The channels are now left as they are (the reference channel is set to 1, as always).
  const Size ref_map = 0; // channel 114, the default reference of iTRAQ 4-plex
  ConsensusMap cm_in;
  ConsensusXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IsobaricNormalizer.consensusXML"), cm_in);
  for (ConsensusFeature& cf : cm_in)
  {
    ConsensusFeature::HandleSetType handles = cf.getFeatures();
    cf.clear();
    for (FeatureHandle handle : handles)
    {
      if (handle.getMapIndex() == ref_map)
      {
        handle.setIntensity(0);
      }
      cf.insert(handle);
    }
  }
  ConsensusMap cm_normalized = cm_in;

  IsobaricNormalizer(&quant_meth).normalize(cm_normalized);

  for (Size i = 0; i < cm_in.size(); ++i)
  {
    auto it_normalized = cm_normalized[i].begin();
    for (auto it = cm_in[i].begin(); it != cm_in[i].end(); ++it, ++it_normalized)
    {
      if (it->getMapIndex() == ref_map)
      {
        TEST_EQUAL(it_normalized->getIntensity(), 1)
      }
      else
      {
        TEST_REAL_SIMILAR(it_normalized->getIntensity(), it->getIntensity())
      }
    }
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/// check the temporary files written above against their XML schema (types without a validator are skipped)
VALIDATE_TMP_FILES

END_TEST
