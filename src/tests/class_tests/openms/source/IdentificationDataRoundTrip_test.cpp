// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/IdentificationDataFile.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/METADATA/ID/IdentificationDataAdapter.h>
#include <OpenMS/test_config.h>
#include <algorithm>
#include <filesystem>
#include <fstream>

// Legacy -> native -> legacy round trips of the golden outputs of TOPP tools. Every identification of an
// output must come back exactly as it was, in the same order, after an import, a native bundle and an export.

using namespace OpenMS;
using ID = IdentificationData;
using Adapter = IdentificationDataAdapter;
namespace fs = std::filesystem;

namespace
{
std::string topp(const std::string& name)
{ return OPENMS_GET_TEST_DATA_PATH("../../../topp/" + name); }

/// The first line in which two text files differ (empty if they are equal), for a readable failure.
std::string firstDifference(const std::string& left, const std::string& right)
{
  std::ifstream a(left), b(right);
  std::string x, y;
  for (Size line = 1;; ++line)
  {
    const bool more_a = static_cast<bool>(std::getline(a, x));
    const bool more_b = static_cast<bool>(std::getline(b, y));
    if (! more_a && ! more_b) return "";
    if (! more_a || ! more_b || x != y) return "line " + std::to_string(line) + ":\n  < " + (more_a ? x : "<end>") + "\n  > " + (more_b ? y : "<end>");
  }
}

/// A new temporary file name (NEW_TMP_FILE names files by source line, so a helper needs a counter).
std::string temporary(const std::string& extension)
{
  static Size count = 0;
  std::string path;
  NEW_TMP_FILE_EXT(path, "_" + std::to_string(++count) + extension)
  return path;
}

/// The meta values in which two objects differ, e.g. "user: 1 (int) != 1 (float)".
std::string metaDifference(const MetaInfoInterface& left, const MetaInfoInterface& right)
{
  std::vector<std::string> keys, other;
  left.getKeys(keys);
  right.getKeys(other);
  keys.insert(keys.end(), other.begin(), other.end());
  std::sort(keys.begin(), keys.end());
  keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
  const auto text = [](const MetaInfoInterface& item, const std::string& key) -> std::string {
    if (! item.metaValueExists(key)) return "<none>";
    const auto& value = item.getMetaValue(key);
    return value.toString() + " (" + DataValue::NamesOfDataType[value.valueType()] + ", unit type " + std::to_string(value.getUnitType())
           + ", unit " + std::to_string(value.getUnit()) + ")";
  };
  std::string result;
  for (const auto& key : keys)
    if (! left.metaValueExists(key) || ! right.metaValueExists(key) || left.getMetaValue(key) != right.getMetaValue(key))
      result += " meta '" + key + "': " + text(left, key) + " != " + text(right, key) + ";";
  return result;
}

/// What differs between two legacy identification lists, field by field (empty if they are equal).
std::string difference(const std::vector<ProteinIdentification>& proteins, const PeptideIdentificationList& peptides,
                       const std::vector<ProteinIdentification>& other_proteins, const PeptideIdentificationList& other_peptides)
{
  std::string result;
  for (Size i = 0; i < std::min(proteins.size(), other_proteins.size()); ++i)
  {
    const auto& a = proteins[i];
    const auto& b = other_proteins[i];
    if (a == b) continue;
    result += "protein run " + std::to_string(i) + ":";
    if (a.getIdentifier() != b.getIdentifier()) result += " identifier " + a.getIdentifier() + " != " + b.getIdentifier() + ";";
    if (a.getSearchEngine() != b.getSearchEngine() || a.getSearchEngineVersion() != b.getSearchEngineVersion()) result += " search engine;";
    if (a.getSearchParameters() != b.getSearchParameters())
      result += " search parameters" + metaDifference(a.getSearchParameters(), b.getSearchParameters()) + ";";
    if (a.getDateTime() != b.getDateTime()) result += " date;";
    if (a.getScoreType() != b.getScoreType() || a.isHigherScoreBetter() != b.isHigherScoreBetter()) result += " score type;";
    if (a.getSignificanceThreshold() != b.getSignificanceThreshold()) result += " threshold;";
    if (a.getProteinGroups() != b.getProteinGroups() || a.getIndistinguishableProteins() != b.getIndistinguishableProteins()) result += " groups;";
    result += metaDifference(a, b);
    if (a.getHits().size() != b.getHits().size()) result += " hit count;";
    for (Size h = 0; h < std::min(a.getHits().size(), b.getHits().size()); ++h)
      if (a.getHits()[h] != b.getHits()[h])
      {
        result += " hit " + a.getHits()[h].getAccession() + ":" + metaDifference(a.getHits()[h], b.getHits()[h]) + ";";
        break;
      }
    result += "\n";
  }
  for (Size i = 0; i < std::min(peptides.size(), other_peptides.size()); ++i)
  {
    const auto& a = peptides[i];
    const auto& b = other_peptides[i];
    if (a == b) continue;
    result += "peptide identification " + std::to_string(i) + ":";
    if (a.getIdentifier() != b.getIdentifier()) result += " identifier " + a.getIdentifier() + " != " + b.getIdentifier() + ";";
    if (a.getScoreType() != b.getScoreType() || a.isHigherScoreBetter() != b.isHigherScoreBetter()) result += " score type;";
    if (a.getSignificanceThreshold() != b.getSignificanceThreshold()) result += " threshold;";
    if (a.getRT() != b.getRT() || a.getMZ() != b.getMZ()) result += " position;";
    result += metaDifference(a, b);
    if (a.getHits().size() != b.getHits().size()) result += " hit count;";
    for (Size h = 0; h < std::min(a.getHits().size(), b.getHits().size()); ++h)
    {
      const auto& x = a.getHits()[h];
      const auto& y = b.getHits()[h];
      if (x == y) continue;
      result += " hit " + std::to_string(h) + ":";
      if (x.getSequence() != y.getSequence()) result += " sequence " + x.getSequence().toString() + " != " + y.getSequence().toString() + ";";
      if (x.getScore() != y.getScore() || x.getRank() != y.getRank() || x.getCharge() != y.getCharge()) result += " score/rank/charge;";
      if (x.getPeptideEvidences() != y.getPeptideEvidences()) result += " evidence;";
      if (x.getPeakAnnotations() != y.getPeakAnnotations()) result += " peak annotations;";
      result += metaDifference(x, y);
      break;
    }
    result += "\n";
    if (result.size() > 2000) break;
  }
  return result;
}

/// Native bundle round trip of @p data, which must not change it.
ID throughBundle(const ID& data)
{
  const auto path = temporary(".idparquet");
  fs::remove_all(path);
  IdentificationDataFile::store(path, data);
  ID reloaded;
  IdentificationDataFile::load(path, reloaded);
  fs::remove_all(path);
  if (! (reloaded == data)) throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The native bundle changed the dataset", path);
  return reloaded;
}

/// idXML text of identifications, as IdXMLFile writes them.
std::string idXML(const std::vector<ProteinIdentification>& proteins, const PeptideIdentificationList& peptides)
{
  const auto path = temporary(".idXML");
  IdXMLFile().store(path, proteins, peptides);
  return path;
}

std::string consensusXML(const ConsensusMap& map)
{
  const auto path = temporary(".consensusXML");
  ConsensusXMLFile().store(path, map);
  return path;
}

/// The identifications of a map as links to the map's own IdentificationData, without legacy copies.
ConsensusMap toNative(const ConsensusMap& map)
{
  auto imported = Adapter::fromConsensusMap(map);
  ConsensusMap native = map;
  native.getProteinIdentifications().clear();
  native.getUnassignedPeptideIdentifications().clear();
  for (auto& feature : native)
    feature.getPeptideIdentifications().clear();
  for (const auto& association : imported.associations)
  {
    if (association.unassigned) continue;
    auto& feature = native[association.feature_path.front()];
    feature.addIDQuery(association.query);
    for (const auto& match : association.matches)
      feature.addIDMatch({association.query.run_uuid, match});
  }
  native.getIdentificationData() = throughBundle(imported.data);
  return native;
}

/// Legacy identifications rebuilt from the native links of @p native.
ConsensusMap toLegacy(const ConsensusMap& native)
{
  auto imported = Adapter::fromConsensusMap(native);
  ConsensusMap legacy = native;
  Adapter::applyToConsensusMap(imported.data, imported.associations, legacy, {}, Adapter::MissingLinkPolicy::REJECT);
  // Only the legacy copies are compared; drop the links again.
  legacy.getIdentificationData().clear();
  for (auto& feature : legacy)
  {
    feature.getIDQueries().clear();
    feature.getIDMatches().clear();
  }
  return legacy;
}
} // namespace

START_TEST(IdentificationDataRoundTrip, "$Id$")

START_SECTION(([EXTRA] idXML outputs of IDMerger, PercolatorAdapter, ProteomicsLFQ and IsobaricWorkflow survive a native round trip))
{
  for (const std::string name : {"IDMerger_1_output.idXML", "IDMerger_2_output.idXML", "IDMerger_3_output.idXML", "IDMerger_4_output.idXML",
                                 "IDMerger_5_output.idXML", "IDMerger_6_output.idXML", "IDMerger_idparquet_output.idXML",
                                 "THIRDPARTY/PercolatorAdapter_1.idXML", "THIRDPARTY/PercolatorAdapter_1_inproc_out.idXML",
                                 "THIRDPARTY/PercolatorAdapter_1_subprocess_out.idXML", "THIRDPARTY/PercolatorAdapter_idparquet_out.idXML",
                                 "IsobaricWorkflow_1_input.idXML", "IsobaricWorkflow_2_input.idXML", "IsobaricWorkflow_duplicate_psm_input.idXML",
                                 "IsobaricWorkflow_strictly_unique_input.idXML", "ProteomicsLFQ_duplicate_psm_1.idXML",
                                 "ProteomicsLFQ_empty_fraction_BSA1_F2.idXML"})
  {
    STATUS(name)
    std::vector<ProteinIdentification> proteins;
    PeptideIdentificationList peptides;
    IdXMLFile().load(topp(name), proteins, peptides);
    const auto native = throughBundle(Adapter::fromLegacy(proteins, peptides));
    const auto exported = Adapter::toLegacy(native); // strict: nothing may be lost
    TEST_EQUAL(exported.proteins.size(), proteins.size())
    TEST_EQUAL(exported.peptides.size(), peptides.size())
    const bool equal = exported.proteins == proteins && exported.peptides == peptides;
    TEST_TRUE(equal)
    if (! equal)
    {
      TEST_EQUAL(firstDifference(idXML(proteins, peptides), idXML(exported.proteins, exported.peptides)), "")
      TEST_EQUAL(difference(proteins, peptides, exported.proteins, exported.peptides), "")
    }
  }
}
END_SECTION

START_SECTION(([EXTRA] consensusXML outputs of ProteomicsLFQ, IsobaricWorkflow and IDMerger survive native links))
{
  for (const std::string name : {"ProteomicsLFQ_1_out.consensusXML", "ProteomicsLFQ_1_subset_out.consensusXML", "ProteomicsLFQ_2_out.consensusXML",
                                 "ProteomicsLFQ_3_out.consensusXML", "ProteomicsLFQ_3_seedRT_out.consensusXML", "ProteomicsLFQ_4_out.consensusXML",
                                 "ProteomicsLFQ_5_out.consensusXML", "ProteomicsLFQ_6_out.consensusXML", "ProteomicsLFQ_7_out.consensusXML",
                                 "ProteomicsLFQ_8_out.consensusXML", "ProteomicsLFQ_9_out.consensusXML", "ProteomicsLFQ_9_seedRT_out.consensusXML",
                                 "IsobaricWorkflow_1_out.consensusXML", "IsobaricWorkflow_strictly_unique_out.consensusXML",
                                 "IDMerger_6_output.consensusXML"})
  {
    STATUS(name)
    ConsensusMap original;
    ConsensusXMLFile().load(topp(name), original);
    const auto restored = toLegacy(toNative(original));
    bool equal = restored.getProteinIdentifications() == original.getProteinIdentifications()
                 && restored.getUnassignedPeptideIdentifications() == original.getUnassignedPeptideIdentifications() && restored.size() == original.size();
    for (Size i = 0; equal && i < original.size(); ++i)
      equal = restored[i].getPeptideIdentifications() == original[i].getPeptideIdentifications();
    TEST_TRUE(equal)
    if (! equal)
    {
      TEST_EQUAL(firstDifference(consensusXML(original), consensusXML(restored)), "")
      TEST_EQUAL(difference(original.getProteinIdentifications(), original.getUnassignedPeptideIdentifications(),
                            restored.getProteinIdentifications(), restored.getUnassignedPeptideIdentifications()), "")
      for (Size i = 0; i < std::min(original.size(), restored.size()); ++i)
        if (original[i].getPeptideIdentifications() != restored[i].getPeptideIdentifications())
        {
          TEST_EQUAL(difference({}, original[i].getPeptideIdentifications(), {}, restored[i].getPeptideIdentifications()), "")
          break;
        }
    }
  }
}
END_SECTION

END_TEST
