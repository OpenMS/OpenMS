// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
//
// Side-by-side regression suite for the in-process Percolator against the
// external `percolator` binary. Sections 2+ are gated on the environment
// variable PERCOLATOR_BINARY; without it, they no-op. Section 1 (PIN stamp
// parity) always runs.

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/ID/Percolator.h>
#include <OpenMS/ANALYSIS/ID/PercolatorTypes.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/PercolatorInfile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>
#include <cstdlib>
#include <iterator>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

using namespace OpenMS;
using namespace std;

///////////////////////////////////////////////////////////////////////////////
// Helpers
///////////////////////////////////////////////////////////////////////////////

namespace
{

/// Reference idXML: same input the TOPP PercolatorAdapter_1 test uses.
String inputIdxml()
{
  return OPENMS_GET_TEST_DATA_PATH("../../../topp/THIRDPARTY/PercolatorAdapter_1.idXML");
}

/// Derive min/max charge from all hits.
void deriveChargeRange(const PeptideIdentificationList& pids,
                       int& min_charge, int& max_charge)
{
  min_charge = 10;
  max_charge = 0;
  for (const auto& pid : pids)
  {
    for (const auto& hit : pid.getHits())
    {
      min_charge = std::min(min_charge, hit.getCharge());
      max_charge = std::max(max_charge, hit.getCharge());
    }
  }
}

/// Build the PIN feature set the adapter would use: standard + any
/// PSMFeatureExtractor-provided extras + Peptide + Proteins trailer.
StringList buildFeatureSet(const vector<ProteinIdentification>& prs,
                           int min_charge, int max_charge)
{
  StringList feature_set = PercolatorInfile::getStandardFeatureSet(min_charge, max_charge);
  if (!prs.empty())
  {
    const auto& sp = prs.front().getSearchParameters();
    if (sp.metaValueExists("extra_features"))
    {
      StringList extras = ListUtils::create<String>(
        sp.getMetaValue("extra_features").toString());
      feature_set.insert(feature_set.end(), extras.begin(), extras.end());
    }
  }
  feature_set.push_back("Peptide");
  feature_set.push_back("Proteins");
  return feature_set;
}

} // namespace

///////////////////////////////////////////////////////////////////////////////

START_TEST(Percolator_subprocess_parity, "$Id$")

///////////////////////////////////////////////////////////////////////////////
// Test 1: PIN file content equals stamped meta values.
//
// This guards the factored helper: PercolatorInfile::preparePin_ internally
// calls stampPinFeaturesOnHits and then serializes the stamped meta values
// to a .pin file. A regression where one of them drifts (e.g. someone adds a
// transformation only on one side) would silently corrupt the in-process
// path's training features.
//
// Approach: pick a realistic idXML (same as the PercolatorAdapter_1 TOPP
// test). Write the .pin via store() AND stamp a parallel copy in memory.
// Parse the .pin back and compare every non-trailing column against the
// stamped hit's meta value.
///////////////////////////////////////////////////////////////////////////////

START_SECTION([EXTRA] PIN file content equals stamped meta values)
{
  vector<ProteinIdentification> prs;
  PeptideIdentificationList pids;
  IdXMLFile().load(inputIdxml(), prs, pids);
  TEST_FALSE(prs.empty())
  TEST_FALSE(pids.empty())

  int min_charge, max_charge;
  deriveChargeRange(pids, min_charge, max_charge);
  TEST_TRUE(min_charge <= max_charge)

  StringList feature_set = buildFeatureSet(prs, min_charge, max_charge);
  const std::string enz = "no_enzyme";

  // (a) Write PIN via store(). This also internally stamps, then serializes.
  const String pin_path = File::getTemporaryFile();
  PercolatorInfile::store(pin_path, pids, feature_set, enz, min_charge, max_charge);

  // (b) Stamp an independent copy so we can read back the meta values
  //     without depending on store()'s internal copy.
  PeptideIdentificationList stamped = pids;
  const auto skipped = PercolatorInfile::stampPinFeaturesOnHits(
    stamped, enz, min_charge, max_charge);

  // (c) Parse the PIN file. Header must equal feature_set verbatim.
  TextFile txt;
  txt.load(pin_path);
  const size_t n_lines = std::distance(txt.begin(), txt.end());
  TEST_TRUE(n_lines >= 2)

  StringList header = ListUtils::create<String>(*txt.begin(), '\t');
  TEST_EQUAL(header.size(), feature_set.size())
  {
    bool header_matches = (header.size() == feature_set.size());
    for (size_t i = 0; header_matches && i < feature_set.size(); ++i)
    {
      if (header[i] != feature_set[i]) header_matches = false;
    }
    TEST_EQUAL(header_matches, true)
  }

  // (d) Build the expected row sequence in the same order preparePin_ writes:
  //     outer loop pids, inner loop hits, skipping any (pid, hit) in `skipped`.
  //     SpecIds are not unique per hit (all hits within a pid share the same
  //     SpecId), so row-by-row alignment is required.
  std::vector<const PeptideHit*> expected_rows;
  for (size_t i = 0; i < stamped.size(); ++i)
  {
    const auto& hits = stamped[i].getHits();
    for (size_t j = 0; j < hits.size(); ++j)
    {
      if (skipped.count({i, j})) continue;
      expected_rows.push_back(&hits[j]);
    }
  }

  const size_t pin_data_rows = n_lines - 1;
  TEST_EQUAL(pin_data_rows, expected_rows.size())

  // (e) Walk rows in parallel. Compare every column except the trailing
  //     Proteins (which embeds tabs and inflates the split count).
  const size_t last_compared = feature_set.size() - 1; // skip Proteins

  size_t compared_rows = 0;
  size_t diffs = 0;
  String first_mismatch;
  size_t row_idx = 0;
  for (auto it = txt.begin() + 1;
       it != txt.end() && row_idx < expected_rows.size();
       ++it, ++row_idx)
  {
    StringList cols = ListUtils::create<String>(*it, '\t');
    if (cols.size() < last_compared)
    {
      diffs++;
      if (first_mismatch.empty()) first_mismatch = "short row: " + *it;
      continue;
    }
    const PeptideHit& hit = *expected_rows[row_idx];

    for (size_t c = 0; c < last_compared; ++c)
    {
      const String& feat = feature_set[c];
      if (!hit.metaValueExists(feat))
      {
        diffs++;
        if (first_mismatch.empty())
          first_mismatch = "row " + String(row_idx) + ": missing meta value " + feat;
        break;
      }
      const String stamped_val = hit.getMetaValue(feat).toString();
      if (cols[c] != stamped_val)
      {
        diffs++;
        if (first_mismatch.empty())
        {
          first_mismatch = "row " + String(row_idx) + " col " + feat
            + ": pin='" + cols[c] + "' stamped='" + stamped_val + "'";
        }
        break;
      }
    }
    compared_rows++;
  }

  if (!first_mismatch.empty())
  {
    TEST_EQUAL(first_mismatch, String())  // makes the diff visible on failure
  }
  TEST_EQUAL(diffs, 0)
  TEST_TRUE(compared_rows > 20)  // sanity: we compared a meaningful number of rows
}
END_SECTION

///////////////////////////////////////////////////////////////////////////////

END_TEST
