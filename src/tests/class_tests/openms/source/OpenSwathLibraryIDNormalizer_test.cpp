// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Exception.h>

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathLibraryIDNormalizer.h>
///////////////////////////

#include <string>
#include <utility>
#include <vector>

using namespace OpenMS;
using namespace std;

START_TEST(OpenSwathLibraryIDNormalizer, "$Id$")

START_SECTION(static void normalizeSourceIDs(OpenSwath::LightTargetedExperiment& exp))
{
  OpenSwath::LightTargetedExperiment exp;

  for (const char* id : {"490", "A", "0", "2"})
  {
    OpenSwath::LightCompound compound;
    compound.id = id;
    exp.compounds.push_back(compound);
  }

  const std::vector<std::pair<std::string, std::string>> transition_ids_and_refs = {
    {"2292", "490"},
    {"transition_A", "A"},
    {"12", "0"},
    {"0", "2"}
  };
  for (const auto& id_and_ref : transition_ids_and_refs)
  {
    OpenSwath::LightTransition transition;
    transition.transition_name = id_and_ref.first;
    transition.peptide_ref = id_and_ref.second;
    exp.transitions.push_back(transition);
  }

  const auto source_ids = OpenSwathLibraryIDNormalizer::normalizeSourceIDs(exp);

  TEST_EQUAL(source_ids.precursor_source_to_canonical.at("0"), "0")
  TEST_EQUAL(source_ids.precursor_source_to_canonical.at("2"), "1")
  TEST_EQUAL(source_ids.precursor_source_to_canonical.at("490"), "2")
  TEST_EQUAL(source_ids.precursor_source_to_canonical.at("A"), "3")
  TEST_EQUAL(source_ids.precursor_canonical_to_source.at("0"), "0")
  TEST_EQUAL(source_ids.precursor_canonical_to_source.at("1"), "2")
  TEST_EQUAL(source_ids.precursor_canonical_to_source.at("2"), "490")
  TEST_EQUAL(source_ids.precursor_canonical_to_source.at("3"), "A")
  TEST_EQUAL(source_ids.transition_canonical_to_source.at("0"), "2292")
  TEST_EQUAL(source_ids.transition_canonical_to_source.at("1"), "transition_A")
  TEST_EQUAL(source_ids.transition_canonical_to_source.at("2"), "12")
  TEST_EQUAL(source_ids.transition_canonical_to_source.at("3"), "0")

  // Source precursor IDs are sorted lexicographically: 0, 2, 490, A.
  TEST_EQUAL(exp.compounds[0].id, "2")
  TEST_EQUAL(exp.compounds[1].id, "3")
  TEST_EQUAL(exp.compounds[2].id, "0")
  TEST_EQUAL(exp.compounds[3].id, "1")

  // Transition IDs follow stable vector order and refs follow the precursor remap.
  TEST_EQUAL(exp.transitions[0].transition_name, "0")
  TEST_EQUAL(exp.transitions[0].peptide_ref, "2")
  TEST_EQUAL(exp.transitions[1].transition_name, "1")
  TEST_EQUAL(exp.transitions[1].peptide_ref, "3")
  TEST_EQUAL(exp.transitions[2].transition_name, "2")
  TEST_EQUAL(exp.transitions[2].peptide_ref, "0")
  TEST_EQUAL(exp.transitions[3].transition_name, "3")
  TEST_EQUAL(exp.transitions[3].peptide_ref, "1")

  // The normalized experiment satisfies the canonical invariant.
  OpenSwathLibraryIDNormalizer::validateCanonicalIDs(exp);
}
END_SECTION

START_SECTION(static void normalizeSourceIDs(OpenSwath::LightTargetedExperiment& exp) -- removes unreferenced compounds without compressing canonical IDs)
{
  OpenSwath::LightTargetedExperiment exp;

  // PEPTIDEA_Extra sorts between the two active IDs but is not referenced by any transition.
  // It reserves canonical ID 1 before being removed from the operational experiment, so
  // PEPTIDEB keeps ID 2 exactly as it would after a source-to-PQP round-trip.
  for (const char* id : {"PEPTIDEA", "PEPTIDEA_Extra", "PEPTIDEB"})
  {
    OpenSwath::LightCompound compound;
    compound.id = id;
    exp.compounds.push_back(compound);
  }

  OpenSwath::LightTransition tr_a;
  tr_a.transition_name = "tr_a";
  tr_a.peptide_ref = "PEPTIDEA";
  exp.transitions.push_back(tr_a);

  OpenSwath::LightTransition tr_b;
  tr_b.transition_name = "tr_b";
  tr_b.peptide_ref = "PEPTIDEB";
  exp.transitions.push_back(tr_b);

  OpenSwathLibraryIDNormalizer::normalizeSourceIDs(exp);

  TEST_EQUAL(exp.compounds.size(), 2)
  TEST_EQUAL(exp.compounds[0].id, "0")
  TEST_EQUAL(exp.compounds[1].id, "2")
  TEST_EQUAL(exp.transitions[0].transition_name, "0")
  TEST_EQUAL(exp.transitions[0].peptide_ref, "0")
  TEST_EQUAL(exp.transitions[1].transition_name, "1")
  TEST_EQUAL(exp.transitions[1].peptide_ref, "2")

  OpenSwathLibraryIDNormalizer::validateCanonicalIDs(exp);
}
END_SECTION

START_SECTION(static void normalizeSourceIDs(OpenSwath::LightTargetedExperiment& exp) -- rebuilds a primed compound-reference cache)
{
  OpenSwath::LightTargetedExperiment exp;

  OpenSwath::LightCompound compound_a;
  compound_a.id = "PEPTIDEA";
  exp.compounds.push_back(compound_a);

  OpenSwath::LightCompound compound_b;
  compound_b.id = "PEPTIDEB";
  exp.compounds.push_back(compound_b);

  OpenSwath::LightTransition tr_a;
  tr_a.transition_name = "tr_a";
  tr_a.peptide_ref = "PEPTIDEA";
  exp.transitions.push_back(tr_a);

  OpenSwath::LightTransition tr_b;
  tr_b.transition_name = "tr_b";
  tr_b.peptide_ref = "PEPTIDEB";
  exp.transitions.push_back(tr_b);

  // Prime LightTargetedExperiment's internal reference lookup before normalization.
  TEST_EQUAL(exp.getCompoundByRef("PEPTIDEA").id, "PEPTIDEA")

  OpenSwathLibraryIDNormalizer::normalizeSourceIDs(exp);

  // A stale source-ID lookup cache used to make this dereference invalid during MS1 extraction.
  TEST_EQUAL(exp.getCompoundByRef("0").id, "0")
  TEST_EQUAL(exp.getCompoundByRef("1").id, "1")
  TEST_EQUAL(exp.transitions[0].peptide_ref, "0")
  TEST_EQUAL(exp.transitions[1].peptide_ref, "1")
}
END_SECTION

START_SECTION(static void normalizeSourceIDs(OpenSwath::LightTargetedExperiment& exp) -- invalid source references)
{
  OpenSwath::LightTargetedExperiment duplicate_compounds;
  OpenSwath::LightCompound c1;
  c1.id = "duplicate";
  duplicate_compounds.compounds.push_back(c1);
  duplicate_compounds.compounds.push_back(c1);
  TEST_EXCEPTION(Exception::InvalidValue,
                 OpenSwathLibraryIDNormalizer::normalizeSourceIDs(duplicate_compounds))

  OpenSwath::LightTargetedExperiment missing_ref;
  OpenSwath::LightCompound c2;
  c2.id = "precursor_A";
  missing_ref.compounds.push_back(c2);
  OpenSwath::LightTransition tr;
  tr.transition_name = "transition_A";
  tr.peptide_ref = "missing_precursor";
  missing_ref.transitions.push_back(tr);
  TEST_EXCEPTION(Exception::InvalidValue,
                 OpenSwathLibraryIDNormalizer::normalizeSourceIDs(missing_ref))
}
END_SECTION

START_SECTION(static SourceIDMapping normalizeSourceIDs(...) -- preserves source-prefix decoy semantics)
{
  OpenSwath::LightTargetedExperiment exp;

  for (const char* id : {"PEPTIDE", "DECOY_PEPTIDE", "REV_PEPTIDE"})
  {
    OpenSwath::LightCompound compound;
    compound.id = id;
    exp.compounds.push_back(compound);
  }

  OpenSwath::LightTransition target;
  target.transition_name = "target_y7";
  target.peptide_ref = "PEPTIDE";
  exp.transitions.push_back(target);

  OpenSwath::LightTransition conventional_decoy;
  conventional_decoy.transition_name = "decoy_y7";
  conventional_decoy.peptide_ref = "DECOY_PEPTIDE";
  exp.transitions.push_back(conventional_decoy);

  OpenSwath::LightTransition custom_decoy;
  custom_decoy.transition_name = "rev_y7";
  custom_decoy.peptide_ref = "REV_PEPTIDE";
  exp.transitions.push_back(custom_decoy);

  OpenSwath::LightTransition custom_transition_decoy;
  custom_transition_decoy.transition_name = "REV_transition_y8";
  custom_transition_decoy.peptide_ref = "PEPTIDE";
  exp.transitions.push_back(custom_transition_decoy);

  const auto source_ids = OpenSwathLibraryIDNormalizer::normalizeSourceIDs(exp);

  // Conventional DECOY*/Decoy*/decoy* semantics are materialized while the
  // source IDs are still available. The custom prefix is materialized through
  // retained provenance before downstream filtering.
  TEST_EQUAL(exp.transitions[0].getDecoy(), false)
  TEST_EQUAL(exp.transitions[1].getDecoy(), true)
  TEST_EQUAL(exp.transitions[2].getDecoy(), false)
  TEST_EQUAL(exp.transitions[3].getDecoy(), false)

  OpenSwathLibraryIDNormalizer::materializeDecoyPrefix(exp, source_ids, "REV_");
  TEST_EQUAL(exp.transitions[2].getDecoy(), true)
  TEST_EQUAL(exp.transitions[3].getDecoy(), true)
  TEST_EQUAL(source_ids.transition_canonical_to_source.at("3"), "REV_transition_y8")

  const std::string canonical_decoy = source_ids.precursor_source_to_canonical.at("DECOY_PEPTIDE");
  TEST_EQUAL(source_ids.precursor_canonical_to_source.at(canonical_decoy), "DECOY_PEPTIDE")
  const std::string canonical_custom_decoy = source_ids.precursor_source_to_canonical.at("REV_PEPTIDE");
  TEST_EQUAL(source_ids.precursor_canonical_to_source.at(canonical_custom_decoy), "REV_PEPTIDE")

  const std::string canonical_target = source_ids.precursor_source_to_canonical.at("PEPTIDE");
  const auto paired_target = OpenSwathLibraryIDNormalizer::canonicalTargetForDecoyPrecursor(
    canonical_decoy, source_ids, "DECOY_");
  TEST_EQUAL(paired_target.has_value(), true)
  TEST_EQUAL(*paired_target, canonical_target)

  const auto paired_custom_target = OpenSwathLibraryIDNormalizer::canonicalTargetForDecoyPrecursor(
    canonical_custom_decoy, source_ids, "REV_");
  TEST_EQUAL(paired_custom_target.has_value(), true)
  TEST_EQUAL(*paired_custom_target, canonical_target)

  TEST_EQUAL(OpenSwathLibraryIDNormalizer::canonicalTargetForDecoyPrecursor(
               canonical_target, source_ids, "DECOY_").has_value(),
             false)
}
END_SECTION

START_SECTION(static bool hasCanonicalIDs(const OpenSwath::LightTargetedExperiment& exp))
{
  OpenSwath::LightTargetedExperiment valid;

  OpenSwath::LightCompound c0;
  c0.id = "0";
  valid.compounds.push_back(c0);

  OpenSwath::LightCompound c7;
  c7.id = "7";
  valid.compounds.push_back(c7);

  OpenSwath::LightTransition t0;
  t0.transition_name = "3";
  t0.peptide_ref = "0";
  valid.transitions.push_back(t0);

  OpenSwath::LightTransition t1;
  t1.transition_name = "100";
  t1.peptide_ref = "7";
  valid.transitions.push_back(t1);

  TEST_EQUAL(OpenSwathLibraryIDNormalizer::hasCanonicalIDs(valid), true)
  TEST_EQUAL(OpenSwathLibraryIDNormalizer::hasCanonicalIDFormat(valid), true)

  OpenSwath::LightTargetedExperiment duplicate_transition = valid;
  duplicate_transition.transitions[1].transition_name = "3";
  TEST_EQUAL(OpenSwathLibraryIDNormalizer::hasCanonicalIDs(duplicate_transition), false)
  TEST_EQUAL(OpenSwathLibraryIDNormalizer::hasCanonicalIDFormat(duplicate_transition), true)

  OpenSwath::LightTargetedExperiment source_style = valid;
  source_style.compounds[0].id = "PEPTIDEA";
  source_style.transitions[0].peptide_ref = "PEPTIDEA";
  TEST_EQUAL(OpenSwathLibraryIDNormalizer::hasCanonicalIDs(source_style), false)
  TEST_EQUAL(OpenSwathLibraryIDNormalizer::hasCanonicalIDFormat(source_style), false)
}
END_SECTION

START_SECTION(static void validateCanonicalIDs(const OpenSwath::LightTargetedExperiment& exp))
{
  // Sparse canonical IDs are valid and must not be renumbered after filtering.
  OpenSwath::LightTargetedExperiment valid;
  for (const char* id : {"0", "7", "15"})
  {
    OpenSwath::LightCompound compound;
    compound.id = id;
    valid.compounds.push_back(compound);
  }

  for (const auto& id_and_ref : std::vector<std::pair<std::string, std::string>>{
         {"3", "0"}, {"100", "15"}, {"902", "7"}})
  {
    OpenSwath::LightTransition transition;
    transition.transition_name = id_and_ref.first;
    transition.peptide_ref = id_and_ref.second;
    valid.transitions.push_back(transition);
  }
  OpenSwathLibraryIDNormalizer::validateCanonicalIDs(valid);

  OpenSwath::LightTargetedExperiment leading_zero = valid;
  leading_zero.compounds[0].id = "007";
  leading_zero.transitions[0].peptide_ref = "007";
  TEST_EXCEPTION(Exception::InvalidValue,
                 OpenSwathLibraryIDNormalizer::validateCanonicalIDs(leading_zero))

  OpenSwath::LightTargetedExperiment negative = valid;
  negative.transitions[0].transition_name = "-1";
  TEST_EXCEPTION(Exception::InvalidValue,
                 OpenSwathLibraryIDNormalizer::validateCanonicalIDs(negative))

  OpenSwath::LightTargetedExperiment duplicate_transition = valid;
  duplicate_transition.transitions[1].transition_name = "3";
  TEST_EXCEPTION(Exception::InvalidValue,
                 OpenSwathLibraryIDNormalizer::validateCanonicalIDs(duplicate_transition))

  OpenSwath::LightTargetedExperiment unknown_ref = valid;
  unknown_ref.transitions[0].peptide_ref = "999";
  TEST_EXCEPTION(Exception::InvalidValue,
                 OpenSwathLibraryIDNormalizer::validateCanonicalIDs(unknown_ref))
}
END_SECTION

END_TEST
