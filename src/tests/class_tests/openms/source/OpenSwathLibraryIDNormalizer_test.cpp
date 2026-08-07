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

  for (const std::string& id : {"490", "A", "0", "2"})
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

  OpenSwathLibraryIDNormalizer::normalizeSourceIDs(exp);

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

START_SECTION(static void validateCanonicalIDs(const OpenSwath::LightTargetedExperiment& exp))
{
  // Sparse canonical IDs are valid and must not be renumbered after filtering.
  OpenSwath::LightTargetedExperiment valid;
  for (const std::string& id : {"0", "7", "15"})
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
