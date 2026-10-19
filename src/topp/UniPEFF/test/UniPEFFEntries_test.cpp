// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include "UniPEFFEntries.h"

using namespace OpenMS;
using namespace OpenMS::UniPEFFEntries;

namespace
{
  UniProtFeature position(const std::string& type, const std::string& description, int pos,
                          const std::string& isoform = "", const std::string& variation = "")
  {
    UniProtFeature f;
    f.type = type;
    f.description = description;
    f.location_sequence = isoform;
    f.has_position = true;
    f.position = pos;
    f.variation = variation;
    return f;
  }

  UniProtFeature range(const std::string& type, const std::string& description, int begin, int end,
                       const std::string& isoform = "", const std::string& variation = "")
  {
    UniProtFeature f;
    f.type = type;
    f.description = description;
    f.location_sequence = isoform;
    f.has_range = true;
    f.begin = begin;
    f.end = end;
    f.variation = variation;
    return f;
  }

  UniProtIsoform isoform(const std::string& id, const std::string& type, const std::string& ref = "")
  {
    return UniProtIsoform{id, id.substr(id.find('-') + 1), type, ref};
  }

  PtmMap ptms()
  {
    PtmMap map;
    map.emplace("Phosphoserine", PtmEntry{"Phosphoserine", "PTM-0253", "MOD:00046", "21"});
    map.emplace("Phosphothreonine", PtmEntry{"Phosphothreonine", "PTM-0254", "MOD:00047", "21"});
    map.emplace("Half cystine", PtmEntry{"Half cystine", "", "MOD:00798", "374"});
    return map;
  }

  /// UniPEFF_5_input_isoform_features.xml: P55555 (30 aa) with a displayed isoform -1, a described isoform -2
  /// (VSP_1: A10 -> ASCHSK, 35 aa), an external isoform -3 and an isoform -4 whose splice variant does not exist
  UniProtEntry isoformFeatures()
  {
    UniProtEntry e;
    e.dataset = "Swiss-Prot";
    e.accession = "P55555";
    e.sequence = "MSEQCKTLDAPGRSTWEVQLKDNACYLRGE";
    e.isoforms = {isoform("P55555-1", "displayed"), isoform("P55555-2", "described", "VSP_1"),
                  isoform("P55555-3", "external"), isoform("P55555-4", "described", "VSP_MISSING")};
    e.features = {position("modified residue", "Phosphoserine", 14),
                  position("modified residue", "Phosphothreonine", 15, "P55555-1"),
                  position("modified residue", "Phosphoserine", 11, "P55555-2"),
                  range("disulfide bond", "", 12, 30, "P55555-2"),
                  position("sequence variant", "In dbSNP.", 33, "P55555-2", "Q"),
                  position("sequence conflict", "In Ref. 1; AAA00001.", 13, "P55555-2", "Y"),
                  position("modified residue", "Phosphothreonine", 7, "P55555-3"),
                  position("modified residue", "Phosphoserine", 2, "P55555-4")};
    UniProtFeature splice = position("splice variant", "In isoform 2.", 10, "", "ASCHSK");
    splice.id = "VSP_1";
    splice.original = "A";
    e.features.push_back(splice);
    return e;
  }

  PreparedEntry prepare(UniProtEntry e, bool record_isoforms, SkippedIsoformCounts& skipped)
  {
    return prepareEntry(std::move(e), ptms(), "", false, true, true, true, record_isoforms, skipped);
  }
}

START_TEST(UniPEFFEntries, "$Id$")

START_SECTION(PreparedEntry prepareEntry(UniProtEntry&&, const PtmMap&, const std::string&, bool, bool, bool, bool, bool, SkippedIsoformCounts&))
{
  // features annotated on a written isoform go to its entry, in its coordinates; those on the displayed isoform stay
  // canonical; those on the external and the unreconstructable isoform are not applied, and a sequence conflict
  // (never written) is not counted
  SkippedIsoformCounts skipped;
  PreparedEntry pe = prepare(isoformFeatures(), true, skipped);
  TEST_EQUAL(pe.prefix, "sp")
  TEST_EQUAL(pe.isoform_scoped_applied, 3)
  TEST_EQUAL(pe.isoform_scoped_skipped, 2)
  TEST_EQUAL(skipped.size(), 1)
  TEST_EQUAL(skipped["external"], 1)
  ABORT_IF(pe.isoforms.size() != 1)
  const IsoformEntry& iso = pe.isoforms[0];
  TEST_EQUAL(iso.accession, "P55555-2")
  TEST_EQUAL(iso.sequence, "MSEQCKTLDASCHSKPGRSTWEVQLKDNACYLRGE")
  // canonical: phosphoserine 14 and the phosphothreonine 15 of the displayed isoform
  ABORT_IF(pe.annotations.mods.size() != 2)
  TEST_EQUAL(pe.annotations.mods[0].position, 14)
  TEST_EQUAL(pe.annotations.mods[1].position, 15)
  TEST_EQUAL(pe.annotations.simple_variants.size(), 0)
  // isoform: phosphoserine 11, the half-cystines 12 and 30 of a labeled disulfide bond, variant R33Q
  ABORT_IF(iso.annotations.mods.size() != 3)
  TEST_EQUAL(iso.annotations.mods[0].position, 11)
  TEST_EQUAL(iso.annotations.mods[1].position, 12)
  TEST_EQUAL(iso.annotations.mods[2].position, 30)
  ABORT_IF(iso.annotations.disulfides.size() != 1)
  TEST_TRUE(iso.annotations.disulfides[0].valid)
  TEST_EQUAL(iso.annotations.disulfides[0].annotation_id, 2)
  ABORT_IF(iso.annotations.simple_variants.size() != 1)
  TEST_EQUAL(iso.annotations.simple_variants[0].position, 33)
  TEST_EQUAL(iso.annotations.simple_variants[0].new_aa, 'Q')

  // without isoform entries (-omit_isoforms), none of the five is applied
  skipped.clear();
  pe = prepare(isoformFeatures(), false, skipped);
  TEST_EQUAL(pe.isoform_scoped_applied, 0)
  TEST_EQUAL(pe.isoform_scoped_skipped, 5)
  TEST_EQUAL(pe.isoforms.size(), 0)
  TEST_EQUAL(skipped.size(), 0)
  TEST_EQUAL(pe.annotations.mods.size(), 2)
}
END_SECTION

START_SECTION(([EXTRA] IsoformScopedTracker counts features and entries))
{
  // UniPEFF_4_input_location_sequence.xml: entries without isoform definitions, so features annotated on their
  // isoforms are not written: four on P11111 next to its canonical features (and a sequence conflict, which is never
  // written and not counted), one on P22222
  UniProtEntry p11111;
  p11111.dataset = "Swiss-Prot";
  p11111.accession = "P11111";
  p11111.sequence = "MSEQNCEMTFQIQRIYNCKDISFEAPNAPCKLTPEYASGW";
  p11111.features = {position("initiator methionine", "Removed", 1),
                     range("chain", "Isoform-scoped feature test protein", 2, 40),
                     position("modified residue", "Phosphoserine", 22),
                     position("modified residue", "Phosphoserine", 20, "P11111-2"),
                     position("modified residue", "Phosphothreonine", 52, "P11111-3"),
                     range("disulfide bond", "Interchain", 6, 18),
                     position("sequence variant", "In dbSNP.", 25, "", "V"),
                     position("sequence variant", "In dbSNP.", 28, "P11111-2", "L"),
                     range("sequence variant", "In dbSNP.", 44, 45, "P11111-3", "E"),
                     position("sequence conflict", "In Ref. 2; AAA00000.", 33, "P11111-2", "K")};
  UniProtEntry p22222;
  p22222.dataset = "Swiss-Prot";
  p22222.accession = "P22222";
  p22222.sequence = "MKLDLDTASS";
  p22222.features = {position("modified residue", "Phosphoserine", 3, "P22222-2")};

  IsoformScopedTracker tracker;
  SkippedIsoformCounts skipped;
  PreparedEntry pe = prepare(p11111, true, skipped);
  TEST_EQUAL(pe.isoform_scoped_skipped, 4)
  TEST_EQUAL(pe.annotations.processed.size(), 2)
  TEST_EQUAL(pe.annotations.mods.size(), 3) // phosphoserine 22 and the half-cystines 6 and 18
  TEST_EQUAL(pe.annotations.simple_variants.size(), 1)
  TEST_EQUAL(pe.annotations.complex_variants.size(), 0)
  tracker.add(pe.isoform_scoped_applied, pe.isoform_scoped_skipped);
  pe = prepare(p22222, true, skipped);
  TEST_EQUAL(pe.isoform_scoped_skipped, 1)
  TEST_EQUAL(pe.annotations.mods.size(), 0)
  tracker.add(pe.isoform_scoped_applied, pe.isoform_scoped_skipped);
  // the run summary: 5 features not applied, in 2 entries
  TEST_EQUAL(tracker.skipped, 5)
  TEST_EQUAL(tracker.skipped_entries, 2)
  TEST_EQUAL(tracker.applied, 0)
  TEST_EQUAL(tracker.applied_entries, 0)

  // the isoform features of UniPEFF_5: 3 written in 1 entry, 2 not applied in 1 entry
  IsoformScopedTracker isoforms;
  pe = prepare(isoformFeatures(), true, skipped);
  isoforms.add(pe.isoform_scoped_applied, pe.isoform_scoped_skipped);
  TEST_EQUAL(isoforms.applied, 3)
  TEST_EQUAL(isoforms.applied_entries, 1)
  TEST_EQUAL(isoforms.skipped, 2)
  TEST_EQUAL(isoforms.skipped_entries, 1)
}
END_SECTION

START_SECTION((bool isResidueCode(char c)))
{
  TEST_TRUE(isResidueCode('A'))
  TEST_TRUE(isResidueCode('q'))
  TEST_TRUE(isResidueCode('*'))
  TEST_FALSE(isResidueCode('1'))
  TEST_FALSE(isResidueCode('-'))
}
END_SECTION

START_SECTION((std::string positionOrUnknown(int p)))
{
  TEST_EQUAL(positionOrUnknown(0), "?")
  TEST_EQUAL(positionOrUnknown(42), "42")
}
END_SECTION

END_TEST
