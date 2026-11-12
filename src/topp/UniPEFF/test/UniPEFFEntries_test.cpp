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

#include <algorithm>

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
    map.emplace("N-acetylmethionine", PtmEntry{"N-acetylmethionine", "PTM-0207", "MOD:00069", "1"});
    map.emplace("N6-acetyllysine", PtmEntry{"N6-acetyllysine", "PTM-0190", "MOD:00064", "1"});
    // made at a terminus of the chain (ptmlist "PP   N-terminal." / "PP   C-terminal.")
    map.emplace("N-acetylalanine", PtmEntry{"N-acetylalanine", "PTM-0199", "MOD:00050", "1", true, false});
    map.emplace("N,N-dimethylproline", PtmEntry{"N,N-dimethylproline", "PTM-0179", "MOD:00075", "529", true, false});
    map.emplace("Alanine amide", PtmEntry{"Alanine amide", "PTM-0057", "MOD:00090", "2", false, true});
    // both termini: a head-to-tail cyclopeptide ("PP   N-terminal-C-terminal.")
    map.emplace("Test cyclopeptide", PtmEntry{"Test cyclopeptide", "PTM-9004", "", "", true, true});
    return map;
  }

  UniProtFeature splice(const std::string& id, int begin, int end, const std::string& original = "",
                        const std::string& variation = "")
  {
    UniProtFeature f = range("splice variant", "", begin, end, "", variation);
    f.id = id;
    f.original = original;
    return f;
  }

  /// P66666 (30 aa) with isoform -2: VSP_A deletes 3-6 (EQCK), VSP_B replaces L20 by LIV (28 aa).
  /// Canonical residues 1-2 keep their positions, 7-19 move by -4, 21-30 by -2.
  UniProtEntry canonicalFeatures()
  {
    UniProtEntry e;
    e.dataset = "Swiss-Prot";
    e.accession = "P66666";
    e.sequence = "MSEQCKTLDAPGRSTWEVQLKDNACYLRGE";
    e.isoforms = {isoform("P66666-1", "displayed"), isoform("P66666-2", "described", "VSP_A VSP_B")};
    e.features = {position("initiator methionine", "Removed", 1),
                  range("chain", "Test protein", 2, 30),
                  range("propeptide", "", 2, 8),                                  // spans VSP_A -> 2-4
                  range("signal peptide", "", 3, 6),                              // deleted
                  position("modified residue", "Phosphoserine", 2),
                  position("modified residue", "Phosphoserine", 14),
                  position("modified residue", "Phosphothreonine", 15),
                  position("modified residue", "Phosphoserine", 20),              // replaced
                  position("modified residue", "Phosphoserine", 0),               // unknown position
                  position("modified residue", "Not in ptmlist", 16),             // no annotation at all
                  range("disulfide bond", "", 5, 25),                             // half-cystine 5 deleted
                  position("sequence variant", "In dbSNP.", 24, "", "V"),
                  range("sequence variant", "In dbSNP.", 18, 21, "", "Q"),       // overlaps VSP_B
                  position("sequence variant", "In dbSNP.", 13, "", "K"),        // also on the isoform
                  position("sequence variant", "In dbSNP.", 9, "P66666-2", "K"),
                  splice("VSP_A", 3, 6),
                  splice("VSP_B", 20, 20, "L", "LIV")};
    return e;
  }

  /// A0AVT1 (UBA6_HUMAN, entry version 163, sequence version 1): canonical sequence, chain, modified residues, natural
  /// variant and the splice variants of isoforms -2 (VSP_023083), -3 (VSP_023086/7), -4 (VSP_023084/5).
  UniProtEntry a0avt1()
  {
    UniProtEntry e;
    e.dataset = "Swiss-Prot";
    e.accession = "A0AVT1";
    e.sequence = "MEGSEPVAAHQGEEASCSSWGTGSTNKNLPIMSTASVEIDDALYSRQRYVLGDTAMQKMAKSHVFLSGMG"
                 "GLGLEIAKNLVLAGIKAVTIHDTEKCQAWDLGTNFFLSEDDVVNKRNRAEAVLKHIAELNPYVHVTSSSV"
                 "PFNETTDLSFLDKYQCVVLTEMKLPLQKKINDFCRSQCPPIKFISADVHGIWSRLFCDFGDEFEVLDTTG"
                 "EEPKEIFISNITQANPGIVTCLENHPHKLETGQFLTFREINGMTGLNGSIQQITVISPFSFSIGDTTELE"
                 "PYLHGGIAVQVKTPKTVFFESLERQLKHPKCLIVDFSNPEAPLEIHTAMLALDQFQEKYSRKPNVGCQQD"
                 "SEELLKLATSISETLEEKPDVNADIVHWLSWTAQGFLSPLAAAVGGVASQEVLKAVTGKFSPLCQWLYLE"
                 "AADIVESLGKPECEEFLPRGDRYDALRACIGDTLCQKLQNLNIFLVGCGAIGCEMLKNFALLGVGTSKEK"
                 "GMITVTDPDLIEKSNLNRQFLFRPHHIQKPKSYTAADATLKINSQIKIDAHLNKVCPTTETIYNDEFYTK"
                 "QDVIITALDNVEARRYVDSRCLANLRPLLDSGTMGTKGHTEVIVPHLTESYNSHRDPPEEEIPFCTLKSF"
                 "PAAIEHTIQWARDKFESSFSHKPSLFNKFWQTYSSAEEVLQKIQSGHSLEGCFQVIKLLSRRPRNWSQCV"
                 "ELARLKFEKYFNHKALQLLHCFPLDIRLKDGSLFWQSPKRPPSPIKFDLNEPLHLSFLQNAAKLYATVYC"
                 "IPFAEEDLSADALLNILSEVKIQEFKPSNKVVQTDETARKPDHVPISSEDERNAIFQLEKAILSNEATKS"
                 "DLQMAVLSFEKDDDHNGHIDFITAASNLRAKMYSIEPADRFKTKRIAGKIIPAIATTTATVSGLVALEMI"
                 "KVTGGYPFEAYKNCFLNLAIPIVVFTETTEVRKTKIRNGISFTIWDRWTVHGKEDFTLLDFINAVKEKYG"
                 "IEPTMVVQGVKMLYVPVMPGHAKRLKLTMHKLVKPTTEKKYVDLTVSFAPDIDGDEDLPGPPVRYYFSHD"
                 "TD";
    e.isoforms = {isoform("A0AVT1-1", "displayed"), isoform("A0AVT1-2", "described", "VSP_023083"),
                  isoform("A0AVT1-3", "described", "VSP_023086 VSP_023087"),
                  isoform("A0AVT1-4", "described", "VSP_023084 VSP_023085")};
    e.features = {range("chain", "Ubiquitin-like modifier-activating enzyme 6", 1, 1052),
                  position("modified residue", "N-acetylmethionine", 1),
                  position("modified residue", "Phosphothreonine", 54),
                  position("modified residue", "Phosphoserine", 301),
                  position("modified residue", "N6-acetyllysine", 544),
                  position("modified residue", "N6-acetyllysine", 729),
                  position("modified residue", "Phosphoserine", 737),
                  splice("VSP_023083", 1, 474),
                  splice("VSP_023084", 321, 340, "APLEIHTAMLALDQFQEKYS", "VNKHFAGLREAAESEMRISE"),
                  splice("VSP_023085", 341, 1052),
                  splice("VSP_023086", 369, 389, "PDVNADIVHWLSWTAQGFLSP", "VTIEIYGCPNICLLIHKCSVY"),
                  splice("VSP_023087", 390, 1052),
                  position("sequence variant", "In dbSNP:rs10010188.", 224, "", "T")};
    return e;
  }

  /// P77777 (30 aa) with isoform -2: VSP_A deletes 3-6 (26 aa), so canonical residues 7-30 move by -4.
  /// Each canonical annotation meets an annotation UniProt makes on the isoform at the same residues:
  /// compatible ones are written side by side, contradicting ones leave the isoform's own.
  UniProtEntry conflictingFeatures()
  {
    const std::string sumo1 = "Glycyl lysine isopeptide (Lys-Gly) (interchain with G-Cter in SUMO1); alternate";
    const std::string sumo2 = "Glycyl lysine isopeptide (Lys-Gly) (interchain with G-Cter in SUMO2)";
    UniProtEntry e;
    e.dataset = "Swiss-Prot";
    e.accession = "P77777";
    e.sequence = "MSEQCKTLDAPGRSTWEVQLKDNACYLRGE";
    e.isoforms = {isoform("P77777-1", "displayed"), isoform("P77777-2", "described", "VSP_A")};
    e.features = {position("initiator methionine", "Removed", 1),                            // -> 1
                  range("chain", "Test protein", 2, 30),                                      // -> 2-26, conflicts
                  range("chain", "Test protein", 2, 20, "P77777-2"),
                  range("disulfide bond", "", 11, 25),                                        // -> 7-21, conflicts
                  range("disulfide bond", "", 7, 14, "P77777-2"),                             // shares C7
                  range("cross-link", "Isoglutamyl lysine isopeptide (Lys-Gln)", 13, 28),     // -> 9-24, conflicts
                  range("cross-link", "Isoglutamyl lysine isopeptide (Lys-Gln)", 22, 24, "P77777-2"),
                  position("cross-link", sumo1, 21),                                          // -> 17, compatible
                  position("cross-link", sumo2, 17, "P77777-2"),
                  position("cross-link", sumo1, 26),                         // -> 22, in the isoform's 22-24: conflicts
                  range("disulfide bond", "", 28, 30),                       // -> 24-26, shares 24 with the cross-link
                  position("modified residue", "Phosphoserine", 14),                          // -> 10, compatible
                  position("modified residue", "Phosphothreonine", 10, "P77777-2"),
                  position("glycosylation site", "N-linked (GlcNAc...) (complex) asparagine", 27),  // -> 23, redundant
                  position("glycosylation site", "N-linked (GlcNAc...) asparagine", 23, "P77777-2"),
                  position("sequence variant", "In dbSNP.", 24, "", "V"),                     // -> 20, compatible
                  position("sequence variant", "In dbSNP.", 20, "P77777-2", "I"),
                  splice("VSP_A", 3, 6)};
    return e;
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

  /// P88888 (30 aa, MAEQCKTLDAPGRSTWEVQLKDNACYLRGE) as in many Swiss-Prot entries: the initiator methionine is removed
  /// and the mature chain 2-30 is N-acetylated at A2. Isoform -2 has an N-terminal extension (M1 -> MPRLLSG, 36 aa),
  /// -3 an alternative first exon (1-5 MAEQC -> MKL, 28 aa), -4 a deletion inside the chain (25-27 CYL, 27 aa).
  UniProtEntry nTerminalFeatures()
  {
    UniProtEntry e;
    e.dataset = "Swiss-Prot";
    e.accession = "P88888";
    e.sequence = "MAEQCKTLDAPGRSTWEVQLKDNACYLRGE";
    e.isoforms = {isoform("P88888-1", "displayed"), isoform("P88888-2", "described", "VSP_EXT"),
                  isoform("P88888-3", "described", "VSP_EXON"), isoform("P88888-4", "described", "VSP_DEL")};
    e.features = {position("initiator methionine", "Removed", 1),
                  range("chain", "Test protein", 2, 30),
                  position("modified residue", "N-acetylalanine", 2),
                  position("modified residue", "Phosphoserine", 14),
                  splice("VSP_EXT", 1, 1, "M", "MPRLLSG"),
                  splice("VSP_EXON", 1, 5, "MAEQC", "MKL"),
                  splice("VSP_DEL", 25, 27, "CYL")};
    return e;
  }

  /// P99999 (30 aa, MAEQCKTLDAPGRSTWEVQLKDNACYLRGE) with signal peptide 1-10, chain 11-24 (N,N-dimethylated at P11,
  /// amidated at A24) and propeptide 25-30. Isoform -2 replaces A10 (A -> AV, 31 aa), so the cleavage site 10|11
  /// differs; -3 replaces 25-26 (CY -> W, 29 aa), so 24|25 differs; -4 deletes 18-19 (VQ, 28 aa), inside the chain.
  UniProtEntry cleavageSiteFeatures()
  {
    UniProtEntry e;
    e.dataset = "Swiss-Prot";
    e.accession = "P99999";
    e.sequence = "MAEQCKTLDAPGRSTWEVQLKDNACYLRGE";
    e.isoforms = {isoform("P99999-1", "displayed"), isoform("P99999-2", "described", "VSP_N"),
                  isoform("P99999-3", "described", "VSP_C"), isoform("P99999-4", "described", "VSP_IN")};
    e.features = {range("signal peptide", "", 1, 10),
                  range("chain", "Test protein", 11, 24),
                  range("propeptide", "", 25, 30),
                  position("modified residue", "N,N-dimethylproline", 11),
                  position("modified residue", "Alanine amide", 24),
                  position("modified residue", "Phosphoserine", 14),
                  splice("VSP_N", 10, 10, "A", "AV"),
                  splice("VSP_C", 25, 26, "CY", "W"),
                  splice("VSP_IN", 18, 19, "VQ")};
    return e;
  }

  PreparedEntry prepare(UniProtEntry e, bool record_isoforms, SkippedIsoformCounts& skipped,
                        bool map_canonical_features = false)
  {
    return prepareEntry(std::move(e), ptms(), "", false, true, true, true, record_isoforms,
                        map_canonical_features, skipped);
  }

  std::vector<int> modPositions(const EntryAnnotations& a)
  {
    std::vector<int> out;
    for (const ModResItem& m : a.mods) out.push_back(m.position);
    return out;
  }
}

START_TEST(UniPEFFEntries, "$Id$")

START_SECTION(PreparedEntry prepareEntry(UniProtEntry&&, const PtmMap&, const std::string&, bool, bool, bool, bool, bool, bool, SkippedIsoformCounts&))
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

START_SECTION(([EXTRA] map_canonical_features: canonical annotations in isoform coordinates))
{
  SkippedIsoformCounts skipped;

  // off (-omit_canonical_features): the isoform carries only the variant UniProt annotates on it
  PreparedEntry pe = prepare(canonicalFeatures(), true, skipped);
  ABORT_IF(pe.isoforms.size() != 1)
  TEST_EQUAL(pe.isoforms[0].sequence, "MSTLDAPGRSTWEVQLIVKDNACYLRGE")
  TEST_EQUAL(pe.isoforms[0].annotations.mods.size(), 0)
  TEST_EQUAL(pe.isoforms[0].annotations.processed.size(), 0)
  TEST_EQUAL(pe.isoforms[0].annotations.simple_variants.size(), 1)
  TEST_EQUAL(pe.canonical_mapped, 0)
  TEST_EQUAL(pe.canonical_unmapped, 0)

  // on: of 13 canonical annotations, 7 are mapped, 5 lie in or cross a splice variant (signal peptide 3-6,
  // phosphoserines 20 and '?', the disulfide bond 5-25, the variant 18-21) and the copy of variant R13K is
  // already annotated on the isoform (as R9K); the feature missing from the ptmlist is no annotation at all
  pe = prepare(canonicalFeatures(), true, skipped, true);
  ABORT_IF(pe.isoforms.size() != 1)
  const EntryAnnotations& iso = pe.isoforms[0].annotations;
  TEST_EQUAL(pe.canonical_mapped, 7)
  TEST_EQUAL(pe.canonical_unmapped, 5)
  TEST_EQUAL(pe.isoform_scoped_applied, 1)
  // the canonical entry is unchanged
  TEST_EQUAL(pe.annotations.processed.size(), 4)
  TEST_EQUAL(pe.annotations.mods.size(), 7) // 5 phosphorylations and 2 half-cystines
  // phosphoserine 2, phosphoserine 14 -> 10, phosphothreonine 15 -> 11; no half-cystines
  TEST_TRUE(modPositions(iso) == std::vector<int>({2, 10, 11}))
  TEST_EQUAL(iso.disulfides.size(), 0)
  // initiator methionine 1, chain 2-30 -> 2-28 (still to the C-terminus), propeptide 2-8 (spanning the deletion) -> 2-4
  ABORT_IF(iso.processed.size() != 3)
  TEST_EQUAL(iso.processed[0].begin, 1)
  TEST_EQUAL(iso.processed[0].end, 1)
  TEST_EQUAL(iso.processed[1].begin, 2)
  TEST_EQUAL(iso.processed[1].end, 28)
  TEST_EQUAL(iso.processed[2].begin, 2)
  TEST_EQUAL(iso.processed[2].end, 4)
  // R9K annotated on the isoform (written once), A24V -> 22
  ABORT_IF(iso.simple_variants.size() != 2)
  TEST_EQUAL(iso.simple_variants[0].position, 9)
  TEST_EQUAL(iso.simple_variants[1].position, 22)
  TEST_EQUAL(iso.simple_variants[1].new_aa, 'V')
  TEST_EQUAL(iso.complex_variants.size(), 0)

  // a copy of an isoform annotation is recognized by what is written: a variant's description is not written,
  // a modification's is (after cleaning, so "; by PKA" does not matter), a position equals a one-residue range
  // at it (A24V -> 22 and the isoform's 22-22 V); a different modification is kept
  UniProtEntry described = canonicalFeatures();
  described.features = {position("sequence variant", "In dbSNP:rs1.", 13, "", "K"),
                        position("sequence variant", "In dbSNP:rs2.", 9, "P66666-2", "K"),
                        position("sequence variant", "In dbSNP.", 24, "", "V"),
                        range("sequence variant", "In dbSNP.", 22, 22, "P66666-2", "V"),
                        position("modified residue", "Phosphoserine", 14),
                        position("modified residue", "Phosphoserine; by PKA", 10, "P66666-2"),
                        position("modified residue", "Phosphothreonine", 15),
                        position("modified residue", "Phosphoserine", 11, "P66666-2"),
                        splice("VSP_A", 3, 6), splice("VSP_B", 20, 20, "L", "LIV")};
  pe = prepare(described, true, skipped, true);
  ABORT_IF(pe.isoforms.size() != 1)
  TEST_EQUAL(pe.isoform_scoped_applied, 4)
  TEST_EQUAL(pe.canonical_mapped, 1)
  TEST_EQUAL(pe.canonical_unmapped, 0)
  TEST_EQUAL(pe.isoforms[0].annotations.simple_variants.size(), 1)
  TEST_EQUAL(pe.isoforms[0].annotations.complex_variants.size(), 1)
  TEST_TRUE(modPositions(pe.isoforms[0].annotations) == std::vector<int>({10, 11, 11}))

  // the disulfide bond survives when both half-cystines do, whatever lies between them
  UniProtEntry bonded = canonicalFeatures();
  bonded.features = {range("disulfide bond", "", 2, 25), splice("VSP_A", 3, 6), splice("VSP_B", 20, 20, "L", "LIV")};
  pe = prepare(bonded, true, skipped, true);
  ABORT_IF(pe.isoforms.size() != 1)
  TEST_TRUE(modPositions(pe.isoforms[0].annotations) == std::vector<int>({2, 23}))
  ABORT_IF(pe.isoforms[0].annotations.disulfides.size() != 1)
  TEST_TRUE(pe.isoforms[0].annotations.disulfides[0].valid)

  // nothing to map onto without isoform entries
  pe = prepare(canonicalFeatures(), false, skipped, true);
  TEST_EQUAL(pe.isoforms.size(), 0)
  TEST_EQUAL(pe.canonical_mapped, 0)
  TEST_EQUAL(pe.canonical_unmapped, 0)
}
END_SECTION

START_SECTION(([EXTRA] map_canonical_features: annotations UniProt makes on the isoform at the same residues))
{
  SkippedIsoformCounts skipped;
  PreparedEntry pe = prepare(conflictingFeatures(), true, skipped, true);
  ABORT_IF(pe.isoforms.size() != 1)
  const IsoformEntry& iso = pe.isoforms[0];
  TEST_EQUAL(iso.sequence, "MSTLDAPGRSTWEVQLKDNACYLRGE")
  TEST_EQUAL(pe.isoform_scoped_applied, 7)
  // of 10 canonical annotations: the initiator methionine, the SUMO1 cross-link at 17, the phosphoserine and the
  // variant are mapped; the chain, the disulfide bond 7-21, the cross-link between two residues, the SUMO1 cross-link
  // at 22 and the disulfide bond 24-26 (both on residues of the isoform's cross-link 22-24) contradict the isoform's
  // own and are not; the glycosylation is already on the isoform (described in less detail) and not written again
  TEST_EQUAL(pe.canonical_mapped, 4)
  TEST_EQUAL(pe.canonical_unmapped, 5)

  // only the isoform's chain 2-20, plus the mapped initiator methionine
  ABORT_IF(iso.annotations.processed.size() != 2)
  TEST_EQUAL(iso.annotations.processed[0].begin, 2)
  TEST_EQUAL(iso.annotations.processed[0].end, 20)
  TEST_EQUAL(iso.annotations.processed[1].begin, 1)
  TEST_EQUAL(iso.annotations.processed[1].end, 1)
  // half-cystines 7 and 14 of the isoform's bond only (none at 21, 24 or 26); phosphothreonine and phosphoserine at 10;
  // SUMO2 and SUMO1 at 17; one glycosylation at 23; the isoform's cross-link 22-24 only (none at 9, no SUMO1 at 22)
  TEST_TRUE(modPositions(iso.annotations) == std::vector<int>({7, 10, 10, 14, 17, 17, 22, 23, 24}))
  ABORT_IF(iso.annotations.disulfides.size() != 1)
  TEST_TRUE(iso.annotations.disulfides[0].valid)
  // R20I annotated on the isoform and the mapped A24V, as two variants of residue 20
  ABORT_IF(iso.annotations.simple_variants.size() != 2)
  TEST_EQUAL(iso.annotations.simple_variants[0].position, 20)
  TEST_EQUAL(iso.annotations.simple_variants[0].new_aa, 'I')
  TEST_EQUAL(iso.annotations.simple_variants[1].position, 20)
  TEST_EQUAL(iso.annotations.simple_variants[1].new_aa, 'V')

  // without the isoform's own annotations, all 10 are mapped
  UniProtEntry plain = conflictingFeatures();
  plain.features.erase(std::remove_if(plain.features.begin(), plain.features.end(),
                                      [](const UniProtFeature& f) { return !f.location_sequence.empty(); }),
                       plain.features.end());
  pe = prepare(plain, true, skipped, true);
  TEST_EQUAL(pe.canonical_mapped, 10)
  TEST_EQUAL(pe.canonical_unmapped, 0)
}
END_SECTION

START_SECTION(([EXTRA] map_canonical_features: cleavage sites and chain-terminal modifications))
{
  SkippedIsoformCounts skipped;

  // initiator methionine removed, N-acetylalanine 2, chain 2-30, phosphoserine 14
  PreparedEntry pe = prepare(nTerminalFeatures(), true, skipped, true);
  ABORT_IF(pe.isoforms.size() != 3)
  TEST_EQUAL(pe.isoforms[0].sequence, "MPRLLSGAEQCKTLDAPGRSTWEVQLKDNACYLRGE")
  TEST_EQUAL(pe.isoforms[1].sequence, "MKLKTLDAPGRSTWEVQLKDNACYLRGE")
  TEST_EQUAL(pe.isoforms[2].sequence, "MAEQCKTLDAPGRSTWEVQLKDNARGE")
  // -2 and -3: the initiator methionine and the N-acetylation are not mapped, the rest is; -4: all four
  TEST_EQUAL(pe.canonical_mapped, 8)
  TEST_EQUAL(pe.canonical_unmapped, 4)
  // N-terminal extension: A2 -> 8 is inside the isoform's chain, so it is not N-acetylated there, and the chain
  // starts at the isoform's N-terminus (not at 8, which would declare 1-7 cleaved)
  const EntryAnnotations& ext = pe.isoforms[0].annotations;
  TEST_TRUE(modPositions(ext) == std::vector<int>({20}))
  ABORT_IF(ext.processed.size() != 1)
  TEST_EQUAL(ext.processed[0].type_name, "mature protein")
  TEST_EQUAL(ext.processed[0].begin, 1)
  TEST_EQUAL(ext.processed[0].end, 36)
  // alternative first exon: the chain starts at the isoform's N-terminus (not at 4, which would declare MKL cleaved)
  const EntryAnnotations& exon = pe.isoforms[1].annotations;
  TEST_TRUE(modPositions(exon) == std::vector<int>({12}))
  ABORT_IF(exon.processed.size() != 1)
  TEST_EQUAL(exon.processed[0].begin, 1)
  TEST_EQUAL(exon.processed[0].end, 28)
  // a splice variant away from the N-terminus changes neither
  const EntryAnnotations& del = pe.isoforms[2].annotations;
  TEST_TRUE(modPositions(del) == std::vector<int>({2, 14}))
  ABORT_IF(del.processed.size() != 2)
  TEST_EQUAL(del.processed[0].type_name, "initiator methionine")
  TEST_EQUAL(del.processed[0].begin, 1)
  TEST_EQUAL(del.processed[0].end, 1)
  TEST_EQUAL(del.processed[1].begin, 2)
  TEST_EQUAL(del.processed[1].end, 27)

  // signal peptide 1-10, chain 11-24 (N,N-dimethylproline 11, alanine amide 24), propeptide 25-30, phosphoserine 14
  pe = prepare(cleavageSiteFeatures(), true, skipped, true);
  ABORT_IF(pe.isoforms.size() != 3)
  TEST_EQUAL(pe.isoforms[0].sequence, "MAEQCKTLDAVPGRSTWEVQLKDNACYLRGE")
  TEST_EQUAL(pe.isoforms[1].sequence, "MAEQCKTLDAPGRSTWEVQLKDNAWLRGE")
  TEST_EQUAL(pe.isoforms[2].sequence, "MAEQCKTLDAPGRSTWELKDNACYLRGE")
  TEST_EQUAL(pe.canonical_mapped, 12)
  TEST_EQUAL(pe.canonical_unmapped, 6)
  // A10 -> AV changes the cleavage site 10|11: neither the signal peptide nor the chain is mapped (they would end and
  // start at a site the isoform does not have), nor the N-terminal N,N-dimethylproline; the propeptide 25-30 -> 26-31
  // and the amidated A24 -> 25 at its cleavage site are
  const EntryAnnotations& n = pe.isoforms[0].annotations;
  TEST_TRUE(modPositions(n) == std::vector<int>({15, 25}))
  ABORT_IF(n.processed.size() != 1)
  TEST_EQUAL(n.processed[0].type_name, "propeptide")
  TEST_EQUAL(n.processed[0].begin, 26)
  TEST_EQUAL(n.processed[0].end, 31)
  // CY25-26 -> W changes the cleavage site 24|25: the chain, the propeptide and the C-terminal amidation are not mapped
  const EntryAnnotations& c = pe.isoforms[1].annotations;
  TEST_TRUE(modPositions(c) == std::vector<int>({11, 14}))
  ABORT_IF(c.processed.size() != 1)
  TEST_EQUAL(c.processed[0].type_name, "signal peptide")
  TEST_EQUAL(c.processed[0].begin, 1)
  TEST_EQUAL(c.processed[0].end, 10)
  // a deletion inside the chain keeps both cleavage sites: everything is mapped
  const EntryAnnotations& in = pe.isoforms[2].annotations;
  TEST_TRUE(modPositions(in) == std::vector<int>({11, 14, 22}))
  ABORT_IF(in.processed.size() != 3)
  TEST_EQUAL(in.processed[0].begin, 1)
  TEST_EQUAL(in.processed[0].end, 10)
  TEST_EQUAL(in.processed[1].begin, 11)
  TEST_EQUAL(in.processed[1].end, 22)
  TEST_EQUAL(in.processed[2].begin, 23)
  TEST_EQUAL(in.processed[2].end, 28)

  // a head-to-tail cyclopeptide links the first and last residue of the chain: like the chain, it is mapped only where
  // both cleavage sites are kept, although its two residues lie outside every splice variant in each isoform. -4 here
  // replaces G12 (G -> GSA, 32 aa) next to the first residue but inside the chain: both are mapped
  UniProtEntry cyclic = cleavageSiteFeatures();
  cyclic.features = {range("chain", "Test protein", 11, 24), range("cross-link", "Test cyclopeptide", 11, 24),
                     splice("VSP_N", 10, 10, "A", "AV"), splice("VSP_C", 25, 26, "CY", "W"), splice("VSP_IN", 12, 12, "G", "GSA")};
  pe = prepare(cyclic, true, skipped, true);
  ABORT_IF(pe.isoforms.size() != 3)
  TEST_EQUAL(pe.isoforms[2].sequence, "MAEQCKTLDAPGSARSTWEVQLKDNACYLRGE")
  TEST_EQUAL(pe.canonical_mapped, 2)
  TEST_EQUAL(pe.canonical_unmapped, 4)
  TEST_EQUAL(pe.isoforms[0].annotations.mods.size(), 0)
  TEST_EQUAL(pe.isoforms[1].annotations.mods.size(), 0)
  TEST_TRUE(modPositions(pe.isoforms[2].annotations) == std::vector<int>({11, 26}))
  ABORT_IF(pe.isoforms[2].annotations.processed.size() != 1)
  TEST_EQUAL(pe.isoforms[2].annotations.processed[0].begin, 11)
  TEST_EQUAL(pe.isoforms[2].annotations.processed[0].end, 26)
}
END_SECTION

START_SECTION(([EXTRA] map_canonical_features: A0AVT1 isoforms))
{
  SkippedIsoformCounts skipped;
  PreparedEntry pe = prepare(a0avt1(), true, skipped, true);
  ABORT_IF(pe.isoforms.size() != 3)
  const IsoformEntry& iso2 = pe.isoforms[0];
  const IsoformEntry& iso3 = pe.isoforms[1];
  const IsoformEntry& iso4 = pe.isoforms[2];
  TEST_EQUAL(iso2.accession, "A0AVT1-2")
  TEST_EQUAL(iso3.accession, "A0AVT1-3")
  TEST_EQUAL(iso4.accession, "A0AVT1-4")
  // the lengths of UniProt's isoform sequences; isoform 2 starts at the internal Met475
  TEST_EQUAL(iso2.sequence.size(), 578)
  TEST_EQUAL(iso3.sequence.size(), 389)
  TEST_EQUAL(iso4.sequence.size(), 340)
  TEST_EQUAL(iso2.sequence.substr(0, 12), "MLKNFALLGVGT")
  TEST_EQUAL(iso3.sequence.substr(377), "NICLLIHKCSVY")
  TEST_EQUAL(iso4.sequence.substr(328), "REAAESEMRISE")

  // isoform 2 (1-474 deleted): the N6-acetyllysines 544 and 729 and phosphoserine 737 move to 70, 255, 263;
  // Met1, Thr54, Ser301 and the variant A224T are not in it
  TEST_TRUE(modPositions(iso2.annotations) == std::vector<int>({70, 255, 263}))
  TEST_EQUAL(iso2.sequence[69], 'K')
  TEST_EQUAL(iso2.sequence[254], 'K')
  TEST_EQUAL(iso2.sequence[262], 'S')
  TEST_EQUAL(iso2.annotations.simple_variants.size(), 0)
  ABORT_IF(iso2.annotations.processed.size() != 1)
  TEST_EQUAL(iso2.annotations.processed[0].begin, 1)
  TEST_EQUAL(iso2.annotations.processed[0].end, 578)

  // isoforms 3 and 4 keep the N-terminal half: Met1, Thr54, Ser301 and A224T; the chain still spans the
  // whole isoform, including the replaced C-terminal residues
  for (const IsoformEntry* iso : {&iso3, &iso4})
  {
    // (ABORT_IF would only leave the loop)
    TEST_TRUE(modPositions(iso->annotations) == std::vector<int>({1, 54, 301}))
    TEST_EQUAL(iso->annotations.simple_variants.size(), 1)
    if (iso->annotations.simple_variants.size() == 1)
    {
      TEST_EQUAL(iso->annotations.simple_variants[0].position, 224)
    }
    TEST_EQUAL(iso->annotations.processed.size(), 1)
    if (iso->annotations.processed.size() == 1)
    {
      TEST_EQUAL(iso->annotations.processed[0].begin, 1)
      TEST_EQUAL(iso->annotations.processed[0].end, static_cast<int>(iso->sequence.size()))
    }
  }

  // 8 canonical annotations per isoform: 4 + 5 + 5 mapped, 4 + 3 + 3 not
  TEST_EQUAL(pe.canonical_mapped, 14)
  TEST_EQUAL(pe.canonical_unmapped, 10)
  IsoformScopedTracker tracker;
  tracker.add(pe.canonical_mapped, pe.canonical_unmapped);
  TEST_EQUAL(tracker.applied_entries, 1)
  TEST_EQUAL(tracker.skipped_entries, 1)
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
