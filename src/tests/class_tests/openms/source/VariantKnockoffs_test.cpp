// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/VariantKnockoffs.h>
///////////////////////////

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <algorithm>
#include <cmath>

using namespace OpenMS;
using namespace std;

namespace
{
  /// eight tryptic peptides GGGGDLGGGK: a D followed by L in each, at 0-based positions 4, 14, ..., 74
  std::string repeats(Size n = 8)
  {
    std::string s;
    for (Size i = 0; i < n; ++i) s += "GGGGDLGGGK";
    return s;
  }

  ProteaseDigestion trypsin()
  {
    ProteaseDigestion d;
    d.setEnzyme("Trypsin");
    return d;
  }

  /// a PSM whose main score is a q-value; its raw score (higher is better) is the meta value "raw"
  PeptideIdentification makeId(double raw, double q, const std::vector<std::string>& accessions)
  {
    PeptideHit hit;
    hit.setSequence(AASequence::fromString("PEPTIDEK"));
    hit.setScore(q);
    hit.setMetaValue("raw", raw);
    hit.setTargetDecoyType(PeptideHit::TargetDecoyType::TARGET);
    std::vector<PeptideEvidence> evidences;
    for (const auto& acc : accessions)
    {
      PeptideEvidence ev;
      ev.setProteinAccession(acc);
      evidences.push_back(ev);
    }
    hit.setPeptideEvidences(evidences);
    PeptideIdentification id;
    id.setHigherScoreBetter(false);
    id.setScoreType("q-value");
    id.insertHit(hit);
    return id;
  }

  double rawScore(const PeptideHit& hit)
  {
    return double(hit.getMetaValue("raw"));
  }

  /// canonical P1; variant D5N (deamidation-like) with two knockoffs; variant D35W (other) with none
  VariantKnockoffs::Index index()
  {
    std::vector<FASTAFile::FASTAEntry> db = {
      {"P1", "", repeats()},
      {"P1_VAR_D5N", "variant_of=P1 variant=D5N", ""},
      {"P1_KO1_D15N_VAR_D5N", "knockoff_of=P1_VAR_D5N", ""},
      {"P1_KO2_D25N_VAR_D5N", "knockoff_of=P1_VAR_D5N", ""},
      {"P1_VAR_D35W", "variant_of=P1 variant=D35W", ""}};
    return VariantKnockoffs::indexDatabase(db);
  }

  const PeptideHit& hitOf(const PeptideIdentificationList& ids, double raw)
  {
    for (const auto& id : ids)
    {
      if (std::fabs(rawScore(id.getHits()[0]) - raw) < 1e-9) return id.getHits()[0];
    }
    throw std::runtime_error("no such hit");
  }
}

START_TEST(VariantKnockoffs, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(static std::string substitutionClass(double delta_mass))
{
  TEST_EQUAL(VariantKnockoffs::substitutionClass(0.984016), "deamidation-like")
  TEST_EQUAL(VariantKnockoffs::substitutionClass(-0.984016), "deamidation-like")
  TEST_EQUAL(VariantKnockoffs::substitutionClass(14.01565), "methyl-like")
  TEST_EQUAL(VariantKnockoffs::substitutionClass(15.994915), "oxidation-like")
  TEST_EQUAL(VariantKnockoffs::substitutionClass(30.0), "other")
}
END_SECTION

START_SECTION(static bool parseEntry(const std::string& accession, const std::string& description, EntryInfo& info))
{
  VariantKnockoffs::EntryInfo info;
  TEST_EQUAL(VariantKnockoffs::parseEntry("sp|P02769|ALBU_BOVIN_VAR_N12D", "variant_of=x tier=sample", info), true)
  TEST_EQUAL(info.knockoff, false)
  TEST_EQUAL(info.family, "sp|P02769|ALBU_BOVIN")
  TEST_EQUAL(info.variant, "N12D")
  TEST_REAL_SIMILAR(info.delta_mass, 0.984016)
  TEST_EQUAL(info.tier, "sample")

  TEST_EQUAL(VariantKnockoffs::parseEntry("P1_KO3_N40D_VAR_N12D", "knockoff_of=P1_VAR_N12D", info), true)
  TEST_EQUAL(info.knockoff, true)
  TEST_EQUAL(info.family, "P1")
  TEST_EQUAL(info.variant, "N12D")
  TEST_EQUAL(info.tier, "")

  TEST_EQUAL(VariantKnockoffs::parseEntry("P1", "", info), false)
  TEST_EQUAL(VariantKnockoffs::parseEntry("P1_VAR_N12", "", info), false)
  TEST_EQUAL(VariantKnockoffs::parseEntry("P1_VAR_B12D", "", info), false) // B is no residue

  TEST_EQUAL(VariantKnockoffs::variantAccession("P1", 'N', 11, 'D'), "P1_VAR_N12D")
  TEST_EQUAL(VariantKnockoffs::knockoffAccession("P1", 'N', 39, 3, 11, 'D'), "P1_KO3_N40D_VAR_N12D")
}
END_SECTION

START_SECTION(static std::vector<Size> knockoffPositions(const std::string& sequence, Size position, const std::vector<Size>& listed, const ProteaseDigestion& digestion, Size per_variant, UInt64 seed))
{
  const std::string seq = repeats();
  const ProteaseDigestion d = trypsin();

  // the same residue, followed by the same residue class, in another tryptic peptide of the same length
  std::vector<Size> qs = VariantKnockoffs::knockoffPositions(seq, 4, {4}, d, 4, 7);
  TEST_EQUAL(qs.size(), 4)
  for (Size q : qs)
  {
    TEST_EQUAL(seq[q], 'D')
    TEST_EQUAL(seq[q + 1], 'L')
    TEST_EQUAL(q >= 14, true)
  }
  std::vector<Size> sorted = qs;
  std::sort(sorted.begin(), sorted.end());
  TEST_EQUAL(std::unique(sorted.begin(), sorted.end()) == sorted.end(), true) // drawn without replacement
  TEST_EQUAL(qs == VariantKnockoffs::knockoffPositions(seq, 4, {4}, d, 4, 7), true) // deterministic
  TEST_EQUAL(VariantKnockoffs::knockoffPositions(seq, 4, {4}, d, 10, 7).size(), 7) // every other D
  TEST_EQUAL(VariantKnockoffs::knockoffPositions(seq, 4, {4}, d, 0, 7).size(), 0)

  // every D is a listed variant
  std::vector<Size> listed;
  for (Size i = 0; i < 8; ++i) listed.push_back(4 + 10 * i);
  TEST_EQUAL(VariantKnockoffs::knockoffPositions(seq, 4, listed, d, 4, 1).size(), 0)
  // every D but the last, which is 10 residues from the nearest
  listed.pop_back();
  qs = VariantKnockoffs::knockoffPositions(seq, 4, listed, d, 4, 1);
  TEST_EQUAL(qs.size(), 1)
  TEST_EQUAL(qs[0], 74)

  // a next residue of another class (D-L at the site, D-E elsewhere)
  std::string other = "GGGGDLGGGK";
  for (int i = 0; i < 3; ++i) other += "GGGGDEGGGK";
  TEST_EQUAL(VariantKnockoffs::knockoffPositions(other, 4, {4}, d, 4, 1).size(), 0)
  // a next residue of the same class (D-L at the site, D-I elsewhere)
  std::string same = "GGGGDLGGGK";
  for (int i = 0; i < 3; ++i) same += "GGGGDIGGGK";
  TEST_EQUAL(VariantKnockoffs::knockoffPositions(same, 4, {4}, d, 4, 1).size(), 3)
  // a peptide too short (5 residues) holds no knockoff
  std::string short_spans = "GGGGDLGGGK";
  for (int i = 0; i < 3; ++i) short_spans += "GDLGK";
  TEST_EQUAL(VariantKnockoffs::knockoffPositions(short_spans, 4, {4}, d, 4, 1).size(), 0)
  // a sequon is matched by a sequon only (N-L-S at the site, N-L-G elsewhere)
  std::string sequons = "GGGGNLSGGK";
  for (int i = 0; i < 3; ++i) sequons += "GGGGNLGGGK";
  TEST_EQUAL(VariantKnockoffs::knockoffPositions(sequons, 4, {4}, d, 4, 1).size(), 0)
}
END_SECTION

START_SECTION(static std::vector<FASTAFile::FASTAEntry> buildDatabase(const std::vector<FASTAFile::FASTAEntry>& db, const std::vector<Variant>& variants, const std::string& enzyme, Size per_variant, DatabaseStats& stats))
{
  const std::string seq = repeats();
  std::vector<FASTAFile::FASTAEntry> db = {{"sp|P00009|KO_HUMAN", "Knockoff test", seq}, {"P2", "", "PEPTIDEK"}};
  std::vector<VariantKnockoffs::Variant> variants(5);
  variants[0] = {"P00009", 4, 'D', 'N', "VAR_000009", "curated"};
  variants[1] = {"P00009", 4, 'D', 'N', "", ""};  // duplicate
  variants[2] = {"P00009", 5, 'D', 'N', "", ""};  // reference mismatch (L)
  variants[3] = {"P99999", 1, 'E', 'Q', "", ""};  // unknown protein
  variants[4] = {"P2", 1, 'E', 'Q', "", ""};      // a protein too short for knockoffs

  VariantKnockoffs::DatabaseStats stats;
  const auto entries = VariantKnockoffs::buildDatabase(db, variants, "Trypsin", 4, stats);
  TEST_EQUAL(stats.variant_entries, 2)
  TEST_EQUAL(stats.knockoff_entries, 4)
  TEST_EQUAL(stats.variants_without_knockoff, 1)
  TEST_EQUAL(stats.duplicates, 1)
  TEST_EQUAL(stats.reference_mismatch, 1)
  TEST_EQUAL(stats.unknown_protein, 1)
  TEST_EQUAL(entries.size(), 6)

  TEST_EQUAL(entries[0].identifier, "sp|P00009|KO_HUMAN_VAR_D5N")
  TEST_EQUAL(entries[0].description, "variant_of=sp|P00009|KO_HUMAN variant=D5N id=VAR_000009 tier=curated")
  TEST_EQUAL(entries[0].sequence[4], 'N')
  TEST_EQUAL(entries[1].identifier, "P2_VAR_E2Q")
  TEST_EQUAL(entries[1].sequence, "PQPTIDEK")
  for (Size i = 2; i < entries.size(); ++i)
  {
    const auto& e = entries[i];
    TEST_EQUAL(e.identifier.rfind("sp|P00009|KO_HUMAN_KO" + std::to_string(i - 1) + "_D", 0), 0)
    TEST_EQUAL(e.description, "knockoff_of=sp|P00009|KO_HUMAN_VAR_D5N tier=curated")
    Size changed = 0;
    for (Size j = 0; j < seq.size(); ++j)
    {
      if (e.sequence[j] != seq[j])
      {
        ++changed;
        TEST_EQUAL(e.sequence[j], 'N')
      }
    }
    TEST_EQUAL(changed, 1)
  }
  // deterministic
  TEST_EQUAL(VariantKnockoffs::buildDatabase(db, variants, "Trypsin", 4, stats)[3].identifier, entries[3].identifier)
  // no knockoffs
  TEST_EQUAL(VariantKnockoffs::buildDatabase(db, variants, "Trypsin", 0, stats).size(), 2)
  TEST_EQUAL(stats.variants_without_knockoff, 0)

  // the index of the database
  std::vector<FASTAFile::FASTAEntry> full = db;
  full.insert(full.end(), entries.begin(), entries.end());
  const VariantKnockoffs::Index idx = VariantKnockoffs::indexDatabase(full);
  TEST_EQUAL(idx.entries.size(), 6)
  TEST_EQUAL(idx.max_knockoffs, 4)
  TEST_EQUAL(idx.knockoffs_per_variant.size(), 1)
  TEST_EQUAL(idx.knockoffs_per_variant.at({"sp|P00009|KO_HUMAN", "D5N"}), 4)
  TEST_EQUAL(idx.entries.at("sp|P00009|KO_HUMAN_VAR_D5N").tier, "curated")
}
END_SECTION

START_SECTION(static std::vector<StratumCount> applyKnockoffQValues(PeptideIdentificationList& ids, const Index& index, const std::function<double(const PeptideHit&)>& rank_of, double cutoff, Size min_tested = 50, bool use_all_hits = false, const std::set<std::string>& decoy_accessions = {}))
{
  // no knockoff entry: nothing changes
  {
    std::vector<FASTAFile::FASTAEntry> db = {{"P1_VAR_D35W", "", ""}};
    PeptideIdentificationList ids;
    ids.push_back(makeId(5.0, 0.001, {"P1_VAR_D35W"}));
    TEST_EQUAL(VariantKnockoffs::applyKnockoffQValues(ids, VariantKnockoffs::indexDatabase(db), rawScore, 0.01).empty(), true)
    TEST_REAL_SIMILAR(ids[0].getHits()[0].getScore(), 0.001)
  }

  // knockoffs above a real PSM raise its q-value and are never accepted
  {
    PeptideIdentificationList ids;
    for (double raw : {8.0, 7.0, 6.0, 5.0}) ids.push_back(makeId(raw, 0.001, {"P1_VAR_D5N"}));
    ids.push_back(makeId(7.5, 0.001, {"P1_KO1_D15N_VAR_D5N"}));
    ids.push_back(makeId(5.5, 0.001, {"P1_KO2_D25N_VAR_D5N"}));
    ids.push_back(makeId(9.0, 0.001, {"P1", "P1_KO1_D15N_VAR_D5N"}));
    const auto counts = VariantKnockoffs::applyKnockoffQValues(ids, index(), rawScore, 0.01);
    // FDP at 8: (1/2 + 0) / 1; at 7: (1/2 + 1/2) / 2; at 6: 1 / 3; at 5: (1/2 + 1) / 4. Running minimum from below.
    TEST_REAL_SIMILAR(hitOf(ids, 5.0).getScore(), 0.375)
    TEST_REAL_SIMILAR(hitOf(ids, 6.0).getScore(), 1.0 / 3.0)
    TEST_REAL_SIMILAR(hitOf(ids, 7.0).getScore(), 1.0 / 3.0)
    TEST_REAL_SIMILAR(hitOf(ids, 8.0).getScore(), 1.0 / 3.0)
    TEST_REAL_SIMILAR(double(hitOf(ids, 8.0).getMetaValue(VariantKnockoffs::KNOCKOFF_Q_VALUE)), 1.0 / 3.0)
    TEST_EQUAL(hitOf(ids, 8.0).getMetaValue(VariantKnockoffs::KNOCKOFF_STRATUM), "deamidation-like")
    TEST_REAL_SIMILAR(hitOf(ids, 7.5).getScore(), 1.0)
    TEST_EQUAL(hitOf(ids, 7.5).metaValueExists(VariantKnockoffs::IS_KNOCKOFF), true)
    TEST_EQUAL(hitOf(ids, 8.0).metaValueExists(VariantKnockoffs::IS_KNOCKOFF), false)
    // canonical PSMs keep their q-value, without the knockoff accession
    const PeptideHit& c = hitOf(ids, 9.0);
    TEST_REAL_SIMILAR(c.getScore(), 0.001)
    TEST_EQUAL(c.getPeptideEvidences().size(), 1)
    TEST_EQUAL(c.getPeptideEvidences()[0].getProteinAccession(), "P1")
    TEST_EQUAL(counts.size(), 1)
    TEST_EQUAL(counts[0].name, "deamidation-like")
    TEST_EQUAL(counts[0].real_before, 4)
    TEST_EQUAL(counts[0].tested_before, 4)
    TEST_EQUAL(counts[0].real_after, 0)
    TEST_REAL_SIMILAR(counts[0].knockoff_weighted, 1.0)
    TEST_REAL_SIMILAR(counts[0].falseShare(), 0.25)
  }

  // real PSMs far above the knockoffs keep a small q-value: (1/2) / 100, the floor of the estimate
  {
    PeptideIdentificationList ids;
    for (int i = 0; i < 100; ++i) ids.push_back(makeId(100.0 - i, 0.001, {"P1_VAR_D5N"}));
    ids.push_back(makeId(-5.0, 0.5, {"P1_KO1_D15N_VAR_D5N"}));
    const auto counts = VariantKnockoffs::applyKnockoffQValues(ids, index(), rawScore, 0.01);
    TEST_REAL_SIMILAR(hitOf(ids, 1.0).getScore(), 0.005)
    TEST_REAL_SIMILAR(hitOf(ids, 100.0).getScore(), 0.005)
    TEST_EQUAL(counts[0].real_after, 100)
    TEST_REAL_SIMILAR(counts[0].knockoff_weighted, 0.0) // the knockoff PSM is not accepted at 0.01
  }

  // a variant without knockoffs is not accepted
  {
    PeptideIdentificationList ids;
    ids.push_back(makeId(5.0, 0.001, {"P1_VAR_D35W"}));
    const auto counts = VariantKnockoffs::applyKnockoffQValues(ids, index(), rawScore, 0.01);
    TEST_REAL_SIMILAR(ids[0].getHits()[0].getScore(), 1.0)
    TEST_EQUAL(counts.size(), 1)
    TEST_EQUAL(counts[0].name, "other")
    TEST_EQUAL(counts[0].real_before, 1)
    TEST_EQUAL(counts[0].tested_before, 0)
    TEST_EQUAL(counts[0].real_after, 0)
    TEST_EQUAL(std::isnan(counts[0].falseShare()), true)
  }

  // a lower-ranked hit of a tested variant is not estimated, so it is not accepted (unless use_all_hits)
  {
    PeptideIdentificationList ids;
    ids.push_back(makeId(50.0, 0.001, {"P1"}));
    PeptideHit second = makeId(40.0, 0.001, {"P1_VAR_D5N"}).getHits()[0];
    ids[0].insertHit(second);
    PeptideIdentificationList all_hits = ids;
    VariantKnockoffs::applyKnockoffQValues(ids, index(), rawScore, 0.01);
    TEST_REAL_SIMILAR(ids[0].getHits()[0].getScore(), 0.001)
    TEST_REAL_SIMILAR(ids[0].getHits()[1].getScore(), 1.0)
    VariantKnockoffs::applyKnockoffQValues(all_hits, index(), rawScore, 0.01, 50, true);
    TEST_EQUAL(all_hits[0].getHits()[1].metaValueExists(VariantKnockoffs::KNOCKOFF_Q_VALUE), true)
  }

  // decoys and PSMs of a canonical protein are left alone
  {
    PeptideIdentificationList ids;
    ids.push_back(makeId(5.0, 0.001, {"P1_VAR_D35W"}));
    ids[0].getHits()[0].setTargetDecoyType(PeptideHit::TargetDecoyType::DECOY);
    ids.push_back(makeId(6.0, 0.001, {"P1", "P1_VAR_D35W"}));
    TEST_EQUAL(VariantKnockoffs::applyKnockoffQValues(ids, index(), rawScore, 0.01).empty(), true)
    TEST_REAL_SIMILAR(hitOf(ids, 5.0).getScore(), 0.001)
    TEST_REAL_SIMILAR(hitOf(ids, 6.0).getScore(), 0.001)
  }

  // a target+decoy hit whose target proteins are all variant entries is a variant PSM: without the decoy accessions
  // it escapes the knockoff q-value
  {
    PeptideIdentificationList ids;
    ids.push_back(makeId(5.0, 0.001, {"DECOY_P9", "P1_VAR_D35W"}));
    ids[0].getHits()[0].setTargetDecoyType(PeptideHit::TargetDecoyType::TARGET_DECOY);
    PeptideIdentificationList with_decoys = ids;
    TEST_EQUAL(VariantKnockoffs::applyKnockoffQValues(ids, index(), rawScore, 0.01).empty(), true)
    TEST_REAL_SIMILAR(ids[0].getHits()[0].getScore(), 0.001)
    const auto counts = VariantKnockoffs::applyKnockoffQValues(with_decoys, index(), rawScore, 0.01, 50, false, {"DECOY_P9"});
    TEST_EQUAL(counts.size(), 1)
    TEST_EQUAL(counts[0].real_before, 1)
    TEST_REAL_SIMILAR(with_decoys[0].getHits()[0].getScore(), 1.0) // a variant without knockoffs
  }

  // prior tiers: a tier with enough tested PSMs gets its own estimate, a small one joins its class's pool
  {
    std::vector<FASTAFile::FASTAEntry> db = {
      {"P1_VAR_D5N", "tier=sample", ""},
      {"P1_KO1_D15N_VAR_D5N", "tier=sample", ""},
      {"P2_VAR_D5N", "tier=catalogue", ""},
      {"P2_KO1_D15N_VAR_D5N", "tier=catalogue", ""},
      {"P3_VAR_D5N", "tier=common", ""},
      {"P3_KO1_D15N_VAR_D5N", "tier=common", ""}};
    const auto idx = VariantKnockoffs::indexDatabase(db);
    PeptideIdentificationList ids;
    // sample: 4 real PSMs, no knockoff PSM; catalogue: 2 real PSMs below a knockoff; common: 1 real PSM
    for (double raw : {20.0, 19.0, 18.0, 17.0}) ids.push_back(makeId(raw, 0.001, {"P1_VAR_D5N"}));
    for (double raw : {10.0, 9.0}) ids.push_back(makeId(raw, 0.001, {"P2_VAR_D5N"}));
    ids.push_back(makeId(11.0, 0.001, {"P2_KO1_D15N_VAR_D5N"}));
    ids.push_back(makeId(30.0, 0.001, {"P3_VAR_D5N"}));
    const auto counts = VariantKnockoffs::applyKnockoffQValues(ids, idx, rawScore, 0.01, 2);
    TEST_EQUAL(counts.size(), 3)
    TEST_EQUAL(counts[0].name, "deamidation-like/catalogue")
    TEST_EQUAL(counts[1].name, "deamidation-like/pooled")
    TEST_EQUAL(counts[2].name, "deamidation-like/sample")
    TEST_REAL_SIMILAR(hitOf(ids, 17.0).getScore(), 0.25)      // (1/1 + 0) / 4
    TEST_REAL_SIMILAR(hitOf(ids, 9.0).getScore(), 1.0)        // (1 + 1) / 2, capped
    TEST_REAL_SIMILAR(hitOf(ids, 30.0).getScore(), 1.0)       // (1 + 0) / 1
    TEST_EQUAL(hitOf(ids, 30.0).getMetaValue(VariantKnockoffs::KNOCKOFF_STRATUM), "deamidation-like/pooled")
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
