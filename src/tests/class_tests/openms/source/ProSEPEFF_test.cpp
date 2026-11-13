// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// $Authors:  $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/PROSE/ProSEPEFF.h>
///////////////////////////

#include <OpenMS/ANALYSIS/ID/PROSE/ProSEAlgorithm.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <algorithm>
#include <fstream>
#include <set>

using namespace OpenMS;
using namespace std;

namespace
{
  // 100 residues; tryptic segments (1-based): MSDER 1-5, EK 6-7, VLGFHQR 8-14, MPNASTICYWDLK 15-27, EGFVR 28-32,
  // THQPSANLDIK 33-43, CMYK 44-47, WTER 48-51, HASGDFLK 52-59, PIVEQNCTMYR 52-70 (KP: no cleavage), GWSADELK 71-78,
  // HPFNQGTICMSYR 79-91, EWDAVLKPH 92-100
  const std::string P01 = "MSDEREKVLGFHQRMPNASTICYWDLKEGFVRTHQPSANLDIKCMYKWTERHASGDFLKPIVEQNCTMYRGWSADELKHPFNQGTICMSYREWDAVLKPH";

  ProSEPEFF::Settings settings(Size max_length, bool segment_decoys = true)
  {
    ProSEPEFF::Settings s;
    s.enzyme = "Trypsin";
    s.peptide_max_length = max_length;
    s.segment_decoys = segment_decoys;
    s.clip_nterm_methionine = true;
    return s;
  }

  // the cleavage sites of Trypsin in @p sequence (segment starts, with 0)
  std::set<Size> sites(const std::string& sequence)
  {
    ProteaseDigestion d;
    d.setEnzyme("Trypsin");
    d.setMissedCleavages(0);
    std::vector<std::pair<Size, Size>> spans;
    d.digestUnmodified(sequence, spans, 1, 0);
    std::set<Size> starts{0};
    for (const auto& span : spans) { starts.insert(span.first); }
    return starts;
  }
}

START_TEST(ProSEPEFF, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((static const std::string& variantMetaValue()))
{
  TEST_STRING_EQUAL(ProSEPEFF::variantMetaValue(), "PEFF_variants")
  TEST_STRING_EQUAL(ProSEPEFF::windowSeparator(), "_VAR_")
}
END_SECTION

START_SECTION((Size addVariantWindows(std::vector<FASTAFile::FASTAEntry>& db, const Settings& settings)))
{
  // G10A (in VLGFHQR), W48C (in WTER), Q13* (stop codon), and annotations that are not searched: S2S (the residue
  // itself), 120A (past the end), K7X (no residue of known mass), a \\VariantComplex; G10A repeated
  std::vector<FASTAFile::FASTAEntry> db = {
    {"sp:P01", "\\PName=Test \\VariantSimple=(10|A)(48|C)(2|S)(13|*)(120|A)(7|X)(10|A) \\VariantComplex=(30|32|AG)", P01},
    {"P02", "plain FASTA", "PEPTIDEKAAAAR"}};
  ProSEPEFF peff;
  TEST_EQUAL(peff.addVariantWindows(db, settings(10)), 3)
  TEST_EQUAL(db.size(), 5)
  TEST_EQUAL(peff.size(), 3)
  TEST_STRING_EQUAL(db[2].identifier, "sp:P01_VAR_G10A")
  TEST_STRING_EQUAL(db[3].identifier, "sp:P01_VAR_W48C")
  TEST_STRING_EQUAL(db[4].identifier, "sp:P01_VAR_Q13*")
  TEST_EQUAL(peff.find("P02") == nullptr, true)

  const std::set<Size> cleavage = sites(P01);
  const Size margin = 10 + 1;
  for (Size i = 2; i < 4; ++i)
  {
    const ProSEPEFF::Window* w = peff.find(db[i].identifier);
    TEST_EQUAL(w != nullptr, true)
    ABORT_IF(w == nullptr)
    TEST_STRING_EQUAL(w->parent, "sp:P01")
    TEST_EQUAL(w->decoy, false)
    TEST_EQUAL(w->length, db[i].sequence.size())
    const Size v = (i == 2 ? 10 : 48) - 1; // 0-based variant position
    // the window is the protein with the variant, from a cleavage site to a cleavage site (or the protein ends)
    std::string expected = P01;
    expected[v] = (i == 2 ? 'A' : 'C');
    TEST_STRING_EQUAL(db[i].sequence, expected.substr(w->shift, w->length))
    TEST_EQUAL(cleavage.count(w->shift), 1)
    TEST_EQUAL(w->shift + w->length == P01.size() || cleavage.count(w->shift + w->length) == 1, true)
    // the region: the variant, the residue before and the two after it
    TEST_EQUAL(w->region.first, v - 1 - w->shift)
    TEST_EQUAL(w->region.second, v + 3 - w->shift)
    // a peptide of at most peptide:max_size residues that overlaps the region does not reach an end of the window
    TEST_EQUAL(w->shift == 0 || w->region.first >= margin, true)
    TEST_EQUAL(w->shift + w->length == P01.size() || w->length - w->region.second >= margin, true)
    // no window starts with M unless it starts the protein
    TEST_EQUAL(w->shift == 0 || db[i].sequence[0] != 'M', true)
  }
  TEST_STRING_EQUAL(peff.find("sp:P01_VAR_G10A")->variant, "G10A")
  TEST_EQUAL(peff.find("sp:P01_VAR_G10A")->shift, 0) // close to the N-terminus

  // the stop codon: the protein up to Q13, without it; the region is its last residue (H12)
  const ProSEPEFF::Window* stop = peff.find("sp:P01_VAR_Q13*");
  TEST_STRING_EQUAL(db[4].sequence, P01.substr(0, 12))
  TEST_EQUAL(stop->shift, 0)
  TEST_EQUAL(stop->region.first, 11)
  TEST_EQUAL(stop->region.second, 12)

  // the status of every annotated variant, in annotation order (repeats once)
  const std::vector<ProSEPEFF::VariantStatus>& status = peff.getVariantStatus();
  TEST_EQUAL(status.size(), 7)
  ABORT_IF(status.size() != 7)
  // the complex variant first (all of an entry's \\VariantComplex before its \\VariantSimple)
  TEST_STRING_EQUAL(status[0].variant, "30-32|AG")
  TEST_EQUAL(status[0].kind == ProSEPEFF::VariantKind::COMPLEX, true)
  TEST_EQUAL(status[0].searched, false)
  TEST_STRING_EQUAL(status[0].reason, "complex_not_supported")
  TEST_EQUAL(status[0].position, 30)
  const std::vector<std::string> variants = {"G10A", "W48C", "S2S", "Q13*", "120A", "K7X"};
  const std::vector<std::string> reasons = {"", "", "same_as_reference", "", "out_of_range", "unknown_residue"};
  for (Size k = 0; k < variants.size(); ++k)
  {
    TEST_STRING_EQUAL(status[k + 1].protein, "sp:P01")
    TEST_STRING_EQUAL(status[k + 1].variant, variants[k])
    TEST_STRING_EQUAL(status[k + 1].reason, reasons[k])
    TEST_EQUAL(status[k + 1].searched, reasons[k].empty())
    TEST_STRING_EQUAL(status[k + 1].window, reasons[k].empty() ? "sp:P01_VAR_" + variants[k] : "")
  }
  TEST_EQUAL(status[4].kind == ProSEPEFF::VariantKind::STOP_GAIN, true)
  TEST_STRING_EQUAL(ProSEPEFF::kindName(status[4].kind), "stop_gain")

  // a stop codon in the middle: P60* (1-based) ends the protein after HASGDFLK (52-59), which the protein itself does not
  // have (KP is not cleaved)
  std::vector<FASTAFile::FASTAEntry> db_stop = {{"sp:P01", "\\VariantSimple=(60|*)", P01}};
  ProSEPEFF peff_stop;
  TEST_EQUAL(peff_stop.addVariantWindows(db_stop, settings(10)), 1)
  const ProSEPEFF::Window* w60 = peff_stop.find("sp:P01_VAR_P60*");
  ABORT_IF(w60 == nullptr)
  TEST_EQUAL(w60->shift + w60->length, 59) // ends at the new C-terminus
  TEST_STRING_EQUAL(db_stop[1].sequence, P01.substr(w60->shift, w60->length))
  TEST_EQUAL(cleavage.count(w60->shift), 1)
  TEST_EQUAL(w60->region.first + w60->shift, 58)
  TEST_EQUAL(w60->region.second + w60->shift, 59)
  TEST_EQUAL(w60->region.first >= margin, true)
  // decoys: the last segment (HASGDFLK, 51-58 0-based) may change
  db_stop.push_back({"DECOY_sp:P01_VAR_P60*", "", ""});
  peff_stop.addDecoys(db_stop, "DECOY_");
  const ProSEPEFF::Window* d60 = peff_stop.find("DECOY_sp:P01_VAR_P60*");
  ABORT_IF(d60 == nullptr)
  TEST_EQUAL(d60->region.first + d60->shift, 51)
  TEST_EQUAL(d60->region.second + d60->shift, 59)

  // peptide:max_size 0 (whole proteins): the truncated protein; decoys per segment: the last segment, in place
  std::vector<FASTAFile::FASTAEntry> db_whole_stop = {{"sp:P01", "\\VariantSimple=(60|*)", P01}, {"DECOY_sp:P01_VAR_P60*", "", ""}};
  ProSEPEFF peff_whole_stop;
  peff_whole_stop.addVariantWindows(db_whole_stop, settings(0));
  peff_whole_stop.addDecoys(db_whole_stop, "DECOY_");
  const ProSEPEFF::Window* ws = peff_whole_stop.find("sp:P01_VAR_P60*");
  const ProSEPEFF::Window* ds = peff_whole_stop.find("DECOY_sp:P01_VAR_P60*");
  ABORT_IF(ws == nullptr || ds == nullptr)
  TEST_EQUAL(ws->shift, 0)
  TEST_EQUAL(ws->length, 59)
  TEST_EQUAL(ds->shift, 0)
  TEST_EQUAL(ds->region.first, 51)
  TEST_EQUAL(ds->region.second, 59)
  // ... decoys reversed whole: the decoy of the truncated protein is the end of the protein's decoy (which keeps its M)
  std::vector<FASTAFile::FASTAEntry> db_whole_mirror = {{"sp:P01", "\\VariantSimple=(60|*)", P01}, {"DECOY_sp:P01_VAR_P60*", "", ""}};
  ProSEPEFF peff_whole_mirror;
  peff_whole_mirror.addVariantWindows(db_whole_mirror, settings(0, false));
  peff_whole_mirror.addDecoys(db_whole_mirror, "DECOY_");
  const ProSEPEFF::Window* dm = peff_whole_mirror.find("DECOY_sp:P01_VAR_P60*");
  ABORT_IF(dm == nullptr)
  TEST_EQUAL(dm->shift, P01.size() - 59)
  TEST_EQUAL(dm->region.first, 0)
  TEST_EQUAL(dm->region.second, 1)

  // the regions for the fragment index
  const auto regions = peff.peptideRegions();
  TEST_EQUAL(regions->size(), 3)
  TEST_EQUAL(regions->at("sp:P01_VAR_W48C") == peff.find("sp:P01_VAR_W48C")->region, true)

  // peptide:max_size 0 (no maximum): the whole protein
  std::vector<FASTAFile::FASTAEntry> db_whole = {{"sp:P01", "\\VariantSimple=(48|C)", P01}};
  ProSEPEFF peff_whole;
  TEST_EQUAL(peff_whole.addVariantWindows(db_whole, settings(0)), 1)
  TEST_EQUAL(peff_whole.find("sp:P01_VAR_W48C")->shift, 0)
  TEST_EQUAL(peff_whole.find("sp:P01_VAR_W48C")->length, P01.size())

  // decoys reverse the whole entry: no cleavage-site alignment, the margin around the region
  std::vector<FASTAFile::FASTAEntry> db_mirror = {{"sp:P01", "\\VariantSimple=(48|C)", P01}};
  ProSEPEFF peff_mirror;
  TEST_EQUAL(peff_mirror.addVariantWindows(db_mirror, settings(10, false)), 1)
  const ProSEPEFF::Window* m = peff_mirror.find("sp:P01_VAR_W48C");
  TEST_EQUAL(m->shift, 47 - 1 - margin) // region starts at 46 (0-based), no M at 35
  TEST_EQUAL(m->length, 4 + 2 * margin)

  // an identifier that exists already gets a number
  std::vector<FASTAFile::FASTAEntry> db_taken = {{"sp:P01", "\\VariantSimple=(48|C)", P01}, {"sp:P01_VAR_W48C", "", "PEPTIDER"}};
  ProSEPEFF peff_taken;
  TEST_EQUAL(peff_taken.addVariantWindows(db_taken, settings(10)), 1)
  TEST_STRING_EQUAL(db_taken.back().identifier, "sp:P01_VAR_W48C_2")
}
END_SECTION

START_SECTION((void addDecoys(const std::vector<FASTAFile::FASTAEntry>& db, const std::string& decoy_prefix)))
{
  std::vector<FASTAFile::FASTAEntry> db = {{"sp:P01", "\\VariantSimple=(48|C)", P01}};
  ProSEPEFF peff;
  peff.addVariantWindows(db, settings(10));
  const ProSEPEFF::Window target = *peff.find("sp:P01_VAR_W48C");
  // as ProSE names generated decoys
  db.push_back({"DECOY_sp:P01", "", P01});
  db.push_back({"DECOY_sp:P01_VAR_W48C", "", db[1].sequence});
  peff.addDecoys(db, "DECOY_");
  TEST_EQUAL(peff.size(), 2)
  const ProSEPEFF::Window* decoy = peff.find("DECOY_sp:P01_VAR_W48C");
  TEST_EQUAL(decoy != nullptr, true)
  ABORT_IF(decoy == nullptr)
  TEST_EQUAL(decoy->decoy, true)
  TEST_STRING_EQUAL(decoy->parent, "DECOY_sp:P01")
  TEST_EQUAL(decoy->shift, target.shift) // decoys per cleavage segment: in place
  // the segments of the changed residues K47, W48C, T49, E50 (1-based): CMYK and WTER, protein positions [43, 51) (0-based)
  TEST_EQUAL(decoy->region.first + decoy->shift, 43)
  TEST_EQUAL(decoy->region.second + decoy->shift, 51)
  TEST_EQUAL(peff.peptideRegions()->count("DECOY_sp:P01_VAR_W48C"), 1)

  // reversed entries: mirrored positions
  std::vector<FASTAFile::FASTAEntry> db_mirror = {{"sp:P01", "\\VariantSimple=(48|C)", P01}};
  ProSEPEFF peff_mirror;
  peff_mirror.addVariantWindows(db_mirror, settings(10, false));
  const ProSEPEFF::Window t = *peff_mirror.find("sp:P01_VAR_W48C");
  db_mirror.push_back({"DECOY_sp:P01_VAR_W48C", "", ""});
  peff_mirror.addDecoys(db_mirror, "DECOY_");
  const ProSEPEFF::Window* d = peff_mirror.find("DECOY_sp:P01_VAR_W48C");
  TEST_EQUAL(d->region.first, t.length - t.region.second)
  TEST_EQUAL(d->region.second, t.length - t.region.first)
  // the decoy protein keeps its initial M in place (peptide:clip_nterm_methionine), the window does not start the protein
  TEST_EQUAL(d->shift, P01.size() - (t.shift + t.length) + 1)
}
END_SECTION

START_SECTION((void mapToParents(std::vector<ProteinIdentification>& protein_ids, PeptideIdentificationList& peptide_ids) const))
{
  std::vector<FASTAFile::FASTAEntry> db = {{"sp:P01", "\\VariantSimple=(48|C)", P01}};
  ProSEPEFF peff;
  peff.addVariantWindows(db, settings(10));
  const ProSEPEFF::Window w = *peff.find("sp:P01_VAR_W48C");
  const Int s = static_cast<Int>(w.shift);

  std::vector<ProteinIdentification> protein_ids(1);
  {
    ProteinHit window_hit;
    window_hit.setAccession("sp:P01_VAR_W48C");
    window_hit.setSequence(db[1].sequence);
    window_hit.setMetaValue("target_decoy", "target");
    protein_ids[0].setHits({window_hit});
  }
  PeptideIdentificationList peptide_ids(1);
  std::vector<PeptideHit> hits(3);
  // 1. CMYKCTER (44-51, 0-based 43-50): only in the window, overlaps the variant
  hits[0].setSequence(AASequence::fromString("CMYKCTER"));
  hits[0].setPeptideEvidences({PeptideEvidence("sp:P01_VAR_W48C", 43 - s, 50 - s, 'K', 'H')});
  hits[0].setMetaValue("protein_references", "unique");
  // 2. THQPSANLDIK (32-42): a peptide of the protein that the window holds outside its region
  hits[1].setSequence(AASequence::fromString("THQPSANLDIK"));
  hits[1].setPeptideEvidences({PeptideEvidence("sp:P01", 32, 42, 'R', 'C'), PeptideEvidence("sp:P01_VAR_W48C", 32 - s, 42 - s, 'R', 'C')});
  hits[1].setMetaValue("protein_references", "non-unique");
  // 3. a peptide in the region that occurs in another protein as well
  hits[2].setSequence(AASequence::fromString("CMYKCTER"));
  hits[2].setPeptideEvidences({PeptideEvidence("sp:P01_VAR_W48C", 43 - s, 50 - s, 'K', 'H'), PeptideEvidence("P99", 0, 7, '[', ']')});
  hits[2].setMetaValue("protein_references", "non-unique");
  peptide_ids[0].setHits(hits);

  peff.mapToParents(protein_ids, peptide_ids);

  const std::vector<PeptideHit>& mapped = peptide_ids[0].getHits();
  TEST_EQUAL(mapped[0].getPeptideEvidences().size(), 1)
  TEST_STRING_EQUAL(mapped[0].getPeptideEvidences()[0].getProteinAccession(), "sp:P01")
  TEST_EQUAL(mapped[0].getPeptideEvidences()[0].getStart(), 43)
  TEST_EQUAL(mapped[0].getPeptideEvidences()[0].getEnd(), 50)
  TEST_EQUAL(mapped[0].getPeptideEvidences()[0].getAABefore(), 'K')
  TEST_STRING_EQUAL(mapped[0].getMetaValue(ProSEPEFF::variantMetaValue()).toString(), "sp:P01:W48C")
  TEST_STRING_EQUAL(mapped[0].getMetaValue("protein_references").toString(), "unique")

  TEST_EQUAL(mapped[1].getPeptideEvidences().size(), 1) // the window's evidence is the protein's
  TEST_STRING_EQUAL(mapped[1].getPeptideEvidences()[0].getProteinAccession(), "sp:P01")
  TEST_EQUAL(mapped[1].metaValueExists(ProSEPEFF::variantMetaValue()), false)
  TEST_STRING_EQUAL(mapped[1].getMetaValue("protein_references").toString(), "unique")

  TEST_EQUAL(mapped[2].getPeptideEvidences().size(), 2)
  TEST_EQUAL(mapped[2].metaValueExists(ProSEPEFF::variantMetaValue()), false) // occurs without the variant

  // the window's protein hit became one of its protein, without the window's residues
  TEST_EQUAL(protein_ids[0].getHits().size(), 1)
  TEST_STRING_EQUAL(protein_ids[0].getHits()[0].getAccession(), "sp:P01")
  TEST_STRING_EQUAL(protein_ids[0].getHits()[0].getSequence(), "")
  TEST_STRING_EQUAL(protein_ids[0].getHits()[0].getMetaValue("target_decoy").toString(), "target")

  // a decoy hit of the decoy window: moved to the decoy protein and labelled as well
  db.push_back({"DECOY_sp:P01_VAR_W48C", "", ""});
  peff.addDecoys(db, "DECOY_");
  const ProSEPEFF::Window d = *peff.find("DECOY_sp:P01_VAR_W48C");
  PeptideIdentificationList decoy_ids(1);
  PeptideHit decoy_hit;
  decoy_hit.setSequence(AASequence::fromString("KYMCRETC"));
  decoy_hit.setPeptideEvidences({PeptideEvidence("DECOY_sp:P01_VAR_W48C", static_cast<Int>(d.region.first), static_cast<Int>(d.region.first) + 7, 'K', 'H')});
  decoy_ids[0].setHits({decoy_hit});
  std::vector<ProteinIdentification> no_proteins;
  peff.mapToParents(no_proteins, decoy_ids);
  const PeptideHit& mapped_decoy = decoy_ids[0].getHits()[0];
  TEST_STRING_EQUAL(mapped_decoy.getPeptideEvidences()[0].getProteinAccession(), "DECOY_sp:P01")
  TEST_EQUAL(mapped_decoy.getPeptideEvidences()[0].getStart(), static_cast<Int>(d.region.first + d.shift))
  TEST_STRING_EQUAL(mapped_decoy.getMetaValue(ProSEPEFF::variantMetaValue()).toString(), "DECOY_sp:P01:W48C")

  // a stop codon at H52: WTER (48-51) ends the truncated protein, which the protein has with H after it. The window's
  // evidence (C-terminal) is the protein's: one evidence, with the protein's residue after it, and no variant
  std::vector<FASTAFile::FASTAEntry> db_stop = {{"sp:P01", "\\VariantSimple=(52|*)", P01}};
  ProSEPEFF peff_stop;
  peff_stop.addVariantWindows(db_stop, settings(10));
  const ProSEPEFF::Window w52 = *peff_stop.find("sp:P01_VAR_H52*");
  const Int s52 = static_cast<Int>(w52.shift);
  PeptideIdentificationList stop_ids(1);
  PeptideHit wter;
  wter.setSequence(AASequence::fromString("WTER"));
  wter.setPeptideEvidences({PeptideEvidence("sp:P01", 47, 50, 'K', 'H'),
                            PeptideEvidence("sp:P01_VAR_H52*", 47 - s52, 50 - s52, 'K', PeptideEvidence::C_TERMINAL_AA)});
  stop_ids[0].setHits({wter});
  peff_stop.mapToParents(no_proteins, stop_ids);
  const PeptideHit& mapped_wter = stop_ids[0].getHits()[0];
  TEST_EQUAL(mapped_wter.getPeptideEvidences().size(), 1)
  TEST_EQUAL(mapped_wter.getPeptideEvidences()[0].getAAAfter(), 'H')
  TEST_EQUAL(mapped_wter.metaValueExists(ProSEPEFF::variantMetaValue()), false)
}
END_SECTION

START_SECTION((std::vector<VariantCounts> countPSMs(const std::vector<const PeptideIdentificationList*>& runs, double q_threshold) const))
{
  std::vector<FASTAFile::FASTAEntry> db = {{"sp:P01", "\\VariantSimple=(10|A)(48|C)(2|S) \\VariantComplex=(30|32|AG)", P01}};
  ProSEPEFF peff;
  peff.addVariantWindows(db, settings(10));
  // status: 30-32|AG (not searched), G10A, W48C, S2S (not searched)
  const auto psm = [](const std::string& sequence, double q, const std::string& variants, const std::string& protein, Int start,
                      bool decoy = false)
  {
    PeptideIdentification id;
    id.setScoreType("q-value");
    id.setHigherScoreBetter(false);
    PeptideHit hit;
    hit.setSequence(AASequence::fromString(sequence));
    hit.setScore(q);
    hit.setMetaValue("target_decoy", decoy ? "decoy" : "target");
    if (!variants.empty()) { hit.setMetaValue(ProSEPEFF::variantMetaValue(), variants); }
    hit.setPeptideEvidences({PeptideEvidence(protein, start, start + static_cast<Int>(sequence.size()) - 1)});
    id.setHits({hit});
    return id;
  };
  PeptideIdentificationList run1, run2;
  run1.push_back(psm("VLAFHQR", 0.001, "sp:P01:G10A", "sp:P01", 7));
  run1.push_back(psm("VLAFHQR", 0.2, "sp:P01:G10A", "sp:P01", 7));             // above the threshold
  run1.push_back(psm("VLAFHQR", 0.0005, "DECOY_sp:P01:G10A", "DECOY_sp:P01", 7, true)); // a decoy: not counted
  run1.push_back(psm("CMYKWTER", 0.004, "", "sp:P01", 43));                   // the reference of W48C (covers position 48)
  run1.push_back(psm("VLAFHQR", 0.002, "", "sp:P01", 7));  // the G10A peptide, also in another entry (no tag): no reference
  run2.push_back(psm("VLAFHQR", 0.003, "sp:P01:G10A", "sp:P01", 7));
  run2.push_back(psm("CMYKWTER", 0.5, "", "sp:P01", 43));                     // above the threshold
  const std::vector<ProSEPEFF::VariantCounts> counts = peff.countPSMs({&run1, &run2}, 0.01);
  const std::vector<ProSEPEFF::VariantStatus>& status = peff.getVariantStatus();
  ABORT_IF(counts.size() != 4 || status.size() != 4)
  TEST_STRING_EQUAL(status[1].variant, "G10A")
  TEST_EQUAL(counts[1].psms, 3)
  TEST_EQUAL(status[1].reference, 'G')
  TEST_EQUAL(counts[1].passing_psms, 2)
  TEST_EQUAL(counts[1].passing_runs, 2)
  TEST_REAL_SIMILAR(counts[1].best_q, 0.001)
  TEST_EQUAL(counts[1].reference_passing_psms, 0)
  TEST_STRING_EQUAL(status[2].variant, "W48C")
  TEST_EQUAL(counts[2].psms, 0)
  TEST_EQUAL(counts[2].reference_passing_psms, 1)
  TEST_EQUAL(counts[0].psms, 0) // not searched

  // the table
  std::string tsv;
  NEW_TMP_FILE(tsv)
  peff.storeVariantStatus(tsv, {&run1, &run2}, 0.01);
  std::ifstream in(tsv);
  std::vector<std::string> lines;
  for (std::string line; std::getline(in, line);) { lines.push_back(line); }
  TEST_EQUAL(lines.size(), 5)
  ABORT_IF(lines.size() != 5)
  TEST_STRING_EQUAL(lines[0], "protein\tposition\tvariant\tkind\tstatus\treason\twindow\tpsms\tpsms_q\truns_q\tbest_q\treference_psms_q")
  TEST_STRING_EQUAL(lines[1], "sp:P01\t30\t30-32|AG\tcomplex\tnot_searched\tcomplex_not_supported\t\t0\t0\t0\t\t0")
  TEST_STRING_EQUAL(lines[2], "sp:P01\t10\tG10A\tsimple\tsearched\t\tsp:P01_VAR_G10A\t3\t2\t2\t0.001\t0")
  TEST_STRING_EQUAL(lines[4], "sp:P01\t2\tS2S\tsimple\tnot_searched\tsame_as_reference\t\t0\t0\t0\t\t0")
  TEST_EXCEPTION(Exception::UnableToCreateFile, peff.storeVariantStatus("/does/not/exist/status.tsv", {&run1}, 0.01))
}
END_SECTION

START_SECTION(([EXTRA] ProSEAlgorithm database:peff_variants identifies a variant peptide))
{
  // G10A turns VLGFHQR into VLAFHQR, P60* ends the protein with HASGDFLK (KP: not cleaved in the protein); the spectra:
  // these variant peptides and peptides of the protein (HASGDFLKPIVEQNCTMYR covers P60)
  const vector<FASTAFile::FASTAEntry> fasta_db = {{"sp:P01", "\\PName=Test \\VariantSimple=(10|A)(48|C)(60|*)", P01}};
  TheoreticalSpectrumGenerator tsg;
  Param tsg_param = tsg.getParameters();
  tsg_param.setValue("add_first_prefix_ion", "true");
  tsg_param.setValue("add_metainfo", "true");
  tsg.setParameters(tsg_param);
  PeakMap spectra;
  double rt = 100.0;
  for (const char* sequence : {"VLAFHQR", "HASGDFLK", "MPNASTIC(Carbamidomethyl)YWDLK", "THQPSANLDIK", "GWSADELK", "EGFVRTHQPSANLDIK",
                               "HASGDFLKPIVEQNC(Carbamidomethyl)TMYR"})
  {
    const AASequence peptide = AASequence::fromString(sequence);
    MSSpectrum spectrum;
    tsg.getSpectrum(spectrum, peptide, 1, 1);
    spectrum.sortByPosition();
    spectrum.setMSLevel(2);
    spectrum.setRT(rt);
    rt += 1.0;
    Precursor precursor;
    precursor.setMZ(peptide.getMZ(2));
    precursor.setCharge(2);
    spectrum.setPrecursors({precursor});
    spectrum.setNativeID("spectrum=" + StringUtils::toStr(spectra.size()));
    spectra.addSpectrum(std::move(spectrum));
  }

  const auto search = [&](const string& peff_variants, vector<ProteinIdentification>& protein_ids, PeptideIdentificationList& peptide_ids,
                          const Param& extra = Param(), std::shared_ptr<const ProSEPEFF>* variants = nullptr)
  {
    ProSEAlgorithm algo;
    Param p = algo.getParameters();
    p.setValue("precursor:mass_tolerance_lower", 10.0);
    p.setValue("precursor:mass_tolerance_upper", 10.0);
    p.setValue("precursor:mass_tolerance_unit", "ppm");
    p.setValue("fragment:mass_tolerance", 20.0);
    p.setValue("fragment:mass_tolerance_unit", "ppm");
    p.setValue("modifications:fixed", vector<string>{"Carbamidomethyl (C)"});
    p.setValue("modifications:variable", vector<string>{});
    p.setValue("decoys", "generate");
    p.setValue("peptide:min_size", 7);
    p.setValue("peptide:max_size", 30);
    p.setValue("peptide:missed_cleavages", 1);
    p.setValue("database:peff_variants", peff_variants);
    p.update(extra, false);
    algo.setParameters(p);
    const ProSEAlgorithm::ExitCodes exit_code = algo.search(spectra, fasta_db, protein_ids, peptide_ids);
    if (variants != nullptr) { *variants = algo.getVariantWindows(); }
    return exit_code;
  };

  // the variant peptide and its evidence on the protein: protein mapping from the index and by PeptideIndexing, the
  // database in chunks (searched chunk by chunk, mapped with PeptideIndexing), decoys with the termini kept
  vector<Param> configurations(4);
  configurations[1].setValue("peptide:protein_mapping", "PeptideIndexing");
  configurations[2].setValue("database:chunk_size", 3);
  configurations[3].setValue("decoy_method", "keep_termini");
  for (const Param& extra : configurations)
  {
    vector<ProteinIdentification> protein_ids;
    PeptideIdentificationList peptide_ids;
    std::shared_ptr<const ProSEPEFF> variants;
    TEST_EQUAL(search("simple", protein_ids, peptide_ids, extra, &variants) == ProSEAlgorithm::ExitCodes::EXECUTION_OK, true)
    bool found_variant = false, found_stop = false;
    for (const PeptideIdentification& id : peptide_ids)
    {
      for (const PeptideHit& hit : id.getHits())
      {
        for (const PeptideEvidence& evidence : hit.getPeptideEvidences())
        {
          TEST_EQUAL(evidence.getProteinAccession().find(ProSEPEFF::windowSeparator()), std::string::npos)
        }
        if (hit.getSequence().toUnmodifiedString() == "HASGDFLK")
        {
          found_stop = true;
          TEST_STRING_EQUAL(hit.getMetaValue(ProSEPEFF::variantMetaValue()).toString(), "sp:P01:P60*")
          ABORT_IF(hit.getPeptideEvidences().size() != 1)
          TEST_STRING_EQUAL(hit.getPeptideEvidences()[0].getProteinAccession(), "sp:P01")
          TEST_EQUAL(hit.getPeptideEvidences()[0].getStart(), 51)
          TEST_EQUAL(hit.getPeptideEvidences()[0].getEnd(), 58)
          TEST_EQUAL(hit.getPeptideEvidences()[0].getAAAfter(), PeptideEvidence::C_TERMINAL_AA) // the new C-terminus
        }
        if (hit.getSequence().toUnmodifiedString() != "VLAFHQR") { continue; }
        found_variant = true;
        TEST_STRING_EQUAL(hit.getMetaValue(ProSEPEFF::variantMetaValue()).toString(), "sp:P01:G10A")
        TEST_EQUAL(hit.getPeptideEvidences().size(), 1)
        ABORT_IF(hit.getPeptideEvidences().size() != 1)
        TEST_STRING_EQUAL(hit.getPeptideEvidences()[0].getProteinAccession(), "sp:P01")
        TEST_EQUAL(hit.getPeptideEvidences()[0].getStart(), 7)
        TEST_EQUAL(hit.getPeptideEvidences()[0].getEnd(), 13)
      }
    }
    TEST_EQUAL(found_variant, true)
    TEST_EQUAL(found_stop, true)

    // the PSMs of each variant (without q-values: all reported PSMs pass)
    ABORT_IF(variants == nullptr)
    const std::vector<ProSEPEFF::VariantStatus>& status = variants->getVariantStatus();
    const std::vector<ProSEPEFF::VariantCounts> counts = variants->countPSMs({&peptide_ids}, 0.01);
    TEST_EQUAL(status.size(), 3)
    TEST_EQUAL(counts.size(), 3)
    ABORT_IF(counts.size() != 3)
    TEST_STRING_EQUAL(status[0].variant, "G10A")
    TEST_EQUAL(counts[0].psms, 1)
    TEST_EQUAL(counts[0].passing_runs, 1)
    TEST_EQUAL(counts[0].reference_passing_psms, 0)
    TEST_STRING_EQUAL(status[1].variant, "W48C")
    TEST_EQUAL(counts[1].psms, 0)
    TEST_EQUAL(counts[1].passing_runs, 0)
    TEST_EQUAL(counts[1].best_q, -1.0)
    TEST_STRING_EQUAL(status[2].variant, "P60*")
    TEST_EQUAL(counts[2].psms, 1)
    TEST_EQUAL(counts[2].reference_passing_psms, 1) // HASGDFLKPIVEQNCTMYR
    ABORT_IF(protein_ids.empty())
    for (const ProteinHit& protein : protein_ids[0].getHits())
    {
      TEST_EQUAL(protein.getAccession().find(ProSEPEFF::windowSeparator()), std::string::npos)
    }
  }

  // without database:peff_variants the variant peptide is not in the database
  {
    vector<ProteinIdentification> protein_ids;
    PeptideIdentificationList peptide_ids;
    TEST_EQUAL(search("none", protein_ids, peptide_ids) == ProSEAlgorithm::ExitCodes::EXECUTION_OK, true)
    for (const PeptideIdentification& id : peptide_ids)
    {
      for (const PeptideHit& hit : id.getHits())
      {
        TEST_EQUAL(hit.getSequence().toUnmodifiedString() == "VLAFHQR", false)
        TEST_EQUAL(hit.metaValueExists(ProSEPEFF::variantMetaValue()), false)
      }
    }
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
