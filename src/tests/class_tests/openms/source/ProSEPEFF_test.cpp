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
    s.specificity = EnzymaticDigestion::SPEC_FULL;
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
  // G10A (in VLGFHQR), W48C (in WTER), and annotations that are skipped: S2S (the residue itself), Q13* (stop), 120A
  // (past the end), K7X (no residue of known mass); G10A repeated
  std::vector<FASTAFile::FASTAEntry> db = {
    {"sp:P01", "\\PName=Test \\VariantSimple=(10|A)(48|C)(2|S)(13|*)(120|A)(7|X)(10|A)", P01},
    {"P02", "plain FASTA", "PEPTIDEKAAAAR"}};
  ProSEPEFF peff;
  TEST_EQUAL(peff.addVariantWindows(db, settings(10)), 2)
  TEST_EQUAL(db.size(), 4)
  TEST_EQUAL(peff.size(), 2)
  TEST_STRING_EQUAL(db[2].identifier, "sp:P01_VAR_G10A")
  TEST_STRING_EQUAL(db[3].identifier, "sp:P01_VAR_W48C")
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

  // the regions for the fragment index
  const auto regions = peff.peptideRegions();
  TEST_EQUAL(regions->size(), 2)
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
}
END_SECTION

START_SECTION(([EXTRA] ProSEAlgorithm database:peff_variants identifies a variant peptide))
{
  // G10A turns VLGFHQR into VLAFHQR; the spectra: the variant peptide and peptides of the protein
  const vector<FASTAFile::FASTAEntry> fasta_db = {{"sp:P01", "\\PName=Test \\VariantSimple=(10|A)(48|C)", P01}};
  TheoreticalSpectrumGenerator tsg;
  Param tsg_param = tsg.getParameters();
  tsg_param.setValue("add_first_prefix_ion", "true");
  tsg_param.setValue("add_metainfo", "true");
  tsg.setParameters(tsg_param);
  PeakMap spectra;
  double rt = 100.0;
  for (const char* sequence : {"VLAFHQR", "MPNASTIC(Carbamidomethyl)YWDLK", "THQPSANLDIK", "GWSADELK", "EGFVRTHQPSANLDIK"})
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
                          const Param& extra = Param())
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
    return algo.search(spectra, fasta_db, protein_ids, peptide_ids);
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
    TEST_EQUAL(search("simple", protein_ids, peptide_ids, extra) == ProSEAlgorithm::ExitCodes::EXECUTION_OK, true)
    bool found_variant = false;
    for (const PeptideIdentification& id : peptide_ids)
    {
      for (const PeptideHit& hit : id.getHits())
      {
        for (const PeptideEvidence& evidence : hit.getPeptideEvidences())
        {
          TEST_EQUAL(evidence.getProteinAccession().find(ProSEPEFF::windowSeparator()), std::string::npos)
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
