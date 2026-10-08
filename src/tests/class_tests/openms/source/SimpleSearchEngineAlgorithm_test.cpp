// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/SimpleSearchEngineAlgorithm.h>
///////////////////////////
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <string>

using namespace OpenMS;
using namespace std;

// exposes the protected static preprocessSpectra_
class SimpleSearchEngineAlgorithm_test : public SimpleSearchEngineAlgorithm
{
public:
  using SimpleSearchEngineAlgorithm::preprocessSpectra_;
};

START_TEST(SimpleSearchEngineAlgorithm, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SimpleSearchEngineAlgorithm* ptr = 0;
SimpleSearchEngineAlgorithm* null_ptr = 0;
START_SECTION(SimpleSearchEngineAlgorithm())
{
	ptr = new SimpleSearchEngineAlgorithm();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~SimpleSearchEngineAlgorithm())
{
	delete ptr;
}
END_SECTION

START_SECTION((ExitCodes search(const std::string &in_mzML, const std::string &in_db, std::vector< ProteinIdentification > &prot_ids, std::vector< PeptideIdentification > &pep_ids) const ))
{
  // tested via tool
  NOT_TESTABLE
}
END_SECTION

START_SECTION(([EXTRA] Stop codons in the database: a trailing one is removed, an inner one does not abort the search))
{
  // Sequences translated from genomes (e.g. SGD's yeast database) end with a stop codon ('*'), and
  // a few contain one. P02's VLGFHQ*R has the precursor mass and fragments of VLGFHQR: it used to
  // be scored, which aborted the search, because AASequence parses '*' as a weightless X.
  // DIVSAGSLYL, the C-terminal peptide of P03, is only searchable without the stop codon.
  // P04 is left without residues once its stop codons are removed, and decoy generation
  // used to crash on it.
  vector<FASTAFile::FASTAEntry> fasta_db = {
    {"P01", "Test", "MSDEREKVLGFHQRMPNASTICYWDLKEGFVRTHQPSANLDIK*"},
    {"P02", "Test", "MSTEKVLGFHQ*RGWSADEK*"},
    {"P03", "Test", "MDSTEKLIHRDIVSAGSLYL*"},
    {"P04", "Test", "**"},
  };

  TheoreticalSpectrumGenerator tsg;
  Param tsg_param = tsg.getParameters();
  tsg_param.setValue("add_first_prefix_ion", "true");
  tsg.setParameters(tsg_param);

  PeakMap spectra;
  for (const std::string seq_str : {"VLGFHQR", "DIVSAGSLYL"})
  {
    const AASequence seq = AASequence::fromString(seq_str);
    MSSpectrum spec;
    tsg.getSpectrum(spec, seq, 1, 1);
    spec.sortByPosition();
    spec.setMSLevel(2);
    spec.setRT(100.0 + spectra.size());
    Precursor prec;
    prec.setMZ(seq.getMZ(2));
    prec.setCharge(2);
    spec.setPrecursors({prec});
    spec.setNativeID("spectrum=" + std::to_string(spectra.size()));
    spectra.addSpectrum(std::move(spec));
  }

  std::string tmp_mzml;
  NEW_TMP_FILE(tmp_mzml)
  tmp_mzml += ".mzML";
  FileHandler().storeExperiment(tmp_mzml, spectra, {FileTypes::MZML});
  std::string tmp_fasta;
  NEW_TMP_FILE(tmp_fasta)
  tmp_fasta += ".fasta";
  FASTAFile().store(tmp_fasta, fasta_db);

  SimpleSearchEngineAlgorithm algo;
  Param p = algo.getParameters();
  p.setValue("precursor:mass_tolerance", 10.0);
  p.setValue("precursor:mass_tolerance_unit", "ppm");
  p.setValue("fragment:mass_tolerance", 20.0);
  p.setValue("fragment:mass_tolerance_unit", "ppm");
  p.setValue("modifications:fixed", vector<string>{});
  p.setValue("modifications:variable", vector<string>{});
  p.setValue("peptide:min_size", 7);
  p.setValue("peptide:max_size", 40);
  p.setValue("peptide:missed_cleavages", 1);
  p.setValue("decoys", "true");
  algo.setParameters(p);

  vector<ProteinIdentification> prot_ids;
  PeptideIdentificationList pep_ids;
  auto ec = algo.search(tmp_mzml, tmp_fasta, prot_ids, pep_ids);
  TEST_EQUAL(ec == SimpleSearchEngineAlgorithm::ExitCodes::EXECUTION_OK, true)

  std::map<std::string, const PeptideHit*> top_hits;
  for (PeptideIdentification& pid : pep_ids)
  {
    pid.sort();
    for (const PeptideHit& hit : pid.getHits())
    {
      TEST_EQUAL(hit.getSequence().toString().find('X'), std::string::npos)
    }
    if (!pid.getHits().empty()) top_hits[pid.getHits()[0].getSequence().toString()] = &pid.getHits()[0];
  }
  TEST_EQUAL(top_hits.size(), 2)
  TEST_EQUAL(top_hits.count("VLGFHQR"), 1)
  ABORT_IF(top_hits.count("DIVSAGSLYL") != 1)
  const std::vector<PeptideEvidence>& evidences = top_hits["DIVSAGSLYL"]->getPeptideEvidences();
  ABORT_IF(evidences.size() != 1)
  TEST_STRING_EQUAL(evidences[0].getProteinAccession(), "P03")
  TEST_EQUAL(evidences[0].getAAAfter(), PeptideEvidence::C_TERMINAL_AA)
}
END_SECTION

START_SECTION(([EXTRA] preprocessSpectra_ keeps a fragment ion that has a small peak one isotope spacing below it))
{
  // b2 of TMTpro-YMATQLLAK-TMTpro (see the same section in ProSEAlgorithm_test): the TMTpro label's
  // isotope impurity puts a 5% peak 1.00335 Da below the ion. That peak must not become the
  // monoisotopic peak of the envelope; previously the ion and its +1 isotope were removed as its
  // isotopes. Regular envelopes are still deisotoped: a 1+ envelope keeps its monoisotopic peak
  // only, and a 2+ envelope is converted to its singly charged m/z.
  PeakMap exp;
  MSSpectrum s;
  s.setMSLevel(2);
  s.setRT(1.0);
  Precursor prec;
  prec.setMZ(600.0);
  prec.setCharge(3);
  s.getPrecursors().push_back(prec);
  const std::vector<std::pair<double, float>> peaks = {
    {450.2500, 0.50f}, {450.7517, 0.20f}, {451.2534, 0.05f},     // 2+ envelope
    {598.3157, 0.04f}, {599.3193, 0.74f}, {600.3250, 0.10f},     // shadow peak, ion, +1 isotope
    {700.4000, 1.00f}, {701.4034, 0.35f}, {702.4067, 0.08f}};    // 1+ envelope
  for (const auto& [mz, intensity] : peaks) s.emplace_back(mz, intensity);
  exp.addSpectrum(s);
  SimpleSearchEngineAlgorithm_test::preprocessSpectra_(exp, 20.0, true);

  auto has_peak = [&exp](double mz)
  {
    return std::any_of(exp[0].begin(), exp[0].end(), [mz](const Peak1D& p) { return std::fabs(p.getMZ() - mz) <= 20e-6 * mz; });
  };
  TEST_EQUAL(has_peak(599.3193), true)   // the fragment ion survives
  TEST_EQUAL(has_peak(700.4000), true)   // 1+ envelope: monoisotopic peak kept ...
  TEST_EQUAL(has_peak(701.4034), false)  // ... isotopes removed
  TEST_EQUAL(has_peak(702.4067), false)
  TEST_EQUAL(has_peak(450.2500), false)  // 2+ envelope converted ...
  TEST_EQUAL(has_peak(450.2500 * 2.0 - Constants::PROTON_MASS_U), true)  // ... to its 1+ m/z
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



