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
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <map>
#include <string>

using namespace OpenMS;
using namespace std;

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
  vector<FASTAFile::FASTAEntry> fasta_db = {
    {"P01", "Test", "MSDEREKVLGFHQRMPNASTICYWDLKEGFVRTHQPSANLDIK*"},
    {"P02", "Test", "MSTEKVLGFHQ*RGWSADEK*"},
    {"P03", "Test", "MDSTEKLIHRDIVSAGSLYL*"},
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


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



