// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>

///////////////////////////
#include <OpenMS/ANALYSIS/XLMS/OpenPepXLAlgorithm.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(OpenPepXLAlgorithm, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

OpenPepXLAlgorithm* ptr = 0;
OpenPepXLAlgorithm* null_ptr = 0;
START_SECTION(OpenPepXLAlgorithm())
{
  ptr = new OpenPepXLAlgorithm();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(virtual ~OpenPepXLAlgorithm())
{
  delete ptr;
}
END_SECTION

START_SECTION(ExitCodes run(PeakMap& unprocessed_spectra, ConsensusMap& cfeatures, std::vector<FASTAFile::FASTAEntry>& fasta_db, std::vector<ProteinIdentification>& protein_ids, std::vector<PeptideIdentification>& peptide_ids, OPXLDataStructs::PreprocessedPairSpectra& preprocessed_pair_spectra, std::vector< std::pair<Size, Size> >& spectrum_pairs, std::vector< std::vector< OPXLDataStructs::CrossLinkSpectrumMatch > >& all_top_csms, PeakMap& spectra))

std::vector<FASTAFile::FASTAEntry> fasta_db;
FASTAFile file;
file.load(OPENMS_GET_TEST_DATA_PATH("OpenPepXL_input.fasta"), fasta_db);

PeakMap unprocessed_spectra;
MzMLFile f;

PeakFileOptions options;
options.clearMSLevels();
options.addMSLevel(2);
f.getOptions() = options;
f.load(OPENMS_GET_TEST_DATA_PATH("OpenPepXL_input.mzML"), unprocessed_spectra);

// load linked features
ConsensusMap cfeatures;
ConsensusXMLFile cf;
cf.load(OPENMS_GET_TEST_DATA_PATH("OpenPepXL_input.consensusXML"), cfeatures);

// initialize solution vectors
vector<ProteinIdentification> protein_ids(1);
vector<PeptideIdentification> peptide_ids;

OPXLDataStructs::PreprocessedPairSpectra preprocessed_pair_spectra(0);
vector< pair<Size, Size> > spectrum_pairs;
vector< vector< OPXLDataStructs::CrossLinkSpectrumMatch > > all_top_csms;
PeakMap spectra;

OpenPepXLAlgorithm search_algorithm;
Param algo_param = search_algorithm.getParameters();
algo_param.setValue("modifications:fixed", std::vector<std::string>{"Carbamidomethyl (C)"});
algo_param.setValue("fragment:mass_tolerance", 0.2, "Fragment mass tolerance");
algo_param.setValue("fragment:mass_tolerance_xlinks", 0.3, "Fragment mass tolerance for cross-link ions");
algo_param.setValue("fragment:mass_tolerance_unit", "Da", "Unit of fragment m");
algo_param.setValue("algorithm:number_top_hits", 5, "Number of top hits reported for each spectrum pair");
search_algorithm.setParameters(algo_param);

// run algorithm
OpenPepXLAlgorithm::ExitCodes exit_code = search_algorithm.run(unprocessed_spectra, cfeatures, fasta_db, protein_ids, peptide_ids, preprocessed_pair_spectra, spectrum_pairs, all_top_csms, spectra);

TEST_EQUAL(exit_code, OpenPepXLAlgorithm::EXECUTION_OK)
TEST_EQUAL(protein_ids.size(), 1)
TEST_EQUAL(peptide_ids.size(), 12)
TEST_EQUAL(spectra.size(), 217)
TEST_EQUAL(spectrum_pairs.size(), 25)
TEST_EQUAL(preprocessed_pair_spectra.spectra_linear_peaks.size(), 25)
TEST_EQUAL(all_top_csms.size(), 12)

for (Size i = 0; i < peptide_ids.size(); i += 10)
{
  auto pep_hits = peptide_ids[i].getHits();
  TEST_EQUAL(pep_hits[0].metaValueExists("xl_chain"), false)
  if (pep_hits[0].getMetaValue("xl_type") == "cross-link")
  {
    TEST_EQUAL(pep_hits[0].metaValueExists("BetaPepEv:pre"), true)
  }
}

END_SECTION

END_TEST
