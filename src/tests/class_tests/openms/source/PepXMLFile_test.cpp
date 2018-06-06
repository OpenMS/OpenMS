// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow, Hendrik Weisser $
// $Authors: Chris Bielow, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/PepXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h> //ONLY used for checking if pepxml transformation produced a reusable id file
///////////////////////////

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/FuzzyStringComparator.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>

using namespace OpenMS;
using namespace std;

START_TEST(PepXMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PepXMLFile * ptr = nullptr;
PepXMLFile* nullPointer = nullptr;
PepXMLFile file;
START_SECTION(PepXMLFile())
ptr = new PepXMLFile();
TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~PepXMLFile())
delete ptr;
END_SECTION

START_SECTION(void load(const String& filename, std::vector<ProteinIdentification>& proteins, std::vector<PeptideIdentification>& peptides, const String& experiment_name, SpectrumMetaDataLookup& lookup))
{
  vector<ProteinIdentification> proteins;
  vector<PeptideIdentification> peptides;
  String pep_file = OPENMS_GET_TEST_DATA_PATH("PepXMLFile_test.pepxml");
  String mz_file = OPENMS_GET_TEST_DATA_PATH("PepXMLFile_test.mzML");
  String exp_name = "PepXMLFile_test";
  PeakMap experiment;
  MzMLFile().load(mz_file, experiment);
  SpectrumMetaDataLookup lookup;
  lookup.readSpectra(experiment.getSpectra());
  file.load(pep_file, proteins, peptides, exp_name, lookup);
  TEST_EQUAL(peptides.size(), 18);
  TEST_EQUAL(proteins.size(), 2);
  PeptideIdentification first = peptides[0];
  TEST_REAL_SIMILAR(first.getRT(), 1.3653);
  TEST_REAL_SIMILAR(first.getMZ(), 538.605);
  // more checks below
}

END_SECTION

START_SECTION(void load(const String& filename, std::vector<ProteinIdentification>& proteins, std::vector<PeptideIdentification>& peptides, const String& experiment_name = ""))
{
  vector<ProteinIdentification> proteins;
  vector<PeptideIdentification> peptides;
  // file contains results from two search runs:
  String filename = OPENMS_GET_TEST_DATA_PATH("PepXMLFile_test.pepxml");
  String exp_name = "PepXMLFile_test";
  file.load(filename, proteins, peptides, exp_name);

  // peptide IDs:
  TEST_EQUAL(peptides.size(), 18);
  PeptideIdentification first = peptides.front(), last = peptides.back();

  bool accu_result = true;   // to avoid spamming TEST_EQUAL's in the "for" loop
  for (Size i = 1; i < 9; ++i) // should be the same for all peptides from the first search run:
  {
    accu_result &= (first.getIdentifier() == peptides[i].getIdentifier());
    accu_result &= (first.getScoreType() == peptides[i].getScoreType());
    accu_result &= (first.isHigherScoreBetter() ==
                    peptides[i].isHigherScoreBetter());
    accu_result &= (first.getSignificanceThreshold() ==
                    peptides[i].getSignificanceThreshold());
  }
  TEST_EQUAL(accu_result, true);

  TEST_REAL_SIMILAR(first.getRT(), 1.3653); // RT of MS2 spectrum
  TEST_REAL_SIMILAR(first.getMZ(), 538.605); // recomputed
  TEST_EQUAL(first.getHits().size(), 1);
  PeptideHit pep_hit = first.getHits()[0];
  TEST_EQUAL(pep_hit.getSequence().toString(), ".(Glu->pyro-Glu)ELNKEMAAEKAKAAAG");
  TEST_EQUAL(pep_hit.getSequence().toUnmodifiedString(), "ELNKEMAAEKAKAAAG");
  TEST_EQUAL(pep_hit.getRank(), 1);

  // no use checking score, because implementation may still change
  TEST_EQUAL(pep_hit.getCharge(), 3);
  vector<PeptideEvidence> pes = pep_hit.getPeptideEvidences();
  TEST_EQUAL(pes.size(), 3);
  TEST_EQUAL(pes[0].getProteinAccession(), "ddb000449223");
  TEST_EQUAL(pes[0].getAABefore(), 'R');
  TEST_EQUAL(pes[0].getAAAfter(), 'E');

  TEST_EQUAL(first.getHits()[0].getSequence().isModified(), true);
  TEST_EQUAL(first.getHits()[0].getSequence().hasNTerminalModification(), true);
  TEST_EQUAL(first.getHits()[0].getSequence().hasCTerminalModification(), false);

  TEST_EQUAL(peptides[1].getHits()[0].getSequence().isModified(), true);
  TEST_EQUAL(peptides[1].getHits()[0].getSequence().hasNTerminalModification(), true);
  TEST_EQUAL(peptides[1].getHits()[0].getSequence().hasCTerminalModification(), false);

  TEST_EQUAL(peptides[5].getHits()[0].getSequence().isModified(), true);
  TEST_EQUAL(peptides[5].getHits()[0].getSequence().hasNTerminalModification(), false);
  TEST_EQUAL(peptides[5].getHits()[0].getSequence().hasCTerminalModification(), false);

  // cursory check of a peptide ID from the second search run:
  pep_hit = last.getHits()[0];
  TEST_EQUAL(pep_hit.getSequence().toString(), "EISPDTTLLDLQNNDISELR");

  // protein ID:
  TEST_EQUAL(proteins.size(), 2);
  TEST_EQUAL(proteins[0].getIdentifier(), first.getIdentifier());
  TEST_EQUAL(proteins[1].getIdentifier(), last.getIdentifier());
  TEST_NOT_EQUAL(proteins[0].getIdentifier(), "");
  TEST_NOT_EQUAL(proteins[1].getIdentifier(), "");
  TEST_NOT_EQUAL(proteins[0].getIdentifier(), proteins[1].getIdentifier());
  TEST_EQUAL(proteins[0].getSearchEngine(), "X! Tandem (k-score)");
  TEST_EQUAL(proteins[1].getSearchEngine(), "SEQUEST");

  vector<ProteinHit> prot_hits = proteins[0].getHits();
  TEST_EQUAL(prot_hits.size(), 20);
  StringList accessions_string;
  for (std::vector<ProteinHit>::iterator it = prot_hits.begin();
       it != prot_hits.end(); ++it)
  {
    accessions_string << it->getAccession();
  }
  // check a sample of the IDs that should be present:
  TEST_EQUAL(ListUtils::contains(accessions_string, "ddb000449223"), true);
  TEST_EQUAL(ListUtils::contains(accessions_string, "ddb000626346"), true);
  TEST_EQUAL(ListUtils::contains(accessions_string, "rev000409159"), true);

  // search parameters:
  ProteinIdentification::SearchParameters params = proteins[0].getSearchParameters();
  TEST_EQUAL(params.db, "./current.fasta");
  TEST_EQUAL(params.mass_type, ProteinIdentification::MONOISOTOPIC);
  TEST_EQUAL(params.digestion_enzyme.getName(), "Trypsin");

  vector<String> fix_mods(params.fixed_modifications), var_mods(params.variable_modifications);
  TEST_EQUAL(fix_mods.size(), 1)
  TEST_EQUAL(var_mods.size(), 5)

  TEST_EQUAL(find(var_mods.begin(), var_mods.end(), "Ammonia-loss (N-term C)") != var_mods.end(), true)
  TEST_EQUAL(find(var_mods.begin(), var_mods.end(), "Glu->pyro-Glu (N-term E)") != var_mods.end(), true)
  TEST_EQUAL(find(var_mods.begin(), var_mods.end(), "Oxidation (M)") != var_mods.end(), true)
  TEST_EQUAL(find(var_mods.begin(), var_mods.end(), "Gln->pyro-Glu (N-term Q)") != var_mods.end(), true)
  TEST_EQUAL(find(var_mods.begin(), var_mods.end(), "M+1") != var_mods.end(), true)

  // wrong "experiment_name" produces an exception:
  TEST_EXCEPTION(Exception::ParseError, file.load(filename, proteins, peptides, "abcxyz"));

  // throw an exception if the pepXML file does not exist:
  TEST_EXCEPTION(Exception::FileNotFound, file.load("this_file_does_not_exist_but_should_be_a_pepXML_file.pepXML", proteins, peptides, exp_name));
}
END_SECTION

START_SECTION([EXTRA] void load(const String& filename, std::vector<ProteinIdentification>& proteins, std::vector<PeptideIdentification>& peptides, const String& experiment_name = ""))
{
  vector<ProteinIdentification> proteins;
  vector<PeptideIdentification> peptides;
  // file contains results from two search runs:
  String filename = OPENMS_GET_TEST_DATA_PATH("PepXMLFile_test_extended.pepxml");
  String exp_name = "PepXMLFile_test";
  file.keepNativeSpectrumName(true);
  file.load(filename, proteins, peptides, exp_name);

  // peptide IDs:
  TEST_EQUAL(peptides.size(), 2);
  PeptideIdentification first = peptides.front(), last = peptides.back();

  TEST_EQUAL(first.getRT(), 1.3653);   // RT of MS2 spectrum
  TEST_REAL_SIMILAR(first.getMZ(), 538.605);   // recomputed
  TEST_EQUAL(first.getHits().size(), 1);

  TEST_EQUAL(last.getRT(), 488.652);   // RT of MS2 spectrum
  TEST_REAL_SIMILAR(last.getMZ(), 585.3166250319);   // recomputed
  TEST_EQUAL(last.getHits().size(), 1);
  TEST_EQUAL(last.metaValueExists("swath_assay"), true);
  TEST_EQUAL(last.metaValueExists("status"), true);
  TEST_EQUAL(last.metaValueExists("pepxml_spectrum_name"), true);
  TEST_EQUAL(last.getExperimentLabel().empty(), false);

  TEST_EQUAL(last.getMetaValue("swath_assay"), "EIVLTQSPGTL2:9");
  TEST_EQUAL(last.getMetaValue("status"), "target");
  TEST_EQUAL(last.getMetaValue("pepxml_spectrum_name"), "hroest_K120718_SM_OGE10_010_IDA.02552.02552.2");
  TEST_EQUAL(last.getExperimentLabel(), "urine");

  PeptideHit pep_hit = last.getHits()[0];
  TEST_EQUAL(pep_hit.getSequence().toString(), "VVITAPGGNDVK");
  TEST_EQUAL(pep_hit.getSequence().toUnmodifiedString(), "VVITAPGGNDVK");
  TEST_EQUAL(pep_hit.getRank(), 1);
  TEST_EQUAL(pep_hit.getCharge(), 2);

  // check the analysis scores
  TEST_EQUAL(pep_hit.getAnalysisResults().size(), 2);

  PeptideHit::PepXMLAnalysisResult a = pep_hit.getAnalysisResults()[0];
  TEST_EQUAL(a.score_type, "peptideprophet");
  TEST_REAL_SIMILAR(a.main_score, 0.0660);

  TEST_EQUAL(a.sub_scores.find("fval") != a.sub_scores.end(), true);
  TEST_EQUAL(a.sub_scores.find("ntt") != a.sub_scores.end(), true);
  TEST_EQUAL(a.sub_scores.find("empir_irt") != a.sub_scores.end(), true);
  TEST_EQUAL(a.sub_scores.find("swath_window") != a.sub_scores.end(), true);

  TEST_REAL_SIMILAR(a.sub_scores.find("fval")->second, 0.7114);
  TEST_REAL_SIMILAR(a.sub_scores.find("ntt")->second, 2);
  TEST_REAL_SIMILAR(a.sub_scores.find("empir_irt")->second, 79.79);
  TEST_REAL_SIMILAR(a.sub_scores.find("swath_window")->second, 9);

  // <analysis_result analysis="peptideprophet">
  //   <peptideprophet_result probability="0.0660" all_ntt_prob="(0.0000,0.0000,0.0660)">
  //     <search_score_summary>
  //       <parameter name="fval" value="0.7114"/>
  //       <parameter name="ntt" value="2"/>
  //       <parameter name="nmc" value="0"/>
  //       <parameter name="massd" value="-0.027"/>
  //       <parameter name="isomassd" value="0"/>

  //       <parameter name="empir_irt" value="79.79"/>
  //       <parameter name="empir_irt_bin" value="53"/>
  //       <parameter name="swath_window" value="9"/>
  //       <parameter name="alt_swath" value="-1"/>

  //     </search_score_summary>
  //   </peptideprophet_result>
  // </analysis_result>

  // <analysis_result analysis="interprophet">
  //   <interprophet_result probability="0.93814" all_ntt_prob="(0,0,0.93814)">
  //     <search_score_summary>
  //       <parameter name="nss" value="0"/>
  //       <parameter name="nrs" value="10.2137"/>
  //       <parameter name="nse" value="0"/>
  //       <parameter name="nsi" value="0.9793"/>
  //       <parameter name="nsm" value="0"/>
  //     </search_score_summary>
  //   </interprophet_result>
  // </analysis_result>

  a = pep_hit.getAnalysisResults()[1];
  TEST_EQUAL(a.score_type, "interprophet");
  TEST_REAL_SIMILAR(a.main_score, 0.93814);

  TEST_EQUAL(a.sub_scores.find("fval") == a.sub_scores.end(), true);
  TEST_EQUAL(a.sub_scores.find("nss") != a.sub_scores.end(), true);
  TEST_REAL_SIMILAR(a.sub_scores.find("nrs")->second, 10.2137);

  // wrong "experiment_name" produces an exception:
  TEST_EXCEPTION(Exception::ParseError, file.load(filename, proteins, peptides, "abcxyz"));

  // throw an exception if the pepXML file does not exist:
  TEST_EXCEPTION(Exception::FileNotFound, file.load("this_file_does_not_exist_but_should_be_a_pepXML_file.pepXML", proteins, peptides, exp_name));
}
END_SECTION

START_SECTION(void store(const String& filename, std::vector<ProteinIdentification>& protein_ids, std::vector<PeptideIdentification>& peptide_ids, const String& mz_file = "", const String& mz_name = "", bool peptideprophet_analyzed = false))
{
  vector<ProteinIdentification> proteins;
  vector<PeptideIdentification> peptides;
  String filename = OPENMS_GET_TEST_DATA_PATH("PepXMLFile_test_store.pepxml");
  PepXMLFile().load(filename, proteins, peptides);

  // Test PeptideProphet-analyzed pepxml.
  String cm_file_out;
  NEW_TMP_FILE(cm_file_out);
  PepXMLFile().store(cm_file_out, proteins, peptides, "", "test", true);

  FuzzyStringComparator fsc;
  fsc.setAcceptableAbsolute(1e-7);
  fsc.setAcceptableRelative(1.0 + 1e-7);
  // fsc.setWhitelist (ListUtils::create<String>("base_name, local_path, <spectrum_query "));
  String filename_out = OPENMS_GET_TEST_DATA_PATH("PepXMLFile_test_out.pepxml");
  TEST_EQUAL(fsc.compareFiles(cm_file_out.c_str(), filename_out.c_str()), true)

  // Test raw_pepxml storage.
  String cm_file_out_1;
  NEW_TMP_FILE(cm_file_out_1);
  PepXMLFile().store(cm_file_out_1, proteins, peptides, "", "test", false);

  FuzzyStringComparator fsc_1;
  fsc_1.setAcceptableAbsolute(1e-7);
  fsc_1.setAcceptableRelative(1.0 + 1e-7);
  // fsc_1.setWhitelist(ListUtils::create<String>("base_name, local_path, <spectrum_query "));
  String filename_out_1 = OPENMS_GET_TEST_DATA_PATH("PepXMLFile_test_out_1.pepxml");
  TEST_EQUAL(fsc_1.compareFiles(cm_file_out_1.c_str(), filename_out_1.c_str()), true)
}
END_SECTION

START_SECTION([EXTRA] void store(const String& filename, std::vector<ProteinIdentification>& protein_ids, std::vector<PeptideIdentification>& peptide_ids, const String& mz_file = "", const String& mz_name = "", bool peptideprophet_analyzed = false))
{
  {  
    vector<ProteinIdentification> proteins;
    vector<PeptideIdentification> peptides;
    // file contains results from two search runs:
    String filename = OPENMS_GET_TEST_DATA_PATH("PepXMLFile_test_extended.pepxml");
    String exp_name = "PepXMLFile_test";
    PepXMLFile file;
    file.keepNativeSpectrumName(true);
    file.load(filename, proteins, peptides, exp_name);

    TEST_EQUAL(peptides.size(), 2);
    PeptideIdentification first = peptides.front(), last = peptides.back();
    TEST_REAL_SIMILAR(first.getMZ(), 538.605);   // recomputed
    TEST_REAL_SIMILAR(last.getMZ(), 585.3166250319);   // recomputed

    // Now try to store the file again ... 
    String cm_file_out;
    NEW_TMP_FILE(cm_file_out);
    file.store(cm_file_out, proteins, peptides, "", exp_name, false); // peptideprophet_analyzed = false is important!

    // And read it back in again
    vector<ProteinIdentification> proteins_new;
    vector<PeptideIdentification> peptides_new;
    file.load(cm_file_out, proteins_new, peptides_new, exp_name);

    TEST_EQUAL(proteins.size(), proteins_new.size())
    TEST_EQUAL(peptides.size(), peptides_new.size())

    // peptide IDs:
    TEST_EQUAL(peptides_new.size(), 2);
    first = peptides_new.front(); last = peptides_new.back();

    TEST_EQUAL(first.getRT(), 1.3653);   // RT of MS2 spectrum
    TEST_REAL_SIMILAR(first.getMZ(), 538.6159248633);   // recomputed
    TEST_EQUAL(first.getHits().size(), 1);

    TEST_EQUAL(last.getRT(), 488.652);   // RT of MS2 spectrum
    TEST_REAL_SIMILAR(last.getMZ(), 585.3304219355);   // recomputed
    TEST_EQUAL(last.getHits().size(), 1);
    PeptideHit pep_hit = last.getHits()[0];
    TEST_EQUAL(pep_hit.getSequence().toString(), "VVITAPGGNDVK");
    TEST_EQUAL(pep_hit.getSequence().toUnmodifiedString(), "VVITAPGGNDVK");
    TEST_EQUAL(pep_hit.getRank(), 1);
    TEST_EQUAL(pep_hit.getCharge(), 2);

    // test extra attributes (correctly read and written)
    TEST_EQUAL(last.metaValueExists("swath_assay"), true);
    TEST_EQUAL(last.metaValueExists("status"), true);
    TEST_EQUAL(last.metaValueExists("pepxml_spectrum_name"), true);
    TEST_EQUAL(last.getExperimentLabel().empty(), false);

    TEST_EQUAL(last.getMetaValue("swath_assay"), "EIVLTQSPGTL2:9");
    TEST_EQUAL(last.getMetaValue("status"), "target");
    TEST_EQUAL(last.getMetaValue("pepxml_spectrum_name"), "hroest_K120718_SM_OGE10_010_IDA.02552.02552.22");
    TEST_EQUAL(last.getMetaValue("pepxml_spectrum_name") == "hroest_K120718_SM_OGE10_010_IDA.02552.02552.22", true);
    TEST_EQUAL(last.getExperimentLabel(), "urine");

    // check the analysis scores
    TEST_EQUAL(pep_hit.getAnalysisResults().size(), 2);

    PeptideHit::PepXMLAnalysisResult a = pep_hit.getAnalysisResults()[0];
    TEST_EQUAL(a.score_type, "peptideprophet");
    TEST_REAL_SIMILAR(a.main_score, 0.0660);

    TEST_EQUAL(a.sub_scores.find("fval") != a.sub_scores.end(), true);
    TEST_EQUAL(a.sub_scores.find("ntt") != a.sub_scores.end(), true);
    TEST_EQUAL(a.sub_scores.find("empir_irt") != a.sub_scores.end(), true);
    TEST_EQUAL(a.sub_scores.find("swath_window") != a.sub_scores.end(), true);

    TEST_REAL_SIMILAR(a.sub_scores.find("fval")->second, 0.7114);
    TEST_REAL_SIMILAR(a.sub_scores.find("ntt")->second, 2);
    TEST_REAL_SIMILAR(a.sub_scores.find("empir_irt")->second, 79.79);
    TEST_REAL_SIMILAR(a.sub_scores.find("swath_window")->second, 9);

  }

  // test keep native spectrum name = false
  {  
    vector<ProteinIdentification> proteins;
    vector<PeptideIdentification> peptides;
    String filename = OPENMS_GET_TEST_DATA_PATH("PepXMLFile_test_extended.pepxml");
    String exp_name = "PepXMLFile_test";
    PepXMLFile file;
    file.keepNativeSpectrumName(false);
    file.load(filename, proteins, peptides, exp_name);

    // Now try to store the file again ... 
    String cm_file_out;
    NEW_TMP_FILE(cm_file_out);
    file.store(cm_file_out, proteins, peptides, "", exp_name, false); // peptideprophet_analyzed = false is important!

    // And read it back in again
    vector<ProteinIdentification> proteins_new;
    vector<PeptideIdentification> peptides_new;
    file.load(cm_file_out, proteins_new, peptides_new, exp_name);

    // peptide IDs:
    PeptideIdentification last = peptides.back();

    // now this should be fales 
    TEST_EQUAL(last.getMetaValue("pepxml_spectrum_name") != "hroest_K120718_SM_OGE10_010_IDA.02552.02552.22", true);
 }
}
END_SECTION

// store PepXML with mzML file information
START_SECTION(void store(const String& filename, std::vector<ProteinIdentification>& protein_ids, std::vector<PeptideIdentification>& peptide_ids, const String& mz_file = "PepXMLFile_test.mzML", const String& mz_name = "", bool peptideprophet_analyzed = false))
{
  vector<ProteinIdentification> proteins;
  vector<PeptideIdentification> peptides;
  String mzML_filename = OPENMS_GET_TEST_DATA_PATH("PepXMLFile_test.mzML");
  String filename = OPENMS_GET_TEST_DATA_PATH("PepXMLFile_test_store.pepxml");
  PepXMLFile().load(filename, proteins, peptides);

  // Test PeptideProphet-analyzed pepxml.
  String cm_file_out;
  NEW_TMP_FILE(cm_file_out);
  PepXMLFile().store(cm_file_out, proteins, peptides, mzML_filename, "test", true);

  FuzzyStringComparator fsc;
  fsc.setAcceptableAbsolute(1e-7);
  fsc.setAcceptableRelative(1.0 + 1e-7);
  // fsc.setWhitelist (ListUtils::create<String>("base_name, local_path, <spectrum_query "));
  String filename_out = OPENMS_GET_TEST_DATA_PATH("PepXMLFile_test_out_mzML.pepxml");
  TEST_EQUAL(fsc.compareFiles(cm_file_out.c_str(), filename_out.c_str()), true)
}
END_SECTION

START_SECTION(void keepNativeSpectrumName(bool keep) )
{
  // tested above in the [EXTRA] store as we store / load once with
  // keepNativeSpectrumName and once without
  NOT_TESTABLE 
}
END_SECTION


START_SECTION(([EXTRA] checking pepxml transformation to reusable identifications))

  // PepXMLFile file; // shadow
  vector<ProteinIdentification> proteins, reread_proteins;
  vector<PeptideIdentification> peptides, reread_peptides;
  String filename = OPENMS_GET_TEST_DATA_PATH("PepXMLFile_test_store.pepxml");
  PepXMLFile().load(filename, proteins, peptides);

  // Test PeptideProphet-analyzed pepxml.
  String cm_file_out;
  NEW_TMP_FILE(cm_file_out);
  IdXMLFile().store(cm_file_out, proteins, peptides);
  IdXMLFile().load(cm_file_out, reread_proteins, reread_peptides);

  ProteinIdentification::SearchParameters params = proteins[0].getSearchParameters();
  ProteinIdentification::SearchParameters reread_params = reread_proteins[0].getSearchParameters();
  TEST_EQUAL(params.db, reread_params.db);
  TEST_EQUAL(params.mass_type, reread_params.mass_type);

  vector<String> fix_mods(params.fixed_modifications), var_mods(params.variable_modifications);
  vector<String> reread_fix_mods(reread_params.fixed_modifications), reread_var_mods(reread_params.variable_modifications);
  TEST_EQUAL(fix_mods.size(), reread_fix_mods.size())
  TEST_EQUAL(var_mods.size(), reread_var_mods.size())

  TEST_EQUAL(find(fix_mods.begin(), fix_mods.end(), reread_fix_mods[0]) != var_mods.end(), true)
  TEST_EQUAL(find(fix_mods.begin(), fix_mods.end(), "Carbamidometyhl (C)") != var_mods.end(), true)

  for (size_t i = 0; i < reread_var_mods.size(); ++i)
  {
    TEST_EQUAL(find(var_mods.begin(), var_mods.end(), reread_var_mods[i]) != var_mods.end(), true)
  }

END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
