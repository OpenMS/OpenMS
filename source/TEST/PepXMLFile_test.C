// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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

///////////////////////////
#include <OpenMS/FORMAT/PepXMLFile.h>
///////////////////////////
#include <OpenMS/CONCEPT/FuzzyStringComparator.h>
#include <OpenMS/FORMAT/MzMLFile.h>

using namespace OpenMS;
using namespace std;

START_TEST(PepXMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PepXMLFile* ptr = 0;
PepXMLFile* nullPointer = 0;
PepXMLFile file;
START_SECTION(PepXMLFile())
	ptr = new PepXMLFile();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~PepXMLFile())
	delete ptr;
END_SECTION

START_SECTION(void load(const String &filename, std::vector<ProteinIdentification>& proteins, std::vector<PeptideIdentification>& peptides, const String& experiment_name, const MSExperiment<>& experiment, bool use_precursor_data = false))
{
	vector<ProteinIdentification> proteins, proteins2;
	vector<PeptideIdentification> peptides, peptides2;
	String pep_file = OPENMS_GET_TEST_DATA_PATH("PepXMLFile_test.pepxml");
	String mz_file = OPENMS_GET_TEST_DATA_PATH("PepXMLFile_test.mzML");
	String exp_name = "PepXMLFile_test";
	MSExperiment<> experiment;
	MzMLFile().load(mz_file, experiment);
	// load with precursor information:
	file.load(pep_file, proteins, peptides, exp_name, experiment, true);
	PeptideIdentification first = peptides[0];
	TEST_EQUAL(first.getMetaValue("RT"), 0.5927);
	TEST_EQUAL(first.getMetaValue("MZ"), 538.605);
	// load without precursor information (tested more thoroughly below):
	file.load(pep_file, proteins2, peptides2);
	TEST_EQUAL(peptides.size(), peptides2.size());
	TEST_EQUAL(proteins.size(), proteins2.size());
	// check that only m/z and RT differ between "load" methods:
	for (Size i = 0; i < peptides.size(); ++i)
	{
		peptides[i].clearMetaInfo();
		peptides2[i].clearMetaInfo();
	}
	TEST_EQUAL(peptides == peptides2, true);
	TEST_EQUAL(proteins == proteins2, true);
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

	bool accu_result = true; // to avoid spamming TEST_EQUAL's in the "for" loop
	for (Size i = 1; i < 9; ++i)
	{ // should be the same for all peptides from the first search run:
		accu_result &= (first.getIdentifier() == peptides[i].getIdentifier());
		accu_result &= (first.getScoreType() == peptides[i].getScoreType());
		accu_result &= (first.isHigherScoreBetter() == 
										peptides[i].isHigherScoreBetter());
		accu_result &= (first.getSignificanceThreshold() ==
										peptides[i].getSignificanceThreshold());
	}
	TEST_EQUAL(accu_result, true);

	TEST_EQUAL(first.getMetaValue("RT"), 1.3653); // RT of MS2 spectrum
	TEST_REAL_SIMILAR(first.getMetaValue("MZ"), 538.605); // recomputed
	TEST_EQUAL(first.getHits().size(), 1);
	PeptideHit pep_hit = first.getHits()[0];
	TEST_EQUAL(pep_hit.getSequence().toString(), "(Glu->pyro-Glu)ELNKEMAAEKAKAAAG");
	TEST_EQUAL(pep_hit.getSequence().toUnmodifiedString(), "ELNKEMAAEKAKAAAG");
	TEST_EQUAL(pep_hit.getRank(), 1);

	// no use checking score, because implementation may still change
	TEST_EQUAL(pep_hit.getCharge(), 3);
	TEST_EQUAL(pep_hit.getProteinAccessions().size(), 3);
	TEST_EQUAL(pep_hit.getProteinAccessions()[0], "ddb000449223");
	TEST_EQUAL(pep_hit.getAABefore(), 'R');
	TEST_EQUAL(pep_hit.getAAAfter(), 'E');

	TEST_EQUAL(first.getHits()[0].getSequence().isModified(), true);
	TEST_EQUAL(first.getHits()[0].getSequence().hasNTerminalModification(), true);
	TEST_EQUAL(first.getHits()[0].getSequence().hasCTerminalModification(), 
false);

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
	StringList accessions;
	for (std::vector<ProteinHit>::iterator it = prot_hits.begin();
			 it != prot_hits.end(); ++it)
	{
		accessions << it->getAccession();
	}
	// check a sample of the IDs that should be present:
	TEST_EQUAL(accessions.contains("ddb000449223"), true);
	TEST_EQUAL(accessions.contains("ddb000626346"), true);
	TEST_EQUAL(accessions.contains("rev000409159"), true);

	// search parameters:
	ProteinIdentification::SearchParameters params = proteins[0].getSearchParameters();
	TEST_EQUAL(params.db, "./current.fasta");
	TEST_EQUAL(params.mass_type, ProteinIdentification::MONOISOTOPIC);
	TEST_EQUAL(params.enzyme, ProteinIdentification::TRYPSIN);

	vector<String> fix_mods(params.fixed_modifications), var_mods(params.variable_modifications);
	TEST_EQUAL(fix_mods.size(), 1)
	TEST_EQUAL(var_mods.size(), 4)

	TEST_EQUAL(find(var_mods.begin(), var_mods.end(), "C-17.0265") != var_mods.end(), true)
	TEST_EQUAL(find(var_mods.begin(), var_mods.end(), "E-18.0106") != var_mods.end(), true)
	TEST_EQUAL(find(var_mods.begin(), var_mods.end(), "M+15.9949") != var_mods.end(), true)
	TEST_EQUAL(find(var_mods.begin(), var_mods.end(), "Q-17.0265") != var_mods.end(), true)

	//TEST_EQUAL(find(var_mods.begin(), var_mods.end(), "Carbamidometyhl (C)") != var_mods.end(), true)	
	//TEST_EQUAL(find(var_mods.begin(), var_mods.end(), "Gln->pyro-Glu (Q)") != var_mods.end(), true)	
	//TEST_EQUAL(find(var_mods.begin(), var_mods.end(), "Glu->pyro-Glu (E)") != var_mods.end(), true)	

	// wrong "experiment_name" produces an exception:
	TEST_EXCEPTION(Exception::ParseError, file.load(filename, proteins, peptides, "abcxyz"));

	// throw an exception if the pepXML file does not exist:
	TEST_EXCEPTION(Exception::FileNotFound, file.load("this_file_does_not_exist_but_should_be_a_pepXML_file.pepXML", proteins, peptides, exp_name));
}
END_SECTION

START_SECTION(void store(const String& filename, std::vector<ProteinIdentification>& protein_ids, std::vector<PeptideIdentification>& peptide_ids))
{
	PepXMLFile file;
	vector<ProteinIdentification> proteins;
	vector<PeptideIdentification> peptides;
	String filename = OPENMS_GET_TEST_DATA_PATH("PepXMLFile_test_store.pepxml");
	file.load(filename, proteins, peptides);
	
	String cm_file_out;
	NEW_TMP_FILE(cm_file_out);
	file.store(cm_file_out, proteins, peptides);
	
	FuzzyStringComparator fsc;
	fsc.setWhitelist (StringList::create("base_name"));
	String filename_out = OPENMS_GET_TEST_DATA_PATH("PepXMLFile_test_out.pepxml");
	TEST_EQUAL(fsc.compareFiles (cm_file_out.c_str(), filename_out.c_str()), true)
}	
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

