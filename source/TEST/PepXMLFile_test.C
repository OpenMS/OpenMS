// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
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
PepXMLFile file;
START_SECTION(PepXMLFile())
	ptr = new PepXMLFile();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION(~PepXMLFile())
	delete ptr;
END_SECTION

START_SECTION(void load(const String &filename, ProteinIdentification &protein, std::vector< PeptideIdentification > &peptides, const String &experiment_name, MSExperiment<> &experiment))
{
	ProteinIdentification protein, protein2;
	std::vector<PeptideIdentification> peptides, peptides2;
	String pep_file = OPENMS_GET_TEST_DATA_PATH("PepXMLFile_test.pepxml");
	String mz_file = OPENMS_GET_TEST_DATA_PATH("PepXMLFile_test.mzML");
	String exp_name = "B09-05181_p";
	MSExperiment<> experiment;
	MzMLFile().load(mz_file, experiment);
	file.load(pep_file, protein, peptides, exp_name, experiment);
	PeptideIdentification first = peptides[0];
	TEST_EQUAL(first.getMetaValue("RT"), 0.5927);
	TEST_EQUAL(first.getMetaValue("MZ"), 538.605);
	// check that only RT and m/z changes compared to the other "load" method:
	file.load(pep_file, protein2, peptides2);
	TEST_EQUAL(peptides.size(), peptides2.size());
	for (Size i = 0; i < peptides.size(); ++i)
	{
		peptides[i].clearMetaInfo();
		peptides2[i].clearMetaInfo();
	}
	TEST_EQUAL(peptides == peptides2, true);
	TEST_EQUAL(protein == protein2, true);
}
END_SECTION

START_SECTION(void load(const String& filename, ProteinIdentification& protein, std::vector<PeptideIdentification>& peptides, const String& experiment_name = ""))
{
	ProteinIdentification protein;
	std::vector<PeptideIdentification> peptides;
	String filename = OPENMS_GET_TEST_DATA_PATH("PepXMLFile_test.pepxml");
	String exp_name = "B09-05181_p";
	file.load(filename, protein, peptides, exp_name);

	// peptide IDs:
	TEST_EQUAL(peptides.size(), 9);
	PeptideIdentification first = peptides[0], last = peptides[8];
	TEST_EQUAL(first.getMetaValue("RT"), 1.3653); // RT of MS2 spectrum
	TEST_REAL_SIMILAR(first.getMetaValue("MZ"), 538.605); // recomputed
	TEST_EQUAL(first.getHits().size(), 1);
	PeptideHit pep_hit = first.getHits()[0];
	TEST_EQUAL(pep_hit.getSequence().toString(), "ELNKEMAAEKAKAAAG");
	TEST_EQUAL(pep_hit.getRank(), 1);
	// no use checking score, because implementation may still change
	TEST_EQUAL(pep_hit.getCharge(), 3);
	TEST_EQUAL(pep_hit.getProteinAccessions().size(), 3);
	TEST_EQUAL(pep_hit.getProteinAccessions()[0], "ddb000449223");
	TEST_EQUAL(pep_hit.getAABefore(), 'R');
	TEST_EQUAL(pep_hit.getAAAfter(), 'E');
	TEST_EQUAL(first.getIdentifier(), last.getIdentifier());
	// TODO: check handling of modifications on 1st/2nd/6th peptide ID

	// protein ID:
	TEST_EQUAL(protein.getIdentifier(), first.getIdentifier());
	TEST_NOT_EQUAL(protein.getIdentifier(), "");
	TEST_EQUAL(protein.getSearchEngine(), "X! Tandem (k-score)");
	std::vector<ProteinHit> prot_hits = protein.getHits();
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
	ProteinIdentification::SearchParameters params = protein.getSearchParameters();
	TEST_EQUAL(params.db, "./current.fasta");
	TEST_EQUAL(params.mass_type, ProteinIdentification::MONOISOTOPIC);
	TEST_EQUAL(params.enzyme, ProteinIdentification::TRYPSIN);
	// TODO: check handling of modification info

	// with the wrong "experiment_name", there are no results:
	file.load(filename, protein, peptides, "test");
	TEST_EQUAL(protein.getSearchParameters() == ProteinIdentification::SearchParameters(), true);
	TEST_EQUAL(protein.getHits().empty(), true);
	TEST_EQUAL(peptides.empty(), true);

	// throw an exception if the pepXML file does not exist:
	NEW_TMP_FILE(filename);
	TEST_EXCEPTION(Exception::FileNotFound,
								 file.load(filename, protein, peptides, exp_name));

}
END_SECTION

START_SECTION(void store(const String& filename, std::vector<ProteinIdentification>& protein_ids, std::vector<PeptideIdentification>& peptide_ids))
{
	PepXMLFile file;
	ProteinIdentification protein;
	std::vector<PeptideIdentification> peptides;
	String filename = OPENMS_GET_TEST_DATA_PATH("PepXMLFile_test.pepxml");
	file.load(filename, protein, peptides);
	
	std::vector<ProteinIdentification> protein_ids;
	protein_ids.push_back(protein);
	
	String cm_file_out;
	NEW_TMP_FILE(cm_file_out);
	String filename_out = OPENMS_GET_TEST_DATA_PATH("PepXMLFile_test_out.pepxml");
	file.store(cm_file_out, protein_ids, peptides);
	
	FuzzyStringComparator fsc;
	fsc.setWhitelist (StringList::create("base_name"));
	TEST_EQUAL(fsc.compareFiles (cm_file_out.c_str(), filename_out.c_str()), true)

}	
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

