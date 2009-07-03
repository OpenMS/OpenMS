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

START_SECTION(void load(const String& filename, ProteinIdentification& protein, std::vector<PeptideIdentification>& peptides, const String& experiment_name))
//TODO Hendrik
ProteinIdentification protein;
std::vector<PeptideIdentification> peptides;
String filename = OPENMS_GET_TEST_DATA_PATH("PepXMLFile_test.pepxml"),
	exp_name = "B08-08318";
file.load(filename, protein, peptides, exp_name);

// peptide IDs:
TEST_EQUAL(peptides.size(), 3);
PeptideIdentification first = peptides[0], last = peptides[2];
TEST_REAL_SIMILAR(first.getMetaValue("RT"), 7.6514);
TEST_REAL_SIMILAR(first.getMetaValue("MZ"), 609.835525);
TEST_EQUAL(first.getHits().size(), 1);
PeptideHit pep_hit = first.getHits()[0];
TEST_EQUAL(pep_hit.getSequence().toString(), "KRGYVPAEGNK");
TEST_EQUAL(pep_hit.getRank(), 1);
// no use checking score, because implementation may still change
TEST_EQUAL(pep_hit.getCharge(), 2);
TEST_EQUAL(pep_hit.getProteinAccessions().size(), 1);
TEST_EQUAL(pep_hit.getProteinAccessions()[0], "ddb000324841");
TEST_EQUAL(pep_hit.getAABefore(), 'K');
TEST_EQUAL(pep_hit.getAAAfter(), 'K');
TEST_EQUAL(first.getIdentifier(), last.getIdentifier());
pep_hit = last.getHits()[0];
TEST_EQUAL(pep_hit.getProteinAccessions().size(), 2);
// TODO: check handling of modifications on 2nd peptide ID ("peptides[1]")

// protein ID:
TEST_EQUAL(protein.getIdentifier(), first.getIdentifier());
TEST_NOT_EQUAL(protein.getIdentifier(), "");
TEST_EQUAL(protein.getSearchEngine(), "X! Tandem");
std::vector<ProteinHit> prot_hits = protein.getHits();
TEST_EQUAL(prot_hits.size(), 4);
StringList accessions;
for (std::vector<ProteinHit>::iterator it = prot_hits.begin();
		 it != prot_hits.end(); ++it)
{
	accessions << it->getAccession();
}
TEST_EQUAL(accessions.contains("ddb000324841"), true);
TEST_EQUAL(accessions.contains("ddb000014895"), true);
TEST_EQUAL(accessions.contains("ddb000000416"), true);
TEST_EQUAL(accessions.contains("ddb000630664"), true);

// search parameters:
ProteinIdentification::SearchParameters params = protein.getSearchParameters();
TEST_EQUAL(params.db, "./current.fasta");
TEST_EQUAL(params.mass_type, ProteinIdentification::MONOISOTOPIC);
TEST_EQUAL(params.enzyme, ProteinIdentification::TRYPSIN);
// TODO: check handling of modification info

END_SECTION

START_SECTION(void store(const String& filename, std::vector<ProteinIdentification>& protein_ids, std::vector<PeptideIdentification>& peptide_ids))
	//TODO Chris
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

