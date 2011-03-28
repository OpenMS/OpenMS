// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Nico Pfeifer $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <string>

#include <OpenMS/FORMAT/MascotXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/DATASTRUCTURES/String.h>

///////////////////////////

START_TEST(PeptideIdentification, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;


DoubleReal peptide_significance_threshold = 42.3;
std::vector<PeptideHit> peptide_hits;
PeptideHit peptide_hit;
ProteinIdentification protein_identification;
vector<PeptideIdentification> identifications; 
MascotXMLFile xml_file;

peptide_hits.push_back(peptide_hit);


PeptideIdentification* ptr = 0;
PeptideIdentification* nullPointer = 0;
START_SECTION((PeptideIdentification()))
	ptr = new PeptideIdentification();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~PeptideIdentification()))
	PeptideIdentification hits;
	delete ptr;
END_SECTION

START_SECTION((PeptideIdentification(const PeptideIdentification& source)))
	PeptideIdentification hits;
	hits.setSignificanceThreshold(peptide_significance_threshold);
	hits.setHits(peptide_hits);
	hits.setMetaValue("label",17);
	hits.setIdentifier("id");
	hits.setScoreType("score_type");
	hits.setHigherScoreBetter(false);
	
	PeptideIdentification hits2(hits);
	
	TEST_EQUAL(hits.getSignificanceThreshold(), hits2.getSignificanceThreshold())
	TEST_EQUAL(hits.getHits().size() == 1, true)
	TEST_EQUAL(*(hits.getHits().begin()) == peptide_hit, true)
	TEST_EQUAL((UInt)hits.getMetaValue("label"),17)
	TEST_EQUAL(hits.getIdentifier(),"id")
	TEST_EQUAL(hits.getScoreType(),"score_type")
	TEST_EQUAL(hits.isHigherScoreBetter(),false)
END_SECTION

START_SECTION((PeptideIdentification& operator=(const PeptideIdentification& source)))
	PeptideIdentification hits;
	hits.setSignificanceThreshold(peptide_significance_threshold);
	hits.setHits(peptide_hits);
	hits.setMetaValue("label",17);
	hits.setIdentifier("id");
	hits.setScoreType("score_type");
	hits.setHigherScoreBetter(false);
	
	PeptideIdentification hits2;
	hits2 = hits;
	
	TEST_EQUAL(hits.getSignificanceThreshold(), hits2.getSignificanceThreshold())
	TEST_EQUAL(hits.getHits().size() == 1, true)
	TEST_EQUAL(*(hits.getHits().begin()) == peptide_hit, true)
	TEST_EQUAL((UInt)hits.getMetaValue("label"),17)
	TEST_EQUAL(hits.getIdentifier(),"id")
	TEST_EQUAL(hits.getScoreType(),"score_type")
	TEST_EQUAL(hits.isHigherScoreBetter(),false)
END_SECTION

START_SECTION((bool operator == (const PeptideIdentification& rhs) const))
	PeptideIdentification search1, search2;
	TEST_EQUAL(search1 == search2, true)
	
	search1.setSignificanceThreshold(peptide_significance_threshold);
	TEST_EQUAL(search1 == search2, false)	
	search1 = search2;

	search2.setMetaValue("label",17);
	TEST_EQUAL(search1 == search2, false)	
	search1 = search2;

	search2.setIdentifier("id");
	TEST_EQUAL(search1 == search2, false)	
	search1 = search2;

	search2.setScoreType("score_type");
	TEST_EQUAL(search1 == search2, false)	
	search1 = search2;

	search2.setHigherScoreBetter(false);
	TEST_EQUAL(search1 == search2, false)	
	search1 = search2;
END_SECTION


START_SECTION((bool operator != (const PeptideIdentification& rhs) const))
	PeptideIdentification search1, search2;
	TEST_EQUAL(search1 != search2, false)	

	search1.setSignificanceThreshold(peptide_significance_threshold);
	TEST_EQUAL(search1 != search2, true)	
	search1 = search2;
	
	//rest does not need to be tested, as it is tested in the operator== test implicitly!
END_SECTION


START_SECTION((DoubleReal getSignificanceThreshold() const))
	PeptideIdentification hits;
	hits.setSignificanceThreshold(peptide_significance_threshold);
	TEST_EQUAL(hits.getSignificanceThreshold(), peptide_significance_threshold)
END_SECTION

START_SECTION((const std::vector<PeptideHit>& getHits() const))
	PeptideIdentification hits;
	hits.insertHit(peptide_hit);
	TEST_EQUAL(hits.getHits().size() == 1, true)
	TEST_EQUAL(hits.getHits()[0] == peptide_hit, true)	
END_SECTION

START_SECTION((void insertHit(const PeptideHit &hit)))
	PeptideIdentification hits;
	hits.insertHit(peptide_hit);
	TEST_EQUAL(hits.getHits().size() == 1, true)
	TEST_EQUAL(*(hits.getHits().begin()) == peptide_hit, true)	
END_SECTION

START_SECTION((void setHits(const std::vector< PeptideHit > &hits)))
	PeptideIdentification hits;
	hits.setHits(peptide_hits);
	TEST_EQUAL(hits.getHits() == peptide_hits, true)
END_SECTION

START_SECTION((void setSignificanceThreshold(DoubleReal value)))
	PeptideIdentification hits;
	hits.setSignificanceThreshold(peptide_significance_threshold);
	TEST_EQUAL(hits.getSignificanceThreshold(), peptide_significance_threshold)
END_SECTION

START_SECTION((String getScoreType() const))
	PeptideIdentification hits;
	TEST_EQUAL(hits.getScoreType(),"")
END_SECTION

START_SECTION((void setScoreType(const String& type)))
	PeptideIdentification hits;
	hits.setScoreType("bla");
	TEST_EQUAL(hits.getScoreType(),"bla")
END_SECTION

START_SECTION((bool isHigherScoreBetter() const))
	PeptideIdentification hits;
	TEST_EQUAL(hits.isHigherScoreBetter(),true)
END_SECTION

START_SECTION((void setHigherScoreBetter(bool value)))
	PeptideIdentification hits;
	hits.setHigherScoreBetter(false);
	TEST_EQUAL(hits.isHigherScoreBetter(),false)
END_SECTION

START_SECTION((const String& getIdentifier() const))
	PeptideIdentification hits;
	TEST_EQUAL(hits.getIdentifier(),"")
END_SECTION

START_SECTION((void setIdentifier(const String& id)))
	PeptideIdentification hits;
	hits.setIdentifier("bla");
	TEST_EQUAL(hits.getIdentifier(),"bla")
END_SECTION

START_SECTION((bool empty() const))
	PeptideIdentification hits;
	TEST_EQUAL(hits.empty(), true)
	hits.setSignificanceThreshold(1);
	TEST_EQUAL(hits.empty(), false)
	hits.setSignificanceThreshold(0);
	hits.insertHit(peptide_hit);
	TEST_EQUAL(hits.empty(), false)
END_SECTION

START_SECTION((void sort()))
	PeptideIdentification id;
	PeptideHit hit;
	hit.setScore(23);
	hit.setSequence("SECONDPROTEIN");
	id.insertHit(hit);
	hit.setScore(45);
	hit.setSequence("FIRSTPROTEIN");
	id.insertHit(hit);
	hit.setScore(7);
	hit.setSequence("THIRDPROTEIN");
	id.insertHit(hit);
	
	//higher score is better
	id.sort();

	TEST_EQUAL(id.getHits()[0].getSequence(), "FIRSTPROTEIN")
	TEST_EQUAL(id.getHits()[1].getSequence(), "SECONDPROTEIN")
	TEST_EQUAL(id.getHits()[2].getSequence(), "THIRDPROTEIN")
	TEST_EQUAL(id.getHits()[0].getScore(), 45)	
	TEST_EQUAL(id.getHits()[1].getScore(), 23)
	TEST_EQUAL(id.getHits()[2].getScore(), 7)

	//lower score is better
	id.setHigherScoreBetter(false);
	id.sort();

	TEST_EQUAL(id.getHits()[0].getSequence(), "THIRDPROTEIN")
	TEST_EQUAL(id.getHits()[1].getSequence(), "SECONDPROTEIN")
	TEST_EQUAL(id.getHits()[2].getSequence(), "FIRSTPROTEIN")
	TEST_EQUAL(id.getHits()[0].getScore(), 7)	
	TEST_EQUAL(id.getHits()[1].getScore(), 23)
	TEST_EQUAL(id.getHits()[2].getScore(), 45)

END_SECTION

START_SECTION((void assignRanks()))
	PeptideIdentification id;
	PeptideHit hit;
	hit.setScore(23);
	hit.setSequence("SECONDPROTEIN");
	id.insertHit(hit);
	hit.setScore(45);
	hit.setSequence("FIRSTPROTEIN");
	id.insertHit(hit);
	hit.setScore(7);
	hit.setSequence("THIRDPROTEIN");
	id.insertHit(hit);

	id.assignRanks();

	TEST_EQUAL(id.getHits()[0].getSequence(), "FIRSTPROTEIN")
	TEST_EQUAL(id.getHits()[1].getSequence(), "SECONDPROTEIN")
	TEST_EQUAL(id.getHits()[2].getSequence(), "THIRDPROTEIN")
	TEST_EQUAL(id.getHits()[0].getRank(), 1)	
	TEST_EQUAL(id.getHits()[1].getRank(), 2)
	TEST_EQUAL(id.getHits()[2].getRank(), 3)
END_SECTION

START_SECTION(void getReferencingHits(const String &protein_accession, std::vector< PeptideHit > &peptide_hits) const)
	PeptideIdentification id;
	PeptideHit hit;
	vector< PeptideHit > peptide_hits;
	
	hit.setScore(23);
	hit.setSequence("FIRSTPROTEIN");
	hit.addProteinAccession("TEST_PROTEIN1");
	id.insertHit(hit);
	
	hit = PeptideHit();
	hit.setScore(10);
	hit.setSequence("SECONDPROTEIN");
	hit.addProteinAccession("TEST_PROTEIN2");
	id.insertHit(hit);

	hit = PeptideHit();
	hit.setScore(11);
	hit.setSequence("THIRDPROTEIN");
	hit.addProteinAccession("TEST_PROTEIN2");
	id.insertHit(hit);

	id.getReferencingHits("TEST_PROTEIN2", peptide_hits);
	TEST_EQUAL(peptide_hits.size(), 2)
	TEST_EQUAL(peptide_hits[0].getSequence(), String("SECONDPROTEIN"))
	TEST_EQUAL(peptide_hits[1].getSequence(), String("THIRDPROTEIN"))
	
END_SECTION

START_SECTION(void getReferencingHits(const std::vector< String > &accessions, std::vector< PeptideHit > &peptide_hits) const)
	PeptideIdentification id;
	PeptideHit hit;
	vector< PeptideHit > peptide_hits;
	vector<String> accessions;
	
	accessions.push_back("TEST_PROTEIN2");
	accessions.push_back("TEST_PROTEIN3");
	
	hit.setScore(23);
	hit.setSequence("FIRSTPROTEIN");
	hit.addProteinAccession("TEST_PROTEIN1");
	id.insertHit(hit);
	
	hit = PeptideHit();
	hit.setScore(10);
	hit.setSequence("SECONDPROTEIN");
	hit.addProteinAccession("TEST_PROTEIN2");
	id.insertHit(hit);

	hit = PeptideHit();
	hit.setScore(11);
	hit.setSequence("THIRDPROTEIN");
	hit.addProteinAccession("TEST_PROTEIN3");
	id.insertHit(hit);

	id.getReferencingHits(accessions, peptide_hits);
	TEST_EQUAL(peptide_hits.size(), 2)
	TEST_EQUAL(peptide_hits[0].getSequence(), String("SECONDPROTEIN"))
	TEST_EQUAL(peptide_hits[1].getSequence(), String("THIRDPROTEIN"))	
END_SECTION

START_SECTION(void getReferencingHits(const std::vector< ProteinHit > &protein_hits, std::vector< PeptideHit > &peptide_hits) const)
	PeptideIdentification id;
	PeptideHit hit;
	vector< PeptideHit > peptide_hits;
	vector<ProteinHit> protein_hits;
	ProteinHit p_hit;
	
	p_hit.setAccession("TEST_PROTEIN2");
	protein_hits.push_back(p_hit);
	p_hit.setAccession("TEST_PROTEIN3");
	protein_hits.push_back(p_hit);
			
	hit.setScore(23);
	hit.setSequence("FIRSTPROTEIN");
	hit.addProteinAccession("TEST_PROTEIN1");
	id.insertHit(hit);
	
	hit = PeptideHit();
	hit.setScore(10);
	hit.setSequence("SECONDPROTEIN");
	hit.addProteinAccession("TEST_PROTEIN2");
	id.insertHit(hit);

	hit = PeptideHit();
	hit.setScore(11);
	hit.setSequence("THIRDPROTEIN");
	hit.addProteinAccession("TEST_PROTEIN3");
	id.insertHit(hit);

	id.getReferencingHits(protein_hits, peptide_hits);
	TEST_EQUAL(peptide_hits.size(), 2)
	TEST_EQUAL(peptide_hits[0].getSequence(), String("SECONDPROTEIN"))
	TEST_EQUAL(peptide_hits[1].getSequence(), String("THIRDPROTEIN"))	
END_SECTION

START_SECTION(void getNonReferencingHits(const String &protein_accession, std::vector< PeptideHit > &peptide_hits) const)
	PeptideIdentification id;
	PeptideHit hit;
	vector< PeptideHit > peptide_hits;
	
	hit.setScore(23);
	hit.setSequence("FIRSTPROTEIN");
	hit.addProteinAccession("TEST_PROTEIN1");
	id.insertHit(hit);
	
	hit = PeptideHit();
	hit.setScore(10);
	hit.setSequence("SECONDPROTEIN");
	hit.addProteinAccession("TEST_PROTEIN2");
	id.insertHit(hit);

	hit = PeptideHit();
	hit.setScore(11);
	hit.setSequence("THIRDPROTEIN");
	hit.addProteinAccession("TEST_PROTEIN2");
	id.insertHit(hit);

	id.getNonReferencingHits("TEST_PROTEIN2", peptide_hits);
	TEST_EQUAL(peptide_hits.size(), 1)
	TEST_EQUAL(peptide_hits[0].getSequence(), String("FIRSTPROTEIN"))
END_SECTION

START_SECTION(void getNonReferencingHits(const std::vector< String > &accessions, std::vector< PeptideHit > &peptide_hits) const)
	PeptideIdentification id;
	PeptideHit hit;
	vector< PeptideHit > peptide_hits;
	vector<String> accessions;
	
	accessions.push_back("TEST_PROTEIN2");
	accessions.push_back("TEST_PROTEIN3");
	
	hit.setScore(23);
	hit.setSequence("FIRSTPROTEIN");
	hit.addProteinAccession("TEST_PROTEIN1");
	id.insertHit(hit);
	
	hit = PeptideHit();
	hit.setScore(10);
	hit.setSequence("SECONDPROTEIN");
	hit.addProteinAccession("TEST_PROTEIN2");
	id.insertHit(hit);

	hit = PeptideHit();
	hit.setScore(11);
	hit.setSequence("THIRDPROTEIN");
	hit.addProteinAccession("TEST_PROTEIN3");
	id.insertHit(hit);

	id.getNonReferencingHits(accessions, peptide_hits);
	TEST_EQUAL(peptide_hits.size(), 1)
	TEST_EQUAL(peptide_hits[0].getSequence(), String("FIRSTPROTEIN"))

END_SECTION

START_SECTION(void getNonReferencingHits(const std::vector< ProteinHit > &protein_hits, std::vector< PeptideHit > &peptide_hits) const)
	PeptideIdentification id;
	PeptideHit hit;
	vector< PeptideHit > peptide_hits;
	vector<ProteinHit> protein_hits;
	ProteinHit p_hit;
	
	p_hit.setAccession("TEST_PROTEIN2");
	protein_hits.push_back(p_hit);
	p_hit.setAccession("TEST_PROTEIN3");
	protein_hits.push_back(p_hit);
			
	hit.setScore(23);
	hit.setSequence("FIRSTPROTEIN");
	hit.addProteinAccession("TEST_PROTEIN1");
	id.insertHit(hit);
	
	hit = PeptideHit();
	hit.setScore(10);
	hit.setSequence("SECONDPROTEIN");
	hit.addProteinAccession("TEST_PROTEIN2");
	id.insertHit(hit);

	hit = PeptideHit();
	hit.setScore(11);
	hit.setSequence("THIRDPROTEIN");
	hit.addProteinAccession("TEST_PROTEIN3");
	id.insertHit(hit);

	id.getNonReferencingHits(protein_hits, peptide_hits);
	TEST_EQUAL(peptide_hits.size(), 1)
	TEST_EQUAL(peptide_hits[0].getSequence(), String("FIRSTPROTEIN"))
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
