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
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

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


double peptide_significance_threshold = 42.3;
std::vector<PeptideHit> peptide_hits;
PeptideHit peptide_hit;
ProteinIdentification protein_identification;
vector<PeptideIdentification> identifications;
MascotXMLFile xml_file;

peptide_hits.push_back(peptide_hit);


PeptideIdentification* ptr = nullptr;
PeptideIdentification* nullPointer = nullptr;
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


START_SECTION((double getRT() const))
	PeptideIdentification pi;
	TEST_EQUAL(pi.hasRT(), false);
	pi.setRT(1024.0);
	TEST_EQUAL(pi.getRT(), 1024.0);
END_SECTION

START_SECTION(void setRT(double mz))
	NOT_TESTABLE // tested above
END_SECTION

START_SECTION(bool hasRT())
	NOT_TESTABLE // tested above
END_SECTION

START_SECTION((double getMZ() const))
	PeptideIdentification pi;
	TEST_EQUAL(pi.hasMZ(), false);
	pi.setMZ(1024.0);
	TEST_EQUAL(pi.getMZ(), 1024.0);
END_SECTION


START_SECTION(bool hasMZ())
  NOT_TESTABLE // tested above
END_SECTION


START_SECTION((double getSignificanceThreshold() const))
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

START_SECTION((void setSignificanceThreshold(double value)))
  PeptideIdentification hits;
  hits.setSignificanceThreshold(peptide_significance_threshold);
  TEST_EQUAL(hits.getSignificanceThreshold(), peptide_significance_threshold)
END_SECTION

START_SECTION((String& getScoreType() const))
	PeptideIdentification hits;
	TEST_EQUAL(hits.getScoreType(), "")
END_SECTION

START_SECTION((void setScoreType(const String& type)))
	PeptideIdentification hits;
	hits.setScoreType("bla");
	TEST_EQUAL(hits.getScoreType(), "bla")
END_SECTION

START_SECTION((bool isHigherScoreBetter() const))
	PeptideIdentification hits;
	TEST_EQUAL(hits.isHigherScoreBetter(), true)
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
  TEST_EQUAL(hits.empty(), true)

  hits.setBaseName("basename");
  TEST_EQUAL(hits.empty(), false)

  hits.setBaseName("");
  TEST_EQUAL(hits.empty(), true)

  hits.insertHit(peptide_hit);
  TEST_EQUAL(hits.empty(), false)
END_SECTION

START_SECTION((void sort()))
  PeptideIdentification id;
  PeptideHit hit;
  hit.setScore(23);
  hit.setSequence(AASequence::fromString("SECONDPROTEIN"));
  id.insertHit(hit);
  hit.setScore(45);
  hit.setSequence(AASequence::fromString("FIRSTPROTEIN"));
  id.insertHit(hit);
  hit.setScore(7);
  hit.setSequence(AASequence::fromString("THIRDPROTEIN"));
  id.insertHit(hit);

  //higher score is better
  id.sort();

  TEST_EQUAL(id.getHits()[0].getSequence(), AASequence::fromString("FIRSTPROTEIN"))
  TEST_EQUAL(id.getHits()[1].getSequence(), AASequence::fromString("SECONDPROTEIN"))
  TEST_EQUAL(id.getHits()[2].getSequence(), AASequence::fromString("THIRDPROTEIN"))
  TEST_EQUAL(id.getHits()[0].getScore(), 45)
  TEST_EQUAL(id.getHits()[1].getScore(), 23)
  TEST_EQUAL(id.getHits()[2].getScore(), 7)

  //lower score is better
  id.setHigherScoreBetter(false);
  id.sort();

  TEST_EQUAL(id.getHits()[0].getSequence(), AASequence::fromString("THIRDPROTEIN"))
  TEST_EQUAL(id.getHits()[1].getSequence(), AASequence::fromString("SECONDPROTEIN"))
  TEST_EQUAL(id.getHits()[2].getSequence(), AASequence::fromString("FIRSTPROTEIN"))
  TEST_EQUAL(id.getHits()[0].getScore(), 7)
  TEST_EQUAL(id.getHits()[1].getScore(), 23)
  TEST_EQUAL(id.getHits()[2].getScore(), 45)

END_SECTION

START_SECTION((void assignRanks()))
  PeptideIdentification id;
  PeptideHit hit;
  hit.setScore(23);
  hit.setSequence(AASequence::fromString("SECONDPROTEIN"));
  id.insertHit(hit);
  hit.setScore(45);
  hit.setSequence(AASequence::fromString("FIRSTPROTEIN"));
  id.insertHit(hit);
  hit.setScore(7);
  hit.setSequence(AASequence::fromString("THIRDPROTEIN"));
  id.insertHit(hit);

  id.assignRanks();

  TEST_EQUAL(id.getHits()[0].getSequence(), AASequence::fromString("FIRSTPROTEIN"))
  TEST_EQUAL(id.getHits()[1].getSequence(), AASequence::fromString("SECONDPROTEIN"))
  TEST_EQUAL(id.getHits()[2].getSequence(), AASequence::fromString("THIRDPROTEIN"))
  TEST_EQUAL(id.getHits()[0].getRank(), 1)
  TEST_EQUAL(id.getHits()[1].getRank(), 2)
  TEST_EQUAL(id.getHits()[2].getRank(), 3)
END_SECTION

START_SECTION(static std::vector<PeptideHit> getReferencingHits(const std::vector<PeptideHit> & , const std::set<String> & accession))
{
  PeptideIdentification id;
  PeptideHit hit;
  vector< PeptideHit > peptide_hits;

  hit.setScore(23);
  hit.setSequence(AASequence::fromString("FIRSTPROTEIN"));
  PeptideEvidence pe;
  pe.setProteinAccession("TEST_PROTEIN1");
  hit.addPeptideEvidence(pe);
  id.insertHit(hit);

  hit = PeptideHit();
  hit.setScore(10);
  hit.setSequence(AASequence::fromString("SECONDPROTEIN"));
  pe.setProteinAccession("TEST_PROTEIN2");
  hit.addPeptideEvidence(pe);
  id.insertHit(hit);

  hit = PeptideHit();
  hit.setScore(11);
  hit.setSequence(AASequence::fromString("THIRDPROTEIN"));
  pe.setProteinAccession("TEST_PROTEIN2");
  hit.addPeptideEvidence(pe);
  id.insertHit(hit);

  set<String> query_accession;
  query_accession.insert("TEST_PROTEIN2");
  peptide_hits = PeptideIdentification::getReferencingHits(id.getHits(), query_accession);
  TEST_EQUAL(peptide_hits.size(), 2)
  TEST_EQUAL(peptide_hits[0].getSequence(), AASequence::fromString("SECONDPROTEIN"))
  TEST_EQUAL(peptide_hits[1].getSequence(), AASequence::fromString("THIRDPROTEIN"))

  query_accession.insert("TEST_PROTEIN3");
  peptide_hits = PeptideIdentification::getReferencingHits(id.getHits(), query_accession);
  TEST_EQUAL(peptide_hits.size(), 2)
  TEST_EQUAL(peptide_hits[0].getSequence(), AASequence::fromString("SECONDPROTEIN"))
  TEST_EQUAL(peptide_hits[1].getSequence(), AASequence::fromString("THIRDPROTEIN"))
}
{
  PeptideIdentification id;
  PeptideHit hit;
  vector< PeptideHit > peptide_hits;

  hit.setScore(23);
  hit.setSequence(AASequence::fromString("FIRSTPROTEIN"));
  PeptideEvidence pe;
  pe.setProteinAccession("TEST_PROTEIN1");
  hit.addPeptideEvidence(pe);
  id.insertHit(hit);

  hit = PeptideHit();
  hit.setScore(10);
  hit.setSequence(AASequence::fromString("SECONDPROTEIN"));
  pe.setProteinAccession("TEST_PROTEIN2");
  hit.addPeptideEvidence(pe);
  id.insertHit(hit);

  hit = PeptideHit();
  hit.setScore(11);
  hit.setSequence(AASequence::fromString("THIRDPROTEIN"));
  pe.setProteinAccession("TEST_PROTEIN3");
  hit.addPeptideEvidence(pe);
  id.insertHit(hit);

  set<String> query_accession;
  query_accession.insert("TEST_PROTEIN2");
  query_accession.insert("TEST_PROTEIN3");
  peptide_hits = PeptideIdentification::getReferencingHits(id.getHits(), query_accession);
  TEST_EQUAL(peptide_hits.size(), 2)
  TEST_EQUAL(peptide_hits[0].getSequence(), AASequence::fromString("SECONDPROTEIN"))
  TEST_EQUAL(peptide_hits[1].getSequence(), AASequence::fromString("THIRDPROTEIN"))
}
END_SECTION

/*
START_SECTION(void getNonReferencingHits(const String &protein_accession, std::vector< PeptideHit > &peptide_hits) const)
  PeptideIdentification id;
  PeptideHit hit;
  vector< PeptideHit > peptide_hits;

  hit.setScore(23);
  hit.setSequence(AASequence::fromString("FIRSTPROTEIN"));
  hit.addProteinAccession("TEST_PROTEIN1");
  id.insertHit(hit);

  hit = PeptideHit();
  hit.setScore(10);
  hit.setSequence(AASequence::fromString("SECONDPROTEIN"));
  hit.addProteinAccession("TEST_PROTEIN2");
  id.insertHit(hit);

  hit = PeptideHit();
  hit.setScore(11);
  hit.setSequence(AASequence::fromString("THIRDPROTEIN"));
  hit.addProteinAccession("TEST_PROTEIN2");
  id.insertHit(hit);

  id.getNonReferencingHits("TEST_PROTEIN2", peptide_hits);
  TEST_EQUAL(peptide_hits.size(), 1)
  TEST_EQUAL(peptide_hits[0].getSequence(), AASequence::fromString("FIRSTPROTEIN"))
END_SECTION

START_SECTION(void getNonReferencingHits(const std::vector< String > &accessions, std::vector< PeptideHit > &peptide_hits) const)
  PeptideIdentification id;
  PeptideHit hit;
  vector< PeptideHit > peptide_hits;
  vector<String> accessions;

  accessions.push_back("TEST_PROTEIN2");
  accessions.push_back("TEST_PROTEIN3");

  hit.setScore(23);
  hit.setSequence(AASequence::fromString("FIRSTPROTEIN"));
  hit.addProteinAccession("TEST_PROTEIN1");
  id.insertHit(hit);

  hit = PeptideHit();
  hit.setScore(10);
  hit.setSequence(AASequence::fromString("SECONDPROTEIN"));
  hit.addProteinAccession("TEST_PROTEIN2");
  id.insertHit(hit);

  hit = PeptideHit();
  hit.setScore(11);
  hit.setSequence(AASequence::fromString("THIRDPROTEIN"));
  hit.addProteinAccession("TEST_PROTEIN3");
  id.insertHit(hit);

  id.getNonReferencingHits(accessions, peptide_hits);
  TEST_EQUAL(peptide_hits.size(), 1)
  TEST_EQUAL(peptide_hits[0].getSequence(), AASequence::fromString("FIRSTPROTEIN"))
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
  hit.setSequence(AASequence::fromString("FIRSTPROTEIN"));
  hit.addProteinAccession("TEST_PROTEIN1");
  id.insertHit(hit);

  hit = PeptideHit();
  hit.setScore(10);
  hit.setSequence(AASequence::fromString("SECONDPROTEIN"));
  hit.addProteinAccession("TEST_PROTEIN2");
  id.insertHit(hit);

  hit = PeptideHit();
  hit.setScore(11);
  hit.setSequence(AASequence::fromString("THIRDPROTEIN"));
  hit.addProteinAccession("TEST_PROTEIN3");
  id.insertHit(hit);

  id.getNonReferencingHits(protein_hits, peptide_hits);
  TEST_EQUAL(peptide_hits.size(), 1)
  TEST_EQUAL(peptide_hits[0].getSequence(), AASequence::fromString("FIRSTPROTEIN"))
END_SECTION
*/
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
