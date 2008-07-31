// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <string>

#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/DATASTRUCTURES/String.h>

///////////////////////////

START_TEST(ProteinHit, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;


Real score = 4.4;
UInt rank = 3;
String sequence = "ARRAY";
String accession = "PROOE34";


ProteinHit* ptr = 0;	
CHECK(ProteinHit())
	ptr = new ProteinHit();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~ProteinHit())
	ptr = new ProteinHit();	
  delete ptr;
RESULT

CHECK((ProteinHit(DoubleReal score, UInt rank, String accession, String sequence)))
	ProteinHit hit(score, rank, accession, sequence);
	TEST_EQUAL(hit.getScore(), score)
	TEST_EQUAL(hit.getRank(), rank)
	TEST_EQUAL(hit.getAccession(), accession)
	TEST_EQUAL(hit.getSequence(), sequence)
RESULT

CHECK(ProteinHit(const ProteinHit& source))
	ProteinHit source;
	source.setScore(score);
	source.setRank(rank);
	source.setAccession(accession);
	source.setSequence(sequence);
	source.setMetaValue("label",17);
  
  ProteinHit hit(source);
	
	TEST_EQUAL(hit.getScore(), score)
	TEST_EQUAL(hit.getRank(), rank)
	TEST_EQUAL(hit.getAccession(), accession)
	TEST_EQUAL(hit.getSequence(), sequence)
	TEST_EQUAL((UInt)hit.getMetaValue("label"),17)
RESULT

CHECK(ProteinHit& operator=(const ProteinHit& source))
	ProteinHit hit;
	ProteinHit hit2(score, rank, accession, sequence);
	hit2.setMetaValue("label",17);
	
	hit = hit2;
	
	TEST_EQUAL(hit.getScore(), score)
	TEST_EQUAL(hit.getRank(), rank)
	TEST_EQUAL(hit.getAccession(), accession)
	TEST_EQUAL(hit.getSequence(), sequence)
	TEST_EQUAL((UInt)hit.getMetaValue("label"),17)
RESULT

CHECK(bool operator == (const ProteinHit& rhs) const)
  ProteinHit hit, hit2;
  TEST_EQUAL(hit==hit2,true);
  
  hit.setScore(score);
  TEST_EQUAL(hit==hit2,false);
	hit = hit2;
  
	hit.setRank(rank);
  TEST_EQUAL(hit==hit2,false);
	hit = hit2;
  
	hit.setAccession(accession);
  TEST_EQUAL(hit==hit2,false);
	hit = hit2;
  
	hit.setSequence(sequence);
  TEST_EQUAL(hit==hit2,false);
	hit = hit2;
	
	hit.setMetaValue("label",17);
  TEST_EQUAL(hit==hit2,false);
	hit = hit2;
RESULT

CHECK(bool operator != (const ProteinHit& rhs) const)
  ProteinHit hit, hit2;
  TEST_EQUAL(hit!=hit2,false);
  
  hit.setScore(score);
  TEST_EQUAL(hit!=hit2,true);
  hit = hit2;
  
  hit.setRank(rank);
  TEST_EQUAL(hit!=hit2,true);
	hit = hit2;
  
  hit.setAccession(accession);
  TEST_EQUAL(hit!=hit2,true);
	hit = hit2;
  
	hit.setSequence(sequence);
  TEST_EQUAL(hit!=hit2,true);
	hit = hit2;

	hit.setMetaValue("label",17);
  TEST_EQUAL(hit!=hit2,true);
	hit = hit2;
RESULT

CHECK(const String& getAccession() const)
	ProteinHit hit(score, rank, accession, sequence);
	TEST_EQUAL(hit.getAccession(), accession)
RESULT

CHECK(const String& getSequence() const)
	ProteinHit hit(score, rank, accession, sequence);
	TEST_EQUAL(hit.getSequence(), sequence)
RESULT

CHECK(Real getScore() const)
	ProteinHit hit(score, rank, accession, sequence);
	TEST_EQUAL(hit.getScore(), score)
RESULT

CHECK(UInt getRank() const)
	ProteinHit hit(score, rank, accession, sequence);
	TEST_EQUAL(hit.getRank(), rank)
RESULT

CHECK(void setRank(UInt newrank))
	ProteinHit hit;
	hit.setRank(rank);
	TEST_EQUAL(hit.getRank(), rank)	
RESULT

CHECK(void setScore(DoubleReal score))
	ProteinHit hit;
	hit.setScore(score);
	TEST_EQUAL(hit.getScore(), score);
RESULT

CHECK(void setSequence(const String& sequence))
	ProteinHit hit;
	hit.setSequence(sequence);
	TEST_EQUAL(hit.getSequence(), sequence)
RESULT

CHECK(void setAccession(const String& accession))
	ProteinHit hit;
	hit.setAccession(accession);
	TEST_EQUAL(hit.getAccession(), accession)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
