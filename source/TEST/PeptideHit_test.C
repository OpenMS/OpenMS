// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/DATASTRUCTURES/String.h>

///////////////////////////

START_TEST(PeptideHit, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

DoubleReal score = 4.4;
uint rank = 3;
AASequence sequence = AASequence("ARRAY");
std::string sequence2 = "  ARRAY  ";
Int charge = 2;

PeptideHit* ptr = 0;
START_SECTION((PeptideHit()))
	ptr = new PeptideHit();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((virtual ~PeptideHit()))
	delete ptr;
END_SECTION

START_SECTION((PeptideHit(DoubleReal score, UInt rank, Int charge, const AASequence &sequence)))
	PeptideHit hit(score, rank, charge, sequence);
	TEST_EQUAL(hit.getScore(), score)
	TEST_EQUAL(hit.getRank(), rank)
	TEST_EQUAL(hit.getCharge(), charge)
	TEST_EQUAL(hit.getSequence(), sequence)
END_SECTION

START_SECTION((PeptideHit& operator=(const PeptideHit& source)))
	PeptideHit hit;
	PeptideHit hit2(score, rank, charge, sequence);
	hit2.setMetaValue("label",17);
	
	hit = hit2;
	
	TEST_EQUAL(hit.getScore(), score)
	TEST_EQUAL(hit.getRank(), rank)
	TEST_EQUAL(hit.getCharge(), charge)
	TEST_EQUAL(hit.getSequence(), sequence)
	TEST_EQUAL((UInt)hit.getMetaValue("label"),17)
END_SECTION

START_SECTION((PeptideHit(const PeptideHit& source)))
	PeptideHit source;
	source.setScore(score);
	source.setRank(rank);
	source.setSequence(sequence);
	source.setMetaValue("label",17);
	
  PeptideHit hit(source);
	
	TEST_EQUAL(hit.getScore(), source.getScore())
	TEST_EQUAL(hit.getRank(), source.getRank())
	TEST_EQUAL(hit.getSequence(), source.getSequence())
	TEST_EQUAL((UInt)hit.getMetaValue("label"),17) 
END_SECTION

START_SECTION((bool operator == (const PeptideHit& rhs) const))
  PeptideHit hit, hit2;
  TEST_EQUAL(hit==hit2,true);

  hit.setScore(score);
  TEST_EQUAL(hit==hit2,false);
	hit=hit2;
	
  hit.setRank(rank);
  TEST_EQUAL(hit==hit2,false);
	hit=hit2;
	
	hit.setSequence(sequence);
  TEST_EQUAL(hit==hit2,false);
	hit=hit2;
	
	hit.setMetaValue("label",17);
  TEST_EQUAL(hit==hit2,false);
	hit=hit2;
END_SECTION

START_SECTION((bool operator != (const PeptideHit& rhs) const))
  PeptideHit hit, hit2;
  TEST_EQUAL(hit!=hit2,false);

  hit.setScore(score);
  TEST_EQUAL(hit!=hit2,true);
	hit=hit2;
	
  hit.setRank(rank);
  TEST_EQUAL(hit!=hit2,true);
	hit=hit2;
	
	hit.setSequence(sequence);
  TEST_EQUAL(hit!=hit2,true);
	hit=hit2;
	
	hit.setMetaValue("label",17);
  TEST_EQUAL(hit!=hit2,true);
	hit=hit2;
END_SECTION

START_SECTION((DoubleReal getScore() const ))
	PeptideHit hit(score, rank, charge, sequence);
	TEST_EQUAL(hit.getScore(), score)
END_SECTION

START_SECTION((UInt getRank() const))
	PeptideHit hit(score, rank, charge, sequence);
	TEST_EQUAL(hit.getRank(), rank)
END_SECTION

START_SECTION((const AASequence& getSequence() const))
	PeptideHit hit(score, rank, charge, sequence);
	TEST_EQUAL(hit.getSequence(), sequence)
END_SECTION

START_SECTION((void setRank(UInt newrank)))
	PeptideHit hit;
	hit.setRank(rank);
	TEST_EQUAL(hit.getRank(), rank)
END_SECTION

START_SECTION((void setScore(DoubleReal score)))
	PeptideHit hit;
	hit.setScore(score);
	TEST_EQUAL(hit.getScore(), score)
END_SECTION

START_SECTION((void setSequence(const AASequence& sequence)))
	PeptideHit hit;
	hit.setSequence(sequence);
	TEST_EQUAL(hit.getSequence(), sequence)
	//hit.setSequence(sequence2);
	// @todo std::string interface?
	TEST_EQUAL(hit.getSequence(), sequence)	
END_SECTION

START_SECTION((void addProteinAccession(const String& accession)))
	String date;
	vector<String> indices;

	date = "2006-12-12 11:59:59";
	PeptideHit hit;

	hit.addProteinAccession("ACC392");
	hit.addProteinAccession("ACD392");
	indices = hit.getProteinAccessions();
	TEST_EQUAL(indices.size(), 2)
	TEST_EQUAL(indices[0] == String("ACC392"), true)
	TEST_EQUAL(indices[1] == String("ACD392"), true)
END_SECTION

START_SECTION((void setProteinAccessions(const std::vector< String > &accessions)))
	vector<String> vec;
	vec.push_back("ACC392");
	vec.push_back("ACD392");
	PeptideHit hit;
	hit.addProteinAccession("ACC392");
	hit.addProteinAccession("ACD392");
	TEST_EQUAL(vec == hit.getProteinAccessions(), true)
END_SECTION

START_SECTION((const std::vector<String>& getProteinAccessions() const))
	PeptideHit hit;
	hit.addProteinAccession("ACC392");
	hit.addProteinAccession("ACD392");
	TEST_EQUAL(hit.getProteinAccessions().size(), 2)
	TEST_EQUAL(hit.getProteinAccessions()[0], "ACC392")
	TEST_EQUAL(hit.getProteinAccessions()[1], "ACD392")
END_SECTION

START_SECTION((Int getCharge() const))
	PeptideHit hit;
	
	hit.setCharge(-43);
	TEST_EQUAL(-43, hit.getCharge())
END_SECTION

START_SECTION((void setCharge(Int charge)))
	PeptideHit hit;
	
	hit.setCharge(-43);
	TEST_EQUAL(-43, hit.getCharge())
END_SECTION

START_SECTION((void setAABefore(char acid)))
	PeptideHit hit;
	
	hit.setAABefore('R');
	TEST_EQUAL(hit.getAABefore(), 'R')
END_SECTION
START_SECTION((char getAABefore() const))
	PeptideHit hit;
	
	hit.setAABefore('R');
	TEST_EQUAL(hit.getAABefore(), 'R')
END_SECTION
START_SECTION((void setAAAfter(char acid)))
	PeptideHit hit;
	
	hit.setAAAfter('R');
	TEST_EQUAL(hit.getAAAfter(), 'R')
END_SECTION
START_SECTION((char getAAAfter() const))
	PeptideHit hit;
	
	hit.setAAAfter('R');
	TEST_EQUAL(hit.getAAAfter(), 'R')
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
