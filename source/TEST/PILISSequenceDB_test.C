// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <iostream>

#include <OpenMS/ANALYSIS/ID/PILISSequenceDB.h>
#include <OpenMS/FORMAT/AnalysisXMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>

///////////////////////////

START_TEST(PILISSequenceDB, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

PILISSequenceDB* ptr = 0;

CHECK((PILISSequenceDB()))
	ptr = new PILISSequenceDB();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~PILISSequenceDB()))
	delete ptr;
RESULT

ptr = new PILISSequenceDB();

CHECK((PILISSequenceDB(const PILISSequenceDB&)))
	PILISSequenceDB copy(*ptr);
	TEST_EQUAL(copy.countPeptides(), ptr->countPeptides())
	TEST_EQUAL(copy.countProteins(), ptr->countProteins())
	vector<PILISSequenceDB::PepStruct> peptides1, peptides2;
	copy.getPeptides(peptides1);
	ptr->getPeptides(peptides2);
	TEST_EQUAL(peptides1.size(), peptides2.size());
	for (Size i = 0; i != peptides1.size(); ++i)
	{
		TEST_EQUAL(peptides1[i].peptide, peptides2[i].peptide)
	}
RESULT

CHECK((PILISSequenceDB& operator = (const PILISSequenceDB& rhs)))
	PILISSequenceDB copy;
	copy = *ptr;

  TEST_EQUAL(copy.countPeptides(), ptr->countPeptides())
  TEST_EQUAL(copy.countProteins(), ptr->countProteins())
  vector<PILISSequenceDB::PepStruct> peptides1, peptides2;
  copy.getPeptides(peptides1);
  ptr->getPeptides(peptides2);
  TEST_EQUAL(peptides1.size(), peptides2.size());
  for (Size i = 0; i != peptides1.size(); ++i)
  {
    TEST_EQUAL(peptides1[i].peptide, peptides2[i].peptide)
  }

RESULT

CHECK((unsigned int countPeptides() const))
	TEST_EQUAL(ptr->countPeptides(), 0)
RESULT

CHECK((unsigned int countProteins() const))
	TEST_EQUAL(ptr->countProteins(), 0)
RESULT

CHECK((void addPeptidesFromFile(const String& filename)))
	ptr->addPeptidesFromFile("data/PILISSequenceDB_sequence.db");
	TEST_EQUAL(ptr->countPeptides(), 542)
RESULT

CHECK((bool has(const String& peptide) const))
	TEST_EQUAL(ptr->has("DFPIANGER"), true)
	TEST_EQUAL(ptr->has("DFPIANGERDFPIANGER"), false)
RESULT

CHECK((void getPeptides(std::vector<PepStruct>& peptides, double range_start = 0, double range_stop = std::numeric_limits<double>::max())))
	vector<PILISSequenceDB::PepStruct> peptides;
	ptr->getPeptides(peptides, 1017.7, 1021.7);
	TEST_EQUAL(peptides.size(), 2)
	
	peptides.clear();
	ptr->getPeptides(peptides, 205, 210);
	TEST_EQUAL(peptides.size(), 0)
RESULT

CHECK((void clearProteins()))
	ptr->clearProteins();
	TEST_EQUAL(ptr->countProteins(), 0)
RESULT

CHECK((void clearPeptides()))
	ptr->clearPeptides();
	TEST_EQUAL(ptr->countPeptides(), 0)
RESULT

CHECK((bool isReplaceXandL() const))
	TEST_EQUAL(ptr->isReplaceXandL(), true);
RESULT

CHECK((double getFactor() const))
	TEST_REAL_EQUAL(ptr->getFactor(), 10.0)
RESULT

CHECK((void digestProteinsTryptic(Size missed_cleavages = 0)))
	ptr->digestProteinsTryptic();
	// TODO
RESULT

CHECK((void setFactor(double factor)))
	ptr->setFactor(200.0);
	TEST_REAL_EQUAL(ptr->getFactor(), 200.0)
RESULT

CHECK((void setReplaceXandL(bool replace = true)))
	ptr->setReplaceXandL(false);
	TEST_EQUAL(ptr->isReplaceXandL(), false)
RESULT

CHECK((void addFASTAFile(const String& filename)))
	TEST_EXCEPTION(Exception::NotImplemented, ptr->addFASTAFile("does_not_exist"))
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
