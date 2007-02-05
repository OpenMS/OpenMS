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

#include <OpenMS/ANALYSIS/ID/PILISIdentification.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>

///////////////////////////

START_TEST(PILISIdentification, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

PILISIdentification* ptr = 0;
PeakSpectrum spec;
DTAFile().load("data/PILISSequenceDB_DFPIANGER_1.dta", spec);
spec.setMSLevel(2);
CHECK(PILISIdentification())
	ptr = new PILISIdentification();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~PILISIdentification())
	delete ptr;
RESULT

ptr = new PILISIdentification();

CHECK(PILISIdentification(const PILISIdentification& source))
	PILISIdentification copy(*ptr);
	TEST_EQUAL(copy.getParameters(), ptr->getParameters())
RESULT

CHECK(PILISIdentification& operator = (const PILISIdentification&))
	PILISIdentification copy;
	copy = *ptr;
	TEST_EQUAL(copy.getParameters(), ptr->getParameters())
RESULT


CHECK(Param& getParameters())
	Param p(ptr->getParameters());
	p.setValue("bla", "blubb");
	ptr->setParameters(p);
	TEST_EQUAL(ptr->getParameters().getValue("bla"), "blubb")
RESULT

CHECK(const Param& getParameters() const)
	Param p(ptr->getParameters());
	TEST_EQUAL(p.getValue("bla"), "blubb")
RESULT

CHECK(void setParameters(const Param& param))
	Param p;
	p.setValue("blubb", "bla");
	ptr->setParameters(p);
RESULT

//CHECK(void resetToDefaultParam())
//	ptr->resetToDefaultParam();
//	PILISIdentification id;
//	TEST_EQUAL(ptr->getParameters(), id.getParameters())
//RESULT

CHECK(void setSequenceDB(PILISSequenceDB* sequence_db))
	PILISSequenceDB* db = new PILISSequenceDB();
	db->addPeptidesFromFile("data/PILISSequenceDB_sequence.db");
	ptr->setSequenceDB(db);
RESULT

CHECK(void setModel(PILISModel* hmm_model))
	PILISModel* model = new PILISModel();
	model->readFromFile("../../data/PILIS/PILIS_default_model.dat");
	ptr->setModel(model);
RESULT

CHECK(void getIdentification(Identification& id, const PeakSpectrum& spectrum))
	Identification id;
	ptr->getIdentification(id, spec);
	TEST_EQUAL(id.getPeptideHits().size(), 3)
	TEST_EQUAL(id.getPeptideHits().begin()->getSequence(), "DFPIANGER")
RESULT

CHECK(void getIdentifications(std::vector<Identification>& ids, const PeakMap& exp))
	vector<Identification> ids;
	PeakMap map;
	map.push_back(spec);
	ptr->getIdentifications(ids, map);
	TEST_EQUAL(ids.size(), map.size())
	TEST_EQUAL(ids.begin()->getPeptideHits().size(), 3)
	TEST_EQUAL(ids.begin()->getPeptideHits().begin()->getSequence(), "DFPIANGER")
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
