// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

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
RichPeakSpectrum spec;
DTAFile().load(OPENMS_GET_TEST_DATA_PATH("PILISSequenceDB_DFPIANGER_1.dta"), spec);
spec.setMSLevel(2);
START_SECTION(PILISIdentification())
	ptr = new PILISIdentification();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION(~PILISIdentification())
	delete ptr;
END_SECTION

ptr = new PILISIdentification();

START_SECTION(PILISIdentification(const PILISIdentification& source))
	PILISIdentification copy(*ptr);
	TEST_EQUAL(copy.getParameters(), ptr->getParameters())
END_SECTION

START_SECTION(PILISIdentification& operator = (const PILISIdentification& source))
	PILISIdentification copy;
	copy = *ptr;
	TEST_EQUAL(copy.getParameters(), ptr->getParameters())
END_SECTION

START_SECTION(void setModel(PILISModel* hmm_model))
	PILISModel* model = new PILISModel();
	model->readFromFile("PILIS/PILIS_default_model.dat");
	ptr->setModel(model);
END_SECTION

START_SECTION(void getIdentification(const std::map<String, UInt>& candidates, PeptideIdentification& id, const RichPeakSpectrum& spectrum))
	map<String, UInt> candidates;
	candidates["DDFPIVIVGNKADIENQR"] = 2;
	candidates["DFPIANGER"] = 1;
	candidates["DFPIADGER"] = 1;
	PeptideIdentification id;
	ptr->getIdentification(candidates, id, spec);
	TEST_EQUAL(id.getHits().size(), 3)
	TEST_EQUAL(id.getHits().begin()->getSequence(), "DFPIANGER")
END_SECTION

START_SECTION(void getIdentifications(const std::vector<std::map<String, UInt> >& candidates, std::vector<PeptideIdentification>& ids, const RichPeakMap& exp))

	map<String, UInt> cand;
	cand["DDFPIVIVGNKADIENQR"] = 2;
	cand["DFPIANGER"] = 1;
	cand["DFPIADGER"] = 1;
	vector<map<String, UInt> > candidates;
	candidates.push_back(cand);

	vector<PeptideIdentification> ids;
	RichPeakMap map;
	map.push_back(spec);
	ptr->getIdentifications(candidates, ids, map);
	TEST_EQUAL(ids.size(), map.size())
	TEST_EQUAL(ids.begin()->getHits().size(), 3)
	TEST_EQUAL(ids.begin()->getHits().begin()->getSequence(), "DFPIANGER")
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
