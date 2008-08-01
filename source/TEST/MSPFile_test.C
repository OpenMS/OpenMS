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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
///////////////////////////

#include <OpenMS/FORMAT/MSPFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(MSPFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MSPFile* ptr = 0;
CHECK((MSPFile()))
	ptr = new MSPFile;
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~MSPFile()))
	delete ptr;
RESULT

CHECK(MSPFile(const MSPFile &rhs))
	MSPFile f1, f2;
	Param p = f1.getParameters();
	p.setValue("instrument", "it");
	f1.setParameters(p);
	TEST_EQUAL(f1.getParameters() == f2.getParameters(), false)
	MSPFile f3(f1);
	TEST_EQUAL(f1.getParameters() == f3.getParameters(), true)
RESULT

CHECK(MSPFile& operator=(const MSPFile &rhs))
	MSPFile f1, f2;
	Param p = f1.getParameters();
	p.setValue("instrument", "it");
	f1.setParameters(p);
	TEST_EQUAL(f1.getParameters() == f2.getParameters(), false)
	f2 = f1;
	TEST_EQUAL(f1.getParameters() == f2.getParameters(), true)
RESULT

CHECK(void load(const String &filename, std::vector< PeptideIdentification > &ids, RichPeakMap &exp))
	MSPFile msp_file;
	vector<PeptideIdentification> ids;
	RichPeakMap exp;
	msp_file.load("data/MSPFile_test.msp", ids, exp);
	TEST_EQUAL(exp.size(), 5)
	TEST_EQUAL(ids.size(), 5)

	Param p(msp_file.getParameters());
	p.setValue("instrument", "qtof");
	msp_file.setParameters(p);
	ids.clear();
	exp.clear();
	msp_file.load("data/MSPFile_test.msp", ids, exp);
	TEST_EQUAL(exp.size(), 2)
	TEST_EQUAL(ids.size(), 2)

	p.setValue("instrument", "it");
	msp_file.setParameters(p);
	ids.clear();
	exp.clear();
	msp_file.load("data/MSPFile_test.msp", ids, exp);
	TEST_EQUAL(exp.size(), 3)
	TEST_EQUAL(ids.size(), 3)
RESULT

CHECK(void store(const String &, const RichPeakMap &) const)
	MSPFile msp_file;
	RichPeakMap exp;
	TEST_EXCEPTION(Exception::NotImplemented, msp_file.store("this_file_will_never_be_created", exp))
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
