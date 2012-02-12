// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/QUANTITATION/ProteinInference.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ProteinInference, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ProteinInference* ptr = 0;
ProteinInference* nullPointer = 0;
START_SECTION(ProteinInference())
{
	ptr = new ProteinInference();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~ProteinInference())
{
	delete ptr;
}
END_SECTION

START_SECTION((ProteinInference(const ProteinInference &cp)))
{
	NOT_TESTABLE
	// has no members - this is useless
}
END_SECTION

START_SECTION((ProteinInference& operator=(const ProteinInference &rhs)))
{
	NOT_TESTABLE
	// has no members - this is useless
}
END_SECTION

START_SECTION((void infer(ConsensusMap &consensus_map, const UInt reference_map)))
{
	
  ConsensusXMLFile cm_file;
	ConsensusMap cm;
	cm_file.load(OPENMS_GET_TEST_DATA_PATH("ProteinInference.consensusXML"),cm);
	
	// delete quantitative info
	for (size_t i=0; i < cm.getProteinIdentifications()[0].getHits().size(); ++i)
	{
		cm.getProteinIdentifications()[0].getHits()[i].clearMetaInfo();
	}
	
	// this should create the quantitation that were in place before deleting them
	ProteinInference inferrer;
	inferrer.infer(cm, 0);
	
	String cm_file_out;// = OPENMS_GET_TEST_DATA_PATH("ItraqQuantifier.consensusXML");
	NEW_TMP_FILE(cm_file_out);
	cm_file.store(cm_file_out,cm);
	
	// TOLERANCE_ABSOLUTE(0.01);
	WHITELIST("<?xml-stylesheet");
	TEST_FILE_SIMILAR(cm_file_out,OPENMS_GET_TEST_DATA_PATH("ProteinInference.consensusXML"));
	
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



