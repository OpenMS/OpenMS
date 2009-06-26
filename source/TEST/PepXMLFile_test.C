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

START_TEST(PepXMLFile, "$Id: PepXMLFile_test.C 4776 2009-03-05 14:14:35Z groepl $")

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

START_SECTION(void load(const String& filename,  std::map<String, std::vector<AASequence> >& peptides))
	//TODO Hendrik
END_SECTION

START_SECTION(void store(const String& filename, std::vector<ProteinIdentification>& protein_ids, std::vector<PeptideIdentification>& peptide_ids))
	//TODO Chris
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

