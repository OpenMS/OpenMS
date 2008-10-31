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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/FORMAT/SemanticValidator.h>
#include <OpenMS/FORMAT/CVMappings.h>
#include <OpenMS/FORMAT/CVMappingFile.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>

///////////////////////////

START_TEST(SemanticValidator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

SemanticValidator* ptr = 0;
CHECK((SemanticValidator()))
	ptr = new SemanticValidator();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((virtual ~SemanticValidator()))
	delete ptr;
RESULT

CHECK(bool validate(const String& filename, const CVMappings& mapping, const ControlledVocabulary& cv))
//	CVMappings mapping;
//	CVMappingFile().load("../../share/OpenMS/MAPPING/ms-mapping.xml",mapping);
//	
//	ControlledVocabulary cv;
//	cv.loadFromOBO("PSI","../../share/OpenMS/CV/psi-ms.obo");
//	cv.loadFromOBO("PATO","../../share/OpenMS/CV/quality.obo");
//	cv.loadFromOBO("UO","../../share/OpenMS/CV/unit.obo");
//	cv.loadFromOBO("brenda","../../share/OpenMS/CV/brenda.obo");
//	cv.loadFromOBO("GO","../../share/OpenMS/CV/goslim_goa.obo");
//	
//	SemanticValidator sv;
//	TEST_EQUAL(sv.validate("data/MzMLFile_1.mzML", mapping, cv),true);
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
