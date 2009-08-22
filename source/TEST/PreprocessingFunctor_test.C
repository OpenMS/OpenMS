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
// $Maintainer: Andreas Bertsch $
// $Authors: Volker Mosthaf, Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Factory.h>

///////////////////////////

#include <OpenMS/FILTERING/TRANSFORMERS/PreprocessingFunctor.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(PreprocessingFunctor, "$Id$")

/////////////////////////////////////////////////////////////

// most of the tests cannot be performed, as this is a pure interface class

START_SECTION((PreprocessingFunctor()))
	NOT_TESTABLE
END_SECTION

START_SECTION((~PreprocessingFunctor()))
	NOT_TESTABLE
END_SECTION

START_SECTION((PreprocessingFunctor(const PreprocessingFunctor& source)))
	NOT_TESTABLE
END_SECTION

START_SECTION((PreprocessingFunctor& operator = (const PreprocessingFunctor& source)))
	NOT_TESTABLE
END_SECTION

START_SECTION((void filterPeakMap(PeakMap& exp)))
	NOT_TESTABLE
END_SECTION

START_SECTION((void filterPeakSpectrum(PeakSpectrum& spectrum)))
	NOT_TESTABLE
END_SECTION

START_SECTION(template<typename SpectrumType> void filterSpectrum(SpectrumType&))
	NOT_TESTABLE
END_SECTION

START_SECTION((static void registerChildren()))
	PreprocessingFunctor* ppf = Factory<PreprocessingFunctor>::create("ThresholdMower");
	TEST_EQUAL(ppf->getName(), "ThresholdMower");
	ppf = Factory<PreprocessingFunctor>::create("WindowMower");
	TEST_EQUAL(ppf->getName(), "WindowMower");
	ppf = Factory<PreprocessingFunctor>::create("Scaler");
	TEST_EQUAL(ppf->getName(), "Scaler");
	ppf = Factory<PreprocessingFunctor>::create("NLargest");
	TEST_EQUAL(ppf->getName(), "NLargest");
	ppf = Factory<PreprocessingFunctor>::create("MarkerMower");
	TEST_EQUAL(ppf->getName(), "MarkerMower");
	ppf = Factory<PreprocessingFunctor>::create("SqrtMower");
	TEST_EQUAL(ppf->getName(), "SqrtMower");
	ppf = Factory<PreprocessingFunctor>::create("Normalizer");
	TEST_EQUAL(ppf->getName(), "Normalizer");
	ppf = Factory<PreprocessingFunctor>::create("ParentPeakMower");
	TEST_EQUAL(ppf->getName(), "ParentPeakMower");	
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
