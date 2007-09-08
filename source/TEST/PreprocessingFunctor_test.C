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

CHECK((PreprocessingFunctor()))
	// nothing to check
RESULT

CHECK((~PreprocessingFunctor()))
	// nothing to check
RESULT

CHECK((PreprocessingFunctor(const PreprocessingFunctor& source)))
	// nothing to check
RESULT

CHECK((PreprocessingFunctor& operator = (const PreprocessingFunctor& source)))
	// nothing to check
RESULT

CHECK((void filterPeakMap(PeakMap& exp)))
	// nothing to check
RESULT

CHECK((void filterPeakSpectrum(PeakSpectrum& spectrum)))
	// nothing to check
RESULT

CHECK(template<typename SpectrumType> void filterSpectrum(SpectrumType&))
	// nothing to check
RESULT

CHECK(static const String getProductName())
	TEST_EQUAL(PreprocessingFunctor::getProductName(), "PreprocessingFunctor")
RESULT

CHECK((static void registerChildren()))
	PreprocessingFunctor* ppf = Factory<PreprocessingFunctor>::create("ThresholdMower");
	ppf = Factory<PreprocessingFunctor>::create("WindowMower");
	ppf = Factory<PreprocessingFunctor>::create("Scaler");
	ppf = Factory<PreprocessingFunctor>::create("NLargest");
	ppf = Factory<PreprocessingFunctor>::create("MarkerMower");
	ppf = Factory<PreprocessingFunctor>::create("SqrtMower");
	ppf = Factory<PreprocessingFunctor>::create("Normalizer");
	ppf = Factory<PreprocessingFunctor>::create("ParentPeakMower");

RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
