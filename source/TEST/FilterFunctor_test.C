// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
///////////////////////////

#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>

///////////////////////////

#include <vector>
#include <iostream>

///////////////////////////
START_TEST(FilterFunctor, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

// pure interface, hardly testable

START_SECTION(FilterFunctor())
	NOT_TESTABLE
END_SECTION

START_SECTION(FilterFunctor(const FilterFunctor& source))
	NOT_TESTABLE
END_SECTION

START_SECTION(FilterFunctor& operator = (const FilterFunctor& source))
	NOT_TESTABLE
END_SECTION

START_SECTION(static void registerChildren())
	FilterFunctor* ff = Factory<FilterFunctor>::create("ComplementFilter");
	TEST_EQUAL(ff->getName(), "ComplementFilter")
	ff = Factory<FilterFunctor>::create("IntensityBalanceFilter");
	TEST_EQUAL(ff->getName(), "IntensityBalanceFilter")	
	ff = Factory<FilterFunctor>::create("NeutralLossDiffFilter");
	TEST_EQUAL(ff->getName(), "NeutralLossDiffFilter")
	ff = Factory<FilterFunctor>::create("IsotopeDiffFilter");
	TEST_EQUAL(ff->getName(), "IsotopeDiffFilter")
	ff = Factory<FilterFunctor>::create("TICFilter");
	TEST_EQUAL(ff->getName(), "TICFilter")
END_SECTION

START_SECTION(template<typename SpectrumType> double apply(SpectrumType&))
	NOT_TESTABLE
END_SECTION

START_SECTION(~FilterFunctor())
	NOT_TESTABLE
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
