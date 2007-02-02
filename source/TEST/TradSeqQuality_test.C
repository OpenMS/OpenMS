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

///////////////////////////

#include <OpenMS/FILTERING/TRANSFORMERS/TradSeqQuality.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterSpectrum.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(TradSeqQuality, "$Id$")

/////////////////////////////////////////////////////////////

TradSeqQuality* e_ptr = 0;
CHECK(TradSeqQuality())
	e_ptr = new TradSeqQuality;
	TEST_NOT_EQUAL(e_ptr, 0)
RESULT

CHECK(~TradSeqQuality())
	delete e_ptr;
RESULT

e_ptr = new TradSeqQuality();

CHECK(TradSeqQuality(const TradSeqQuality& source))
	TradSeqQuality copy(*e_ptr);
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
RESULT

CHECK(TradSeqQuality& operator = (const TradSeqQuality& source))
	TradSeqQuality copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
RESULT

CHECK(double operator () (const ClusterSpectrum& spec))
	DTAFile dta_file;
	PeakSpectrum spec;
	dta_file.load("data/Transformers_tests.dta", spec);

	ClusterSpectrum cs(spec, 1, 1);
	
	double filter = (*e_ptr)(cs);
	TEST_REAL_EQUAL(filter, -1000)
RESULT

CHECK(static FilterFunctor* create())
	FilterFunctor* ff = TradSeqQuality::create();
	TradSeqQuality filter;
	TEST_EQUAL(ff->getParameters(), filter.getParameters())
	TEST_EQUAL(ff->getName(), filter.getName())
RESULT

CHECK(static const String getName())
	TEST_EQUAL(e_ptr->getName(), "TradSeqQuality")
RESULT

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
