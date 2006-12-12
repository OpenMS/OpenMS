// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/COMPARISON/SPECTRA/PeakSpectrumCompareFunctor.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(PeakSpectrumCompareFunctor, "$Id:$")

/////////////////////////////////////////////////////////////

PeakSpectrumCompareFunctor* e_ptr = 0;
CHECK(PeakSpectrumCompareFunctor())
	e_ptr = new PeakSpectrumCompareFunctor;
	TEST_NOT_EQUAL(e_ptr, 0)
RESULT

CHECK(~PeakSpectrumCompareFunctor())
	delete e_ptr;
RESULT

e_ptr = new PeakSpectrumCompareFunctor();

CHECK(PeakSpectrumCompareFunctor(const PeakSpectrumCompareFunctor& source))
	PeakSpectrumCompareFunctor copy(*e_ptr);
	TEST_EQUAL(*e_ptr == copy, true)
RESULT

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
