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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/KERNEL/StandardTypes.h>

///////////////////////////

START_TEST(StandardTypes, "$Id: StandardTypes_test.C,v 1.5 2006/04/05 10:22:39 marc_sturm Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;


// RawData1D* drd1d_ptr = 0;
// CHECK(RawData1D())
// 	drd1d_ptr = new RawData1D;
// 	TEST_NOT_EQUAL(drd1d_ptr, 0)
// RESULT

// CHECK(~RawData1D())
// 	delete drd1d_ptr;
// RESULT

CHECK(StandardTypes)



RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


