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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/KERNEL/DimensionDescription.h>

///////////////////////////

template < unsigned int DIMENSION >
struct GetDim
{
  unsigned int getDim ()
  { return DIMENSION; }
};
template < unsigned int DIMENSION >
unsigned int getDim ()
{ return DIMENSION; }

///////////////////////////

START_TEST(DimensionDescription<LCMS_Tag>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

typedef DimensionDescription<LCMS_Tag> DDLCMS;

DDLCMS * ddlcms_ptr = 0;

STATUS("DDLCMS means DimensionDescription<LCMS_Tag>")

CHECK(DDLCMS())
  ddlcms_ptr = new DDLCMS;
  TEST_NOT_EQUAL(ddlcms_ptr,0)
RESULT

CHECK(~DDLCMS())
  delete ddlcms_ptr;
RESULT

CHECK(DDLCMS::DimensionId)
  DDLCMS::DimensionId dim;
  dim = DDLCMS::RT;
  TEST_EQUAL(dim,0)
  dim = DDLCMS::MZ;
  TEST_EQUAL(dim,1)
  dim = DDLCMS::DIMENSION;
  TEST_EQUAL(dim,2)
RESULT

CHECK("Using enums RT and MZ als template arguments")
  TEST_EQUAL(GetDim<DDLCMS::RT>().getDim(),DDLCMS::RT)
  TEST_EQUAL(getDim<DDLCMS::RT>(),DDLCMS::RT)
  TEST_EQUAL(GetDim<DDLCMS::MZ>().getDim(),DDLCMS::MZ)
  TEST_EQUAL(getDim<DDLCMS::MZ>(),DDLCMS::MZ)
  TEST_EQUAL(GetDim<DDLCMS::DIMENSION>().getDim(),DDLCMS::DIMENSION)
  TEST_EQUAL(getDim<DDLCMS::DIMENSION>(),DDLCMS::DIMENSION)
RESULT

  // more convenient
  enum { RT=DDLCMS::RT, MZ=DDLCMS::MZ };

STATUS("From now on, RT means DDLCMS::RT, etc.")

CHECK(DDLCMS::dimension_name_short[])
  TEST_STRING_EQUAL(DDLCMS::dimension_name_short[RT],"RT")
  TEST_STRING_EQUAL(DDLCMS::dimension_name_short[MZ],"MZ")
RESULT

CHECK(DDLCMS::dimension_name_full[])
  TEST_STRING_EQUAL(DDLCMS::dimension_name_full[RT],"retention time")
  TEST_STRING_EQUAL(DDLCMS::dimension_name_full[MZ],"mass-to-charge")
RESULT

CHECK(DDLCMS::dimension_unit_short[])
  TEST_STRING_EQUAL(DDLCMS::dimension_unit_short[RT],"sec")
  TEST_STRING_EQUAL(DDLCMS::dimension_unit_short[MZ],"Th")
RESULT

CHECK(DDLCMS::dimension_unit_full[])
  TEST_STRING_EQUAL(DDLCMS::dimension_unit_full[RT],"Seconds")
  TEST_STRING_EQUAL(DDLCMS::dimension_unit_full[MZ],"Thomson")
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
