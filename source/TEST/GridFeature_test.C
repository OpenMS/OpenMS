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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/GridFeature.h>

///////////////////////////

START_TEST(GridFeature, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

GridFeature* gf_ptr = 0;
GridFeature* gf_nullPointer = 0;

START_SECTION((GridFeature(const BaseFeature& feature, Size map_index, Size feature_index)))
{
	BaseFeature bf;
	gf_ptr = new GridFeature(bf, 0, 0);
  TEST_NOT_EQUAL(gf_ptr, gf_nullPointer);
}
END_SECTION

START_SECTION((~GridFeature()))
{
	delete gf_ptr;
}
END_SECTION

START_SECTION((const BaseFeature& getFeature() const))
{
	BaseFeature bf;
	bf.setRT(1.1);
	bf.setMZ(2.2);
	bf.setCharge(3);
	const BaseFeature bf_const(bf);
	GridFeature gf(bf_const, 0, 0);
	TEST_EQUAL(gf.getFeature() == bf_const, true);
}
END_SECTION

START_SECTION((Size getMapIndex() const))
{
	BaseFeature bf;
	GridFeature gf(bf, 123, 0);
	TEST_EQUAL(gf.getMapIndex(), 123);
}
END_SECTION

START_SECTION((Size getFeatureIndex() const))
{
	BaseFeature bf;
	GridFeature gf(bf, 0, 123);
	TEST_EQUAL(gf.getFeatureIndex(), 123);
}
END_SECTION

START_SECTION((Int getID() const))
{
	BaseFeature bf;
	GridFeature gf(bf, 0, 123);
	TEST_EQUAL(gf.getID(), 123);
}
END_SECTION

START_SECTION((const std::set<AASequence>& getAnnotations() const))
{
	BaseFeature bf;
	GridFeature gf(bf, 0, 0);
	TEST_EQUAL(gf.getAnnotations().size(), 0);
	bf.getPeptideIdentifications().resize(2);
	PeptideHit hit;
	hit.setSequence("AAA");
	bf.getPeptideIdentifications()[0].insertHit(hit);
	hit.setSequence("CCC");
	bf.getPeptideIdentifications()[1].insertHit(hit);
	GridFeature gf2(bf, 0, 0);
	TEST_EQUAL(gf2.getAnnotations().size(), 2);
	TEST_EQUAL(*(gf2.getAnnotations().begin()), "AAA");
	TEST_EQUAL(*(gf2.getAnnotations().rbegin()), "CCC");
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
