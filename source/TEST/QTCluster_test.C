// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/DATASTRUCTURES/QTCluster.h>

///////////////////////////

START_TEST(QTCluster, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

QTCluster* qtc_ptr = 0;
QTCluster* qtc_nullPointer = 0;

BaseFeature bf;
bf.setRT(1.1);
bf.setMZ(2.2);
bf.setCharge(3);
bf.getPeptideIdentifications().resize(2);
PeptideHit hit;
hit.setSequence("AAA");
bf.getPeptideIdentifications()[0].insertHit(hit);
hit.setSequence("CCC");
bf.getPeptideIdentifications()[1].insertHit(hit);
GridFeature gf(bf, 123, 456);

START_SECTION((QTCluster(GridFeature* center_point, Size num_maps, DoubleReal max_distance, bool use_IDs)))
{
	qtc_ptr = new QTCluster(&gf, 2, 11.1, false);
  TEST_NOT_EQUAL(qtc_ptr, qtc_nullPointer);
}
END_SECTION

START_SECTION((~QTCluster()))
{
	delete qtc_ptr;
}
END_SECTION

QTCluster cluster(&gf, 2, 11.1, true);

START_SECTION((DoubleReal getCenterRT() const))
{
	TEST_EQUAL(cluster.getCenterRT(), 1.1);
}
END_SECTION

START_SECTION((DoubleReal getCenterMZ() const))
{
	TEST_EQUAL(cluster.getCenterMZ(), 2.2);
}
END_SECTION

START_SECTION((Size size() const))
{
	TEST_EQUAL(cluster.size(), 1);
}
END_SECTION

GridFeature gf2(bf, 789, 1011);

START_SECTION((void add(GridFeature* element, DoubleReal distance)))
{
	cluster.add(&gf2, 3.3);
	TEST_EQUAL(cluster.size(), 2);
}
END_SECTION

START_SECTION((bool operator<(QTCluster& cluster)))
{
	QTCluster cluster2(&gf, 2, 11.1, false);
	TEST_EQUAL(cluster2 < cluster, true);
}
END_SECTION

START_SECTION((void getElements(std::map<Size, GridFeature*>& elements)))
{
	map<Size, GridFeature*> elements;
	cluster.getElements(elements);
	TEST_EQUAL(elements.size(), 2);
	TEST_EQUAL(elements[123], &gf);
	TEST_EQUAL(elements[789], &gf2);
}
END_SECTION

START_SECTION((bool update(const std::map<Size, GridFeature*>& removed)))
{
	map<Size, GridFeature*> removed;
	removed[789] = &gf2;
	TEST_EQUAL(cluster.update(removed), true);
	TEST_EQUAL(cluster.size(), 1);
	removed[123] = &gf;
	// removing the center invalidates the cluster:
	TEST_EQUAL(cluster.update(removed), false);
}
END_SECTION

START_SECTION((DoubleReal getQuality()))
{
	TEST_EQUAL(cluster.getQuality(), 0.0);
	cluster.add(&gf2, 3.3);
	TEST_EQUAL(cluster.getQuality(), (11.1 - 3.3) / 11.1);
}
END_SECTION

START_SECTION((const set<AASequence>& getAnnotations()))
{
	TEST_EQUAL(cluster.getAnnotations().size(), 2);
	TEST_EQUAL(*(cluster.getAnnotations().begin()), "AAA");
	TEST_EQUAL(*(cluster.getAnnotations().rbegin()), "CCC");
	QTCluster cluster2(&gf, 2, 11.1, false);
	TEST_EQUAL(cluster2.getAnnotations().empty(), true);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
