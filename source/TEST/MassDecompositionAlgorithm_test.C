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
// $Maintainer: Sandro Andreotti $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/MassDecompositionAlgorithm.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MassDecompositionAlgorithm, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MassDecompositionAlgorithm* ptr = 0;
MassDecompositionAlgorithm* nullPointer = 0;
START_SECTION(MassDecompositionAlgorithm())
{
	ptr = new MassDecompositionAlgorithm();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(virtual ~MassDecompositionAlgorithm())
{
	delete ptr;
}
END_SECTION

START_SECTION((void getDecompositions(std::vector<MassDecomposition>& decomps, DoubleReal weight)))
{
  vector<MassDecomposition> decomps;
	DoubleReal mass = AASequence("DFPIANGER").getMonoWeight(Residue::Internal);
	cerr << mass << endl;

	MassDecompositionAlgorithm mda;
	Param p(mda.getParameters());
	p.setValue("tolerance", 0.0001);
	mda.setParameters(p);

	mda.getDecompositions(decomps, mass);
	TEST_EQUAL(decomps.size(), 842)

	p.setValue("tolerance", 0.001);
	mda.setParameters(p);
	decomps.clear();
	mda.getDecompositions(decomps, mass);
	TEST_EQUAL(decomps.size(), 911);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



