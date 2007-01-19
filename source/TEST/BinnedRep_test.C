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

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/COMPARISON/CLUSTERING/BinnedRep.h>

///////////////////////////

#include <OpenMS/KERNEL/DSpectrum.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterSpectrum.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <iostream>

///////////////////////////
START_TEST(BinnedRep, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

BinnedRep* brp = 0;

double binsize = 0.5;
uint binspread = 2;

DTAFile dtafile;
PeakSpectrum spec;
dtafile.load("data/Transformers_tests.dta",spec);

CHECK(BinnedRep())
	brp = new BinnedRep();
	TEST_NOT_EQUAL(brp, 0)
RESULT

CHECK(BinnedRep(const double, const uint = 0))
	BinnedRep br(2.5, 1);
	TEST_REAL_EQUAL(br.getBinSize(), 2.5)
	TEST_EQUAL(br.getBinSpread(), 1)
RESULT

delete brp;
brp = new BinnedRep(spec, binsize, binspread);

CHECK(BinnedRep(const PeakSpectrum& spec, double binsize = 1.0, uint binspread = 0))
	BinnedRep br(spec);
	TEST_REAL_EQUAL(br.getBinSize(), 1.0)
	TEST_REAL_EQUAL(br.getBinSpread(), 0.0)
	BinnedRep br2(spec, 2.0, 1);
	TEST_REAL_EQUAL(br2.getBinSize(), 2.0)
	TEST_EQUAL(br2.getBinSpread(), 1)
RESULT

CHECK(BinnedRep& operator = (const BinnedRep& source))
  BinnedRep br = *brp;
  BinnedRep::const_iterator brit1 = brp->begin();
  BinnedRep::const_iterator brit2 = br.begin();
  while (brit2 != br.end())
  {
    TEST_REAL_EQUAL(*brit1++,*brit2++);
  }
  bool br2end = ( brit2 != br.end() );
  TEST_EQUAL(br2end,0)
RESULT

CHECK(BinnedRep(const BinnedRep& source))
  BinnedRep br(*brp);
  BinnedRep::const_iterator brit1 = brp->begin();
  BinnedRep::const_iterator brit2 = br.begin();
  while (brit2 != br.end())
  {
    TEST_REAL_EQUAL(*brit1++,*brit2++);
  }
  bool br2end = ( brit2 != br.end() );
  TEST_EQUAL(br2end,0)
RESULT

CHECK(uint id() const)
	TEST_EQUAL(brp->id(), 0)
RESULT

CHECK(double getBinSize() const)
	TEST_REAL_EQUAL(brp->getBinSize(), binsize)
RESULT

CHECK(double max() const)
	TEST_REAL_EQUAL(brp->max(), 1251.0)
RESULT

CHECK(double min() const)
	TEST_REAL_EQUAL(brp->min(), 104.0)
RESULT

CHECK(double getRetention() const)
	TEST_REAL_EQUAL(brp->getRetention(), -1.0)
RESULT

CHECK(unsigned int getBinSpread() const)
	TEST_EQUAL(brp->getBinSpread(), binspread)
RESULT

CHECK(double getParentmz() const)
	PRECISION(0.001)
	TEST_REAL_EQUAL(brp->getParentmz(), 370.385)
RESULT

CHECK(unsigned int getPrecursorPeakCharge() const)
	TEST_EQUAL(brp->getPrecursorPeakCharge(), 2)
RESULT

CHECK(String str() const)
	//TEST_EQUAL(brp->str(), "DASDAS")
	// TODO
RESULT

CHECK(intensity operator[] (int) const)
	
RESULT

CHECK((friend void operator << (BinnedRep& bin_rep, const PeakSpectrum& spec)))
	
RESULT

CHECK(void normalize())

RESULT

CHECK(unsigned int size() const)
	TEST_EQUAL(brp->size(), 2297)
RESULT

CHECK(ConstIterator begin() const)
	// TODO
RESULT

CHECK(ConstIterator end() const)
	// TODO
RESULT

CHECK(Iterator begin())
	// TODO
RESULT

CHECK(Iterator end())
	// TODO
RESULT

CHECK(~BinnedRep())
  delete brp;
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
