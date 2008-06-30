// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Mathias Walzer$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrum.h>
#include <OpenMS/FORMAT/DTAFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(BinnedSpectrum, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

BinnedSpectrum* ptr = 0;
CHECK(BinnedSpectrum())
{
	ptr = new BinnedSpectrum();
	TEST_NOT_EQUAL(ptr, 0)
  TEST_EXCEPTION(BinnedSpectrum::NoSpectrumIntegrated,ptr->setBinning();)
}
RESULT

CHECK(~BinnedSpectrum())
{
	delete ptr;
}
RESULT

  BinnedSpectrum* bs1;
  DTAFile dtafile;
  PeakSpectrum s1;
  DTAFile().load("data/PILISSequenceDB_DFPIANGER_1.dta", s1);

CHECK((BinnedSpectrum(Real size, UInt spread, PeakSpectrum ps)))
{
  bs1 = new BinnedSpectrum(1.5,2,s1);
  TEST_NOT_EQUAL(bs1,0)
}
RESULT

CHECK((BinnedSpectrum(const BinnedSpectrum &source)))
{
  BinnedSpectrum copy(*bs1);
  TEST_EQUAL(copy.getName(), bs1->getName());
  TEST_EQUAL(copy.getBinSize(), bs1->getBinSize());
  TEST_EQUAL((UInt)copy.getPrecursorPeak().getPosition()[0],(UInt)bs1->getPrecursorPeak().getPosition()[0]);
}
RESULT

CHECK((BinnedSpectrum& operator=(const BinnedSpectrum &source)))
{
  BinnedSpectrum copy(*bs1);
  bs1 = new BinnedSpectrum(1.5,2,s1);
  TEST_EQUAL(copy.getName(), bs1->getName());
  TEST_EQUAL(copy.getBinSize(), bs1->getBinSize());
  TEST_EQUAL((UInt)copy.getPrecursorPeak().getPosition()[0],(UInt)bs1->getPrecursorPeak().getPosition()[0]);
}
RESULT

CHECK((BinnedSpectrum& operator=(const PeakSpectrum &source)))
{
  bs1 = new BinnedSpectrum();
  *bs1 = s1;
  TEST_EQUAL(bs1->getPrecursorPeak().getPosition()[0],s1.getPrecursorPeak().getPosition()[0]);
  bs1->setBinSize(1.5);
  bs1->setBinSpread(2);
}
RESULT

CHECK((bool operator==(const BinnedSpectrum &rhs) const ))
{
	BinnedSpectrum copy = *bs1;
	TEST_EQUAL((*bs1==copy),true)
}
RESULT

CHECK((bool operator!=(const BinnedSpectrum &rhs) const ))
{
	BinnedSpectrum copy = *bs1;
	TEST_EQUAL((*bs1!=copy),false)
}
RESULT

CHECK((bool operator==(const PeakSpectrum &rhs) const ))
{
	TEST_EQUAL((*bs1==s1),true)
}
RESULT

CHECK((bool operator!=(const PeakSpectrum &rhs) const ))
{
	BinnedSpectrum copy = *bs1;
	TEST_EQUAL((*bs1!=s1),false)
}
RESULT

CHECK((double getBinSize() const ))
{
	TEST_EQUAL(bs1->getBinSize(),1.5)
}
RESULT

CHECK((UInt getBinSpread() const ))
{
	TEST_EQUAL(bs1->getBinSpread(),2)
}
RESULT

CHECK((UInt getBinNumber() const ))
{
	TEST_EQUAL(bs1->getBinNumber(),659)
}
RESULT

CHECK((UInt getFilledBinNumber() const ))
{
	TEST_EQUAL(bs1->getFilledBinNumber(),347)
}
RESULT

CHECK((const SparseVector<Real>& getBins() const  throw (NoSpectrumIntegrated)))
{
	TEST_EQUAL(bs1->getBins().at(658),501645)
}
RESULT

CHECK((SparseVector<Real>& getBins() throw (NoSpectrumIntegrated)))
{
	TEST_EQUAL(bs1->getBins().at(658),501645)
}
RESULT

CHECK((const_bin_iterator begin() const ))
{
	UInt c(0);
	for (BinnedSpectrum::const_bin_iterator it1 = bs1->begin(); it1 != bs1->end(); ++it1)
	{
		++c;
	}
	TEST_EQUAL(bs1->getBinNumber(),c)
}
RESULT

CHECK((const_bin_iterator end() const ))
{
	NOT_TESTABLE
	//tested above
}
RESULT

CHECK((bin_iterator begin()))
{
	UInt c(0);
	for (BinnedSpectrum::bin_iterator it1 = bs1->begin(); it1 != bs1->end(); ++it1)
	{
		++c;
	}
	TEST_EQUAL(bs1->getBinNumber(),c)
}
RESULT

CHECK((bin_iterator end()))
{
  NOT_TESTABLE
	//tested above
}
RESULT

CHECK((void setBinSize(double s) throw (NoSpectrumIntegrated)))
{
	TEST_EQUAL(bs1->getBinSize(),1.5)
}
RESULT

CHECK((void setBinSpread(UInt s) throw (NoSpectrumIntegrated)))
{
	TEST_EQUAL(bs1->getBinSpread(),2)
}
RESULT

CHECK((void setBinning() throw (NoSpectrumIntegrated)))
{
	NOT_TESTABLE
	//tested within another test
}
RESULT

CHECK((bool checkCompliance(const BinnedSpectrum &bs) const ))
{
	TEST_EQUAL(bs1->checkCompliance(BinnedSpectrum()),false)
}
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



