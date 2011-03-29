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
// $Maintainer: Mathias Walzer$
// $Authors: $
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
BinnedSpectrum* nullPointer = 0;
START_SECTION(BinnedSpectrum())
{
	ptr = new BinnedSpectrum();
	TEST_NOT_EQUAL(ptr, nullPointer)
  TEST_EXCEPTION(BinnedSpectrum::NoSpectrumIntegrated,ptr->setBinning();)
}
END_SECTION

START_SECTION(~BinnedSpectrum())
{
	delete ptr;
}
END_SECTION

  BinnedSpectrum* bs1;
  DTAFile dtafile;
  PeakSpectrum s1;
  DTAFile().load(OPENMS_GET_TEST_DATA_PATH("PILISSequenceDB_DFPIANGER_1.dta"), s1);

START_SECTION((BinnedSpectrum(Real size, UInt spread, PeakSpectrum ps)))
{
  bs1 = new BinnedSpectrum(1.5,2,s1);
  TEST_NOT_EQUAL(bs1,nullPointer)
}
END_SECTION

START_SECTION((BinnedSpectrum(const BinnedSpectrum &source)))
{
  BinnedSpectrum copy(*bs1);
  TEST_EQUAL(copy.getName(), bs1->getName());
  TEST_EQUAL(copy.getBinSize(), bs1->getBinSize());
  TEST_EQUAL((UInt)copy.getPrecursors()[0].getMZ(),(UInt)bs1->getPrecursors()[0].getMZ());
}
END_SECTION

START_SECTION((BinnedSpectrum& operator=(const BinnedSpectrum &source)))
{
  BinnedSpectrum copy(*bs1);
  bs1 = new BinnedSpectrum(1.5,2,s1);
  TEST_EQUAL(copy.getName(), bs1->getName());
  TEST_EQUAL(copy.getBinSize(), bs1->getBinSize());
  TEST_EQUAL((UInt)copy.getPrecursors()[0].getMZ(),(UInt)bs1->getPrecursors()[0].getMZ());
}
END_SECTION

START_SECTION((BinnedSpectrum& operator=(const PeakSpectrum &source)))
{
  bs1 = new BinnedSpectrum();
  *bs1 = s1;
  TEST_EQUAL(bs1->getPrecursors()[0].getMZ(),s1.getPrecursors()[0].getMZ());
  bs1->setBinSize(1.5);
  bs1->setBinSpread(2);
}
END_SECTION

START_SECTION((bool operator==(const BinnedSpectrum &rhs) const ))
{
	BinnedSpectrum copy = *bs1;
	TEST_EQUAL((*bs1==copy),true)
}
END_SECTION

START_SECTION((bool operator!=(const BinnedSpectrum &rhs) const ))
{
	BinnedSpectrum copy = *bs1;
	TEST_EQUAL((*bs1!=copy),false)
}
END_SECTION

START_SECTION((bool operator==(const PeakSpectrum &rhs) const ))
{
	TEST_EQUAL((*bs1==s1),true)
}
END_SECTION

START_SECTION((bool operator!=(const PeakSpectrum &rhs) const ))
{
	BinnedSpectrum copy = *bs1;
	TEST_EQUAL((*bs1!=s1),false)
}
END_SECTION

START_SECTION((double getBinSize() const ))
{
	TEST_EQUAL(bs1->getBinSize(),1.5)
}
END_SECTION

START_SECTION((UInt getBinSpread() const ))
{
	TEST_EQUAL(bs1->getBinSpread(),2)
}
END_SECTION

START_SECTION((UInt getBinNumber() const ))
{
	TEST_EQUAL(bs1->getBinNumber(),659)
}
END_SECTION

START_SECTION((UInt getFilledBinNumber() const ))
{
	TEST_EQUAL(bs1->getFilledBinNumber(),347)
}
END_SECTION

START_SECTION((const SparseVector<Real>& getBins() const))
{
	TEST_EQUAL(bs1->getBins().at(658),501645)
}
END_SECTION

START_SECTION((SparseVector<Real>& getBins()))
{
	TEST_EQUAL(bs1->getBins().at(658),501645)
}
END_SECTION

START_SECTION((const_bin_iterator begin() const ))
{
	UInt c(0);
	for (BinnedSpectrum::const_bin_iterator it1 = bs1->begin(); it1 != bs1->end(); ++it1)
	{
		++c;
	}
	TEST_EQUAL(bs1->getBinNumber(),c)
}
END_SECTION

START_SECTION((const_bin_iterator end() const ))
{
	NOT_TESTABLE
	//tested above
}
END_SECTION

START_SECTION((bin_iterator begin()))
{
	UInt c(0);
	for (BinnedSpectrum::bin_iterator it1 = bs1->begin(); it1 != bs1->end(); ++it1)
	{
		++c;
	}
	TEST_EQUAL(bs1->getBinNumber(),c)
}
END_SECTION

START_SECTION((bin_iterator end()))
{
  NOT_TESTABLE
	//tested above
}
END_SECTION

START_SECTION((void setBinSize(double s)))
{
	TEST_EQUAL(bs1->getBinSize(),1.5)
}
END_SECTION

START_SECTION((void setBinSpread(UInt s)))
{
	TEST_EQUAL(bs1->getBinSpread(),2)
}
END_SECTION

START_SECTION((void setBinning()))
{
	NOT_TESTABLE
	//tested within another test
}
END_SECTION

START_SECTION((bool checkCompliance(const BinnedSpectrum &bs) const ))
{
	TEST_EQUAL(bs1->checkCompliance(BinnedSpectrum()),false)
}
END_SECTION

START_SECTION(([BinnedSpectrum::NoSpectrumIntegrated] NoSpectrumIntegrated(const char *file, int line, const char *function, const char *message="BinnedSpectrum hasn't got a PeakSpectrum to base on yet")))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION(([BinnedSpectrum::NoSpectrumIntegrated] virtual ~NoSpectrumIntegrated()))
{
  NOT_TESTABLE
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



