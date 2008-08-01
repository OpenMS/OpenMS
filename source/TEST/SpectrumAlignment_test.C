// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <iostream>

#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>
#include <OpenMS/FORMAT/DTAFile.h>

///////////////////////////

START_TEST(SpectrumAlignment, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

SpectrumAlignment* ptr = 0;

CHECK(SpectrumAlignment())
	ptr = new SpectrumAlignment();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(virtual ~SpectrumAlignment())
	delete ptr;
RESULT

ptr = new SpectrumAlignment();

CHECK(SpectrumAlignment(const SpectrumAlignment &source))
  SpectrumAlignment sas1;
  Param p(sas1.getParameters());
  p.setValue("tolerance", 0.2);
  sas1.setParameters(p);

  SpectrumAlignment sas2(sas1);

  TEST_EQUAL(sas1.getName(), sas2.getName())
  TEST_EQUAL(sas1.getParameters(), sas2.getParameters())
RESULT

CHECK(SpectrumAlignment& operator=(const SpectrumAlignment &source))
  SpectrumAlignment sas1;
  Param p(sas1.getParameters());
  p.setValue("tolerance", 0.2);
  sas1.setParameters(p);

  SpectrumAlignment sas2;

	sas2 = sas1;

  TEST_EQUAL(sas1.getName(), sas2.getName())
  TEST_EQUAL(sas1.getParameters(), sas2.getParameters())

RESULT
		    
CHECK(template <typename SpectrumType> void getSpectrumAlignment(std::vector< std::pair< UInt, UInt > > &alignment, const SpectrumType &s1, const SpectrumType &s2) const)
	PeakSpectrum s1, s2;
  DTAFile().load("data/PILISSequenceDB_DFPIANGER_1.dta", s1);
  DTAFile().load("data/PILISSequenceDB_DFPIANGER_1.dta", s2);

  PRECISION(0.01)

	SpectrumAlignment sas1;
	vector<pair<UInt, UInt > > alignment;
	sas1.getSpectrumAlignment(alignment, s1, s2);
	
	for (vector<pair<UInt, UInt > >::const_iterator it = alignment.begin(); it != alignment.end(); ++it)
	{
		cerr << it->first << " " << it->second << endl;
	}


  TEST_REAL_EQUAL(alignment.size(), s1.size())

  s2.getContainer().resize(100);

	alignment.clear();
  sas1.getSpectrumAlignment(alignment, s1, s2);

  TEST_REAL_EQUAL(alignment.size(), 100)

RESULT

ptr = new SpectrumAlignment();

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
