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

BinnedRep* brp;

double binsize = 0.5;
uint binspread = 2;

DTAFile dtafile;
MSSpectrum< DPeak<1> > spec;
dtafile.load("data/spectrum.dta",spec);

CHECK(BinnedRep::BinnedRep())
  brp = new BinnedRep(binsize,binspread);
  TEST_NOT_EQUAL(brp, 0)
RESULT

CHECK(operator<<(const BinnedRep&, const MSSpectrum< DPeak<1> >&))
  *brp << spec;
  TEST_REAL_EQUAL(spec.getRetentionTime(),brp->getRetention())
  TEST_REAL_EQUAL(spec.getPrecursorPeak().getCharge(),brp->getPrecursorPeakCharge())
  TEST_REAL_EQUAL(spec.getPrecursorPeak().getPosition()[0],brp->getParentmz())
  TEST_REAL_EQUAL(brp->getBinSize(),binsize)
  TEST_EQUAL(brp->getBinSpread(),binspread)
RESULT

CHECK(BinnedRep::operator=(const BinnedRep&))
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

CHECK(BinnedRep::BinnedRep(const BinnedRep&))
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

CHECK(BinnedRep::~BinnedRep())
  delete brp;
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
