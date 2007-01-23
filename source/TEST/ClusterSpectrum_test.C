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
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterSpectrum.h>
#include <OpenMS/METADATA/Identification.h>
///////////////////////////

#include <OpenMS/COMPARISON/CLUSTERING/ClusterSpectrum.h>

///////////////////////////

///////////////////////////
START_TEST(ClusterSpectrum, "$Id ClusterSpectrum_test.C,v 1.3 2005/02/21 20:00:59 fukuryu Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

ClusterSpectrum* cspec;

DTAFile dtafile;

MSSpectrum< DPeak<1> >* spec = new MSSpectrum< DPeak<1> >();
MSSpectrum< DPeak<1> >* spec2 = new MSSpectrum< DPeak<1> >();
dtafile.load("data/Transformers_tests.dta",*spec);
dtafile.load("data/Transformers_tests.dta",*spec2);

Identification dbs;
dbs.insertPeptideHit(PeptideHit(27.0,"Mascot",1, 1, "RRYA"));
spec->getIdentifications().push_back(dbs);

CHECK(ClusterSpectrum::ClusterSpectrum())
  cspec = new ClusterSpectrum();
  TEST_NOT_EQUAL(cspec,0)
RESULT

CHECK(ClusterSpectrum::~ClusterSpectrum())
  delete cspec;
RESULT

CHECK(ClusterSpectrum::ClusterSpectrum(MSSpectrum< DPeak<1> >*))
  cspec = new ClusterSpectrum(spec,0.5,2);
  TEST_NOT_EQUAL(cspec,0)
RESULT

CHECK(ClusterSpectrum::ClusterSpectrum(const ClusterSpectrum& source))
  ClusterSpectrum* cspec2 = new ClusterSpectrum(*cspec);
  delete cspec;
  cspec2->getSpec();
  cspec2->getBinrep();
  cspec = cspec2;
RESULT

CHECK(ClusterSpectrum::ClusterSpectrum::operator=(const ClusterSpectrum& source))
  ClusterSpectrum* cspec2 = new ClusterSpectrum();
  *cspec2 = *cspec;
  delete cspec;
  cspec2->getSpec();
  cspec2->getBinrep();
  cspec = cspec2;
RESULT

CHECK(ClusterSpectrum::getParentionCharge())
  TEST_EQUAL(cspec->getParentionCharge(),UnsignedInt(spec2->getPrecursorPeak().getCharge()))
RESULT

CHECK(ClusterSpectrum::getParentMass())
  TEST_EQUAL(cspec->getParentMass(),spec2->getPrecursorPeak().getPosition()[0])
RESULT

CHECK(ClusterSpectrum::getTophit())
	std::cout << "Tophit: " << cspec->getTophit().getSequence() << std::endl;
  TEST_EQUAL(cspec->getTophit().getSequence(),"RRYA")
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
