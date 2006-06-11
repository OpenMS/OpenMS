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
// $Id: PeakMarker_test.C,v 1.8 2006/04/05 11:18:25 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterSpectrum.h>
///////////////////////////

#include <OpenMS/COMPARISON/CLUSTERING/ClusterFactory.h>
#include <OpenMS/FILTERING/TRANSFORMERS/PeakMarker.h>

///////////////////////////

///////////////////////////
START_TEST(PeakMarker, "$Id: PeakMarker_test.C,v 1.8 2006/04/05 11:18:25 marc_sturm Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

PeakMarker* pmp;

ClusterFactory* factoryp = ClusterFactory::instance();

DTAFile dtafile;

MSSpectrum< DPeak<1> >* spec = new MSSpectrum< DPeak<1> >();
dtafile.load("data/spectrum.dta",*spec);
ClusterSpectrum cspec(spec);

vector<String> catalogue = factoryp->catalogue("PeakMarker");

// go through all registered FilterFunctors and check if they accept a spectrum 
// and return something
for ( vector<String>::const_iterator cvit = catalogue.begin();
    cvit != catalogue.end(); ++cvit )
{
  CHECK()
    STATUS("ClusterFactory::create("+*cvit+")")
    pmp = dynamic_cast<PeakMarker*>(factoryp->create(*cvit));
    TEST_NOT_EQUAL(pmp, 0)
  RESULT

  CHECK(PeakMarker::operator())
    STATUS(*cvit+"::operator()")
    map<double,bool> result = (*pmp)(cspec.spec());
    TEST_NOT_EQUAL(result.size(),0)
  RESULT
  
  CHECK(PeakMarker::~PeakMarker())
    delete pmp;
  RESULT
}

factoryp->destroy();

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
