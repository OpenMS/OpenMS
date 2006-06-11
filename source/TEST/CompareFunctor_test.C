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
// $Id: CompareFunctor_test.C,v 1.12 2006/06/09 23:47:35 nicopfeifer Exp $
// $Author: nicopfeifer $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterSpectrum.h>
#include <OpenMS/METADATA/Identification.h>
///////////////////////////

#include <OpenMS/COMPARISON/CLUSTERING/ClusterFactory.h>
#include <OpenMS/COMPARISON/SPECTRA/CompareFunctor.h>

///////////////////////////

#include <vector>
#include <iostream>

///////////////////////////
START_TEST(CompareFunctor, "$Id: CompareFunctor_test.C,v 1.12 2006/06/09 23:47:35 nicopfeifer Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

CompareFunctor* cfp;

ClusterFactory* factoryp = ClusterFactory::instance();

DTAFile dtafile;

MSSpectrum< DPeak<1> >* spec = new MSSpectrum< DPeak<1> >();
MSSpectrum< DPeak<1> >* spec2 = new MSSpectrum< DPeak<1> >();
MSSpectrum< DPeak<1> >* spec3 = new MSSpectrum< DPeak<1> >();
dtafile.load("data/spectrum.dta",*spec);
dtafile.load("data/spectrum2.dta",*spec2);
dtafile.load("data/spectrum2.dta",*spec3);

Identification dbs;
dbs.insertPeptideHit(PeptideHit(27.0,"Mascot",1,"RRYA"));
spec->getIdentification().push_back(dbs);
spec2->getIdentification().push_back(dbs);
ClusterSpectrum cspec(spec,0,0.5,2);
ClusterSpectrum cspec2(spec2,0,0.5,2);
ClusterSpectrum cspec3(spec3,0,1,1);

vector<String> catalogue = factoryp->catalogue("CompareFunctor");

// go through all registered FilterFunctors and check if they accept a spectrum 
// and return something

for ( vector<String>::const_iterator cvit = catalogue.begin();
    cvit != catalogue.end(); ++cvit )
{
  CHECK()
    STATUS("ClusterFactory::create("+*cvit+")")
    cfp = dynamic_cast<CompareFunctor*>(factoryp->create(*cvit));
    TEST_NOT_EQUAL(cfp, 0)
  RESULT


  CHECK(CompareFunctor::operator())
    STATUS(*cvit+"::operator()")
    double result = (*cfp)(cspec,cspec2);
    bool ispositive = result >= 0;
    TEST_EQUAL(ispositive, true)
  RESULT


  CHECK(CompareFunctor::WrongRepresentation)
    if ( cfp->usebins() )
    {
      TEST_EXCEPTION(ClusterSpectrum::WrongRepresentation,(*cfp)(cspec,cspec3))
    }
  RESULT
  
  delete cfp;

}

factoryp->destroy();

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
