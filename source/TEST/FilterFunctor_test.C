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
// $Id: FilterFunctor_test.C,v 1.13 2006/06/09 23:47:35 nicopfeifer Exp $
// $Author: nicopfeifer $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/Identification.h>
///////////////////////////

#include <OpenMS/COMPARISON/CLUSTERING/ClusterFactory.h>
#include <OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>

///////////////////////////

#include <vector>
#include <iostream>

///////////////////////////
START_TEST(FilterFunctor, "$Id: FilterFunctor_test.C,v 1.13 2006/06/09 23:47:35 nicopfeifer Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

FilterFunctor* ffp;

ClusterFactory* factoryp = ClusterFactory::instance();

DTAFile dtafile;

MSSpectrum<DPeak<1> >* spec = new MSSpectrum<DPeak<1> >();
dtafile.load("data/spectrum.dta",*spec);
Identification dbs;
dbs.insertPeptideHit(PeptideHit(27.0,"Mascot",1,"RRYA"));
spec->getIdentification().push_back(dbs);
//ClusterSpectrum cspec(spec);
//const ClusterSpectrum& ccspec = cspec;

vector<String> catalogue = factoryp->catalogue("FilterFunctor");

// go through all registered FilterFunctors and check if they accept a spectrum 
// and return something
for (vector<String>::const_iterator cvit = catalogue.begin(); cvit != catalogue.end(); ++cvit)
{
  CHECK()
    STATUS("ClusterFactory::create("+*cvit+")")
    ffp = dynamic_cast<FilterFunctor*>(factoryp->create(*cvit));
    TEST_NOT_EQUAL(ffp, 0)
  RESULT

  CHECK(FilterFunctor::operator())
    STATUS(*cvit+"::operator()")
    //double result = ffp->apply(*spec);
    //TEST_NOT_EQUAL(result, 0)
  RESULT

  delete ffp;
}

factoryp->destroy();

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
