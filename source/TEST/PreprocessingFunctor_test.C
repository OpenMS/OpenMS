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
// $Id: PreprocessingFunctor_test.C,v 1.8 2006/04/05 11:18:25 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/DTAFile.h>
///////////////////////////

#include <OpenMS/COMPARISON/CLUSTERING/ClusterFactory.h>
#include <OpenMS/FILTERING/TRANSFORMERS/PreprocessingFunctor.h>

///////////////////////////

#include <vector>
#include <iostream>

///////////////////////////
START_TEST(PreprocessingFunctor, "$Id: PreprocessingFunctor_test.C,v 1.8 2006/04/05 11:18:25 marc_sturm Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

PreprocessingFunctor* mfp;

ClusterFactory* factoryp = ClusterFactory::instance();

DTAFile dtafile;

MSSpectrum< DPeak<1> > spec;
dtafile.load("data/spectrum.dta",spec);

vector<String> catalogue = factoryp->catalogue("PreprocessingFunctor");

// go through all registered PreprocessingFunctors and check if they accept a spectrum 
for ( vector<String>::const_iterator cvit = catalogue.begin();
    cvit != catalogue.end(); ++cvit )
{
  CHECK()
    STATUS("ClusterFactory::create("+*cvit+")")
    mfp = dynamic_cast<PreprocessingFunctor*>(factoryp->create(*cvit));
    TEST_NOT_EQUAL(mfp, 0)
  RESULT

  CHECK(PreprocessingFunctor::operator())
    STATUS(*cvit+"::operator()")
    (*mfp)(spec);
  RESULT

  delete mfp;
}

factoryp->destroy();

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
