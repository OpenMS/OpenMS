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

#include <OpenMS/COMPARISON/CLUSTERING/ClusterFactory.h>
#include <OpenMS/COMPARISON/CLUSTERING/AnalysisFunctor.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterExperiment.h>

///////////////////////////


///////////////////////////
START_TEST(AnalysisFunctor, "$Id: AnalysisFunctor_test.C,v 1.4 2006/03/28 12:53:13 marc_sturm Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

AnalysisFunctor* afp;

ClusterExperiment cexp;
cexp.load("data/clusterexperiment.xml");

ClusterFactory* factoryp = ClusterFactory::instance();

vector<String> catalogue = factoryp->catalogue("AnalysisFunctor");

// go through all registered FilterFunctors and check if they accept a spectrum 
// and return something
// todo
for ( vector<String>::const_iterator cvit = catalogue.begin();
    cvit != catalogue.end(); ++cvit )
{
  CHECK()
    STATUS( *cvit + "::create() + " + *cvit + "::" + *cvit + "()" )
    afp = dynamic_cast<AnalysisFunctor*>(factoryp->create(*cvit));
    TEST_NOT_EQUAL(afp, 0)
  RESULT

  CHECK()
    STATUS(*cvit + "::operator()")
    if ( afp->needsDBAdapter() )
    {
      STATUS(*cvit + "cannot be tested without DataBase")
    }
    else
    {
      if ( afp->needsClusterRun() )
      {
        afp->setClusterRun(&cexp[0]);
      }
      (*afp)(cexp[1].getClustering());
    }
  RESULT
  
  CHECK()
    STATUS(*cvit + "::~" + *cvit + "()")
    delete afp;
  RESULT
}

factoryp->destroy();

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
