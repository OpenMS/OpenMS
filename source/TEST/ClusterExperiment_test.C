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
// $Id: ClusterExperiment_test.C,v 1.9 2006/04/05 11:18:25 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:  $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/COMPARISON/CLUSTERING/ClusterExperiment.h>

///////////////////////////

#include <OpenMS/COMPARISON/CLUSTERING/ClusterFactory.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/KERNEL/DSpectrum.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterSpectrum.h>
#include <OpenMS/COMPARISON/SPECTRA/CompareFunctor.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterFunctor.h>
#include <OpenMS/FILTERING/TRANSFORMERS/MowerFunctor.h>
#include <OpenMS/COMPARISON/CLUSTERING/AnalysisFunctor.h>

///////////////////////////
START_TEST(ClusterExperiment, "$Id: ClusterExperiment_test.C,v 1.9 2006/04/05 11:18:25 marc_sturm Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

ClusterFactory* factoryp = ClusterFactory::instance();

DTAFile dtafile;

MSSpectrum< DPeak<1> >* spec = new MSSpectrum< DPeak<1> >();
MSSpectrum< DPeak<1> >* spec2 = new MSSpectrum< DPeak<1> >();

dtafile.load("data/spectrum.dta",*spec);
dtafile.load("data/spectrum2.dta",*spec2);

ClusterExperiment* cexp = 0;
ClusterExperiment::ClusterRun* crp = 0;

CHECK(ClusteExperiment::ClusterExperiment())
  cexp = new ClusterExperiment();
  TEST_NOT_EQUAL(cexp,0)
RESULT

CHECK(ClusterExperiment::createrun())
  cexp->createrun();
  TEST_EQUAL(cexp->size(),1)
RESULT

CHECK(ClusterExperiment::setBinSize())
  cexp->setBinSize(1.234);
  TEST_REAL_EQUAL((*cexp)[0].getBinSize(),1.234);
RESULT

CHECK(ClusterExperiment::setBinSpread())
  cexp->setBinSpread(8);
  TEST_EQUAL((*cexp)[0].getBinSpread(),8);
RESULT

CHECK(ClusterExperiment::setNorm())
  cexp->setNorm(geometric);
  TEST_EQUAL((*cexp)[0].getNorm(),geometric)
RESULT

CHECK(ClusterExperiment::setSimFunc())
  cexp->setSimFunc(dynamic_cast<CompareFunctor*>(factoryp->create("BinnedRepSpectrumContrastAngle")));
  TEST_EQUAL((*cexp)[0].getSimFunc()->getName(), "BinnedRepSpectrumContrastAngle")
RESULT

CHECK(ClusterExperiment::setClusterFunc())
  cexp->setClusterFunc(dynamic_cast<ClusterFunctor*>(factoryp->create("LinkageCluster")));
  TEST_EQUAL((*cexp)[0].getClusterFunc()->getName(),"LinkageCluster")
RESULT

CHECK(ClusterExperiment::addMower())
  cexp->addMower(dynamic_cast<MowerFunctor*>(factoryp->create("ParentPeakMower")));
  cexp->addMower(dynamic_cast<MowerFunctor*>(factoryp->create("Normalizer")));
  TEST_EQUAL((*cexp)[0].getPreprocessqueue()[0]->getName(),"ParentPeakMower")
  TEST_EQUAL((*cexp)[0].getPreprocessqueue()[1]->getName(),"Normalizer")
RESULT

CHECK(ClusterExperiment::addAnalysisFunctor())
  cexp->addAnalysisFunctor(dynamic_cast<AnalysisFunctor*>(factoryp->create("ClusterCompareFunctor")));
  cexp->addAnalysisFunctor(dynamic_cast<AnalysisFunctor*>(factoryp->create("DistanceAnalyzer")));
  TEST_EQUAL((*cexp)[0][0].name(),"ClusterCompareFunctor")
  TEST_EQUAL((*cexp)[0][1].name(),"DistanceAnalyzer")
RESULT

CHECK(ClusterExperiment::ClusterRun::ClusterRun())
  crp = new ClusterExperiment::ClusterRun();
  TEST_NOT_EQUAL(crp,0)
RESULT

CHECK(ClusterExperiment::ClusterRun::ClusterRun(const ClusterExperiment::ClusterRun&))
  ClusterSpectrum cspec(*spec,0,(*cexp)[0].getBinSize(),(*cexp)[0].getBinSpread());
  ClusterSpectrum cspec2(*spec2,0,(*cexp)[0].getBinSize(),(*cexp)[0].getBinSpread());
  double sim1 = (*cexp)[0].similarity(cspec,cspec2);
  ClusterExperiment::ClusterRun cr((*cexp)[0]);
  ClusterSpectrum cspec3(*spec,0,cr.getBinSize(),cr.getBinSpread());
  ClusterSpectrum cspec4(*spec2,0,cr.getBinSize(),cr.getBinSpread());
  double sim2 = cr.similarity(cspec3,cspec4);
  TEST_REAL_EQUAL(sim1,sim2)
RESULT

CHECK(ClusterExperiment::ClusterRun::operator=(const ClusterExperiment::ClusterRun&))
  ClusterSpectrum cspec(*spec,0,(*cexp)[0].getBinSize(),(*cexp)[0].getBinSpread());
  ClusterSpectrum cspec2(*spec2,0,(*cexp)[0].getBinSize(),(*cexp)[0].getBinSpread());
  double sim1 = (*cexp)[0].similarity(cspec,cspec2);
  ClusterExperiment::ClusterRun cr;
  cr = (*cexp)[0];
  ClusterSpectrum cspec3(*spec,0,cr.getBinSize(),cr.getBinSpread());
  ClusterSpectrum cspec4(*spec2,0,cr.getBinSize(),cr.getBinSpread());
  double sim2 = cr.similarity(cspec3,cspec4);
  TEST_REAL_EQUAL(sim1,sim2)
RESULT

CHECK(ClusterExperiment::ClusterRun::~ClusterRun())
  delete crp;
RESULT

CHECK(ClusterExperiment::ClusterRun::ClusterRun(const ClusterRun&))
RESULT

CHECK(ClusteExperiment::~ClusterExperiment())
  delete cexp;
RESULT

factoryp->destroy();

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
