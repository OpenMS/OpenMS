// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  Copyright The OpenMS team, Eberhard Karls University Tübingen,
//  ETH Zürich and FU Berlin 2001-2012.
//  This software is released under a BSD license. For a full list of
//  authors, refer to the file AUTHORS. For full licensing conditions
//  refer to the file LICENSE.
//
// --------------------------------------------------------------------------
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger, Hannes Roest, Witold Wolski $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/TraMLFile.h>

#include <boost/assign/std/vector.hpp>
#include <boost/assign/list_of.hpp>

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/MRMFeatureAccessOpenMS.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MRMFeatureAccessOpenMS, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

//FeatureOpenMS
{
FeatureOpenMS* ptr = 0;
FeatureOpenMS* nullPointer = 0;

START_SECTION(FeatureOpenMS())
{
  Feature f;
  ptr = new FeatureOpenMS(f);
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~FeatureOpenMS())
{
  delete ptr;
}
END_SECTION
}

//MRMFeatureOpenMS
{
MRMFeatureOpenMS* ptr = 0;
MRMFeatureOpenMS* nullPointer = 0;

START_SECTION(MRMFeatureOpenMS())
{
  MRMFeature f;
  ptr = new MRMFeatureOpenMS(f);
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~MRMFeatureOpenMS())
{
  delete ptr;
}
END_SECTION
}

//TransitionGroupOpenMS
{
TransitionGroupOpenMS<MSSpectrum, Peak1D, ReactionMonitoringTransition>* ptr = 0;
TransitionGroupOpenMS<MSSpectrum, Peak1D, ReactionMonitoringTransition>* nullPointer = 0;

START_SECTION(TransitionGroupOpenMS())
{
  MRMTransitionGroup<MSSpectrum, Peak1D, ReactionMonitoringTransition> trgroup;
  ptr = new TransitionGroupOpenMS<MSSpectrum, Peak1D, ReactionMonitoringTransition>(trgroup);
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~TransitionGroupOpenMS())
{
  delete ptr;
}
END_SECTION
}

//SignalToNoiseOpenMS
{
SignalToNoiseOpenMS<Peak1D>* ptr = 0;
SignalToNoiseOpenMS<Peak1D>* nullPointer = 0;

START_SECTION(SignalToNoiseOpenMS())
{
  OpenMS::MSSpectrum<> chromat;
  ptr = new SignalToNoiseOpenMS<Peak1D>(chromat, 1.0, 3.0);
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~SignalToNoiseOpenMS())
{
  delete ptr;
}
END_SECTION
}
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



