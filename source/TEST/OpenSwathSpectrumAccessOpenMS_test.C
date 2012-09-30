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
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SpectrumAccessOpenMS, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SpectrumAccessOpenMS* ptr = 0;
SpectrumAccessOpenMS* nullPointer = 0;

START_SECTION(SpectrumAccessOpenMS())
{
  MSExperiment<> exp;
  ptr = new SpectrumAccessOpenMS(exp);
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~SpectrumAccessOpenMS())
{
  delete ptr;
}
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



