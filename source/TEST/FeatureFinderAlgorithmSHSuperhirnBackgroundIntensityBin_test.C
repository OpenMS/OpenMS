// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Florian Zeller$
// $Authors: Florian Zeller$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ms_peak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/BackgroundIntensityBin.h>

///////////////////////////

START_TEST(BackgroundIntensityBin, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

BackgroundIntensityBin* ptr;
START_SECTION((BackgroundIntensityBin(double, double)))
	ptr = new BackgroundIntensityBin(300, 12);
	TEST_NOT_EQUAL(ptr,0)
END_SECTION

START_SECTION((~BackgroundIntensityBin()))
	delete ptr;
END_SECTION

START_SECTION((checkBelonging(ms_peak*)))
  ptr = new BackgroundIntensityBin(300, 12);
  ms_peak *p = new ms_peak();
  TEST_EQUAL(ptr->checkBelonging(p), false);
  delete p;

  ms_peak* p2 = new ms_peak(1, 300, 100);  // (int IN_scan, double IN_mass, float IN_intens)
  p2->set_retention_time(12);
  TEST_EQUAL(ptr->checkBelonging(p2), true);
  delete p2;
END_SECTION

START_SECTION((addIntensity( double )))
  ptr = new BackgroundIntensityBin(300, 12);
  TEST_EQUAL(ptr->getIntensityMap()->size(), 0);
  ptr->addIntensity(100);
  TEST_EQUAL(ptr->getIntensityMap()->size(), 1);
END_SECTION

START_SECTION((addMSPeak( ms_peak* )))
  ptr = new BackgroundIntensityBin(300, 12);
  ms_peak* p = new ms_peak(1, 300, 100);  // (int IN_scan, double IN_mass, float IN_intens)
  TEST_EQUAL(ptr->getIntensityMap()->size(), 0);
  ptr->addMSPeak(p);
  TEST_EQUAL(ptr->getIntensityMap()->size(), 1);
  delete p;
END_SECTION

START_SECTION((processIntensities()))
  ptr = new BackgroundIntensityBin(300, 12);
  ptr->processIntensities();
  TEST_REAL_SIMILAR(ptr->getMean(), 0);
END_SECTION

START_SECTION((getIntensityHist()))
  ptr = new BackgroundIntensityBin(300, 12);
  TEST_NOT_EQUAL(ptr->getIntensityHist(), 0);
END_SECTION

START_SECTION((getMean()))
  ptr = new BackgroundIntensityBin(300, 12);
  ptr->processIntensities();
  //TEST_EQUAL(ptr->getMean(), 0)
  TEST_REAL_SIMILAR(ptr->getMean(), 0);
END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
