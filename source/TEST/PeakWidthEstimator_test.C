// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:expandtab
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2011 -- Bastian Blank
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Bastian Blank $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/PeakWidthEstimator.h>

#include <OpenMS/FORMAT/MzDataFile.h>

using namespace OpenMS;

START_TEST(PeakWidthEstimator, "$Id$")

MSExperiment<> input;
MzDataFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_orbitrap.mzData"), input);

START_SECTION(static void estimateSpectrumFWHM(const MSSpectrum<> &, std::set<boost::tuple<DoubleReal, DoubleReal, DoubleReal> > &))
{
  typedef std::set<boost::tuple<DoubleReal, DoubleReal, DoubleReal> > Fwhm;
  Fwhm fwhm;
  PeakWidthEstimator::estimateSpectrumFWHM(input[0], fwhm);
  TEST_EQUAL(fwhm.size(), 151);
  Fwhm::const_reverse_iterator it = fwhm.rbegin();
  TEST_REAL_SIMILAR(it->get<0>(), 202394.);
  TEST_REAL_SIMILAR(it->get<1>(), 591.358);
  TEST_REAL_SIMILAR(it->get<2>(), .010647);
}
END_SECTION

START_SECTION(static Result estimateFWHM(const MSExperiment<> &))
{
  PeakWidthEstimator::Result r(PeakWidthEstimator::estimateFWHM(input));
  TEST_REAL_SIMILAR(r.c0, -14.15849);
  TEST_REAL_SIMILAR(r.c1, 1.50632);
}
END_SECTION

END_TEST

