// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse, Holger Plattfaut $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/FILTERING/DATAREDUCTION/SILACFilter.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SILACFiltering.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/PeakWidthEstimator.h>

using namespace OpenMS;
using namespace std;

START_TEST(SILACFiltering, "$Id$")

MSExperiment<> input;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("SILACFiltering_test.mzML"), input);
const PeakWidthEstimator::Result peak_width(PeakWidthEstimator::estimateFWHM(input));

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

std::vector<DoubleReal> mass_separations;
mass_separations.push_back(8.0142);
SILACFiltering filtering(input, peak_width, 0);
SILACFilter filter(mass_separations, 2, 2, 3, 0, .9, false);

START_SECTION((SILACFiltering(MSExperiment< Peak1D > &exp, const PeakWidthEstimator::Result &, const DoubleReal intensity_cutoff, const String debug_filebase_="")))
{
  TEST_EQUAL(filtering.filters_.size(), 0);
  TEST_EQUAL(filtering.blacklist.size(), 0);
}
END_SECTION

START_SECTION((void addFilter(SILACFilter &filter)))
{
  filtering.addFilter(filter);
  TEST_EQUAL(filtering.filters_.size(), 1);
}
END_SECTION

START_SECTION((void filterDataPoints()))
{
  filtering.filterDataPoints();
  SILACFiltering::Filters::iterator filter_it = filtering.filters_.begin();

  std::vector<SILACPattern> &p = filter_it->getElements();
  TEST_EQUAL(p.size(), 3);
  TEST_REAL_SIMILAR(p[0].rt, 830);
  TEST_REAL_SIMILAR(p[0].mz, 670.84);
  TEST_REAL_SIMILAR(p[1].rt, 830);
  TEST_REAL_SIMILAR(p[1].mz, 670.84);
  TEST_REAL_SIMILAR(p[2].rt, 833);
  TEST_REAL_SIMILAR(p[2].mz, 670.84);
}
END_SECTION

START_SECTION([SILACFiltering::SpectrumInterpolation] SpectrumInterpolation(const MSSpectrum<> &, const SILACFiltering &))
  SILACFiltering::SpectrumInterpolation si(input[0], filtering);
END_SECTION

START_SECTION([SILACFiltering::SpectrumInterpolation] ~SpectrumInterpolation())
  SILACFiltering::SpectrumInterpolation si(input[0], filtering);
END_SECTION

START_SECTION([SILACFiltering::SpectrumInterpolation] DoubleReal operator()(DoubleReal mz) const)
  SILACFiltering::SpectrumInterpolation si(input[0], filtering);
  TEST_REAL_SIMILAR(si(670.5), 0)
  TEST_REAL_SIMILAR(si(671.1), 0)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

