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
// $Maintainer: Stephan Aiche$
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/TraceFitter.h>
#include <OpenMS/KERNEL/Peak1D.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(TraceFitter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TraceFitter<Peak1D>* ptr = 0;
START_SECTION(TraceFitter())
{
	ptr = new TraceFitter<Peak1D>();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~TraceFitter())
{
	delete ptr;
}
END_SECTION

START_SECTION((TraceFitter(const TraceFitter &source)))
{
  NOT_TESTABLE
  // has no public members to check if copy has same proberties
}
END_SECTION

START_SECTION((virtual TraceFitter& operator=(const TraceFitter &source)))
{
  NOT_TESTABLE
  // has no public members to check if copy has same proberties
}
END_SECTION

TraceFitter<Peak1D> trace_fitter;
START_SECTION((virtual void fit(FeatureFinderAlgorithmPickedHelperStructs::MassTraces< PeakType > &)))
{
  FeatureFinderAlgorithmPickedHelperStructs::MassTraces<Peak1D> m;
  TEST_EXCEPTION(Exception::NotImplemented, trace_fitter.fit(m))
}
END_SECTION

START_SECTION((virtual DoubleReal getLowerRTBound() const ))
{
  TEST_EXCEPTION(Exception::NotImplemented, trace_fitter.getLowerRTBound())
}
END_SECTION

START_SECTION((virtual DoubleReal getUpperRTBound() const ))
{
  TEST_EXCEPTION(Exception::NotImplemented, trace_fitter.getUpperRTBound())
}
END_SECTION

START_SECTION((virtual DoubleReal getHeight() const ))
{
  TEST_EXCEPTION(Exception::NotImplemented, trace_fitter.getHeight())
}
END_SECTION

START_SECTION((virtual DoubleReal getCenter() const ))
{
  TEST_EXCEPTION(Exception::NotImplemented, trace_fitter.getCenter())
}
END_SECTION

START_SECTION((DoubleReal computeTheoretical(const FeatureFinderAlgorithmPickedHelperStructs::MassTrace< PeakType > &, Size )))
{
  FeatureFinderAlgorithmPickedHelperStructs::MassTrace<Peak1D> mt;
  Size i = 0;
  TEST_EXCEPTION(Exception::NotImplemented, trace_fitter.computeTheoretical(mt,i))
}
END_SECTION

START_SECTION((bool checkMinimalRTSpan(const std::pair< DoubleReal, DoubleReal > &, const DoubleReal)))
{
  std::pair<DoubleReal, DoubleReal> p;
  DoubleReal x = 0.0;
  TEST_EXCEPTION(Exception::NotImplemented, trace_fitter.checkMinimalRTSpan(p,x))
}
END_SECTION

START_SECTION((virtual bool checkMaximalRTSpan(const DoubleReal)))
{
  DoubleReal x = 0.0;
  TEST_EXCEPTION(Exception::NotImplemented, trace_fitter.checkMaximalRTSpan(x))
}
END_SECTION

START_SECTION((virtual DoubleReal getFeatureIntensityContribution()))
{
  TEST_EXCEPTION(Exception::NotImplemented, trace_fitter.getFeatureIntensityContribution())
}
END_SECTION

START_SECTION((virtual String getGnuplotFormula(FeatureFinderAlgorithmPickedHelperStructs::MassTrace<PeakType> const & , const char , const DoubleReal, const DoubleReal)))
{
  FeatureFinderAlgorithmPickedHelperStructs::MassTrace<Peak1D> mt;
  DoubleReal shift = 0.0;
  DoubleReal baseline = 0.0;
  char f = 'f';
  TEST_EXCEPTION(Exception::NotImplemented, trace_fitter.getGnuplotFormula(mt, f, baseline, shift))
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



