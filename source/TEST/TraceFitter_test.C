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

// dummy implementation for the test
template <class PeakType>
class DerivedTraceFitter
    : public TraceFitter<PeakType>
{

public:

    void fit(FeatureFinderAlgorithmPickedHelperStructs::MassTraces<PeakType>&)
    {
        throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
    }

    DoubleReal getLowerRTBound() const
    {
        throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
    }

    DoubleReal getUpperRTBound() const
    {
        throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
    }

    DoubleReal getHeight() const
    {
        throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
    }

    DoubleReal getCenter() const
    {
        throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
    }

    DoubleReal getFWHM() const
    {
        throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
    }

    DoubleReal computeTheoretical(const FeatureFinderAlgorithmPickedHelperStructs::MassTrace<PeakType>&, Size)
    {
        throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
    }

    bool checkMinimalRTSpan(const std::pair<double, double>&, const DoubleReal)
    {
        throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
    }

    bool checkMaximalRTSpan(const DoubleReal)
    {
        throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
    }

    DoubleReal getFeatureIntensityContribution()
    {
        throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
    }

    String getGnuplotFormula(FeatureFinderAlgorithmPickedHelperStructs::MassTrace<PeakType> const &, const char, const DoubleReal, const DoubleReal)
    {
        throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
    }

    void printState_(SignedSize, gsl_multifit_fdfsolver*)
    {
        throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
    }

    void getOptimizedParameters_(gsl_multifit_fdfsolver*)
    {
        throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
    }

};

START_TEST(TraceFitter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TraceFitter<Peak1D>* ptr = 0;
TraceFitter<Peak1D>* nullPointer = 0;
START_SECTION(TraceFitter())
{
    ptr = new DerivedTraceFitter<Peak1D>();
	TEST_NOT_EQUAL(ptr, nullPointer)
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

DerivedTraceFitter<Peak1D> trace_fitter;
START_SECTION((virtual void fit(FeatureFinderAlgorithmPickedHelperStructs::MassTraces< PeakType > &traces)=0))
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

START_SECTION((virtual DoubleReal computeTheoretical(const FeatureFinderAlgorithmPickedHelperStructs::MassTrace< PeakType > &trace, Size k)=0))
{
  FeatureFinderAlgorithmPickedHelperStructs::MassTrace<Peak1D> mt;
  Size i = 0;
  TEST_EXCEPTION(Exception::NotImplemented, trace_fitter.computeTheoretical(mt,i))
}
END_SECTION

START_SECTION((virtual bool checkMinimalRTSpan(const std::pair< DoubleReal, DoubleReal > &rt_bounds, const DoubleReal min_rt_span)=0))
{
  std::pair<DoubleReal, DoubleReal> p;
  DoubleReal x = 0.0;
  TEST_EXCEPTION(Exception::NotImplemented, trace_fitter.checkMinimalRTSpan(p,x))
}
END_SECTION

START_SECTION((virtual bool checkMaximalRTSpan(const DoubleReal max_rt_span)=0))
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

START_SECTION((virtual String getGnuplotFormula(FeatureFinderAlgorithmPickedHelperStructs::MassTrace< PeakType > const &trace, const char function_name, const DoubleReal baseline, const DoubleReal rt_shift)=0))
{
  FeatureFinderAlgorithmPickedHelperStructs::MassTrace<Peak1D> mt;
  DoubleReal shift = 0.0;
  DoubleReal baseline = 0.0;
  char f = 'f';
  TEST_EXCEPTION(Exception::NotImplemented, trace_fitter.getGnuplotFormula(mt, f, baseline, shift))
}
END_SECTION

START_SECTION((virtual DoubleReal getFWHM() const))
{
  TEST_EXCEPTION(Exception::NotImplemented, trace_fitter.getFWHM())
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



