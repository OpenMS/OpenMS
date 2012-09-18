// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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



