// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/TraceFitter.h>
///////////////////////////

#include <OpenMS/KERNEL/Peak1D.h>

using namespace OpenMS;
using namespace std;

// dummy implementation for the test
class DerivedTraceFitter
    : public TraceFitter
{

public:

    void fit(FeatureFinderAlgorithmPickedHelperStructs::MassTraces&) override
    {
        throw Exception::NotImplemented(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION);
    }

    double getLowerRTBound() const override
    {
        throw Exception::NotImplemented(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION);
    }

    double getUpperRTBound() const override
    {
        throw Exception::NotImplemented(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION);
    }

    double getHeight() const override
    {
        throw Exception::NotImplemented(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION);
    }

    double getCenter() const override
    {
        throw Exception::NotImplemented(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION);
    }

    double getFWHM() const override
    {
        throw Exception::NotImplemented(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION);
    }

    double getValue(double /* rt */) const override
    {
        throw Exception::NotImplemented(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION);
    }

    bool checkMinimalRTSpan(const std::pair<double, double>&, const double) override
    {
        throw Exception::NotImplemented(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION);
    }

    bool checkMaximalRTSpan(const double) override
    {
        throw Exception::NotImplemented(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION);
    }

    double getArea() override
    {
        throw Exception::NotImplemented(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION);
    }

    String getGnuplotFormula(const FeatureFinderAlgorithmPickedHelperStructs::MassTrace&, const char, const double, const double) override
    {
        throw Exception::NotImplemented(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION);
    }

    void getOptimizedParameters_(const Eigen::VectorXd&) override
    {
        throw Exception::NotImplemented(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION);
    }

};

START_TEST(TraceFitter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TraceFitter* ptr = nullptr;
TraceFitter* nullPointer = nullptr;
START_SECTION(TraceFitter())
{
  ptr = new DerivedTraceFitter();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~TraceFitter())
{
	delete ptr;
}
END_SECTION

START_SECTION((TraceFitter(const TraceFitter& source)))
{
  NOT_TESTABLE
  // has no public members to check if copy has same proberties
}
END_SECTION

START_SECTION((virtual TraceFitter& operator=(const TraceFitter& source)))
{
  NOT_TESTABLE
  // has no public members to check if copy has same proberties
}
END_SECTION

DerivedTraceFitter trace_fitter;
START_SECTION((virtual void fit(FeatureFinderAlgorithmPickedHelperStructs::MassTraces& traces)=0))
{
  FeatureFinderAlgorithmPickedHelperStructs::MassTraces m;
  TEST_EXCEPTION(Exception::NotImplemented, trace_fitter.fit(m))
}
END_SECTION

START_SECTION((virtual double getLowerRTBound() const ))
{
  TEST_EXCEPTION(Exception::NotImplemented, trace_fitter.getLowerRTBound())
}
END_SECTION

START_SECTION((virtual double getUpperRTBound() const ))
{
  TEST_EXCEPTION(Exception::NotImplemented, trace_fitter.getUpperRTBound())
}
END_SECTION

START_SECTION((virtual double getHeight() const ))
{
  TEST_EXCEPTION(Exception::NotImplemented, trace_fitter.getHeight())
}
END_SECTION

START_SECTION((virtual double getCenter() const ))
{
  TEST_EXCEPTION(Exception::NotImplemented, trace_fitter.getCenter())
}
END_SECTION

START_SECTION((virtual double getValue(double rt) const ))
{
  TEST_EXCEPTION(Exception::NotImplemented, trace_fitter.getValue(0.0))
}
END_SECTION

START_SECTION((double computeTheoretical(const FeatureFinderAlgorithmPickedHelperStructs::MassTrace& trace, Size k)))
{
  FeatureFinderAlgorithmPickedHelperStructs::MassTrace mt;
  Peak1D p;
  mt.peaks.push_back(make_pair(1.0, &p));
  TEST_EXCEPTION(Exception::NotImplemented, trace_fitter.computeTheoretical(mt, 0))
}
END_SECTION

START_SECTION((virtual bool checkMinimalRTSpan(const std::pair<double, double>& rt_bounds, const double min_rt_span)=0))
{
  std::pair<double, double> p;
  double x = 0.0;
  TEST_EXCEPTION(Exception::NotImplemented, trace_fitter.checkMinimalRTSpan(p,x))
}
END_SECTION

START_SECTION((virtual bool checkMaximalRTSpan(const double max_rt_span)=0))
{
  double x = 0.0;
  TEST_EXCEPTION(Exception::NotImplemented, trace_fitter.checkMaximalRTSpan(x))
}
END_SECTION

START_SECTION((virtual double getArea()))
{
  TEST_EXCEPTION(Exception::NotImplemented, trace_fitter.getArea())
}
END_SECTION

START_SECTION((virtual String getGnuplotFormula(const FeatureFinderAlgorithmPickedHelperStructs::MassTrace& trace, const char function_name, const double baseline, const double rt_shift)=0))
{
  FeatureFinderAlgorithmPickedHelperStructs::MassTrace mt;
  double shift = 0.0;
  double baseline = 0.0;
  char f = 'f';
  TEST_EXCEPTION(Exception::NotImplemented, trace_fitter.getGnuplotFormula(mt, f, baseline, shift))
}
END_SECTION

START_SECTION((virtual double getFWHM() const))
{
  TEST_EXCEPTION(Exception::NotImplemented, trace_fitter.getFWHM())
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



