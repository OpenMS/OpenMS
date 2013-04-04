// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Mathias Walzer$
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrumCompareFunctor.h>
#include <OpenMS/CONCEPT/Factory.h>

using namespace OpenMS;
using namespace std;

START_TEST(BinnedSpectrumCompareFunctor, "$Id$")

/////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

START_SECTION(BinnedSpectrumCompareFunctor())
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION(~BinnedSpectrumCompareFunctor())
{
  NOT_TESTABLE
}
END_SECTION

//interface class is not testable

START_SECTION((BinnedSpectrumCompareFunctor(const BinnedSpectrumCompareFunctor &source)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((BinnedSpectrumCompareFunctor& operator=(const BinnedSpectrumCompareFunctor &source)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual double operator()(const BinnedSpectrum &spec1, const BinnedSpectrum &spec2) const =0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual double operator()(const BinnedSpectrum &spec) const =0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((static void registerChildren()))
{
  BinnedSpectrumCompareFunctor* c1 = Factory<BinnedSpectrumCompareFunctor>::create("BinnedSharedPeakCount");
  TEST_EQUAL(c1->getName(), "BinnedSharedPeakCount")
  c1 = Factory<BinnedSpectrumCompareFunctor>::create("BinnedSpectralContrastAngle");
  TEST_EQUAL(c1->getName(), "BinnedSpectralContrastAngle")
  c1 = Factory<BinnedSpectrumCompareFunctor>::create("BinnedSumAgreeingIntensities");
  TEST_EQUAL(c1->getName(), "BinnedSumAgreeingIntensities")
}
END_SECTION

START_SECTION((static const String getProductName()))
{
	TEST_EQUAL(BinnedSpectrumCompareFunctor::getProductName(), "BinnedSpectrumCompareFunctor")
}
END_SECTION

START_SECTION(([BinnedSpectrumCompareFunctor::IncompatibleBinning] IncompatibleBinning(const char *file, int line, const char *function, const char *message="compared spectra have different settings in binsize and/or binspread")))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION(([BinnedSpectrumCompareFunctor::IncompatibleBinning] virtual ~IncompatibleBinning()))
{
  NOT_TESTABLE
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



