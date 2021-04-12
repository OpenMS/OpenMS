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
// $Maintainer: Mathias Walzer$
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/COMPARISON/SPECTRA/BinnedSumAgreeingIntensities.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrum.h>
#include <OpenMS/FORMAT/DTAFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(BinnedSumAgreeingIntensities, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

BinnedSumAgreeingIntensities* ptr = nullptr;
BinnedSumAgreeingIntensities* nullPointer = nullptr;
START_SECTION(BinnedSumAgreeingIntensities())
{
  ptr = new BinnedSumAgreeingIntensities();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~BinnedSumAgreeingIntensities())
{
  delete ptr;
}
END_SECTION

ptr = new BinnedSumAgreeingIntensities();

START_SECTION((BinnedSumAgreeingIntensities(const BinnedSumAgreeingIntensities &source)))
{
  BinnedSumAgreeingIntensities copy(*ptr);
  TEST_EQUAL(copy.getName(), ptr->getName());
  TEST_EQUAL(copy.getParameters(), ptr->getParameters());
}
END_SECTION

START_SECTION((BinnedSumAgreeingIntensities& operator=(const BinnedSumAgreeingIntensities &source)))
{
  BinnedSumAgreeingIntensities copy;
  copy = *ptr;
  TEST_EQUAL(copy.getName(), ptr->getName());
  TEST_EQUAL(copy.getParameters(), ptr->getParameters());
}
END_SECTION

START_SECTION((double operator()(const BinnedSpectrum &spec1, const BinnedSpectrum &spec2) const))
{
  PeakSpectrum s1, s2;
  DTAFile().load(OPENMS_GET_TEST_DATA_PATH("PILISSequenceDB_DFPIANGER_1.dta"), s1);
  DTAFile().load(OPENMS_GET_TEST_DATA_PATH("PILISSequenceDB_DFPIANGER_1.dta"), s2);
  s2.pop_back();
  BinnedSpectrum bs1(s1, 1.5, false, 2, 0.0 );
  BinnedSpectrum bs2(s2, 1.5, false, 2, 0.0);

  double score = (*ptr)(bs1, bs2);
  TEST_REAL_SIMILAR(score, 0.99707)

  score = (*ptr)(bs1, bs1);
  TEST_REAL_SIMILAR(score, 1.0)
}
END_SECTION

START_SECTION((double operator()(const BinnedSpectrum &spec) const ))
{
  PeakSpectrum s1;
  DTAFile().load(OPENMS_GET_TEST_DATA_PATH("PILISSequenceDB_DFPIANGER_1.dta"), s1);
  BinnedSpectrum bs1(s1, 1.5, false, 2, BinnedSpectrum::DEFAULT_BIN_OFFSET_LOWRES);
  double score = (*ptr)(bs1);
  TEST_REAL_SIMILAR(score, 1);
}
END_SECTION

START_SECTION((static BinnedSpectrumCompareFunctor* create()))
{
  BinnedSpectrumCompareFunctor* bsf = BinnedSumAgreeingIntensities::create();
  BinnedSumAgreeingIntensities bsp;
  TEST_EQUAL(bsf->getParameters(), bsp.getParameters())
  TEST_EQUAL(bsf->getName(), bsp.getName())
}
END_SECTION

START_SECTION((static const String getProductName()))
{
  TEST_EQUAL(ptr->getProductName(), "BinnedSumAgreeingIntensities")
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



