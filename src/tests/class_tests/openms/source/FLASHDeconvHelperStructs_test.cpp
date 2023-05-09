// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Maintainer: Jihyung Kim$
// $Authors: Jihyung Kim$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

typedef FLASHDeconvHelperStructs::LogMzPeak LogMzPeak;
typedef FLASHDeconvHelperStructs::PrecalculatedAveragine PrecalculatedAveragine;

START_TEST(FLASHDeconvHelperStructs, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FLASHDeconvHelperStructs* ptr = 0;
FLASHDeconvHelperStructs* null_ptr = 0;
START_SECTION(FLASHDeconvHelperStructs())
{
  ptr = new FLASHDeconvHelperStructs();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~FLASHDeconvHelperStructs())
{
  delete ptr;
}
END_SECTION


START_SECTION((static double getLogMz(const double mz, const bool positive)))
{
  double mz = 1300;
  double tmp_lmz1 = OpenMS::FLASHDeconvHelperStructs::getLogMz(mz, true);
  double tmp_lmz2 = OpenMS::FLASHDeconvHelperStructs::getLogMz(mz, false);
  TOLERANCE_ABSOLUTE(0.1);
  TEST_REAL_SIMILAR(tmp_lmz1, 7.169344415063863);
  TEST_REAL_SIMILAR(tmp_lmz2, 7.170119121465);
}
END_SECTION

START_SECTION((static double getChargeMass(const bool positive)))
{
  double temp_pos = OpenMS::FLASHDeconvHelperStructs::getChargeMass(true);
  double temp_neg = OpenMS::FLASHDeconvHelperStructs::getChargeMass(false);
  TEST_REAL_SIMILAR(temp_pos, Constants::PROTON_MASS_U);
  TEST_REAL_SIMILAR(temp_neg, -Constants::PROTON_MASS_U);
}
END_SECTION


/// testing LogMzPeak
START_SECTION(([FLASHDeconvHelperStructs::LogMzPeak] LogMzPeak()=default))
{
  LogMzPeak* lmp_ptr = new LogMzPeak();
  LogMzPeak* lmp_null_ptr = 0;

  TEST_NOT_EQUAL(lmp_ptr, lmp_null_ptr);
  delete lmp_ptr;
}
END_SECTION

// test data
Peak1D tmp_p1;
tmp_p1.setIntensity(443505.625);
tmp_p1.setMZ(1125.5118055019082);
START_SECTION(([FLASHDeconvHelperStructs::LogMzPeak] LogMzPeak(const Peak1D &peak, const bool positive)))
{
  LogMzPeak tmp_peak(tmp_p1, true);
  TEST_REAL_SIMILAR(tmp_peak.mz, 1125.5118055019082);
  TEST_REAL_SIMILAR(tmp_peak.intensity, 443505.625);
  TEST_REAL_SIMILAR(tmp_peak.logMz, 7.0250977989903145);
}
END_SECTION

LogMzPeak test_peak(tmp_p1, true);
START_SECTION(([FLASHDeconvHelperStructs::LogMzPeak] LogMzPeak(const LogMzPeak &)))
{
  LogMzPeak tmp_p(test_peak);
  TEST_REAL_SIMILAR(test_peak.mz, tmp_p.mz);
  TEST_REAL_SIMILAR(test_peak.intensity, tmp_p.intensity);
  TEST_REAL_SIMILAR(test_peak.logMz, tmp_p.logMz);
}
END_SECTION

START_SECTION(([FLASHDeconvHelperStructs::LogMzPeak] double getUnchargedMass()))
{
  test_peak.abs_charge = 2;
  TEST_REAL_SIMILAR(test_peak.getUnchargedMass(), 2249.0090580702745);
}
END_SECTION

LogMzPeak test_peak2(test_peak);
test_peak2.logMz = 8.;
START_SECTION(([FLASHDeconvHelperStructs::LogMzPeak] bool operator<(const LogMzPeak &a) const))
{
  bool is_p2_larger = test_peak < test_peak2;
  TEST_EQUAL(is_p2_larger, true);
}
END_SECTION

START_SECTION(([FLASHDeconvHelperStructs::LogMzPeak] bool operator>(const LogMzPeak &a) const))
{
  bool is_p2_larger = test_peak2 > test_peak;
  TEST_EQUAL(is_p2_larger, true);
}
END_SECTION

START_SECTION(([FLASHDeconvHelperStructs::LogMzPeak] bool operator==(const LogMzPeak &other) const))
{
  test_peak2 = test_peak;
  bool are_two_ps_same = test_peak2 == test_peak;
  TEST_EQUAL(are_two_ps_same, true);
}
END_SECTION
///


/// testing PrecalculatedAveragine
START_SECTION(([FLASHDeconvHelperStructs::PrecalculatedAveragine] PrecalculatedAveragine()))
{
  PrecalculatedAveragine* p_avg_ptr = new PrecalculatedAveragine();
  PrecalculatedAveragine* p_avg_null_ptr = 0;

  TEST_NOT_EQUAL(p_avg_ptr, p_avg_null_ptr)
  delete p_avg_ptr;
}
END_SECTION

// test data
CoarseIsotopePatternGenerator generator = CoarseIsotopePatternGenerator();
PrecalculatedAveragine p_avg_test;
START_SECTION(([FLASHDeconvHelperStructs::PrecalculatedAveragine] PrecalculatedAveragine(const double min_mass, const double max_mass, const double delta, CoarseIsotopePatternGenerator& generator, const bool use_RNA_averagine)))
{
  p_avg_test = PrecalculatedAveragine(50, 100, 25, generator, false);
  Size temp_a_idx = p_avg_test.getApexIndex(75);
  double temp_m_diff = p_avg_test.getAverageMassDelta(75);
  TEST_EQUAL(temp_a_idx, 0);
  TOLERANCE_ABSOLUTE(0.3);
  TEST_REAL_SIMILAR(temp_m_diff, 0.04);
}
END_SECTION

START_SECTION(([FLASHDeconvHelperStructs::PrecalculatedAveragine] IsotopeDistribution get(const double mass) const))
{
  IsotopeDistribution tmp_iso = p_avg_test.get(60);
  TOLERANCE_ABSOLUTE(2);
  TEST_REAL_SIMILAR(tmp_iso.getMin(), 53.);
  TEST_REAL_SIMILAR(tmp_iso.getMax(), 55.);
}
END_SECTION

START_SECTION(([FLASHDeconvHelperStructs::PrecalculatedAveragine] void setMaxIsotopeIndex(const int index)))
{
  p_avg_test.setMaxIsotopeIndex(4);
  TEST_EQUAL(p_avg_test.getMaxIsotopeIndex(), 4);
}
END_SECTION

START_SECTION(([FLASHDeconvHelperStructs::PrecalculatedAveragine] int getMaxIsotopeIndex() const))
{
  int tmp_max_idx = p_avg_test.getMaxIsotopeIndex();
  TEST_EQUAL(tmp_max_idx, 4);
}
END_SECTION

START_SECTION(([FLASHDeconvHelperStructs::PrecalculatedAveragine] Size getLeftCountFromApex(const double mass) const))
{
  Size tmp_left = p_avg_test.getLeftCountFromApex(75);
  TEST_EQUAL(tmp_left, 2);
}
END_SECTION

START_SECTION(([FLASHDeconvHelperStructs::PrecalculatedAveragine] Size getRightCountFromApex(const double mass) const))
{
  Size temp_right = p_avg_test.getRightCountFromApex(75);
  TEST_EQUAL(temp_right, 2);
}
END_SECTION

START_SECTION(([FLASHDeconvHelperStructs::PrecalculatedAveragine] Size getApexIndex(const double mass) const))
{
  Size tmp_apex = p_avg_test.getApexIndex(75);
  TEST_EQUAL(tmp_apex, 0);
}
END_SECTION

START_SECTION(([FLASHDeconvHelperStructs::PrecalculatedAveragine] double getAverageMassDelta(const double mass) const))
{
  double tmp_m_delta = p_avg_test.getAverageMassDelta(50);
  TOLERANCE_ABSOLUTE(0.1);
  TEST_REAL_SIMILAR(tmp_m_delta, 0.025);
}
END_SECTION

START_SECTION(([FLASHDeconvHelperStructs::PrecalculatedAveragine] double getMostAbundantMassDelta(const double mass) const))
{
  double tmp_m_delta = p_avg_test.getMostAbundantMassDelta(1000);
  TOLERANCE_ABSOLUTE(0.1);
  TEST_REAL_SIMILAR(tmp_m_delta, 0);
}
END_SECTION

START_SECTION(([FLASHDeconvHelperStructs::PrecalculatedAveragine] Size getLastIndex(const double mass) const))
{
  double last_index = p_avg_test.getLastIndex(50);
  TEST_EQUAL(last_index, 2);
}
END_SECTION
///


/// testing TopPicItem part is skipped

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST