// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Authors: Timo Sachsenberg$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/DTAFile.h>

///////////////////////////
#include <OpenMS/ANALYSIS/RNPXL/PScore.h>
///////////////////////////
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>

using namespace OpenMS;
using namespace std;

START_TEST(PScore, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PScore* ptr = nullptr;
PScore* null_ptr = nullptr;
START_SECTION(PScore())
{
	ptr = new PScore();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~PScore())
{
	delete ptr;
}
END_SECTION

START_SECTION((static std::vector<Size> calculateIntensityRankInMZWindow(const std::vector< double > &mz, const std::vector< double > &intensities, double mz_window)))
{
  std::vector< double > mz;
  std::vector< double > intensities;

  // simple increasing sequence
  for (Size m = 0; m < 100; ++m)
  {
    mz.push_back(m);
    intensities.push_back(m);
  }

  // test window size
  std::vector<Size> ranks = PScore::calculateIntensityRankInMZWindow(mz, intensities, 9.9);
  TEST_EQUAL(ranks.size(), mz.size());

  for (Size i = 0; i != ranks.size() - 4; ++i)
  {
    TEST_EQUAL(ranks[i], 4);
  }

  ranks = PScore::calculateIntensityRankInMZWindow(mz, intensities, 10.1);
  TEST_EQUAL(ranks.size(), mz.size());

  for (Size i = 0; i != ranks.size() - 5; ++i)
  {
    TEST_EQUAL(ranks[i], 5);
  }
}
END_SECTION

START_SECTION((static std::vector<std::vector<Size> > calculateRankMap(const PeakMap &peak_map, double mz_window=100)))
{
  // Convenience function. Calculations tested via calculateIntensityRankInMZWindow
}
END_SECTION

START_SECTION((static std::map<Size, PeakSpectrum> calculatePeakLevelSpectra(const PeakSpectrum &spec, const std::vector< Size > &ranks, Size min_level=2, Size max_level=10)))
{
  DTAFile dta_file;
  PeakSpectrum spec;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("PScore_test.dta"), spec);
  vector<double> mz, intensities;
  for (Size i = 0; i != spec.size(); ++i)
  {
    mz.push_back(spec[i].getMZ());
    intensities.push_back(spec[i].getIntensity());
  }

  std::vector<Size> ranks = PScore::calculateIntensityRankInMZWindow(mz, intensities, 100.0);
  TEST_EQUAL(ranks.size(), spec.size())

  std::map<Size, PeakSpectrum > pls = PScore::calculatePeakLevelSpectra(spec, ranks, 0, 1);
  TEST_EQUAL(pls.size(), 2)

  // top intensity peaks in +- 50 Th neighborhood
  TEST_REAL_SIMILAR(pls[0][0].getMZ(), 169.65);
  TEST_REAL_SIMILAR(pls[0][1].getMZ(), 231.51);
  TEST_REAL_SIMILAR(pls[0][2].getMZ(), 362.22);
  TEST_REAL_SIMILAR(pls[0][3].getMZ(), 508.47);
  TEST_REAL_SIMILAR(pls[0][4].getMZ(), 579.61);
  TEST_REAL_SIMILAR(pls[0][5].getMZ(), 629.66);
  TEST_REAL_SIMILAR(pls[0][6].getMZ(), 712.18);

  // top two intensity peaks in +- 50 Th neighborhood
  TEST_REAL_SIMILAR(pls[1][0].getMZ(), 149.93);
  TEST_REAL_SIMILAR(pls[1][1].getMZ(), 169.65);
  TEST_REAL_SIMILAR(pls[1][2].getMZ(), 231.51);
  TEST_REAL_SIMILAR(pls[1][3].getMZ(), 263.88);
  TEST_REAL_SIMILAR(pls[1][4].getMZ(), 318.38);
  TEST_REAL_SIMILAR(pls[1][5].getMZ(), 362.22);
  TEST_REAL_SIMILAR(pls[1][6].getMZ(), 389.84);
  TEST_REAL_SIMILAR(pls[1][7].getMZ(), 489.86);
  TEST_REAL_SIMILAR(pls[1][8].getMZ(), 508.47);
  TEST_REAL_SIMILAR(pls[1][9].getMZ(), 562.72);
  TEST_REAL_SIMILAR(pls[1][10].getMZ(), 579.61);
  TEST_REAL_SIMILAR(pls[1][11].getMZ(), 629.66);
  TEST_REAL_SIMILAR(pls[1][12].getMZ(), 712.18);
}
END_SECTION

START_SECTION((static double computePScore(double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, const std::map< Size, PeakSpectrum > &peak_level_spectra, const std::vector< PeakSpectrum > &theo_spectra, double mz_window=100.0)))
{
  // Convenience function. Calculations tested via computePScore
}
END_SECTION

START_SECTION((static double computePScore(double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, const std::map< Size, PeakSpectrum > &peak_level_spectra, const PeakSpectrum &theo_spectrum, double mz_window=100.0)))
{
  DTAFile dta_file;
  PeakSpectrum spec;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("PScore_test.dta"), spec);
  vector<double> mz, intensities;
  for (Size i = 0; i != spec.size(); ++i)
  {
    mz.push_back(spec[i].getMZ());
    intensities.push_back(spec[i].getIntensity());
  }

  PeakSpectrum theo_spec;
  for (Size i = 0; i != spec.size(); ++i)
  {
    Peak1D p;
    p.setMZ(spec[i].getMZ());
    p.setIntensity(spec[i].getIntensity());
    theo_spec.push_back(p);
  }

  std::vector<Size> ranks = PScore::calculateIntensityRankInMZWindow(mz, intensities, 100.0);
  std::map<Size, PeakSpectrum > pls = PScore::calculatePeakLevelSpectra(spec, ranks, 0, 0);

  double pscore_all_match_top_1 = PScore::computePScore(0.1, true, pls, theo_spec);
  pls = PScore::calculatePeakLevelSpectra(spec, ranks, 0, 1);
  double pscore_all_match_top_2 = PScore::computePScore(0.1, true, pls, theo_spec);

  TEST_REAL_SIMILAR(pscore_all_match_top_1, 83.867454)
  TEST_REAL_SIMILAR(pscore_all_match_top_2, 154.682242)

  AASequence peptide = AASequence::fromString("IFSQVGK");
  TheoreticalSpectrumGenerator tg;
  Param param(tg.getParameters());
  param.setValue("add_first_prefix_ion", "true");
  tg.setParameters(param);
  spec.clear(true);
  tg.getSpectrum(spec, peptide, 1, 1);
  TEST_EQUAL(spec.size(), 12)

  mz.clear();
  intensities.clear();

  for (Size i = 0; i != spec.size(); ++i)
  {
    mz.push_back(spec[i].getMZ());
    intensities.push_back(spec[i].getIntensity());
  }

  ranks = PScore::calculateIntensityRankInMZWindow(mz, intensities, 100.0);
  pls = PScore::calculatePeakLevelSpectra(spec, ranks, 0, 0);
  double all_match = PScore::computePScore(0.1, true, pls, spec);
  TEST_REAL_SIMILAR(all_match, 240)
}
END_SECTION

START_SECTION((static double massCorrectionTerm(double mass)))
{
  // Not tested
}
END_SECTION

START_SECTION((static double cleavageCorrectionTerm(Size cleavages, bool consecutive_cleavage)))
{
  // Not tested
}
END_SECTION

START_SECTION((static double modificationCorrectionTerm(Size modifications)))
{
  // Not tested
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



