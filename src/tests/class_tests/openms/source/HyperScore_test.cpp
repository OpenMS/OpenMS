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

///////////////////////////
#include <OpenMS/ANALYSIS/RNPXL/HyperScore.h>
///////////////////////////

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>

using namespace OpenMS;
using namespace std;

START_TEST(HyperScore, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

HyperScore* ptr = nullptr;
HyperScore* null_ptr = nullptr;

TheoreticalSpectrumGenerator tsg;
Param param = tsg.getParameters();
param.setValue("add_metainfo", "true");
tsg.setParameters(param);

START_SECTION(HyperScore())
{
  ptr = new HyperScore();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~HyperScore())
{
  delete ptr;
}
END_SECTION

START_SECTION((static double compute(double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, const PeakSpectrum &exp_spectrum, const RichPeakSpectrum &theo_spectrum)))
{
  PeakSpectrum exp_spectrum;
  PeakSpectrum theo_spectrum;

  AASequence peptide = AASequence::fromString("PEPTIDE");
  
  // empty spectrum
  tsg.getSpectrum(theo_spectrum, peptide, 1, 1);
  TEST_REAL_SIMILAR(HyperScore::compute(0.1, false, exp_spectrum, theo_spectrum), 0.0);

  // full match, 11 identical masses, identical intensities (=1)
  tsg.getSpectrum(exp_spectrum, peptide, 1, 1);
  TEST_REAL_SIMILAR(HyperScore::compute(0.1, false, exp_spectrum, theo_spectrum), 13.8516496);
  TEST_REAL_SIMILAR(HyperScore::compute(10, true, exp_spectrum, theo_spectrum), 13.8516496);

  exp_spectrum.clear(true);
  theo_spectrum.clear(true);

  // no match
  tsg.getSpectrum(exp_spectrum, peptide, 1, 3);
  tsg.getSpectrum(theo_spectrum, AASequence::fromString("YYYYYY"), 1, 3);
  TEST_REAL_SIMILAR(HyperScore::compute(1e-5, false, exp_spectrum, theo_spectrum), 0.0);
  
  exp_spectrum.clear(true);
  theo_spectrum.clear(true);

  // full match, 33 identical masses, identical intensities (=1)
  tsg.getSpectrum(exp_spectrum, peptide, 1, 3);
  tsg.getSpectrum(theo_spectrum, peptide, 1, 3);
  TEST_REAL_SIMILAR(HyperScore::compute(0.1, false, exp_spectrum, theo_spectrum), 67.8210771);
  TEST_REAL_SIMILAR(HyperScore::compute(10, true, exp_spectrum, theo_spectrum), 67.8210771);

  // full match if ppm tolerance and partial match for Da tolerance
  for (Size i = 0; i < theo_spectrum.size(); ++i)
  {
    double mz = pow( theo_spectrum[i].getMZ(), 2);
    exp_spectrum[i].setMZ(mz);
    theo_spectrum[i].setMZ(mz + 9 * 1e-6 * mz); // +9 ppm error
  }

  TEST_REAL_SIMILAR(HyperScore::compute(0.1, false, exp_spectrum, theo_spectrum), 3.401197);
  TEST_REAL_SIMILAR(HyperScore::compute(10, true, exp_spectrum, theo_spectrum), 67.8210771);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

