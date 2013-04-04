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

