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
// $Maintainer: Timo Sachsenberg $
// $Authors: Bastian Blank $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/PeakWidthEstimator.h>

#include <OpenMS/FORMAT/MzMLFile.h>

using namespace OpenMS;

START_TEST(PeakWidthEstimator, "$Id$")

MSExperiment<> input;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_orbitrap.mzML"), input);

TOLERANCE_RELATIVE(1.001);

START_SECTION(static void estimateSpectrumFWHM(const MSSpectrum<> &, std::set<boost::tuple<DoubleReal, DoubleReal, DoubleReal> > &))
{
  typedef std::set<boost::tuple<DoubleReal, DoubleReal, DoubleReal> > Fwhm;
  Fwhm fwhm;
  PeakWidthEstimator::estimateSpectrumFWHM(input[0], fwhm);
  TEST_EQUAL(fwhm.size(), 155);
  Fwhm::const_reverse_iterator it = fwhm.rbegin();
  TEST_REAL_SIMILAR(it->get<0>(), 202394.);
  TEST_REAL_SIMILAR(it->get<1>(), 591.358);
  TEST_REAL_SIMILAR(it->get<2>(), .010647);
}
END_SECTION

START_SECTION(static Result estimateFWHM(const MSExperiment<> &))
{
  PeakWidthEstimator::Result r(PeakWidthEstimator::estimateFWHM(input));
  TEST_REAL_SIMILAR(r.c0, -14.15849);
  TEST_REAL_SIMILAR(r.c1, 1.50632);
}
END_SECTION

END_TEST

