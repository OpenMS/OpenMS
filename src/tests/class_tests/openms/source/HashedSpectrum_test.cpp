// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/COMPARISON/SPECTRA/HashedSpectrum.h>

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>

using namespace OpenMS;

START_TEST(MultiplexFiltering, "$Id$")

// read data
MSExperiment<Peak1D> exp;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MultiplexFiltering.mzML"), exp);
exp.updateRanges();
MSExperiment<Peak1D>::Iterator it_rt = exp.begin();    // RT = 1595.25192 sec

// set parameters
double mz_bin = 1000.0;
double mz_tolerance = 5.0;
bool mz_unit_ppm = true;

START_SECTION(HashedSpectrum(const MSSpectrum<Peak1D>& raw_spectrum, const double mz_bin, const bool mz_unit_ppm))
    HashedSpectrum spectrum(*it_rt, mz_bin, mz_unit_ppm);
    TEST_REAL_SIMILAR(spectrum.getMzBin(), 1000.0);
END_SECTION

START_SECTION(double getMzBin() const)
    HashedSpectrum spectrum(*it_rt, mz_bin, mz_unit_ppm);
    TEST_REAL_SIMILAR(spectrum.getMzBin(), 1000.0);
END_SECTION

START_SECTION(bool getMzUnitPpm() const)
    HashedSpectrum spectrum(*it_rt, mz_bin, mz_unit_ppm);
    TEST_EQUAL(spectrum.getMzUnitPpm(), true);
END_SECTION

// There is a peak within the tolerance.
START_SECTION(MSSpectrum<Peak1D>::ConstIterator findNearest(const double mz, const double mz_tolerance, const bool mz_unit_ppm) const)
    HashedSpectrum spectrum(*it_rt, mz_bin, mz_unit_ppm);
    TEST_REAL_SIMILAR(spectrum.findNearest(751.41, mz_tolerance, mz_unit_ppm)->getMZ(), 751.408386230469);
END_SECTION

// There is no peak within the tolerance.
START_SECTION(MSSpectrum<Peak1D>::ConstIterator findNearest(const double mz, const double mz_tolerance, const bool mz_unit_ppm) const)
    HashedSpectrum spectrum(*it_rt, mz_bin, mz_unit_ppm);
    TEST_EQUAL(spectrum.findNearest(822.0, mz_tolerance, mz_unit_ppm) == it_rt->end(), true);
END_SECTION

// m/z outside the range of the spectrum
START_SECTION(MSSpectrum<Peak1D>::ConstIterator findNearest(const double mz, const double mz_tolerance, const bool mz_unit_ppm) const)
    HashedSpectrum spectrum(*it_rt, mz_bin, mz_unit_ppm);
    TEST_EQUAL(spectrum.findNearest(200.0, mz_tolerance, mz_unit_ppm) == it_rt->end(), true);
END_SECTION

// m/z outside the range of the spectrum
START_SECTION(MSSpectrum<Peak1D>::ConstIterator findNearest(const double mz, const double mz_tolerance, const bool mz_unit_ppm) const)
    HashedSpectrum spectrum(*it_rt, mz_bin, mz_unit_ppm);
    TEST_EQUAL(spectrum.findNearest(3000.0, mz_tolerance, mz_unit_ppm) == it_rt->end(), true);
END_SECTION

END_TEST
