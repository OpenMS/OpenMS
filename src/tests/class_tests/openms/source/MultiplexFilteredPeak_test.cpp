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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteredPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexSatelliteCentroided.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexSatelliteProfile.h>

using namespace OpenMS;

START_TEST(MultiplexFilteredPeak, "$Id$")

MultiplexFilteredPeak* nullPointer = nullptr;
MultiplexFilteredPeak* ptr;

START_SECTION(MultiplexFilteredPeak())
    MultiplexFilteredPeak peak(654.32, 2345.67, 24, 42);
    TEST_EQUAL(peak.getMZidx(), 24);
    ptr = new MultiplexFilteredPeak(654.32, 2345.67, 24, 42);
    TEST_NOT_EQUAL(ptr, nullPointer);
    delete ptr;
END_SECTION

MultiplexFilteredPeak peak(654.32, 2345.67, 24, 42);
MultiplexSatelliteCentroided satellite_centroided(26, 44);
MultiplexSatelliteProfile satellite_profile(2346.67, 655.32, 1000.0);
peak.addSatellite(25, 43, 3);
peak.addSatellite(satellite_centroided, 3);
peak.addSatelliteProfile(2347.67, 656.32, 1010.0, 4);
peak.addSatelliteProfile(satellite_profile, 4);
size_t n;

START_SECTION(double getMZ())
  TEST_REAL_SIMILAR(peak.getMZ(), 654.32);
END_SECTION

START_SECTION(double getRT())
  TEST_REAL_SIMILAR(peak.getRT(), 2345.67);
END_SECTION

START_SECTION(size_t getMZidx())
  TEST_EQUAL(peak.getMZidx(), 24);
END_SECTION

START_SECTION(size_t getRTidx())
  TEST_EQUAL(peak.getRTidx(), 42);
END_SECTION

START_SECTION(void addSatellite(size_t rt_idx, size_t mz_idx, size_t pattern_idx))
  n = peak.getSatellites().size();
  peak.addSatellite(25, 43, 3);
  TEST_EQUAL(peak.getSatellites().size(), n+1);
END_SECTION

START_SECTION(void addSatellite(const MultiplexSatelliteCentroided& satellite, size_t pattern_idx))
  n = peak.getSatellites().size();
  MultiplexSatelliteCentroided satellite_centroided_temp(27, 45);
  peak.addSatellite(satellite_centroided_temp, 3);
  TEST_EQUAL(peak.getSatellites().size(), n+1);
END_SECTION

START_SECTION(void addSatelliteProfile(double rt, double mz, double intensity, size_t pattern_idx))
  n = peak.getSatellitesProfile().size();
  peak.addSatelliteProfile(2348.67, 657.32, 1020.0, 5);
  TEST_EQUAL(peak.getSatellitesProfile().size(), n+1);
END_SECTION

START_SECTION(void addSatelliteProfile(const MultiplexSatelliteProfile& satellite, size_t pattern_idx))
  n = peak.getSatellitesProfile().size();
  MultiplexSatelliteProfile satellite_profile_temp(2349.67, 658.32, 1030.0);
  peak.addSatelliteProfile(satellite_profile_temp, 6);
  TEST_EQUAL(peak.getSatellitesProfile().size(), n+1);
END_SECTION

START_SECTION(getSatellites())
  TEST_EQUAL(peak.getSatellites().size(), 4);
END_SECTION

START_SECTION(getSatellitesProfile())
  TEST_EQUAL(peak.getSatellitesProfile().size(), 4);
END_SECTION

START_SECTION(size_t size())
  TEST_EQUAL(peak.size(), 4);
END_SECTION

START_SECTION(size_t sizeProfile())
  TEST_EQUAL(peak.sizeProfile(), 4);
END_SECTION

END_TEST
