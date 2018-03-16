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
// $Maintainer: Florian Zeller$
// $Authors: Florian Zeller$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include "test_config.h"

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/CentroidPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MSPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2Fragment.h>

///////////////////////////

START_TEST(MS2Fragment, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

MS2Fragment* ptr = nullptr;
MS2Fragment* nullPtr = nullptr;

START_SECTION((MS2Fragment()))
  ptr = new MS2Fragment();
  TEST_NOT_EQUAL(ptr,nullPtr)
END_SECTION

START_SECTION((~MS2Fragment()))
  delete ptr;
END_SECTION

ptr = new MS2Fragment();

START_SECTION(MS2Fragment(double iPrecursorMZ, int iPrecursorCHRG, double iTR, int iScan, int iZ, double iFragmentMZ, double iIntensityArea,
                int iScanStart, int iScanEnd, double iTrStart, double iTrEnd))
  MS2Fragment ms2fragment = MS2Fragment(400.25, 2, 1800, 500, 1, 500.25, 1000, 
                                        499, 501, 1750, 1850);
END_SECTION

START_SECTION(MS2Fragment(double iPrecursorMZ, int iPrecursorCHRG, double iTR, int iScan, int iZ, double iFragmentMZ, double iIntensityArea))
  MS2Fragment ms2fragment = MS2Fragment(400.25, 2, 1800, 500, 1, 500.25, 1000);
END_SECTION

START_SECTION(MS2Fragment(const MS2Fragment & tmp))
  MS2Fragment ms2fragment = MS2Fragment(400.25, 2, 1800, 500, 1, 500.25, 1000);
  MS2Fragment other = MS2Fragment(ms2fragment);
END_SECTION

START_SECTION(MS2Fragment(MS2Fragment::operator=(const MS2Fragment & tmp)))
  MS2Fragment ms2fragment = MS2Fragment(400.25, 2, 1800, 500, 1, 500.25, 1000);
  MS2Fragment other = ms2fragment;
END_SECTION

START_SECTION(void MS2Fragment::show_info())
  MS2Fragment ms2fragment = MS2Fragment(400.25, 2, 1800, 500, 1, 500.25, 1000);
  ms2fragment.show_info();
END_SECTION

START_SECTION(( void setPrecursorMZ(double iMZ)))
  MS2Fragment ms2fragment = MS2Fragment();
  ms2fragment.setPrecursorMZ(500.0);
  TEST_REAL_SIMILAR(ms2fragment.getPrecursorMZ(), 500.0);
END_SECTION

START_SECTION((double getFragmentPeakArea()))
  MS2Fragment ms2fragment = MS2Fragment();
  ms2fragment.setFragmentPeakArea(500.0);
  TEST_REAL_SIMILAR(ms2fragment.getFragmentPeakArea(), 500.0);
END_SECTION

START_SECTION((double getFragmentMz()))
  MS2Fragment ms2fragment = MS2Fragment();
  ms2fragment.setFragmentMz(500.0);
  TEST_REAL_SIMILAR(ms2fragment.getFragmentMz(), 500.0);
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
