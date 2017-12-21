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
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2ConsensusSpectrum.h>

///////////////////////////

START_TEST(MS2ConsensusSpectrum, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

MS2ConsensusSpectrum* ptr = nullptr;
MS2ConsensusSpectrum* nullPtr = nullptr;

START_SECTION((MS2ConsensusSpectrum()))
	ptr = new MS2ConsensusSpectrum();
	TEST_NOT_EQUAL(ptr,nullPtr)
END_SECTION

START_SECTION((~MS2ConsensusSpectrum()))
	delete ptr;
END_SECTION

ptr = new MS2ConsensusSpectrum();

START_SECTION( MS2ConsensusSpectrum(double iPrecursorMZ, double iTR, int iChrg, int iApexScan) )
  MS2ConsensusSpectrum constr_test = MS2ConsensusSpectrum(400.25, 1800, 2, 5);
END_SECTION

START_SECTION( double getPrecursorMZ())
  MS2ConsensusSpectrum constr_test = MS2ConsensusSpectrum(400.25, 1800, 2, 5);
  TEST_REAL_SIMILAR(constr_test.getPrecursorMZ(), 400.25)
END_SECTION

START_SECTION( double getTR())
  MS2ConsensusSpectrum constr_test = MS2ConsensusSpectrum(400.25, 1800, 2, 5);
  TEST_REAL_SIMILAR(constr_test.getTR(), 1800.0)
END_SECTION

START_SECTION( double getStartTR())
  MS2ConsensusSpectrum constr_test = MS2ConsensusSpectrum(400.25, 1800, 2, 5);
  TEST_REAL_SIMILAR(constr_test.getStartTR(), 1800.0)
END_SECTION

START_SECTION( double getEndTR())
  MS2ConsensusSpectrum constr_test = MS2ConsensusSpectrum(400.25, 1800, 2, 5);
  TEST_REAL_SIMILAR(constr_test.getEndTR(), 1800.0)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
