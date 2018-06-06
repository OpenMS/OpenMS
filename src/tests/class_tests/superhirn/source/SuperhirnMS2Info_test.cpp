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
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2Info.h>

///////////////////////////

START_TEST(MS2Info, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

MS2Info* ptr = nullptr;
MS2Info* nullPtr = nullptr;

START_SECTION((MS2Info()))
	ptr = new MS2Info();
	TEST_NOT_EQUAL(ptr,nullPtr)
END_SECTION

START_SECTION((~MS2Info()))
	delete ptr;
END_SECTION


ptr = new MS2Info();

START_SECTION( bool check_MODIFICATION() )
  MS2Info ms2info = MS2Info();
  TEST_EQUAL(ms2info.check_MODIFICATION(), false);
END_SECTION

START_SECTION( MS2Info(string IN_AC, string IN_SQ, float IN_PEP) )
  MS2Info ms2info = MS2Info("test1", "test2", 0.95);
END_SECTION

START_SECTION(MS2Info(string IN_AC, string IN_SQ, int IN_CHRG, float IN_PEP))
  MS2Info ms2info = MS2Info("test1", "test2", 2, 0.95);
END_SECTION

START_SECTION(MS2Info(string IN_AC, string IN_SQ, float IN_PEP, int IN_CHRG, int IN_SCAN))
  MS2Info ms2info = MS2Info("test1", "test2", 0.95, 2, 500);
END_SECTION

START_SECTION(MS2Info(int IN_ID))
  MS2Info ms2info = MS2Info(42);
END_SECTION

START_SECTION(MS2Info(const MS2Info & tmp))
  MS2Info ms2info = MS2Info("test1", "test2", 0.95, 2, 500);
  MS2Info other = MS2Info(ms2info);
END_SECTION

START_SECTION(MS2Info::operator=(const MS2Info & tmp))
  MS2Info ms2info = MS2Info("test1", "test2", 0.95, 2, 500);
  MS2Info other = ms2info;
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
