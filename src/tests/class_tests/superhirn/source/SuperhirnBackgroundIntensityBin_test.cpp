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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MSPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/BackgroundIntensityBin.h>

///////////////////////////

START_TEST(BackgroundIntensityBin, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

BackgroundIntensityBin* ptr = nullptr;
BackgroundIntensityBin* nullPtr = nullptr;
START_SECTION((BackgroundIntensityBin(double, double)))
	ptr = new BackgroundIntensityBin(300, 12);
	TEST_NOT_EQUAL(ptr,nullPtr)
END_SECTION

START_SECTION((~BackgroundIntensityBin()))
	delete ptr;
END_SECTION

START_SECTION((bool checkBelonging(MSPeak *  )))
  ptr = new BackgroundIntensityBin(300, 12);
  MSPeak *p = new MSPeak();
  TEST_EQUAL(ptr->checkBelonging(p), false);
  delete p;

  MSPeak* p2 = new MSPeak(1, 300, 100);  // (int IN_scan, double IN_mass, float IN_intens)
  p2->set_retention_time(12);
  TEST_EQUAL(ptr->checkBelonging(p2), true);
  delete p2;
END_SECTION

START_SECTION((void addIntensity( double )))
  ptr = new BackgroundIntensityBin(300, 12);
  TEST_EQUAL(ptr->getIntensityMap()->size(), 0);
  ptr->addIntensity(100);
  TEST_EQUAL(ptr->getIntensityMap()->size(), 1);
END_SECTION

START_SECTION((void addMSPeak( MSPeak* )))
  ptr = new BackgroundIntensityBin(300, 12);
  MSPeak* p = new MSPeak(1, 300, 100);  // (int IN_scan, double IN_mass, float IN_intens)
  TEST_EQUAL(ptr->getIntensityMap()->size(), 0);
  ptr->addMSPeak(p);
  TEST_EQUAL(ptr->getIntensityMap()->size(), 1);
  delete p;
END_SECTION

START_SECTION((void processIntensities()))
  ptr = new BackgroundIntensityBin(300, 12);
  ptr->processIntensities();
  TEST_REAL_SIMILAR(ptr->getMean(), 0);
END_SECTION

START_SECTION((std::map<double, double> * getIntensityHist()))
  ptr = new BackgroundIntensityBin(300, 12);
  std::map<double, double>* intensityHistNullPtr = nullptr;
  TEST_NOT_EQUAL(ptr->getIntensityHist(), intensityHistNullPtr);
END_SECTION

START_SECTION((double getMean()))
  ptr = new BackgroundIntensityBin(300, 12);
  ptr->processIntensities();
  //TEST_EQUAL(ptr->getMean(), 0)
  TEST_REAL_SIMILAR(ptr->getMean(), 0);
END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
