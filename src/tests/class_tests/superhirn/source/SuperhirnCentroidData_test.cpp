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
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/RawData.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/CentroidData.h>

///////////////////////////

START_TEST(CentroidData, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

CentroidData* ptr = nullptr;
CentroidData* nullPtr = nullptr;

vector<double>* centroidMasses = new vector<double>(); // Centroided masses
vector<double>* centroidIntens = new vector<double>(); // Centroided intensities
boost::shared_ptr<RawData> raw(new RawData(*centroidMasses, *centroidIntens));


START_SECTION((CentroidData()))
	ptr = new CentroidData(1, raw, true);
	TEST_NOT_EQUAL(ptr,nullPtr)
END_SECTION

START_SECTION((~CentroidData()))
	delete ptr;
END_SECTION

ptr = new CentroidData(1, raw, true);


START_SECTION((setAndGet()))
ptr = new CentroidData(1, raw, true);

centroidMasses->push_back(1.0);
centroidIntens->push_back(2.0);

ptr->set(*centroidMasses, *centroidIntens);

list<CentroidPeak>* centroidPeaks = new list<CentroidPeak>();
ptr->get(*centroidPeaks);
TEST_EQUAL(centroidPeaks->size(), 1);
END_SECTION

START_SECTION(TODO)

// CentroidPeak
// 			void show_info();
// 			void subtractIntensity(double);

// DeconvPeak
// std::vector<CentroidPeak> getIsotopicPeaks();
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
