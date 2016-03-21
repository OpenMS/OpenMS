// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Witold Wolski $
// $Authors: Witold Wolski $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/OPENSWATH/DIAPrescoring.h>
#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/MockObjects.h"

using namespace std;
using namespace OpenMS;
using namespace OpenSwath;


START_TEST(DiaPrescore2, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DiaPrescore* ptr = 0;
DiaPrescore* nullPointer = 0;

START_SECTION(DiaPrescore())
{
  ptr = new DiaPrescore();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~DiaPrescore())
{
  delete ptr;
}
END_SECTION


START_SECTION ( testscorefunction)
{

  OpenSwath::LightTransition mock_tr1;
  mock_tr1.product_mz = 500.;
  mock_tr1.product_charge = 1;
  mock_tr1.transition_name = "group1";
  mock_tr1.library_intensity = 5.;

  OpenSwath::LightTransition mock_tr2;
  mock_tr2.product_mz = 600.;
  mock_tr2.product_charge = 1;
  mock_tr2.transition_name = "group2";
  mock_tr2.library_intensity = 5.;

  OpenSwath::SpectrumPtr sptr = (OpenSwath::SpectrumPtr)(new OpenSwath::Spectrum);
  std::vector<OpenSwath::BinaryDataArrayPtr> binaryDataArrayPtrs;
  OpenSwath::BinaryDataArrayPtr data1(new OpenSwath::BinaryDataArray);
  OpenSwath::BinaryDataArrayPtr data2(new OpenSwath::BinaryDataArray);

  static const double arr1[] = {
    /*
    10, 20, 50, 100, 50, 20, 10, // peak at 499
    3, 7, 15, 30, 15, 7, 3,      // peak at 500
    1, 3, 9, 15, 9, 3, 1,        // peak at 501
    3, 9, 3,                     // peak at 502
    */

    10, 20, 50, 100, 50, 20, 10, // peak at 600
    3, 7, 15, 30, 15, 7, 3,      // peak at 601
    1, 3, 9, 15, 9, 3, 1,        // peak at 602
    3, 9, 3                      // peak at 603
  };
  std::vector<double> intensity (arr1, arr1 + sizeof(arr1) / sizeof(double) );
  static const double arr2[] = {
    /*
    498.97, 498.98, 498.99, 499.0, 499.01, 499.02, 499.03,
    499.97, 499.98, 499.99, 500.0, 500.01, 500.02, 500.03,
    500.97, 500.98, 500.99, 501.0, 501.01, 501.02, 501.03,
    501.99, 502.0, 502.01,
    */
    599.97, 599.98, 599.99, 600.0, 600.01, 600.02, 600.03,
    600.97, 600.98, 600.99, 601.0, 601.01, 601.02, 601.03,
    601.97, 601.98, 601.99, 602.0, 602.01, 602.02, 602.03,
    602.99, 603.0, 603.01
  };
  std::vector<double> mz (arr2, arr2 + sizeof(arr2) / sizeof(double) );
  data1->data = mz;
  data2->data = intensity;
  sptr->setMZArray( data1);
  sptr->setIntensityArray( data2);

  std::vector<OpenSwath::LightTransition> transitions;
  transitions.push_back(mock_tr1);
  transitions.push_back(mock_tr2);

  DiaPrescore diaprescore(0.05);
  double manhattan = 0., dotprod = 0.;
  diaprescore.score(sptr, transitions , dotprod, manhattan);
  //std::cout << "dotprod : " << dotprod << std::endl;
  //std::cout << "manhattan : " << manhattan << std::endl;
  // >> exp = [240, 74, 39, 15, 0]
  // >> theo = [1, 0.325757771553019, 0.0678711748364005, 0.0105918703087134, 0.00134955223787482]
  // >> from scipy.stats.stats import pearsonr
  // >> pearsonr(exp, theo)
  // (0.99463189043051314, 0.00047175434098498532)
  //
  TEST_REAL_SIMILAR(dotprod, 0.644473950768828)
  TEST_REAL_SIMILAR(manhattan, 1.00199893289589)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

