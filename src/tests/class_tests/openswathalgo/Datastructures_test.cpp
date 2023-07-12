// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include "OpenMS/OPENSWATHALGO/OpenSwathAlgoConfig.h"
#include "OpenMS/OPENSWATHALGO/DATAACCESS/DataStructures.h"

#include <OpenMS/CONCEPT/ClassTest.h>
using namespace OpenMS;
using namespace std;
using namespace OpenSwath;

///////////////////////////

START_TEST(DataStructures, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(OSSpectrum_empty)
{
  OSSpectrum s;

  TEST_EQUAL (s.getMZArray() == nullptr, false)
  TEST_EQUAL (s.getIntensityArray() == nullptr, false)
  TEST_EQUAL (s.getDriftTimeArray() == nullptr, true)

  TEST_EQUAL (s.getMZArray()->data.size(), 0)
  TEST_EQUAL (s.getIntensityArray()->data.size(), 0)
}
END_SECTION

START_SECTION(OSSpectrum_data)
{
  OSSpectrum s;

  BinaryDataArrayPtr mz(new BinaryDataArray);
  mz->data.push_back(1.5);
  BinaryDataArrayPtr inten(new BinaryDataArray);
  inten->data.push_back(100.1);
  BinaryDataArrayPtr im(new BinaryDataArray);
  im->data.push_back(300.1);
  im->description = "Ion Mobility"; // old format

  s.setMZArray(mz);
  s.setIntensityArray(inten);
  s.getDataArrays().push_back(im);

  TEST_EQUAL (s.getMZArray() == nullptr, false)
  TEST_EQUAL (s.getIntensityArray() == nullptr, false)
  TEST_EQUAL (s.getDriftTimeArray() == nullptr, false)

  TEST_EQUAL (s.getMZArray()->data.size(), 1)
  TEST_EQUAL (s.getIntensityArray()->data.size(), 1)
  TEST_EQUAL (s.getDriftTimeArray()->data.size(), 1)

  TEST_REAL_SIMILAR (s.getMZArray()->data[0], 1.5)
  TEST_REAL_SIMILAR (s.getIntensityArray()->data[0], 100.1)
  TEST_REAL_SIMILAR (s.getDriftTimeArray()->data[0], 300.1)
}
END_SECTION

START_SECTION(OSSpectrum_data_2)
{
  OSSpectrum s;

  BinaryDataArrayPtr mz(new BinaryDataArray);
  mz->data.push_back(1.5);
  BinaryDataArrayPtr inten(new BinaryDataArray);
  inten->data.push_back(100.1);
  BinaryDataArrayPtr im(new BinaryDataArray);
  im->data.push_back(300.1);
  im->description = "Ion Mobility (MS:1002476)"; // new format

  s.setMZArray(mz);
  s.setIntensityArray(inten);
  s.getDataArrays().push_back(im);

  TEST_EQUAL (s.getMZArray() == nullptr, false)
  TEST_EQUAL (s.getIntensityArray() == nullptr, false)
  TEST_EQUAL (s.getDriftTimeArray() == nullptr, false)

  TEST_EQUAL (s.getMZArray()->data.size(), 1)
  TEST_EQUAL (s.getIntensityArray()->data.size(), 1)
  TEST_EQUAL (s.getDriftTimeArray()->data.size(), 1)

  TEST_REAL_SIMILAR (s.getMZArray()->data[0], 1.5)
  TEST_REAL_SIMILAR (s.getIntensityArray()->data[0], 100.1)
  TEST_REAL_SIMILAR (s.getDriftTimeArray()->data[0], 300.1)

  s.getDataArrays().back()->description = "";
  TEST_EQUAL (s.getDriftTimeArray() == nullptr, true)
  s.getDataArrays().back()->description = "Ion Mobility (blah)";
  TEST_EQUAL (s.getDriftTimeArray() == nullptr, false)
  s.getDataArrays().back()->description = "Ion mOBILITY (blah)"; // wrong
  TEST_EQUAL (s.getDriftTimeArray() == nullptr, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
