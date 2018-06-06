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
// $Maintainer: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/CalibrationData.h>
///////////////////////////

#include <OpenMS/MATH/MISC/MathFunctions.h>

using namespace OpenMS;
using namespace std;

START_TEST(Adduct, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

CalibrationData* ptr = nullptr;
CalibrationData* nullPointer = nullptr;
START_SECTION(CalibrationData())
{
	ptr = new CalibrationData();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~CalibrationData())
{
	delete ptr;
}
END_SECTION


CalibrationData cd;
for (Size i = 0; i < 10; ++i)
{
  cd.insertCalibrationPoint(100.100 + i, 200.200 + i, 128.5 + i, 200.0 + i, 1, 66);
  cd.insertCalibrationPoint(120.100 + i + 0.5, 400.200 + i, 128.5 + i, 200.0 + i, 1, 77);
}

START_SECTION(CalDataType::CoordinateType getMZ(Size i) const)
  TEST_REAL_SIMILAR(cd.getMZ(0), 200.200 + 0)
  TEST_REAL_SIMILAR(cd.getMZ(3), 400.200 + 1)
END_SECTION

START_SECTION(CalDataType::CoordinateType getRT(Size i) const)
  TEST_REAL_SIMILAR(cd.getRT(0), 100.100 + 0)
  TEST_REAL_SIMILAR(cd.getRT(3), 120.100 + 1 + 0.5)
END_SECTION
      
START_SECTION(CalDataType::CoordinateType getIntensity(Size i) const)
  TEST_REAL_SIMILAR(cd.getIntensity(0), 128.5 + 0)
  TEST_REAL_SIMILAR(cd.getIntensity(3), 128.5 + 1)
END_SECTION

START_SECTION(const_iterator begin() const)
  TEST_EQUAL(cd.size(), Size(cd.end() - cd.begin()))
END_SECTION
      
START_SECTION(const_iterator end() const)
  TEST_EQUAL(cd.size(), Size(cd.end() - cd.begin()))
END_SECTION

START_SECTION(Size size() const)
  TEST_EQUAL(cd.size(), Size(cd.end() - cd.begin()))
END_SECTION


START_SECTION(bool empty() const)
  TEST_EQUAL(cd.empty(), false)
  TEST_EQUAL(CalibrationData().empty(), true)
END_SECTION

START_SECTION(void clear())
  CalibrationData cd2(cd);
  TEST_EQUAL(cd2.empty(), false)
  cd2.clear();
  TEST_EQUAL(cd2.empty(), true)
END_SECTION
      
START_SECTION(void setUsePPM(bool usePPM))
  cd.setUsePPM(false);
  TEST_EQUAL(cd.usePPM(), false)
  cd.setUsePPM(true);
  TEST_EQUAL(cd.usePPM(), true)
END_SECTION

      
START_SECTION(bool usePPM() const)
  NOT_TESTABLE // tested above
END_SECTION

      
START_SECTION(void insertCalibrationPoint(CalDataType::CoordinateType rt, CalDataType::CoordinateType mz_obs, CalDataType::IntensityType intensity, CalDataType::CoordinateType mz_ref, double weight, int group = -1))
  NOT_TESTABLE // tested above
END_SECTION

      
START_SECTION(Size getNrOfGroups() const)
  TEST_EQUAL(cd.getNrOfGroups(), 2);
  CalibrationData cd2;
  TEST_EQUAL(cd2.getNrOfGroups(), 0);
END_SECTION

      
START_SECTION(CalDataType::CoordinateType getError(Size i) const)
  TEST_REAL_SIMILAR(cd.getError(0), Math::getPPM(200.200 + 0, 200.0))
  TEST_REAL_SIMILAR(cd.getError(3), Math::getPPM(400.200 + 1, 200.0 + 1))
  cd.setUsePPM(false); // use absolute error
  TEST_REAL_SIMILAR(cd.getError(0), (200.200 + 0 - 200.0))
  TEST_REAL_SIMILAR(cd.getError(3), (400.200 + 1 - (200.0 + 1)))
  cd.setUsePPM(true); // reset
END_SECTION

      
START_SECTION(CalDataType::CoordinateType getRefMZ(Size i) const)
  TEST_REAL_SIMILAR(cd.getRefMZ(0), 200.0 + 0)
  TEST_REAL_SIMILAR(cd.getRefMZ(3), 200.0 + 1)
END_SECTION

      
START_SECTION(CalDataType::CoordinateType getWeight(Size i) const)
  TEST_REAL_SIMILAR(cd.getWeight(0), 1.0)
  TEST_REAL_SIMILAR(cd.getWeight(3), 1.0)
END_SECTION

      
START_SECTION(int getGroup(Size i) const)
  TEST_EQUAL(cd.getGroup(0), 66)
  TEST_EQUAL(cd.getGroup(3), 77)
END_SECTION

      
START_SECTION(static StringList getMetaValues())
  TEST_EQUAL(ListUtils::concatenate(CalibrationData::getMetaValues(), ","), "mz_ref,ppm_error,weight")
END_SECTION

      
START_SECTION(CalibrationData median(double rt_left, double rt_right) const)
  CalibrationData m = cd.median(0, 1e6);
  TEST_EQUAL(m.size(), 2); // two medians (of two groups)
  TEST_REAL_SIMILAR(m.getMZ(0),  200.200 + 9.0/2)
  TEST_REAL_SIMILAR(m.getMZ(1),  400.200 + 9.0/2)
END_SECTION

      
START_SECTION(void sortByRT())
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



