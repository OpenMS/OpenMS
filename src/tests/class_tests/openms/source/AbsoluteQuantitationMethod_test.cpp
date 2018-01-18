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
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitationMethod.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(AbsoluteQuantitationMethod, "$Id$")

/////////////////////////////////////////////////////////////

AbsoluteQuantitationMethod* ptr = nullptr;
AbsoluteQuantitationMethod* nullPointer = nullptr;

START_SECTION(AbsoluteQuantitationMethod())
{
  ptr = new AbsoluteQuantitationMethod();
  TEST_NOT_EQUAL(ptr, nullPointer);
}
END_SECTION

START_SECTION(~AbsoluteQuantitationMethod())
{
  delete ptr;
}
END_SECTION

START_SECTION(all setters and getters)
{
  AbsoluteQuantitationMethod aqm;

  aqm.setComponentName("component");
  aqm.setFeatureName("feature");
  aqm.setISName("IS");
  aqm.setLLOD(1.2);
  aqm.setULOD(3.4);
  aqm.setLLOQ(5.6);
  aqm.setULOQ(7.8);
  aqm.setNPoints(9);
  aqm.setCorrelationCoefficient(0.44);
  aqm.setConcentrationUnits("uM");
  aqm.setTransformationModel("TransformationModelLinear");
  Param params1;
  params1.setValue("slope", 1);
  aqm.setTransformationModelParams(params1);

  TEST_EQUAL(aqm.getComponentName(), "component")
  TEST_EQUAL(aqm.getFeatureName(), "feature")
  TEST_EQUAL(aqm.getISName(), "IS")
  TEST_REAL_SIMILAR(aqm.getLLOD(), 1.2)
  TEST_REAL_SIMILAR(aqm.getULOD(), 3.4)
  TEST_REAL_SIMILAR(aqm.getLLOQ(), 5.6)
  TEST_REAL_SIMILAR(aqm.getULOQ(), 7.8)
  TEST_EQUAL(aqm.getNPoints(), 9)
  TEST_REAL_SIMILAR(aqm.getCorrelationCoefficient(), 0.44)
  TEST_EQUAL(aqm.getConcentrationUnits(), "uM")
  TEST_EQUAL(aqm.getTransformationModel(), "TransformationModelLinear")
  Param params2 = aqm.getTransformationModelParams();
  TEST_EQUAL(params2.getValue("slope"), 1)
}
END_SECTION

START_SECTION(bool checkLOD(const double value) const)
{
  AbsoluteQuantitationMethod aqm;
  const double value = 2.0;
  aqm.setLLOD(0.0);
  aqm.setULOD(4.0);
  TEST_EQUAL(aqm.checkLOD(value), true);
  aqm.setULOD(1.0);
  TEST_EQUAL(aqm.checkLOD(value), false);
  aqm.setLLOD(3.0);
  aqm.setULOD(4.0);
  TEST_EQUAL(aqm.checkLOD(value), false);
}
END_SECTION

START_SECTION(bool checkLOQ(const double value) const)
{
  AbsoluteQuantitationMethod aqm;
  const double value = 2.0;
  aqm.setLLOQ(0.0);
  aqm.setULOQ(4.0);
  TEST_EQUAL(aqm.checkLOQ(value), true);
  aqm.setULOQ(1.0);
  TEST_EQUAL(aqm.checkLOQ(value), false);
  aqm.setLLOQ(3.0);
  aqm.setULOQ(4.0);
  TEST_EQUAL(aqm.checkLOQ(value), false);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST