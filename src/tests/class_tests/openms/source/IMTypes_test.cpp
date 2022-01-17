// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/IONMOBILITY/IMDataConverter.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(MSRunIMSplitter, "$Id$")

/////////////////////////////////////////////////////////////

IMTypes* e_ptr = nullptr;
IMTypes* e_nullPointer = nullptr;

START_SECTION((IMTypes()))
	e_ptr = new IMTypes;
  TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION((~IMTypes()))
	delete e_ptr;
END_SECTION

START_SECTION((DriftTimeUnit toDriftTimeUnit(const String& dtu_string)))
  TEST_EQUAL(toDriftTimeUnit("<NONE>") == DriftTimeUnit::NONE, true)
  for (size_t i = 0; i < (size_t)DriftTimeUnit::SIZE_OF_DRIFTTIMEUNIT; ++i)
  {
    TEST_EQUAL((size_t)toDriftTimeUnit(NamesOfDriftTimeUnit[i]), i)
  }
  TEST_EXCEPTION(Exception::InvalidValue, toDriftTimeUnit("haha"));
END_SECTION

START_SECTION(const String& toString(const DriftTimeUnit value))
  TEST_EQUAL(toString(DriftTimeUnit::NONE), "<NONE>")
  for (size_t i = 0; i < (size_t)DriftTimeUnit::SIZE_OF_DRIFTTIMEUNIT; ++i)
  {
    TEST_EQUAL(toString(DriftTimeUnit(i)), NamesOfDriftTimeUnit[i])
  }
  TEST_EXCEPTION(Exception::InvalidValue, toString(DriftTimeUnit::SIZE_OF_DRIFTTIMEUNIT));
END_SECTION


START_SECTION((IMFormat toIMFormat(const String& IM_format)))
  TEST_EQUAL(toIMFormat("none") == IMFormat::NONE, true)
  for (size_t i = 0; i < (size_t) IMFormat::SIZE_OF_IMFORMAT; ++i)
  {
    TEST_EQUAL((size_t) toIMFormat(NamesOfIMFormat[i]), i)
  }
  TEST_EXCEPTION(Exception::InvalidValue, toIMFormat("haha"));
END_SECTION

START_SECTION(const String& toString(const IMFormat value))
  TEST_EQUAL(toString(IMFormat::NONE), "none")
  for (size_t i = 0; i < (size_t)IMFormat::SIZE_OF_IMFORMAT; ++i)
  {
    TEST_EQUAL(toString(IMFormat(i)), NamesOfIMFormat[i])
  }
  TEST_EXCEPTION(Exception::InvalidValue, toString(IMFormat::SIZE_OF_IMFORMAT));

END_SECTION


// single IM value for whole spec
const MSSpectrum IMwithDrift = [&]() { 
  MSSpectrum spec; 
  spec.setDriftTime(123.4);
  spec.setDriftTimeUnit(DriftTimeUnit::VSSC);
  return spec;
}();

// convert to IM-Frame with float meta-data array
const MSSpectrum IMwithFDA = [&]() {
  MSExperiment exp;
  exp.addSpectrum(IMwithDrift);
  auto single = IMDataConverter::collapseFramesToSingle(exp);
  return single[0];
}();

START_SECTION(static IMFormat determineIMFormat(const MSExperiment& exp))

  TEST_EQUAL(IMTypes::determineIMFormat(MSExperiment()) == IMFormat::NONE, true)

  {
    MSExperiment exp;
    exp.addSpectrum(MSSpectrum());
    exp.addSpectrum(MSSpectrum());
    TEST_EQUAL(IMTypes::determineIMFormat(exp) == IMFormat::NONE, true)
  }
  
  {
    MSExperiment exp;
    exp.addSpectrum(MSSpectrum());
    exp.addSpectrum(IMwithDrift);
    TEST_EQUAL(IMTypes::determineIMFormat(exp) == IMFormat::MULTIPLE_SPECTRA, true)
  }

  {
    MSExperiment exp;
    exp.addSpectrum(MSSpectrum());
    exp.addSpectrum(IMwithFDA);
    TEST_EQUAL(IMTypes::determineIMFormat(exp) == IMFormat::CONCATENATED, true)
  }

  {
    MSExperiment exp;
    exp.addSpectrum(IMwithDrift);
    exp.addSpectrum(IMwithFDA);
    TEST_EQUAL(IMTypes::determineIMFormat(exp) == IMFormat::MIXED, true)
  }

  {
    // set both ... invalid!
    auto IMwithFDA2 = IMwithFDA;
    IMwithFDA2.setDriftTime(123.4);
    MSExperiment exp;
    exp.addSpectrum(IMwithDrift);
    exp.addSpectrum(IMwithFDA);
    exp.addSpectrum(IMwithFDA2);
    TEST_EXCEPTION(Exception::InvalidValue, IMTypes::determineIMFormat(exp))
  }

END_SECTION

START_SECTION(static IMFormat determineIMFormat(const MSSpectrum& spec))
   TEST_EQUAL(IMTypes::determineIMFormat(MSSpectrum()) == IMFormat::NONE, true)
   
   // single IM value for whole spec
   TEST_EQUAL(IMTypes::determineIMFormat(IMwithDrift) == IMFormat::MULTIPLE_SPECTRA, true)

   // convert to IM-Frame with float meta-data array
   TEST_EQUAL(IMTypes::determineIMFormat(IMwithFDA) == IMFormat::CONCATENATED, true)

   // set both ... invalid!
   auto IMwithFDA2 = IMwithFDA;
   IMwithFDA2.setDriftTime(123.4);
   TEST_EXCEPTION(Exception::InvalidValue, IMTypes::determineIMFormat(IMwithFDA2))

END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
