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
// $Maintainer: Timo Sachsenberg $
// $Authors: Sandro Andreotti $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <iostream>

#include <OpenMS/CHEMISTRY/SvmTheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/FORMAT/MzMLFile.h>

///////////////////////////

START_TEST(SvmTheoreticalSpectrumGenerator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

SvmTheoreticalSpectrumGenerator* ptr = nullptr;
SvmTheoreticalSpectrumGenerator* nullPointer = nullptr;

START_SECTION(SvmTheoreticalSpectrumGenerator())
  ptr = new SvmTheoreticalSpectrumGenerator();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(SvmTheoreticalSpectrumGenerator(const SvmTheoreticalSpectrumGenerator& source))
  SvmTheoreticalSpectrumGenerator copy(*ptr);
  TEST_EQUAL(copy.getParameters(), ptr->getParameters())
END_SECTION

START_SECTION(~SvmTheoreticalSpectrumGenerator())
  delete ptr;
END_SECTION

ptr = new SvmTheoreticalSpectrumGenerator();
AASequence peptide = AASequence::fromString("IFSQVGK");

START_SECTION(SvmTheoreticalSpectrumGenerator& operator = (const SvmTheoreticalSpectrumGenerator& tsg))
  SvmTheoreticalSpectrumGenerator copy;
  copy = *ptr;
  TEST_EQUAL(copy.getParameters(), ptr->getParameters())
END_SECTION


START_SECTION(void simulate(PeakSpectrum &spectrum, const AASequence &peptide, boost::random::mt19937_64&rng, Size precursor_charge))
  // init rng
  boost::random::mt19937_64 rnd_gen (0);
  PeakSpectrum spec;

  Param p = ptr->getDefaults();
  p.setValue ("hide_losses", "true");
  p.setValue ("add_metainfo", "true");
  ptr->setParameters (p);

  ptr->load();
  ptr->simulate(spec, peptide, rnd_gen, 1);

  PeakMap exp;
  MzMLFile mz_file;

#if OPENMS_BOOST_VERSION_MINOR < 56
  mz_file.load(OPENMS_GET_TEST_DATA_PATH("SvmTheoreticalSpectrumGenerator_test.mzML"),exp);
  TEST_EQUAL(spec.size(), 7);
#else
  mz_file.load(OPENMS_GET_TEST_DATA_PATH("SvmTheoreticalSpectrumGenerator_test_boost58.mzML"),exp);
  TEST_EQUAL(spec.size(), 8);
  // the extra peak:
  TEST_EQUAL(spec.getStringDataArrays()[0][2], "YIon  0++") // TODO: ion_nr is always zero, its actually y4++
  TEST_EQUAL(spec.getIntegerDataArrays()[0][2], 2)
#endif

  TEST_EQUAL(exp.size(), 1);
  if(exp.size())
  {
    TEST_EQUAL(spec.size(), exp[0].size());
    Size min_size = min(spec.size(), exp[0].size());

    for(Size i = 0; i<min_size; ++i)
    {
      TEST_REAL_SIMILAR(spec[i].getPosition()[0],(exp[0][i]).getPosition()[0]);
      TEST_REAL_SIMILAR(spec[i].getIntensity(),(exp[0][i]).getIntensity());
      }
  }
END_SECTION
delete ptr;

START_SECTION(void load())
//This method is already used(and therefore tested) in the simulation test
NOT_TESTABLE
END_SECTION

START_SECTION(const std::vector<IonType>& getIonTypes())
//This method is already used(and therefore tested) in the simulation test
NOT_TESTABLE
END_SECTION

SvmTheoreticalSpectrumGenerator::IonType* ptr_t = nullptr;
SvmTheoreticalSpectrumGenerator::IonType* nullPointer_t = nullptr;
START_SECTION([SvmTheoreticalSpectrumGenerator::IonType] IonType())
  ptr_t = new SvmTheoreticalSpectrumGenerator::IonType();
  TEST_NOT_EQUAL(ptr_t, nullPointer_t)
  delete ptr_t;
END_SECTION

START_SECTION([SvmTheoreticalSpectrumGenerator::IonType] IonType(Residue::ResidueType residue, EmpiricalFormula loss=EmpiricalFormula(), Int charge=1))
  SvmTheoreticalSpectrumGenerator::IonType type(Residue::BIon, EmpiricalFormula(""), 2);
  TEST_EQUAL(type.residue, Residue::BIon)
  TEST_EQUAL(type.loss, EmpiricalFormula(""))
  TEST_EQUAL(type.charge, 2)
END_SECTION

START_SECTION([SvmTheoreticalSpectrumGenerator::IonType] IonType(const IonType &rhs))
  SvmTheoreticalSpectrumGenerator::IonType type(Residue::BIon, EmpiricalFormula(""), 2);
  SvmTheoreticalSpectrumGenerator::IonType copy(type);
  TEST_EQUAL(type.residue, copy.residue)
  TEST_EQUAL(type.charge, copy.charge)
  TEST_EQUAL(type.loss, copy.loss)
END_SECTION

START_SECTION([SvmTheoreticalSpectrumGenerator::IonType] IonType& operator=(const IonType &rhs))
  SvmTheoreticalSpectrumGenerator::IonType type(Residue::BIon, EmpiricalFormula(""), 2);
  SvmTheoreticalSpectrumGenerator::IonType copy;
  copy=type;
  TEST_EQUAL(type.residue, copy.residue)
  TEST_EQUAL(type.charge, copy.charge)
  TEST_EQUAL(type.loss, copy.loss)
END_SECTION

START_SECTION([SvmTheoreticalSpectrumGenerator::IonType] bool operator<(const IonType &rhs) const)
  SvmTheoreticalSpectrumGenerator::IonType type(Residue::BIon, EmpiricalFormula(""), 2);
  SvmTheoreticalSpectrumGenerator::IonType type2(Residue::YIon, EmpiricalFormula(""), 2);
  TEST_EQUAL(type<type2, true)
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST

