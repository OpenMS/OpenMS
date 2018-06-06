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

#include <OpenMS/CHEMISTRY/SvmTheoreticalSpectrumGeneratorSet.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/FORMAT/MzMLFile.h>

///////////////////////////

START_TEST(SvmTheoreticalSpectrumGeneratorSet, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

SvmTheoreticalSpectrumGeneratorSet* ptr = nullptr;
SvmTheoreticalSpectrumGeneratorSet* nullPointer = nullptr;

START_SECTION(SvmTheoreticalSpectrumGeneratorSet())
  ptr = new SvmTheoreticalSpectrumGeneratorSet();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(SvmTheoreticalSpectrumGeneratorSet(const SvmTheoreticalSpectrumGeneratorSet& source))
  NOT_TESTABLE //is tested in getSupportedCharges test
END_SECTION

START_SECTION(SvmTheoreticalSpectrumGeneratorSet& operator =(const SvmTheoreticalSpectrumGeneratorSet& tsg))
  NOT_TESTABLE //is tested in getSupportedCharges test
END_SECTION

START_SECTION(~SvmTheoreticalSpectrumGeneratorSet())
  delete ptr;
END_SECTION

SvmTheoreticalSpectrumGeneratorSet gen_set;

START_SECTION(void load(String))
    gen_set.load("examples/simulation/SvmModelSet.model");
    NOT_TESTABLE //is implicitly tested by the following two tests
END_SECTION

START_SECTION(void getSupportedCharges(std::set<Size>&charges))

    std::set<Size>charges;
    gen_set.getSupportedCharges(charges);
    TEST_EQUAL(charges.size(), 3)
    TEST_EQUAL(*(charges.begin()) , 1)
    TEST_EQUAL(*(--charges.end()), 3)

    charges.clear ();
    SvmTheoreticalSpectrumGeneratorSet gen_set_copy(gen_set);
    gen_set_copy.getSupportedCharges(charges);
    TEST_EQUAL(charges.size(), 3)
    TEST_EQUAL(*(charges.begin()) , 1)
    TEST_EQUAL(*(--charges.end()), 3)

    charges.clear ();
    SvmTheoreticalSpectrumGeneratorSet gen_set_assign;
    gen_set_assign = gen_set;
    gen_set_assign.getSupportedCharges(charges);
    TEST_EQUAL(charges.size(), 3)
    TEST_EQUAL(*(charges.begin()) , 1)
    TEST_EQUAL(*(--charges.end()), 3)
END_SECTION

START_SECTION(SvmTheoreticalSpectrumGenerator & getSvmModel(Size))
    NOT_TESTABLE
END_SECTION

START_SECTION(void simulate(PeakSpectrum & spectrum, const AASequence & peptide, boost::random::mt19937_64& rng, Size precursor_charge))

    PeakMap exp;
    boost::random::mt19937_64 rnd_gen (0);
    PeakSpectrum spec;
    AASequence peptide = AASequence::fromString("IFSQVGK");

    Param p = gen_set.getSvmModel(2).getDefaults();
    p.setValue("hide_losses", "true");
    gen_set.getSvmModel(2).setParameters(p);

    gen_set.simulate(spec, peptide, rnd_gen, 2);

#if OPENMS_BOOST_VERSION_MINOR < 56
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("SvmTheoreticalSpectrumGenerator_test.mzML"),exp);
#else
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("SvmTheoreticalSpectrumGenerator_test_boost58.mzML"),exp);
#endif

    if(exp.size())
    {
      TEST_EQUAL(spec.size(), exp[0].size());
      Size min_size = min(spec.size(), exp[0].size());

      for(Size i = 0; i<min_size; ++i)
      {
        TEST_REAL_SIMILAR(spec[i].getPosition()[0],(exp[0][i]).getPosition()[0])
        TEST_EQUAL(spec[i].getIntensity(),(exp[0][i]).getIntensity())
        }
    }
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST

