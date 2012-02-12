// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Sandro Andreotti $
// $Authors: Sandro Andreotti $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <iostream>

#include <OpenMS/CHEMISTRY/SvmTheoreticalSpectrumGeneratorSet.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/FORMAT/MzMLFile.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

///////////////////////////

START_TEST(SvmTheoreticalSpectrumGeneratorSet, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

SvmTheoreticalSpectrumGeneratorSet* ptr = 0;
SvmTheoreticalSpectrumGeneratorSet* nullPointer = 0;

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

START_SECTION(void simulate(RichPeakSpectrum &spectrum, const AASequence &peptide, const gsl_rng *rng, Size precursor_charge))

    RichPeakMap exp;
    gsl_rng* rnd_gen = gsl_rng_alloc (gsl_rng_taus);
    gsl_rng_set(rnd_gen, 0);
    RichPeakSpectrum spec;
    AASequence peptide("IFSQVGK");

    Param p = gen_set.getSvmModel(2).getDefaults();
    p.setValue("hide_losses", "true");
    gen_set.getSvmModel(2).setParameters(p);

    gen_set.simulate(spec, peptide,rnd_gen,2);
    gsl_rng_free(rnd_gen);

    MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("SvmTheoreticalSpectrumGenerator_test.mzML"),exp);
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

