// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/CHEMISTRY/SvmTheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/FORMAT/MzMLFile.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

///////////////////////////

START_TEST(SvmTheoreticalSpectrumGenerator, "$Id: SvmTheoreticalSpectrumGenerator_test.C 7353 2010-09-03 21:18:31Z aiche $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

SvmTheoreticalSpectrumGenerator* ptr = 0;
SvmTheoreticalSpectrumGenerator* nullPointer = 0;

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
AASequence peptide("IFSQVGK");

START_SECTION(SvmTheoreticalSpectrumGenerator& operator = (const SvmTheoreticalSpectrumGenerator& tsg))
  SvmTheoreticalSpectrumGenerator copy;
  copy = *ptr;
  TEST_EQUAL(copy.getParameters(), ptr->getParameters())
END_SECTION


START_SECTION(void simulate(RichPeakSpectrum &spectrum, const AASequence &peptide, const gsl_rng *rng, Size precursor_charge))
  // init rng
  gsl_rng* rnd_gen = gsl_rng_alloc (gsl_rng_taus);
  gsl_rng_set(rnd_gen, 0);
  RichPeakSpectrum spec;

  Param p = ptr->getDefaults();
  p.setValue ("hide_losses", "true");
  ptr->setParameters (p);

  ptr->load();
  ptr->simulate(spec, peptide,rnd_gen,1);
  gsl_rng_free(rnd_gen);

  MSExperiment<RichPeak1D>exp;
  //MSExperiment<RichPeak1D>exp2;
  //exp2.push_back(spec);
  MzMLFile mz_file;
  //MzMLFile().store(OPENMS_GET_TEST_DATA_PATH("SvmTheoreticalSpectrumGenerator_test.mzML"),exp2);

  mz_file.load(OPENMS_GET_TEST_DATA_PATH("SvmTheoreticalSpectrumGenerator_test.mzML"),exp);

  TEST_EQUAL(exp.size(), 1);
  if(exp.size())
  {
    TEST_EQUAL(spec.size(), exp[0].size());
    Size min_size = min(spec.size(), exp[0].size());

    for(Size i = 0; i<min_size; ++i)
    {
      TEST_REAL_SIMILAR(spec[i].getPosition()[0],(exp[0][i]).getPosition()[0]);
      TEST_EQUAL(spec[i].getIntensity(),(exp[0][i]).getIntensity());
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

SvmTheoreticalSpectrumGenerator::IonType* ptr_t = 0;
SvmTheoreticalSpectrumGenerator::IonType* nullPointer_t = 0;
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


SvmTheoreticalSpectrumGenerator::SvmModel* ptr_m = 0;
SvmTheoreticalSpectrumGenerator::SvmModel* nullPointer_m = 0;
START_SECTION([SvmTheoreticalSpectrumGenerator::SvmModel] SvmModel())
    ptr_m = new SvmTheoreticalSpectrumGenerator::SvmModel;
    TEST_NOT_EQUAL(ptr_m, nullPointer_m)
END_SECTION

START_SECTION([SvmTheoreticalSpectrumGenerator::SvmModel] ~SvmModel())
    delete ptr_m;
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST

