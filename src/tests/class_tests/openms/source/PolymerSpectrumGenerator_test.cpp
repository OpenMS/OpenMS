// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: OpenMS Team $
// $Authors: OpenMS Team $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/PolymerSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/NucleicAcidSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <iostream>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PolymerSpectrumGenerator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PolymerSpectrumGenerator<TheoreticalSpectrumGenerator>* ptr_tsg = nullptr;
PolymerSpectrumGenerator<TheoreticalSpectrumGenerator>* nullPointer_tsg = nullptr;

PolymerSpectrumGenerator<NucleicAcidSpectrumGenerator>* ptr_nasg = nullptr;
PolymerSpectrumGenerator<NucleicAcidSpectrumGenerator>* nullPointer_nasg = nullptr;

START_SECTION(PolymerSpectrumGenerator(const String& name))
{
  ptr_tsg = new TheoreticalSpectrumGenerator();
  ptr_nasg = new NucleicAcidSpectrumGenerator();
  TEST_NOT_EQUAL(ptr_tsg, nullptr);
  TEST_NOT_EQUAL(ptr_nasg, nullptr);
}
END_SECTION

START_SECTION(Derived& derived())
{
  TheoreticalSpectrumGenerator& tsg_ref = ptr_tsg->derived();
  TEST_EQUAL(typeid(tsg_ref).name(), typeid(TheoreticalSpectrumGenerator).name());
  
  NucleicAcidSpectrumGenerator& nasg_ref = ptr_nasg->derived();
  TEST_EQUAL(typeid(nasg_ref).name(), typeid(NucleicAcidSpectrumGenerator).name());
}
END_SECTION

START_SECTION(const Derived& derived() const)
{
  const TheoreticalSpectrumGenerator& tsg_ref = ptr_tsg->derived();
  TEST_EQUAL(typeid(tsg_ref).name(), typeid(TheoreticalSpectrumGenerator).name());
  
  const NucleicAcidSpectrumGenerator& nasg_ref = ptr_nasg->derived();
  TEST_EQUAL(typeid(nasg_ref).name(), typeid(NucleicAcidSpectrumGenerator).name());
}
END_SECTION

START_SECTION(String getSequenceType() const)
{
  TEST_EQUAL(ptr_tsg->getSequenceType(), "AASequence");
  TEST_EQUAL(ptr_nasg->getSequenceType(), "NASequence");
}
END_SECTION

START_SECTION(template<typename SeqType> bool supportsSequenceType() const)
{
  // Test TheoreticalSpectrumGenerator supports AASequence
  TEST_EQUAL(ptr_tsg->supportsSequenceType<AASequence>(), true);
  TEST_EQUAL(ptr_tsg->supportsSequenceType<NASequence>(), false);
  
  // Test NucleicAcidSpectrumGenerator supports NASequence  
  TEST_EQUAL(ptr_nasg->supportsSequenceType<AASequence>(), false);
  TEST_EQUAL(ptr_nasg->supportsSequenceType<NASequence>(), true);
}
END_SECTION

START_SECTION(generateSpectrum integration test)
{
  // Test that the generateSpectrum method works through CRTP interface
  AASequence peptide = AASequence::fromString("PEPTIDE");
  MSSpectrum spectrum;
  
  ptr_tsg->generateSpectrum(spectrum, peptide, 1, 2);
  
  // Should have generated some peaks
  TEST_EQUAL(spectrum.size() > 0, true);
  TEST_EQUAL(spectrum.getMSLevel(), 2);
  
  // Test with NASequence
  NASequence oligo = NASequence::fromString("AUCG");
  MSSpectrum na_spectrum;
  
  ptr_nasg->generateSpectrum(na_spectrum, oligo, -1, -2);
  
  // Should have generated some peaks
  TEST_EQUAL(na_spectrum.size() > 0, true);
}
END_SECTION

START_SECTION(CRTP inheritance verification)
{
  // Verify that both generators properly inherit from PolymerSpectrumGenerator
  TheoreticalSpectrumGenerator tsg;
  NucleicAcidSpectrumGenerator nasg;
  
  // Should be able to use base class interface
  PolymerSpectrumGenerator<TheoreticalSpectrumGenerator>& tsg_base = tsg;
  PolymerSpectrumGenerator<NucleicAcidSpectrumGenerator>& nasg_base = nasg;
  
  TEST_EQUAL(tsg_base.getSequenceType(), "AASequence");
  TEST_EQUAL(nasg_base.getSequenceType(), "NASequence");
  
  // Test derived() returns correct type
  TEST_EQUAL(&tsg_base.derived(), &tsg);
  TEST_EQUAL(&nasg_base.derived(), &nasg);
}
END_SECTION

START_SECTION(parameter handling through base class)
{
  // Test that parameter handling works through base class
  TheoreticalSpectrumGenerator tsg;
  PolymerSpectrumGenerator<TheoreticalSpectrumGenerator>& base = tsg;
  
  // Should be able to access DefaultParamHandler interface
  Param params = base.getParameters();
  TEST_EQUAL(params.empty(), false); // Should have some parameters
  
  // Test setting parameters
  params.setValue("add_b_ions", "true");
  base.setParameters(params);
  
  Param retrieved = base.getParameters();
  TEST_EQUAL(retrieved.getValue("add_b_ions").toString(), "true");
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

delete ptr_tsg;
delete ptr_nasg;

END_TEST