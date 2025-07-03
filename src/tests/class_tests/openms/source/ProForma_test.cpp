// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Ayesha Feroz $
// $Authors: Ayesha Feroz, Samuel Wein $
// --------------------------------------------------------------------------
//
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/CHEMISTRY/ProForma.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

using namespace OpenMS;
using namespace std;

START_TEST(ProForma, "$Id$")

ProForma* ptr = nullptr;

START_SECTION(ProForma(const AASequence& seq))
  ptr = new ProForma(AASequence::fromString("ACDEFGHIK"));
  TEST_NOT_EQUAL(ptr, nullptr);
END_SECTION

START_SECTION(~ProForma())
  delete ptr;
END_SECTION


String proforma_str = "ACDEFGHIK";
AASequence seq = AASequence::fromString(proforma_str);
ProForma proforma(seq);

START_SECTION(AASequence FromProFormaString(const string& proforma_str))
  proforma_str = "A[Phospho]CDEFGHIK";
  proforma.fromProFormaString(proforma_str);
  TEST_EQUAL(proforma.toProFormaString(), proforma_str);

  //parseMassShiftModification
  proforma_str = "A[+15.99]CDEFGHIK";
  proforma.fromProFormaString(proforma_str);
  TEST_EQUAL(proforma.toProFormaString(), proforma_str);

  //parseNTerminalModification
  proforma_str = "[Acetyl]-ACDEFGHIK";
  proforma.fromProFormaString(proforma_str);
  TEST_EQUAL(proforma.toProFormaString(), proforma_str);

  //parseCTerminalModification
  proforma_str = "ACDEFGHIK-[Amidation]";
  proforma.fromProFormaString(proforma_str);
  TEST_EQUAL(proforma.toProFormaString(), proforma_str);

  //parseAmbiguousModification
  proforma_str = "A[+15.99]?CDEFGHIK";
  proforma.fromProFormaString(proforma_str);
  TEST_EQUAL(proforma.toProFormaString(), proforma_str);

  //invalidModificationFormat
  proforma_str = "A[PhosphoCDEFGHIK"; // Missing closing bracket
  TEST_EXCEPTION(Exception::ParseError, proforma.fromProFormaString(proforma_str));

  //test_unsupportedCV
  proforma_str = "A[XYZ:12345]CDEFGHIK"; // Unsupported CV
  TEST_EXCEPTION(Exception::IllegalArgument, proforma.fromProFormaString(proforma_str));

  //validProFormaPSIMOD
  proforma_str = "EM[MOD:00719]EVEES[MOD:00046]PEK";
  proforma.fromProFormaString(proforma_str);
  TEST_EQUAL(proforma.toProFormaString(), proforma_str);

  //validProFormaRESID
  proforma_str = "EM[RESID:AA0581]EVEES[RESID:AA0037]PEK";
  proforma.fromProFormaString(proforma_str);
  TEST_EQUAL(proforma.toProFormaString(), proforma_str);

  //validProFormaUNIMOD
  proforma_str = "EM[UNIMOD:35]EVEES[UNIMOD:56]PEK";
  proforma.fromProFormaString(proforma_str);
  TEST_EQUAL(proforma.toProFormaString(), proforma_str);

  //invalidProFormaPSIMODSyntax
  proforma_str = "EM[M:00719]EVEES[M:00046]PEK";
  TEST_EXCEPTION(Exception::IllegalArgument, proforma.fromProFormaString(proforma_str));

  //invalidProFormaUNIMODSyntax
  proforma_str = "EM[U:35]EVEES[U:56]PEK";
  TEST_EXCEPTION(Exception::IllegalArgument, proforma.fromProFormaString(proforma_str));

  //invalidProFormaRESIDSyntax
  TEST_EXCEPTION(Exception::IllegalArgument, proforma.fromProFormaString("EM[R:AA0581]EVEES[R:AA0037]PEK"));

  //parseRangeModification
  proforma_str = "A(CDE)[+12.5]FGH[+15.99]IK";
  proforma.fromProFormaString(proforma_str);
  TEST_EQUAL(proforma.toProFormaString(), proforma_str);

END_SECTION

//removeModification
START_SECTION(void removeModification(size_t position))  
  proforma_str = "A[Phospho]CDEFGHIK";
  proforma.fromProFormaString(proforma_str);
  proforma.removeModification(1);
  TEST_EQUAL(proforma.toProFormaString(), "ACDEFGHIK");
  TEST_EXCEPTION(Exception::InvalidValue, proforma.removeModification(1));
  TEST_EXCEPTION(Exception::OutOfRange, proforma.removeModification(9999));
  
END_SECTION

END_TEST