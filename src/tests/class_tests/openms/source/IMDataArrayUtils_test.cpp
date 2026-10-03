// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/IONMOBILITY/IMDataArrayUtils.h>
#include <OpenMS/IONMOBILITY/IMDataConverter.h>
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>

#include <set>
#include <string>

using namespace OpenMS;

namespace
{
  // unit of a PSI-MS term, chosen from its units as the vocabulary lookup did before the table
  DriftTimeUnit unitFromOntology(const ControlledVocabulary::CVTerm& term)
  {
    if (term.units.contains("MS:1002814")) return DriftTimeUnit::VSSC;        // volt-second per square centimeter
    if (term.units.contains("UO:0000028")) return DriftTimeUnit::MILLISECOND; // millisecond
    if (term.units.contains("UO:0000324")) return DriftTimeUnit::CCS;         // square angstrom
    return DriftTimeUnit::NONE;
  }

  std::string joined(const std::set<std::string>& values)
  {
    std::string result;
    for (const std::string& value : values)
    {
      result += (result.empty() ? "" : ", ") + value;
    }
    return result;
  }
}

START_TEST(IMDataArrayUtils, "$Id$")

START_SECTION((unit round trips and compatibility entry points))
  for (const auto expected : {DriftTimeUnit::MILLISECOND, DriftTimeUnit::VSSC})
  {
    DataArrays::FloatDataArray direct, legacy;
    direct.push_back(1.25f);
    IMDataArrayUtils::setIMUnit(direct, expected);
    IMDataConverter::setIMUnit(legacy, expected);
    TEST_EQUAL(direct.getName(), legacy.getName())
    TEST_REAL_SIMILAR(direct[0], 1.25f)
    DriftTimeUnit actual = DriftTimeUnit::NONE;
    TEST_TRUE(IMDataArrayUtils::getIMUnit(direct, actual))
    TEST_TRUE(actual == expected)
    actual = DriftTimeUnit::NONE;
    TEST_TRUE(IMDataConverter::getIMUnit(direct, actual))
    TEST_TRUE(actual == expected)
  }
  DataArrays::FloatDataArray array;
  IMDataArrayUtils::setIMUnit(array, DriftTimeUnit::MILLISECOND);
  TEST_STRING_EQUAL(array.getName(), "mean ion mobility array")
  IMDataArrayUtils::setIMUnit(array, DriftTimeUnit::VSSC);
  TEST_STRING_EQUAL(array.getName(), "raw inverse reduced ion mobility array")
  TEST_EXCEPTION(Exception::InvalidValue, IMDataArrayUtils::setIMUnit(array, DriftTimeUnit::NONE))
  TEST_EXCEPTION(Exception::InvalidValue, IMDataArrayUtils::setIMUnit(array, DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE))
  TEST_EXCEPTION(Exception::InvalidValue, IMDataArrayUtils::setIMUnit(array, DriftTimeUnit::CCS))
  TEST_STRING_EQUAL(array.getName(), "raw inverse reduced ion mobility array")
END_SECTION

START_SECTION((ontology and vendor names retain their units))
  DataArrays::FloatDataArray array;
  DriftTimeUnit unit = DriftTimeUnit::NONE;
  array.setName("mean inverse reduced ion mobility array");
  TEST_TRUE(IMDataArrayUtils::getIMUnit(array, unit))
  TEST_TRUE(unit == DriftTimeUnit::VSSC)
  array.setName(Constants::UserParam::ION_MOBILITY);
  TEST_TRUE(IMDataArrayUtils::getIMUnit(array, unit))
  TEST_TRUE(unit == DriftTimeUnit::MILLISECOND)
  array.setName(Constants::UserParam::INVERSE_REDUCED_ION_MOBILITY);
  TEST_TRUE(IMDataArrayUtils::getIMUnit(array, unit))
  TEST_TRUE(unit == DriftTimeUnit::VSSC)
  array.setName(std::string(Constants::UserParam::ION_MOBILITY) + " MS:1002954");
  TEST_TRUE(IMDataArrayUtils::getIMUnit(array, unit))
  TEST_TRUE(unit == DriftTimeUnit::CCS)
  array.setName("unrelated auxiliary data");
  TEST_FALSE(IMDataArrayUtils::getIMUnit(array, unit))
  TEST_TRUE(unit == DriftTimeUnit::CCS)
END_SECTION

START_SECTION((every PSI-MS ion mobility array name gets the unit of its term))
  for (const auto& term : IMDataArrayUtils::getArrayTerms())
  {
    DataArrays::FloatDataArray array;
    array.setName(std::string(term.name));
    DriftTimeUnit unit = DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE; // never a result
    TEST_TRUE(IMDataArrayUtils::getIMUnit(array, unit))
    TEST_EQUAL(driftTimeUnitToString(unit), driftTimeUnitToString(term.unit))
  }
END_SECTION

START_SECTION((CV names match exactly, before the vendor fallbacks))
  DataArrays::FloatDataArray array;
  DriftTimeUnit unit = DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE;
  // neither CV names (case and whitespace matter) nor vendor names: not recognized, unit unchanged
  for (const std::string name : {"MEAN ION MOBILITY ARRAY", "mean ion mobility array ", "ion mobility array", "ion mobility", ""})
  {
    array.setName(name);
    TEST_FALSE(IMDataArrayUtils::getIMUnit(array, unit))
    TEST_TRUE(unit == DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE)
  }
  // a CV name followed by text is no CV name, but the vendor prefix still applies
  array.setName("mean inverse reduced ion mobility array (1/K0)");
  TEST_TRUE(IMDataArrayUtils::getIMUnit(array, unit))
  TEST_TRUE(unit == DriftTimeUnit::VSSC)
  array.setName("inverse reduced ion mobility array");
  TEST_TRUE(IMDataArrayUtils::getIMUnit(array, unit))
  TEST_TRUE(unit == DriftTimeUnit::VSSC)
  // 'Ion Mobility' names: VSSC for MS:1002815 or MS:1003006, then CCS for MS:1002954, else milliseconds
  array.setName(std::string(Constants::UserParam::ION_MOBILITY) + " MS:1002815");
  TEST_TRUE(IMDataArrayUtils::getIMUnit(array, unit))
  TEST_TRUE(unit == DriftTimeUnit::VSSC)
  array.setName(std::string(Constants::UserParam::ION_MOBILITY) + " MS:1003006 MS:1002954");
  TEST_TRUE(IMDataArrayUtils::getIMUnit(array, unit))
  TEST_TRUE(unit == DriftTimeUnit::VSSC)
  array.setName(std::string(Constants::UserParam::ION_MOBILITY) + " array");
  TEST_TRUE(IMDataArrayUtils::getIMUnit(array, unit))
  TEST_TRUE(unit == DriftTimeUnit::MILLISECOND)
END_SECTION

// The only use of the CV files in this test. It runs last, so the sections above run without them loaded.
START_SECTION((static std::span<const ArrayTerm> getArrayTerms()))
  // The compiled-in table must hold exactly the child terms of 'MS:1002893 ! ion mobility array' in
  // the bundled PSI-MS ontology, with their names and the units a lookup in the ontology yields.
  // If this fails after an update of psi-ms.obo, update the table in IMDataArrayUtils.cpp.
  const ControlledVocabulary& cv = ControlledVocabulary::getPSIMSCV();
  std::set<std::string> ontology_terms;
  cv.getAllChildTerms(ontology_terms, "MS:1002893");
  std::set<std::string> table_terms;
  for (const auto& term : IMDataArrayUtils::getArrayTerms())
  {
    table_terms.insert(std::string(term.accession));
  }
  TEST_EQUAL(table_terms.size(), IMDataArrayUtils::getArrayTerms().size()) // no accession twice
  TEST_EQUAL(joined(table_terms), joined(ontology_terms))

  for (const auto& term : IMDataArrayUtils::getArrayTerms())
  {
    const std::string accession(term.accession);
    const std::string name(term.name);
    if (!cv.exists(accession))
    {
      continue; // reported by the comparison of the accessions above
    }
    const ControlledVocabulary::CVTerm& cv_term = cv.getTerm(accession);
    TEST_STRING_EQUAL(name, cv_term.name)
    TEST_EQUAL(driftTimeUnitToString(term.unit), driftTimeUnitToString(unitFromOntology(cv_term)))
    // the name lookup of the vocabulary resolves the name to this term (no other term has the name)
    const ControlledVocabulary::CVTerm* by_name = cv.checkAndGetTermByName(name);
    TEST_EQUAL(by_name == nullptr ? std::string("no term") : by_name->id, accession)
  }
END_SECTION

END_TEST
