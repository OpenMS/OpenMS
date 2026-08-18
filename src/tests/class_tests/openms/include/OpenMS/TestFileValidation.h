// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Clemens Groepl, Chris Bielow, Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

// XML schema validation of the temporary files a test creates.
//
// This lives in the class-test project rather than in the test framework
// (OpenMSTestFramework) on purpose: validation needs the FORMAT layer, while the
// framework must stay usable by tests of libraries that do not have it. It is
// header-only, so a test that wants its output validated just includes it and
// calls #VALIDATE_TMP_FILES before END_TEST -- nothing to link, nothing to register.

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Types.h>            // Size
#include <OpenMS/DATASTRUCTURES/TypeAliases.h> // StringList

#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <iostream>
#include <string>
#include <vector>

namespace OpenMS::Internal::ClassTest
{
  /**
    @brief Validates the given files against their XML schema (if one applies).

    Files whose type is not one of the validatable XML formats are skipped.

    @param[in] file_names The files to check (usually TEST::tmp_file_list)
    @return true if all applicable files passed validation

    @ingroup ClassTest
  */
  inline bool validateTmpFiles(const std::vector<std::string>& file_names)
  {
    std::cout << "checking (created temporary files)...\n";
    bool passed_all = true;
    for (Size i = 0; i < file_names.size(); ++i)
    {
      if (!File::exists(file_names[i]))
      {
        continue;
      }
      FileTypes::Type type = FileHandler::getType(file_names[i]);
      bool passed_single = true;
      bool skipped = false;
      switch (type)
      {
      case FileTypes::MZML:
      {
        if (!MzMLFile().isValid(file_names[i]))
        {
          std::cout << " - Error: mzML file does not validate against XML schema '" << file_names[i] << "'\n";
          passed_single = false;
        }
        StringList errors, warnings;
        if (!MzMLFile().isSemanticallyValid(file_names[i], errors, warnings))
        {
          std::cout << " - Error: mzML file semantically invalid '" << file_names[i] << "'\n";
          for (Size j = 0; j < errors.size(); ++j)
          {
            std::cout << "Error - " << errors[j] << '\n';
          }
          passed_single = false;
        }
      }
      break;

      case FileTypes::MZDATA:
        if (!MzDataFile().isValid(file_names[i], std::cerr))
        {
          std::cout << " - Error: Invalid mzData file '" << file_names[i] << "'\n";
          passed_single = false;
        }
        break;

      case FileTypes::MZXML:
        if (!MzXMLFile().isValid(file_names[i], std::cerr))
        {
          std::cout << " - Error: Invalid mzXML file '" << file_names[i] << "'\n";
          passed_single = false;
        }
        break;

      case FileTypes::FEATUREXML:
        if (!FeatureXMLFile().isValid(file_names[i], std::cerr))
        {
          std::cout << " - Error: Invalid featureXML file '" << file_names[i] << "'\n";
          passed_single = false;
        }
        break;

      case FileTypes::IDXML:
        if (!IdXMLFile().isValid(file_names[i], std::cerr))
        {
          std::cout << " - Error: Invalid idXML file '" << file_names[i] << "'\n";
          passed_single = false;
        }
        break;

      case FileTypes::CONSENSUSXML:
        if (!ConsensusXMLFile().isValid(file_names[i], std::cerr))
        {
          std::cout << " - Error: Invalid consensusXML file '" << file_names[i] << "'\n";
          passed_single = false;
        }
        break;

      case FileTypes::INI:
        if (!ParamXMLFile().isValid(file_names[i], std::cerr))
        {
          std::cout << " - Error: Invalid Param file '" << file_names[i] << "'\n";
          passed_single = false;
        }
        break;

      case FileTypes::TRANSFORMATIONXML:
        if (!TransformationXMLFile().isValid(file_names[i], std::cerr))
        {
          passed_single = false;
        }
        break;

      default:
        skipped = true;
        break;
      }

      if (skipped)
      {
        std::cout << " +  skipped file '" << file_names[i] << "' (type: " << FileTypes::typeToName(type) << ")\n";
      }
      else if (passed_single)
      {
        std::cout << " +  valid file '" << file_names[i] << "' (type: " << FileTypes::typeToName(type) << ")\n";
      }
      else
      {
        passed_all = false;
        std::cout << " -  invalid file '" << file_names[i] << "' (type: " << FileTypes::typeToName(type) << ")\n";
      }
    }

    std::cout << (passed_all ? ": passed" : ": failed") << std::endl << '\n';
    return passed_all;
  }
} // namespace OpenMS::Internal::ClassTest

/**
  @brief Validate every file created with #NEW_TMP_FILE against its XML schema.

  Place this right before #END_TEST in tests that write mzML, mzXML, mzData,
  featureXML, idXML, consensusXML, transformationXML or INI files. Files of any
  other type are skipped, so it is always safe to call.

  Requires <OpenMS/TestFileValidation.h>. The test fails if any file is invalid.

  @ingroup ClassTest
  @hideinitializer
*/
#define VALIDATE_TMP_FILES                                                                \
  {                                                                                       \
    ++TEST::test_count;                                                                   \
    TEST::test_line = __LINE__;                                                           \
    TEST::this_test = OpenMS::Internal::ClassTest::validateTmpFiles(TEST::tmp_file_list);  \
    TEST::test = TEST::test && TEST::this_test;                                           \
    TEST::all_tests = TEST::all_tests && TEST::this_test;                                 \
    if (!TEST::this_test)                                                                 \
    {                                                                                     \
      stdcout << " -  line " << TEST::test_line << ": VALIDATE_TMP_FILES\n";              \
      TEST::failed_lines_list.push_back(TEST::test_line);                                 \
    }                                                                                     \
  }

/**
  @brief Validate a single file against its XML schema, at the point of the call.

  Same checks as #VALIDATE_TMP_FILES, but for one file and reported at this line,
  which makes a validation failure easier to attribute than a check at the end of
  the test. Also usable for files that were not created via #NEW_TMP_FILE.

  @ingroup ClassTest
  @hideinitializer
*/
#define VALIDATE_FILE(filename)                                                           \
  {                                                                                       \
    ++TEST::test_count;                                                                   \
    TEST::test_line = __LINE__;                                                           \
    TEST::this_test = OpenMS::Internal::ClassTest::validateTmpFiles(                      \
                        std::vector<std::string> {(filename)});                           \
    TEST::test = TEST::test && TEST::this_test;                                           \
    TEST::all_tests = TEST::all_tests && TEST::this_test;                                 \
    if (!TEST::this_test)                                                                 \
    {                                                                                     \
      stdcout << " -  line " << TEST::test_line << ": VALIDATE_FILE(" # filename ")\n";   \
      TEST::failed_lines_list.push_back(TEST::test_line);                                 \
    }                                                                                     \
  }
