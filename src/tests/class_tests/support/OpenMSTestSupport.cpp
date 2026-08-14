// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Clemens Groepl, Chris Bielow, Timo Sachsenberg $
// --------------------------------------------------------------------------

// This translation unit wires libOpenMS-specific behavior into the otherwise
// library-agnostic class-test framework (OpenMSTestFramework):
//  - #START_TEST seeds the UniqueIdGenerator, so that unique ids are deterministic
//  - #END_TEST validates the XML files created via #NEW_TMP_FILE against their schema
//
// It is linked (as an OBJECT library, see CMakeLists.txt) into the class tests of
// libOpenMS and libOpenMS_GUI; the registration below runs from a static initializer,
// i.e. before main() of the test program. Tests of libraries that do not depend on
// libOpenMS (e.g. OpenSwathAlgo) do not link this library and get no validation,
// which is fine since they cannot produce OpenMS XML files.

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
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

using namespace OpenMS;

namespace
{
  /// Validates the given files against the XML schema (if available); returns true if all files passed
  bool validateTmpFiles(const std::vector<std::string>& file_names)
  {
    std::cout << "checking (created temporary files)...\n";
    bool passed_all = true;
    for (Size i = 0; i < file_names.size(); ++i)
    {
      if (File::exists(file_names[i]))
      {
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
          if (!MzMLFile().isSemanticallyValid(file_names[i], errors,
                                              warnings))
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
        //output for single file
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
    }
    //output for all files
    if (passed_all)
    {
      std::cout << ": passed" << std::endl << '\n';
    }
    else
    {
      std::cout << ": failed" << std::endl << '\n';
    }
    return passed_all;
  }

  /// Runs at the start of every test program (via START_TEST), before any test code
  void initTest()
  {
    // fixed seed, so that tests which store unique ids produce reproducible output
    UniqueIdGenerator::setSeed(2453440375);
  }

  // registration happens at static-initialization time, i.e. before main() of the test
  struct TestSupportRegistrar
  {
    TestSupportRegistrar()
    {
      Internal::ClassTest::setTestInitHook(&initTest);
      Internal::ClassTest::setTmpFileValidator(&validateTmpFiles);
    }
  };
  [[maybe_unused]] const TestSupportRegistrar registrar_instance;
}
