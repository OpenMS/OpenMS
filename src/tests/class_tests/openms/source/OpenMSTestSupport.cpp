// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

// Registers libOpenMS-specific behavior with the standard-library-only class-test
// framework. This file is compiled into every openms / openms_gui class-test
// executable (see the add_executable() calls in the CMakeLists.txt of the test
// projects); tests of libraries that do not link libOpenMS (e.g. OpenSwathAlgo)
// simply do not include it. It is NOT a test itself.
//
// Two registrations happen at static initialization (i.e. before main(), and
// therefore before any test code runs):
//
//  1. UniqueIdGenerator is seeded with a fixed value so that unique IDs written
//     into test output files are reproducible against the reference files.
//     (Previously done inside the framework at START_TEST; static init is even
//     earlier.) Do NOT remove this file from a test target: ID-bearing
//     reference-file comparisons (e.g. featureXML) would become
//     nondeterministic.
//
//  2. An exception translator (see ClassTest.h) so that unexpected OpenMS
//     exceptions in failure reports show their name and origin. Without it,
//     reports fall back to std::exception::what() -- degraded output, never
//     wrong results.

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>

#include <cstring>
#include <sstream>
#include <string>

namespace
{
  bool describeOpenMSException(std::string& description)
  {
    try
    {
      throw; // rethrow the exception currently being handled
    }
    catch (const OpenMS::Exception::BaseException& e)
    {
      std::ostringstream os;
      os << "OpenMS exception of type '" << e.getName() << "'";
      if ((e.getLine() > 0) && std::strcmp(e.getFile(), "") != 0)
      {
        os << " thrown in line " << e.getLine() << " of file '" << e.getFile()
           << "' in function '" << e.getFunction() << "'";
      }
      os << " - Message: " << e.what();
      description = os.str();
      return true;
    }
    catch (...)
    {
      return false; // not ours; let the framework (or another translator) describe it
    }
  }

  struct OpenMSTestSupport
  {
    OpenMSTestSupport()
    {
      // Fixed seed so unique ids stored in test output are reproducible.
      OpenMS::UniqueIdGenerator::setSeed(2453440375);
      OpenMS::Internal::ClassTest::registerExceptionTranslator(&describeOpenMSException);
    }
  };

  const OpenMSTestSupport openms_test_support_init;
}
