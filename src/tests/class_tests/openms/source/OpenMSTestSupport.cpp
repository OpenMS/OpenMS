// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

// libOpenMS-specific setup for the std-only class-test framework. Compiled into
// every openms/openms_gui test executable (not tests that don't link libOpenMS,
// e.g. OpenSwathAlgo). Registers, at static init:
//  1. a fixed UniqueIdGenerator seed, so IDs in output files match the reference
//     files -- do NOT drop this file from a target or ID comparisons (e.g.
//     featureXML) become nondeterministic;
//  2. an OpenMS exception translator, so failure reports name unexpected OpenMS
//     exceptions (without it: std::exception::what() -- degraded, not wrong).

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
      // Key TEST_PRECONDITION/POSTCONDITION_VIOLATED on the same OPENMS_ASSERTIONS
      // (from <OpenMS/config.h>) that Macros.h keys the library's checks on, so the
      // two cannot diverge across build type or generator.
#ifdef OPENMS_ASSERTIONS
      OpenMS::Internal::ClassTest::setPreconditionTestsEnabled(true);
#endif
    }
  };

  const OpenMSTestSupport openms_test_support_init;
}
