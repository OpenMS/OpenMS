// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include "OpenDIA.h"

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/LogStream.h>

int main()
{
  // Both calls invoke the same compiled orchestration that the executable uses.
  const char* help_args[] = {"OpenDIA", "--help"};
  if (OpenDIA::run(2, help_args) != OpenMS::TOPPBase::EXECUTION_OK)
  {
    OPENMS_LOG_ERROR << "OpenDIA help did not succeed\n";
    return 1;
  }

  // A second invocation must reject an unknown option before scientific work.
  const char* malformed_args[] = {"OpenDIA", "-wp01_unknown_option"};
  if (OpenDIA::run(2, malformed_args) != OpenMS::TOPPBase::ILLEGAL_PARAMETERS)
  {
    OPENMS_LOG_ERROR << "OpenDIA did not reject an unknown option\n";
    return 1;
  }
  return 0;
}
