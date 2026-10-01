// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

// A TOPP-style tool of an external project. It exercises the installed
// OpenMS::OpenMS_CLI target: the headers under OpenMS/APPLICATIONS/, the generated
// OpenMS_CLIConfig.h with the OPENMS_CLI_DLLAPI export macro, and linking TOPPBase.

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <iostream>

using namespace OpenMS;

/// A minimal TOPP tool: prints the value of its only option.
class TOPPExternalTool : public TOPPBase
{
public:
  /// Not a registered TOPP tool, so the ToolHandler registry check is skipped.
  TOPPExternalTool() :
    TOPPBase("TestExternalCodeCLI", "Tool of an external project built against the installed OpenMS_CLI library.", {}, false)
  {
  }

protected:
  /// Declares the single string option '-greeting'.
  void registerOptionsAndFlags_() override
  {
    registerStringOption_("greeting", "<text>", "hello", "Text to print", false);
  }

  /// Prints the greeting and reports success.
  ExitCodes main_(int, const char**) override
  {
    std::cout << getStringOption_("greeting") << std::endl;
    return EXECUTION_OK;
  }
};

/// Runs the tool through TOPPBase::main(), like every TOPP tool does.
int main(int argc, const char** argv)
{
  TOPPExternalTool tool;
  return tool.main(argc, argv);
}
