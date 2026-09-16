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

class TOPPExternalTool : public TOPPBase
{
public:
  TOPPExternalTool() :
    TOPPBase("TestExternalCodeCLI", "Tool of an external project built against the installed OpenMS_CLI library.", false)
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerStringOption_("greeting", "<text>", "hello", "Text to print", false);
  }

  ExitCodes main_(int, const char**) override
  {
    std::cout << getStringOption_("greeting") << std::endl;
    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  TOPPExternalTool tool;
  return tool.main(argc, argv);
}
