//! [doxygen_snippet_Template]
// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Maintainer $
// $Authors: Author1, Author2 $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_NewTool

    @brief Template for a new Tool

    This tool can be used for scientific stuff.

    And more scientific applications.

    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_NewTool.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_NewTool.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPNewTool :
  public TOPPBase
{
public:
  TOPPNewTool() :
    TOPPBase("NewTool", "Template for Tool creation", false)
  {

  }

protected:

  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_()
  {

  }


  // the main_ function is called after all parameters are read
  ExitCodes main_(int, const char **)
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------


    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------


    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------


    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    return ExitCodes::EXECUTION_OK;
  }

};


// the actual main function needed to create an executable
int main(int argc, const char ** argv)
{
  TOPPNewTool tool;
  return tool.main(argc, argv);
}
/// @endcond

//! [doxygen_snippet_Template]
