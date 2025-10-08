// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/DigestionEnzymeProtein.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_OpenMSDatabasesInfo OpenMSDatabasesInfo

@brief Information about OpenMS' internal databases

This util prints the content of OpenMS' enzyme and modification databases to a TSV file.
<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_OpenMSDatabasesInfo.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_OpenMSDatabasesInfo.html
 */

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class OpenMSDatabasesInfo final :
    public TOPPBase
{
public:
  OpenMSDatabasesInfo() :
      TOPPBase("OpenMSDatabasesInfo", "Prints the content of OpenMS' enzyme and modification databases to TSV")
  {
  }

protected:
  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_() override
  {
    // Output CSV file
    registerOutputFile_("enzymes_out", "<out>", "", "Currently supported enzymes as TSV", true, false);
    setValidFormats_("enzymes_out", ListUtils::create<String>("tsv"));
    registerOutputFile_("mods_out", "<out>", "", "Currently supported modifications as TSV", true, false);
    setValidFormats_("mods_out", ListUtils::create<String>("tsv"));
  }

  // the main_ function is called after all parameters are read
  ExitCodes main_(int, const char**) final 
  {
    auto* enz_db = ProteaseDB::getInstance();
    enz_db->writeTSV(getStringOption_("enzymes_out"));

    auto* mod_db = ModificationsDB::getInstance();
    mod_db->writeTSV(getStringOption_("mods_out"));
    
    return EXECUTION_OK;
  }
};

// the actual main function needed to create an executable
int main(int argc, const char **argv)
{
  OpenMSDatabasesInfo tool;
  return tool.main(argc, argv);
}


/// @endcond
