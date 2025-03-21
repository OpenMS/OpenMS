// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/OMSFile.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_JSONExporter JSONExporter

@brief Converts .oms (SQLite) files to JSON

<CENTER>
<table>
<tr>
<th ALIGN = "center"> potential predecessor tools </td>
<td VALIGN="middle" ROWSPAN=2> &rarr; JSONExporter &rarr;</td>
<th ALIGN = "center"> potential successor tools </td>
</tr>
<tr>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFileConverter </td>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> external tools</td>
</tr>
</table>
</CENTER>

...

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_JSONExporter.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_JSONExporter.html
 */

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class JSONExporter: public TOPPBase
{
  public:
  JSONExporter(): TOPPBase("JSONExporter", "Exports .oms (SQLite) files in JSON format")
  {
  }

  protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file");
    setValidFormats_("in", {"oms"});
    registerOutputFile_("out", "<file>", "", "Output file");
    setValidFormats_("out", {"json"});
  }

  /// Main function
  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    String in = getStringOption_("in");
    String out = getStringOption_("out");

    OMSFile oms(log_type_);
    oms.exportToJSON(in, out);

    return EXECUTION_OK;
  }
};


int main(int argc, const char** argv)
{
  JSONExporter tool;
  return tool.main(argc, argv);
}

/// @endcond
