// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Author: Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/QcMLFile.h>

#include <QByteArray>
#include <QFile>
#include <QString>
#include <QFileInfo>

//~ #include <QIODevice>
#include <fstream>
#include <vector>
#include <map>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_QCShrinker QCShrinker

@brief This application is used to remove extra verbose table attachments from a qcML file that are not needed anymore, e.g. for a final report.

<CENTER>
  <table>
    <tr>
    <th ALIGN = "center"> pot. predecessor tools </td>
    <td VALIGN="middle" ROWSPAN=3> &rarr; QCShrinker &rarr;</td>
    <th ALIGN = "center"> pot. successor tools </td>
    </tr>
    <tr>
    <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_QCMerger </td>
    </tr>
    <tr>
    <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> ... </td>
    </tr>
  </table>
</CENTER>

If there is a lot of verbose or deprecated information in the given qcml file at @p in that can be purged.

- @p qp_accessions A list of cv accessions that should be removed. If empty, the usual suspects will be removed.
- @p run the file that defined the run under which the qp for the attachment is aggregated as MZML file. The file is only used to extract the run name from the file name;
- @p name if no file for the run was given (or if the target qp is contained in a set), at least a name of the target run/set containing the the qp for the attachment has to be given;

Output is in qcML format (see parameter @p out) which can be viewed directly in a modern browser (chromium, firefox, safari).

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_QCShrinker.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_QCShrinker.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class TOPPQCShrinker :
  public TOPPBase
{
public:
  TOPPQCShrinker() :
    TOPPBase("QCShrinker", "Remove unneeded or verbose table attachments from a qcml file.", 
    true, {{ "Walzer M, Pernas LE, Nasso S, Bittremieux W, Nahnsen S, Kelchtermans P,  Martens, L", "qcML: An Exchange Format for Quality Control Metrics from Mass Spectrometry Experiments", "Molecular & Cellular Proteomics 2014; 13(8)" , "10.1074/mcp.M113.035907"}})
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input qcml file");
    setValidFormats_("in", ListUtils::create<String>("qcML"));
    //~ registerFlag_("tables", "Remove all tables. (Of all runs and sets if these are not given with parameter name or run.)");
    registerStringList_("qp_accessions", "<names>", StringList(), "A list of cv accessions that should be removed. If empty, the usual suspects will be removed!", false);
    registerStringOption_("name", "<string>", "", "The name of the target run or set that contains the requested quality parameter.", false);
    registerInputFile_("run", "<file>", "", "The file from which the name of the target run that contains the requested quality parameter is taken. This overrides the name parameter!", false);
    setValidFormats_("run", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "Output extended/reduced qcML file");
    setValidFormats_("out", ListUtils::create<String>("qcML"));
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String in                   = getStringOption_("in");
    String out                  = getStringOption_("out");
    String target_run           = getStringOption_("name");
    String target_file          = getStringOption_("run");
    StringList qp_accs          = getStringList_("qp_accessions");

    //-------------------------------------------------------------
    // reading input
    //------------------------------------------------------------
    if (!target_file.empty())
    {
      target_run = QFileInfo(QString::fromStdString(target_file)).baseName();
    }

    //~ !getFlag_("tables")

    QcMLFile qcmlfile;
    qcmlfile.load(in);

    if (qp_accs.empty())
    {
      qp_accs.push_back("QC:0000044");
      qp_accs.push_back("QC:0000047");
      qp_accs.push_back("QC:0000022");
      qp_accs.push_back("QC:0000038");
      qp_accs.push_back("QC:0000049");
    }

    //TODO care for QualityParameter s
    if (target_run.empty())
    {
      for (Size i = 0; i < qp_accs.size(); ++i)
      {
        qcmlfile.removeAllAttachments(qp_accs[i]);
      }
    }
    else
    {
      for (Size i = 0; i < qp_accs.size(); ++i)
      {
        qcmlfile.removeAttachment(target_run,qp_accs[i]);
      }
    }

    qcmlfile.store(out);
    return EXECUTION_OK;
  }

};
int main(int argc, const char** argv)
{
  TOPPQCShrinker tool;
  return tool.main(argc, argv);
}

/// @endcond
