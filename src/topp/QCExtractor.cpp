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
#include <iostream>
#include <fstream>
#include <vector>
#include <map>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_QCExtractor QCExtractor

@brief Extracts a table attachment of a given quality parameter from a qcML file as tabular (text) format.

<CENTER>
  <table>
    <tr>
    <th ALIGN = "center"> pot. predecessor tools </td>
    <td VALIGN="middle" ROWSPAN=3> &rarr; QCExtractor &rarr;</td>
    <th ALIGN = "center"> pot. successor tools </td>
    </tr>
    <tr>
    <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_QCEmbedder </td>
    </tr>
    <tr>
    <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_QCShrinker </td>
    </tr>
  </table>
</CENTER>

If there is a table attached to a given qp that is needed as a single file, e.g. for easy input to plotting software, this can be extracted to a tabular (text) format.

- @p qp defines the qp name to which the table is attached;
- @p run the file that defined the run under which the qp for the attachment is aggregated as mzML file. The file is only used to extract the run name from the file name.
- @p name if no file for the run was given (or if the target qp is contained in a set), at least a name of the target run/set containing the the qp for the attachment has to be given.
- @p set/run if the target qp is contained in a set, this has to be set here;

Output is in csv format (see parameter @p out_csv) which can be easily parsed by many programs. 

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_QCExtractor.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_QCExtractor.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class TOPPQCExtractor :
  public TOPPBase
{
public:
  TOPPQCExtractor():
    TOPPBase("QCExtractor", "Extracts a table attachment to a given qc parameter.", 
    true, {{ "Walzer M, Pernas LE, Nasso S, Bittremieux W, Nahnsen S, Kelchtermans P,  Martens, L", "qcML: An Exchange Format for Quality Control Metrics from Mass Spectrometry Experiments", "Molecular & Cellular Proteomics 2014; 13(8)" , "10.1074/mcp.M113.035907"}})

  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input qcml file");
    setValidFormats_("in", ListUtils::create<String>("qcML"));
    registerStringOption_("qp", "<string>", "", "Target attachment qp.");
    registerInputFile_("run", "<file>", "", "The file that defined the run under which the qp for the attachment is aggregated as mzML file. The file is only used to extract the run name from the file name.", false);
    setValidFormats_("run", ListUtils::create<String>("mzML"));
    registerStringOption_("name", "<string>", "", "If no file for the run was given (or if the target qp is contained in a set), at least a name of the target run/set containing the the qp for the attachment has to be given.", false);
    registerOutputFile_("out_csv", "<file>", "", "Output csv formatted table.");
    setValidFormats_("out_csv", ListUtils::create<String>("csv"));
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String in                    = getStringOption_("in");
    String csv                  = getStringOption_("out_csv");
    String target_qp        = getStringOption_("qp");
    String target_run       = getStringOption_("name");
    String target_file        = getStringOption_("run");

    //-------------------------------------------------------------
    // reading input
    //------------------------------------------------------------
    if (!target_file.empty())
    {
      target_run = QFileInfo(QString::fromStdString(target_file)).baseName();
    }

    QcMLFile qcmlfile;
    qcmlfile.load(in);

    if (target_run.empty())
    {
      //~ check if only one run in file
      std::vector<String> nas;
      qcmlfile.getRunNames(nas);
      if (nas.size() == 1)
      {
        target_run = nas.front();
      }
      else
      {
        cerr << "Error: You have to give at least one of the following parameter (in ascending precedence): name, run. Aborting!" << endl;
        return ILLEGAL_PARAMETERS;
      }
    }

    String csv_str = "";
    if (target_qp == "set id")
    {
      if (qcmlfile.existsSet(target_run,true))
      {
        csv_str = qcmlfile.exportIDstats(target_run);
      }
      else
      {
        cerr << "Error: You have to specify a existing set for this qp. " << target_run << " seems not to exist. Aborting!" << endl;
        return ILLEGAL_PARAMETERS;
      }
    }
    else
    {
      //TODO warn when target_run is empty or not present in qcml
      csv_str = qcmlfile.exportAttachment(target_run, target_qp);
    }

    ofstream fout(csv.c_str());
    fout << csv_str << endl;
    fout.close();
    //~ qcmlfile.store(out);

    return EXECUTION_OK;
    //~ TODO export table containing all given qp
  }

};
int main(int argc, const char** argv)
{
  TOPPQCExtractor tool;
  return tool.main(argc, argv);
}

/// @endcond
