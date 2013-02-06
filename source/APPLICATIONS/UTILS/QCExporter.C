// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
    @page UTILS_QCExporter QCExporter

    @brief This application is used to provide data export from raw, id and feature data files generated via TOPP pipelines. It is intended to provide tables that can be read into R where QC metrics will be calculated.

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_QCExporter.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_QCExporter.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class TOPPQCExporter :
  public TOPPBase
{
public:
  TOPPQCExporter() :
    TOPPBase("QCExporter", "produces qcml files", false)
  {
  }

protected:
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Input qcml file");
    setValidFormats_("in", StringList::create("qcML"));
    registerStringOption_("qp", "<choice>", "", "Target attachment table.");
    setValidStrings_("qp", StringList::create("precursor tables,charge tables,total ion current tables,delta ppm tables,feature tables,set id"));
    registerStringOption_("name", "<string>", "", "The name of the target run or set that contains the requested quality parameter.", false);
    registerInputFile_("run", "<file>", "", "The file from which the name of the target run that contains the requested quality parameter is taken. This overrides the name parameter!", false);
    setValidFormats_("run", StringList::create("mzML"));
    registerOutputFile_("out_csv", "<file>", "", "Output csv formated quality parameter or extended qcML file");
    setValidFormats_("out_csv", StringList::create("csv"));
  }

  ExitCodes main_(int, const char**)
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String in                   = getStringOption_("in");
    String csv                  = getStringOption_("out_csv");
    String target_qp            = getStringOption_("qp");
    String target_run           = getStringOption_("name");
    String target_file          = getStringOption_("run");

    //-------------------------------------------------------------
    // reading input
    //------------------------------------------------------------
    if (target_file != "")
    {
      target_run = QFileInfo(QString::fromStdString(target_file)).baseName();
    }

    QcMLFile qcmlfile;
    qcmlfile.load(in);

    if (target_run == "")
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
      if (qcmlfile.existsSet(target_run))
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
  }

};
int main(int argc, const char** argv)
{
  TOPPQCExporter tool;
  return tool.main(argc, argv);
}

/// @endcond
