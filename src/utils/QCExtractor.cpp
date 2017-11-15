// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
    @page UTILS_QCExtractor QCExtractor

    @brief Extracts a table attachment of a given quality parameter from a qcML file as tabular (text) format.

    <CENTER>
      <table>
        <tr>
        <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
        <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ QCExtractor \f$ \longrightarrow \f$</td>
        <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref UTILS_QCEmbedder </td>
        </tr>
        <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref UTILS_QCShrinker </td>
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
    @verbinclude UTILS_QCExtractor.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_QCExtractor.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class TOPPQCExtractor :
  public TOPPBase
{
public:
  TOPPQCExtractor():
    TOPPBase("QCExtractor", "Extracts a table attachment to a given qc parameter.", false, {{ "Walzer M, Pernas LE, Nasso S, Bittremieux W, Nahnsen S, Kelchtermans P,  Martens, L", "qcML: An Exchange Format for Quality Control Metrics from Mass Spectrometry Experiments", "Molecular & Cellular Proteomics 2014; 13(8)" , "10.1074/mcp.M113.035907"}})

  {
  }

protected:
  void registerOptionsAndFlags_()
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

  ExitCodes main_(int, const char**)
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
