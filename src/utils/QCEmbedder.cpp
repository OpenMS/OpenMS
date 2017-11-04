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
#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>

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
    @page UTILS_QCEmbedder QCEmbedder

    @brief This application is used to embed tables or plots generated externally as attachments to existing quality parameters in qcML files.

    <CENTER>
      <table>
        <tr>
        <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
        <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ QCEmbedder \f$ \longrightarrow \f$</td>
        <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref UTILS_QCExporter </td>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref UTILS_QCMerger </td>
        </tr>
        <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_XTandemAdapter </td>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref UTILS_QCShrinker </td>
        </tr>
      </table>
    </CENTER>

    If there is additional data from external tools to a certain quality parameter (qp) in the qcML file at @p in, it can be attached in tabluar (csv) format or as png image file.
    If no corresponding quality parameter is present an empty value one will be generated with the name of "default set name"/"default mzML file".

    - @p qp_att_acc defines the qp cv accession of the qp to which the table/image is attached.
    - @p cv_acc defines the cv accession of the attachment.
    - @p run the file that defined the run under which the qp for the attachment is aggregated as mzML file. The file is only used to extract the run name from the file name.
    - @p name if no file for the run was given (or if the target qp is contained in a set), at least a name of the target run/set containing the the qp for the attachment has to be given.

    - @p plot if a plot image is to be attached to a qp, this has to be specified here.
    - @p table if a table is to be attached to a qp, this has to be specified here.

    Output is in qcML format (see parameter @p out) which can be viewed directly in a modern browser (chromium, firefox, safari).

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_QCEmbedder.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_QCEmbedder.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class TOPPQCEmbedder :
  public TOPPBase
{
public:
  TOPPQCEmbedder() :
    TOPPBase("QCEmbedder", "Attaches a table or an image to a given qc parameter.", false, {{ "Walzer M, Pernas LE, Nasso S, Bittremieux W, Nahnsen S, Kelchtermans P,  Martens, L", "qcML: An Exchange Format for Quality Control Metrics from Mass Spectrometry Experiments", "Molecular & Cellular Proteomics 2014; 13(8)" , "10.1074/mcp.M113.035907"}})
  {
  }

protected:
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Input qcml file", false);
    setValidFormats_("in", ListUtils::create<String>("qcML"));
    registerStringOption_("qp_att_acc", "<string>", "", "Defines the qp cv accession of the qp to which the table/image is attached.", false);
    registerStringOption_("cv_acc", "<string>", "", "Defines the cv accession of the attachment.");
    registerInputFile_("run", "<file>", "", "The file that defined the run under which the qp for the attachment is aggregated as mzML file. The file is only used to extract the run name from the file name.", false);
    setValidFormats_("run", ListUtils::create<String>("mzML"));
    registerStringOption_("name", "<String>", "", "If no file for the run was given (or if the target qp is contained in a set), at least a name of the target run/set containing the the qp for the attachment has to be given.", false);
    registerInputFile_("plot", "<file>", "", "If a plot image is to be attached to a qp, this has to be specified here.", false);
    setValidFormats_("plot", ListUtils::create<String>("PNG"));
    registerInputFile_("table", "<file>", "", "If a table is to be attached to a qp, this has to be specified here.", false);
    setValidFormats_("table", ListUtils::create<String>("csv"));
    registerOutputFile_("out", "<file>", "", "Output extended qcML file");
    setValidFormats_("out", ListUtils::create<String>("qcML"));
  }

  ExitCodes main_(int, const char**)
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String in                   = getStringOption_("in");
    String out                 = getStringOption_("out");
    String target_qp       = getStringOption_("qp_att_acc");
    String target_acc      = getStringOption_("cv_acc");
    String target_run      = getStringOption_("name");
    String target_file       = getStringOption_("run");
    String plot_file          = getStringOption_("plot");
    String tab                 = getStringOption_("table");

    //-------------------------------------------------------------
    // fetch vocabularies
    //------------------------------------------------------------
    ControlledVocabulary cv;
    cv.loadFromOBO("PSI-MS", File::find("/CV/psi-ms.obo"));
    cv.loadFromOBO("QC", File::find("/CV/qc-cv.obo"));

    //-------------------------------------------------------------
    // reading input
    //------------------------------------------------------------
    if (target_file != "")
    {
      target_run = QFileInfo(QString::fromStdString(target_file)).baseName();
    }

    QcMLFile qcmlfile;
    if (in != "")
    {
      qcmlfile.load(in);
    }

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

    QFile f(plot_file.c_str());
    String plot_b64;
    if (f.open(QIODevice::ReadOnly))
    {
      QByteArray ba = f.readAll();
      f.close();
      plot_b64 = String(QString(ba.toBase64()));
    }

    QcMLFile::Attachment at;
    at.cvAcc = target_acc;
    at.id = String(UniqueIdGenerator::getUniqueId());
    at.cvRef = "QC"; //TODO assign right cv reference

    if (plot_b64 != "" || tab != "")
    {
      if (plot_b64 != "")
      {
        try
        {
          const ControlledVocabulary::CVTerm& term = cv.getTerm(target_acc);
          at.name = term.name; ///< Name
          //~ at.unitRef; //TODO MIME type
          //~ at.unitAcc;
        }
        catch (...)
        {
          cerr << "Error: You have to give the accession of a existing cv term. Aborting!" << endl;
          return ILLEGAL_PARAMETERS;
        }
        at.binary = plot_b64;
      }
      else if (tab != "")
      {
        try
        {
          const ControlledVocabulary::CVTerm& term = cv.getTerm(target_acc);
          at.name = term.name; ///< Name
          //~ at.unitRef; //TODO MIME type
          //~ at.unitAcc;
        }
        catch (...)
        {
          cerr << "Error: You have to give the accession of a existing cv term. Aborting!" << endl;
          return ILLEGAL_PARAMETERS;
        }

        CsvFile csv_file(tab);
        if (csv_file.rowCount() > 1)
        {
          StringList li;
          csv_file.getRow(0, li);
          for (Size i = 0; i < li.size(); ++i)
          {
            at.colTypes.push_back(li[i]);
          }
          for (UInt i = 1; i < csv_file.rowCount(); ++i)
          {
            StringList li;
            std::vector<String> v;
            csv_file.getRow(i, li);
            //TODO throw error if li.size() != at.colTypes.size()
            for (Size i = 0; i < li.size(); ++i)
            {
              v.push_back(li[i]);
            }
            at.tableRows.push_back(v);
          }
        }
      }
      else
      {
        cerr << "Error: Nothing valid to attach. Aborting!" << endl;
        return ILLEGAL_PARAMETERS;
      }

      std::vector<String> ids;
      qcmlfile.existsRunQualityParameter(target_run, target_qp, ids);
      if (!ids.empty())
      {
        at.qualityRef = ids.front();
        qcmlfile.addRunAttachment(target_run, at);
      }
      else
      {
        qcmlfile.existsSetQualityParameter(target_run, target_qp, ids);
        if (!ids.empty())
        {
          at.qualityRef = ids.front();
          qcmlfile.addSetAttachment(target_run, at);
        }
        else
        {
          cerr << "Error: You have to give the accession of a existing cv term to attacht to. Aborting!" << endl;
          return ILLEGAL_PARAMETERS;
        }
      }
    }
    qcmlfile.store(out);
    return EXECUTION_OK;
  }

};
int main(int argc, const char** argv)
{
  TOPPQCEmbedder tool;
  return tool.main(argc, argv);
}

/// @endcond
