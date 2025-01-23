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
// TODO this is currently needed for attachments
#include <OpenMS/FORMAT/QcMLFile.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>

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
@page TOPP_QCEmbedder QCEmbedder

@brief This application is used to embed tables or plots generated externally as attachments to existing quality parameters in qcML files.

<CENTER>
  <table>
    <tr>
    <th ALIGN = "center"> pot. predecessor tools </td>
    <td VALIGN="middle" ROWSPAN=3> &rarr; QCEmbedder &rarr;</td>
    <th ALIGN = "center"> pot. successor tools </td>
    </tr>
    <tr>
    <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_QCExporter </td>
    <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_QCMerger </td>
    </tr>
    <tr>
    <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_XTandemAdapter </td>
    <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_QCShrinker </td>
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
@verbinclude TOPP_QCEmbedder.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_QCEmbedder.html

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
  void registerOptionsAndFlags_() override
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

  ExitCodes main_(int, const char**) override
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
    cv.loadFromOBO("QC", File::find("/CV/qc-cv-legacy.obo"));


    //-------------------------------------------------------------
    // reading input
    //------------------------------------------------------------
    if (!target_file.empty())
    {
      target_run = QFileInfo(QString::fromStdString(target_file)).baseName();
    }

    QcMLFile qcmlfile;
    if (!in.empty())
    {
      qcmlfile.load(in);
    }

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

    if (!plot_b64.empty() || !tab.empty())
    {
      if (!plot_b64.empty())
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
      else if (!tab.empty())
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
