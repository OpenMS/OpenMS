// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
    @page UTILS_QCEmbedder QCEmbedder

    @brief This application is used embed tables or pictures generated externally as attachments to existing quality parameters in the targeted run/set meant to have attachments. If no quality parameter is present an empty value one will be generated with the name of "default set name"/"default mzML file".

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
    TOPPBase("QCEmbedder", "produces qcml files", false)
  {
  }

protected:
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Input qcml file",false);
    setValidFormats_("in", StringList::create("qcML"));
    registerStringOption_("qp", "<string>", "", "Target attachment table.");
    registerStringOption_("qp_acc", "<string>", "", "The accession number of the given qp, only needed if qp is not yet contained in the run/set.", false);
    registerStringOption_("name", "<String>", "", "The name of the target run or set that contains the requested quality parameter.", false);
    registerInputFile_("run", "<file>", "", "The file from which the name of the target run that contains the requested quality parameter is taken. This overrides the name parameter!", false);
    setValidFormats_("run", StringList::create("mzML"));
    registerInputFile_("plot", "<file>", "", "Plot file to be added to target quality parameter. (Plot file generated from csv output.)", false);
    setValidFormats_("plot", StringList::create("PNG"));
    registerInputFile_("table", "<file>", "", "Table file that will be added as attachment to the given qc.", false);
    setValidFormats_("table", StringList::create("csv"));
    registerOutputFile_("out", "<file>", "", "Output extended/reduced qcML file");
    setValidFormats_("out", StringList::create("qcML"));
		registerStringOption_("set/run", "<choice>", "", "If no in file, decides where to attach.",false);
    setValidStrings_("set/run", StringList::create("set,run"));

  }

  ExitCodes main_(int, const char**)
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String in                   = getStringOption_("in");
    String out                  = getStringOption_("out");
    String target_qp            = getStringOption_("qp");
    String target_run           = getStringOption_("name");
    String target_file          = getStringOption_("run");
    String plot_file            = getStringOption_("plot");
    String target_acc           = getStringOption_("qp_acc");
    String tab                  = getStringOption_("table");
		String setrun            = getStringOption_("set/run");
		
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

    if (plot_b64 != "" || tab != "")
    {
      if (plot_b64 != "")
      {
        QcMLFile::Attachment at;
        at.name = target_qp;
        //~ at.unitRef; //TODO MIME type
        //~ at.unitAcc;
        at.binary = plot_b64;
        at.cvRef = "QC";
        at.cvAcc = "QC:xxxxxxxx"; //TODO determine cv! btw create cv for plots?!

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
            QcMLFile::QualityParameter qp;
            if (target_acc != "" && target_qp != "")
            {
              qp.name = target_qp; ///< Name
              qp.id = target_run + "_" + target_acc; ///< Identifier
              qp.cvRef = "QC"; ///< cv reference
              qp.cvAcc = target_acc;
              qp.value = target_run;
              qcmlfile.addRunQualityParameter(target_run, qp);
              //TODO check if the qp are in the obo as soon as there is one

              at.qualityRef = qp.id;
              qcmlfile.addRunAttachment(target_run, at);
            }
            else
            {
              cerr << "Error: You have to specify a correct cv with accession and name. Aborting!" << endl;
              return ILLEGAL_PARAMETERS;
            }
          }
        }
      }
      if (tab != "")
      {
        QcMLFile::Attachment at;
        CsvFile csv_file(tab);
        if (csv_file.size()>1)
        {
          at.name = target_qp;
          //~ at.unitRef; //TODO MIME type
          //~ at.unitAcc;
          at.cvRef = "QC";
          at.cvAcc = "QC:xxxxxxxx"; //TODO determine cv! btw create cv for plots?!

          StringList li;
          csv_file.getRow(0, li);
          for (Size i = 0; i < li.size(); ++i)
          {
            at.colTypes.push_back(li[i]);
          }
          for (UInt i = 1; i < csv_file.size(); ++i)
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
						//if exists set/run TODO
            QcMLFile::QualityParameter qp;
            if (target_acc != "" && target_qp != "")
            {
							QcMLFile::QualityParameter def;
              qp.name = target_qp; ///< Name
              qp.id = target_run + "_" + target_acc; ///< Identifier
              qp.cvRef = "QC"; ///< cv reference
              qp.cvAcc = target_acc;
              qp.value = target_run;
              
              //TODO check if the qp are in the obo as soon as there is one

              at.qualityRef = qp.id;
							if (qcmlfile.existsSet(target_run) || setrun == "set") //TODO default name-qp  if file created new (no set/run exists)
							{
								if (in == "")
								{
									QcMLFile::QualityParameter def;
									def.name = "set name"; ///< Name
									def.id = "default set name"; ///< Identifier
									def.cvRef = "QC"; ///< cv reference
									def.cvAcc = "QC:0000058";
									def.value = "default set name";
									qcmlfile.addSetQualityParameter(target_run, def);
								}
								qcmlfile.addSetQualityParameter(target_run, qp);
								qcmlfile.addSetAttachment(target_run, at);
							}
							else
							{
								if (in == "")
								{
									def.name = "mzML file"; ///< Name
									def.id = "default mzML file"; ///< Identifier
									def.cvRef = "MS"; ///< cv reference
									def.cvAcc = "MS:1000584";
									def.value = "default mzML file";
									qcmlfile.addSetQualityParameter(target_run, def);
								}
								qcmlfile.addRunQualityParameter(target_run, qp);
								qcmlfile.addRunAttachment(target_run, at);
							}
            }
            else
            {
              cerr << "Error: You have to specify a correct cv with accession and name. Aborting!" << endl;
              return ILLEGAL_PARAMETERS;
            }
          }
        }
      }
      qcmlfile.store(out);
    }

    return EXECUTION_OK;
  }

};
int main(int argc, const char** argv)
{
  TOPPQCEmbedder tool;
  return tool.main(argc, argv);
}

/// @endcond
