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
// $Author: Mathias Walzer, Sven Nahnsen $
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

    @brief This application is used to provide data export from raw, id and feature data files generated via TOPP pipelines. It is intended to provide tables that can be read into R where QC metrics will be calculated.

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
    registerInputFile_("in", "<file>", "", "Input qcml file");
    setValidFormats_("in", StringList::create("qcML"));
    registerStringOption_("qp", "<choice>", "", "Target QualityParameter. ('ticpoints': total ion current stats, 'precursorpoints': stats about the precursors from the rawfile, 'pepidpoints': identifications from the rawfile, 'featurepoints': stats about features from the rawfile, 'consensusstats': stats about consensus features from the rawfile)\n");
    setValidStrings_("qp", StringList::create("ticpoints,precursorpoints,pepidpoints,featurepoints,consensuspoints"));
    registerStringOption_("run_name", "<string>", "", "The name of the target msrun file of the respective quality parameter.");
    registerInputFile_("plot", "<file>", "", "Plot file to be added to target quality parameter. (Plot file generated from csv output.)");
    setValidFormats_("plot", StringList::create("PNG"));
    registerOutputFile_("out", "<file>", "", "Output extended/reduced qcML file");
    setValidFormats_("out",StringList::create("qcML"));
  }

  ExitCodes main_(int, const char **)
  {
    String plot_file = "";
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String in                   = getStringOption_("in");
    String out                  = getStringOption_("out");
    String target_qp            = getStringOption_("qp");
    String target_run           = getStringOption_("msrun_filename");
    plot_file                   = getStringOption_("plot");
    //-------------------------------------------------------------
    // reading input
    //------------------------------------------------------------

    QcMLFile qcmlfile;
    qcmlfile.load(in);


    QFile f(plot_file.c_str());
    String plot_b64;
    if (f.open(QIODevice::ReadOnly))
    {
      QByteArray ba = f.readAll();
      f.close();
      plot_b64 = String(QString(ba.toBase64()));
    }
    if (plot_b64 != "")
    {
      QcMLFile::Attachment at;
      at.name = target_qp + "_plot";
      at.cvRef = "QC";
      at.cvAcc = "QC:xxxxxxxx";
      //~ at.unitRef; //TODO MIME type
      //~ at.unitAcc;
      at.binary = plot_b64;
      at.qualityRef = qcmlfile.existsQualityParameter(target_run, target_qp);
      qcmlfile.addAttachment(target_run, at);
      qcmlfile.store(out);
    }

    return EXECUTION_OK;
  }

};
int main(int argc, const char ** argv)
{
  TOPPQCEmbedder tool;
  return tool.main(argc, argv);
}

/// @endcond
