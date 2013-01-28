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
    @page UTILS_QCShrinker QCShrinker

    @brief This application is used to provide data export from raw, id and feature data files generated via TOPP pipelines. It is intended to provide tables that can be read into R where QC metrics will be calculated.

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_QCShrinker.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_QCShrinker.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class TOPPQCShrinker :
  public TOPPBase
{
public:
  TOPPQCShrinker() :
    TOPPBase("QCShrinker", "produces qcml files", false)
  {
  }

protected:
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Input qcml file");
    setValidFormats_("in", StringList::create("qcML"));
    registerStringOption_("qp", "<choice>", "", "Target attachment table.");
    setValidStrings_("qp", StringList::create("precursor tables,charge tables,total ion current tables,delta ppm tables,feature tables"));
    registerInputFile_("name", "<file>", "", "The name of the target run or set that contains the requested quality parameter.");
    registerStringOption_("at", "<string>", "", "If given, only those attachments are being removed.",false);
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
    String target_run           = getStringOption_("name");
    String target_at            = getStringOption_("at");
    //-------------------------------------------------------------
    // reading input
    //------------------------------------------------------------

    if (target_run.hasSuffix(".mzML"))
    {
      target_run = QFileInfo(QString::fromStdString(target_run)).baseName();
    }

    QcMLFile qcmlfile;
    qcmlfile.load(in);

    std::vector<String> ids;
    qcmlfile.existsRunQualityParameter(target_run, target_qp, ids); //TODO this only works if the attachments are referencing the qp - okay for now
    qcmlfile.removeAttachment(target_run, ids, target_at);
    if ( target_at != "" )
    {
      qcmlfile.removeQualityParameter(target_run, ids);
    }

    qcmlfile.store(out);
    return EXECUTION_OK;
  }

};
int main(int argc, const char ** argv)
{
  TOPPQCShrinker tool;
  return tool.main(argc, argv);
}

/// @endcond
