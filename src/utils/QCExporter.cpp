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

    @brief Will extract several quality parameter from several run/sets from a qcML file into a tabular (text) format - counterpart to QCImporter.

    <CENTER>
      <table>
        <tr>
        <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
        <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ QCExporter \f$ \longrightarrow \f$</td>
        <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> ? </td>
        </tr>
        <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref UTILS_QCEmbedder </td>
        </tr>
      </table>
    </CENTER>

    The data contained as values of the qp of a qcML file at @p in can be exported in tabluar (csv) format.

    - @p names The name of the target runs or sets to be exported from. If empty, from all will be exported.
    - @p mapping The mapping of the exported table's headers to the according qp cvs. The first row is considered containing the headers as for the exported the table. The second row is considered the according qp cv accessions of the qp to be exported.

    Output is in csv format (see parameter @p out_csv) which can be easily viewed/parsed by many programs.

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
    TOPPBase("QCExporter", "Will extract several qp from several run/sets in a tabular format.", false, {{ "Walzer M, Pernas LE, Nasso S, Bittremieux W, Nahnsen S, Kelchtermans P,  Martens, L", "qcML: An Exchange Format for Quality Control Metrics from Mass Spectrometry Experiments", "Molecular & Cellular Proteomics 2014; 13(8)" , "10.1074/mcp.M113.035907"}})
  {
  }

protected:
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Input qcml file");
    setValidFormats_("in", ListUtils::create<String>("qcML"));
    registerStringList_("names", "<names>", StringList(), "The name of the target runs or sets to be exported from. If empty, from all will be exported.", false);
    registerInputFile_("mapping", "<file>", "", "The mapping of the exported table's headers to the according qp cvs. The first row is considered containing the headers as for the exported the table. The second row is considered the according qp cv accessions of the qp to be exported.", true);
    setValidFormats_("mapping", ListUtils::create<String>("csv"));
    registerOutputFile_("out_csv", "<file>", "", "Output csv formatted quality parameter.");
    setValidFormats_("out_csv", ListUtils::create<String>("csv"));
  }

  ExitCodes main_(int, const char**)
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String in = getStringOption_("in");
    String csv = getStringOption_("out_csv");
    StringList names = getStringList_("names");
    String mappi = getStringOption_("mapping");

    ControlledVocabulary cv;
    cv.loadFromOBO("PSI-MS", File::find("/CV/psi-ms.obo"));
    cv.loadFromOBO("QC", File::find("/CV/qc-cv.obo"));
    //-------------------------------------------------------------
    // reading input
    //------------------------------------------------------------
    QcMLFile qcmlfile;
    qcmlfile.load(in);

    if (mappi != "")
    {
      CsvFile map_file(mappi);

      if (map_file.rowCount() < 2) //assumed that first row is the header of table and second row is the according qc
      {
        cerr << "Error: You have to give a mapping of your table (first row is the header of table and second row is the according qc). Aborting!" << endl;
        return ILLEGAL_PARAMETERS;
      }
      StringList header, according;
      map_file.getRow(0, header);
      map_file.getRow(1, according);
      if (header.size() != according.size())
      {
        cerr << "Error: You have to give a mapping of your table (first row is the header of table and second row is the according qc). Aborting!" << endl;
        return ILLEGAL_PARAMETERS;
      }
      //~ std::map<String,String> mapping;
      //~ std::transform( header.begin(), header.end(), according.begin(), std::inserter(mapping, mapping.end() ), std::make_pair<String,String> );
      //~ Size runset_col;
      for (Size i = 0; i < according.size(); ++i)
      {
        if (!cv.exists(according[i]))
        {
          try
          {
            const ControlledVocabulary::CVTerm& term = cv.getTermByName(according[i]);
            header[i] = term.name;
            according[i] = term.id;
          }
          catch (...)
          {
            cerr << "Error: You have to specify a correct cv with accession or name in col " << String(i) << ". Aborting!" << endl;
            return ILLEGAL_PARAMETERS;
          }
        }
        //~ else
        //~ {
        //~ const ControlledVocabulary::CVTerm& term = cv.getTerm(according[i]);
        //~ header[i] = term.name; //TODO what if custom headers are needed?!
        //~ }
        //~ if (header[i] == "raw file name")
        //~ {
        //~ runset_col = i;
        //~ }
      }

      if (names.size() < 1)
      {
        std::vector<String> ns;
        qcmlfile.getRunIDs(ns); //n.b. names are ids
        names = StringList(ns); //TODO also  sets
      }

      String csv_str = ListUtils::concatenate(header, ",");
      csv_str += '\n';
      for (Size i = 0; i < names.size(); ++i)
      {
        //~ if (qcmlfile.existsRun(names[i]))
        //~ {
        csv_str += qcmlfile.exportQPs(names[i], according);
        csv_str += '\n';
        //~ }
        //~ else if (qcmlfile.existsSet(names[i]))
        //~ {
        //~ csv_str += qcmlfile.exportSetQP(names[i],according);
        //~ }
        //~ else
        //~ {
        //~ cerr << "Error: You have to specify a existing set for this qp. " << names[i] << " seems not to exist. Aborting!" << endl;
        //~ return ILLEGAL_PARAMETERS;
        //~ }
      }

      ofstream fout(csv.c_str());
      fout << csv_str << endl;
      fout.close();
      //~ qcmlfile.store(out);

      //~ return EXECUTION_OK;
      //~ TODO export table containing all given qp
    }
    return EXECUTION_OK;
  }

};
int main(int argc, const char** argv)
{
  TOPPQCExporter tool;
  return tool.main(argc, argv);
}

/// @endcond
