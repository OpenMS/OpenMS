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
    @page UTILS_QCImporter QCImporter

    @brief Will import several quality parameter from a tabular (text) format into a qcML file - counterpart to QCExporter.

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

    If there is additional data from external tools in tabular format containing additional quality parameter (qp) to runs or sets, or even new runs, these can be imported into the qcML file. For an example see the examples in the share directory.

    - @p table The table containing the additional qp values in the columns. First row is considered containing the header. The target run or set names/ids are indicated by column "raw data file", so each row after the header will contain the values of qps for that run.
    - @p mapping The mapping of the table header to the according qp cvs, also in csv format. The first row is considered containing the headers as in the table. The second row is considered the according qp cv accessions.

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_QCImporter.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_QCImporter.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class TOPPQCImporter :
  public TOPPBase
{
public:
  TOPPQCImporter() :
    TOPPBase("QCImporter", "Imports tables with quality control parameters into qcml files.", false, {{ "Walzer M, Pernas LE, Nasso S, Bittremieux W, Nahnsen S, Kelchtermans P,  Martens, L", "qcML: An Exchange Format for Quality Control Metrics from Mass Spectrometry Experiments", "Molecular & Cellular Proteomics 2014; 13(8)" , "10.1074/mcp.M113.035907"}})
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input qcml file", false);
    setValidFormats_("in", ListUtils::create<String>("qcML"));
    registerInputFile_("table", "<file>", "", "The table containing the additional qp values in the columns. First row is considered containing the header. The target run or set names/ids are indicated by column \"raw data file\", so each row after the header will contain the values of qps for that run. (csv without \"!)", true);
    setValidFormats_("table", ListUtils::create<String>("csv"));
    registerInputFile_("mapping", "<file>", "", "The mapping of the table header to the according qp cvs, also in csv format. The first row is considered containing the headers as in the table. The second row is considered the according qp cv accessions. (csv without \"!)", true);
    setValidFormats_("mapping", ListUtils::create<String>("csv"));
    registerOutputFile_("out", "<file>", "", "Output extended qcML file", true);
    setValidFormats_("out", ListUtils::create<String>("qcML"));
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    String mappi = getStringOption_("mapping");
    String tab = getStringOption_("table");

    ControlledVocabulary cv;
    cv.loadFromOBO("PSI-MS", File::find("/CV/psi-ms.obo"));
    cv.loadFromOBO("QC", File::find("/CV/qc-cv.obo"));

    //-------------------------------------------------------------
    // reading input
    //------------------------------------------------------------
    QcMLFile qcmlfile;
    if (in != "")
    {
      qcmlfile.load(in);
    }

    if (mappi != "" && tab != "")
    {
      CsvFile csv_file(tab);
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
      int runset_col = -1;
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
            //~ cerr << "Header was: "<< header[i] << " , according value was: " << according[i] << endl;
            return ILLEGAL_PARAMETERS;
          }
        }
        else
        {
          const ControlledVocabulary::CVTerm& term = cv.getTerm(according[i]);
          header[i] = term.name;
        }
        if (header[i] == "raw data file") //TODO add set name as possibility!
        {
          runset_col = i;
        }
      }
      if (runset_col < 0)
      {
        cerr << "Error: You have to give a mapping of your table - rows to runs/sets. Aborting!" << endl;
        return ILLEGAL_PARAMETERS;
      }

      if (csv_file.rowCount() > 1)
      {
        StringList li;
        for (Size i = 1; i < csv_file.rowCount(); ++i)
        {
          StringList li;
          csv_file.getRow(i, li);
          if (li.size() < according.size())
          {
            cerr << "Error: You have to give a correct mapping of your table - row " << String(i + 1) << " is too short. Aborting!" << endl;
            return ILLEGAL_PARAMETERS;
          }

          std::vector<QcMLFile::QualityParameter> qps;
          String id;
          bool set = false;
          for (Size j = 0; j < li.size(); ++j)
          {
            if (j == static_cast<Size>(runset_col))
            {
              if (qcmlfile.existsRun(li[j])) //TODO this only works for real run IDs
              {
                id = li[j];
              }
              else if (qcmlfile.existsSet(li[j])) //TODO this only works for real set IDs
              {
                id = li[j];
                set = true;
              }
              else
              {
                id = li[j];
                qcmlfile.registerRun(id, id);
                //TODO warn that if this was supposed to be a set - now it is not!
              }
            }
            QcMLFile::QualityParameter def;
            def.name = header[j]; ///< Name
            def.id = String(UniqueIdGenerator::getUniqueId());
            def.cvRef = "QC"; ///< cv reference ('full name')
            def.cvAcc = according[j];
            def.value = li[j];
            qps.push_back(def);
          }
          if (id != "")
          {
            for (std::vector<QcMLFile::QualityParameter>::const_iterator qit = qps.begin(); qit != qps.end(); ++qit)
            {
              if (!set)
              {
                qcmlfile.addRunQualityParameter(id, *qit);
              }
              else
              {
                qcmlfile.addSetQualityParameter(id, *qit);
              }
            }
          }
        }
      }
    }
    qcmlfile.store(out);
    return EXECUTION_OK;
  }

};
int main(int argc, const char** argv)
{
  TOPPQCImporter tool;
  return tool.main(argc, argv);
}

/// @endcond
