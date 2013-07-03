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
    @page UTILS_QCImporter QCImporter

    @brief This application is used embed tables or pictures generated externally as attachments to existing quality parameters in the targeted run/set meant to have attachments. If no quality parameter is present an empty value one will be generated with the name of "default set name"/"default mzML file".

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
    TOPPBase("QCImporter", "produces qcml files", false)
  {
  }

protected:
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Input qcml file",false);
    setValidFormats_("in", StringList::create("qcML"));
    registerInputFile_("table", "<file>", "", "Table file that will be imported into the given qc file .", true);
    setValidFormats_("table", StringList::create("csv"));
    registerInputFile_("mapping", "<file>", "", "Mapping table of which column in the import will be represented as which qc.", true);
    setValidFormats_("mapping", StringList::create("csv"));
    registerOutputFile_("out", "<file>", "", "Output extended/reduced qcML file",true);
    setValidFormats_("out", StringList::create("qcML"));
  }

  ExitCodes main_(int, const char**)
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String in                   = getStringOption_("in");
    String out                 = getStringOption_("out");
    String mappi          	= getStringOption_("mapping");
    String tab                 = getStringOption_("table");
    
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
      
      if (map_file.size()<2) //assumed that first row is the header of table and second row is the according qc
      {
        cerr << "Error: You have to give a mapping of your table (first row is the header of table and second row is the according qc). Aborting!" << endl;
        return ILLEGAL_PARAMETERS;
      }
      StringList header,according;
      map_file.getRow(0, header);
      map_file.getRow(1, according);
      
      if (header.size() != according.size())
      {
        cerr << "Error: You have to give a mapping of your table (first row is the header of table and second row is the according qc). Aborting!" << endl;
        return ILLEGAL_PARAMETERS;
      }
      //~ std::map<String,String> mapping;
      //~ std::transform( header.begin(), header.end(), according.begin(), std::inserter(mapping, mapping.end() ), std::make_pair<String,String> );
      Size runset_col;
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
            cerr << "Error: You have to specify a correct cv with accession or name in col "<< String(i) <<". Aborting!" << endl;
            return ILLEGAL_PARAMETERS;
          }
        }
        else
        {
          const ControlledVocabulary::CVTerm& term = cv.getTerm(according[i]);
          header[i] = term.name;
        }
        if (header[i] == "raw file name") //TODO add set name as possibility!
        {
          runset_col = i;
        }
      }

      if (csv_file.size()>1)
      {
        StringList li;
        for (Size i = 0; i < csv_file.size(); ++i)
        {
          StringList li;
          csv_file.getRow(i, li);
          if (li.size() < according.size())
          {
            cerr << "Error: You have to give a correct mapping of your table - row " << String(i) <<" is too short. Aborting!" << endl;
            return ILLEGAL_PARAMETERS;
          }
          
          for (Size j = 0; j < li.size(); ++j)
          {
            if (j==runset_col)
            {
              continue;
            }
            QcMLFile::QualityParameter def;
            def.name = header[i]; ///< Name
            def.id = "default"; ///TODO !!!
            def.cvRef = "QC"; ///< cv reference ('full name')
            def.cvAcc = according[j];
            def.value = li[j];
            
            if (qcmlfile.existsRun(header[runset_col]))
            {
              qcmlfile.addRunQualityParameter(header[runset_col], def);
            }
            else if (qcmlfile.existsSet(header[runset_col]))
            {
              qcmlfile.addSetQualityParameter(header[runset_col], def);
            }
            else
            {
              cerr << "Error: You have to give a existing run or set - row " << String(i) <<" has none. Aborting!" << endl;
              return ILLEGAL_PARAMETERS;
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
