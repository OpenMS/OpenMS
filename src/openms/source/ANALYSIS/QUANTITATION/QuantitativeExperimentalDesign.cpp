// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Immanuel Luhn$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/QUANTITATION/QuantitativeExperimentalDesign.h>

#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>

#include <QDir>

using std::map;
using std::vector;
using std::set;
using std::pair;
using std::endl;

//using namespace std;

namespace OpenMS
{
  QuantitativeExperimentalDesign::QuantitativeExperimentalDesign() :
    DefaultParamHandler("QuantitativeExperimentDesign")
  {
    defaults_.setValue("designer:experiment", "ExperimentalSetting", "Identifier for the experimental design.");
    defaults_.setValue("designer:file", "File", "Identifier for the file name.");

    defaults_.setValue("designer:separator", "tab", "Separator, which should be used to split a row into columns");
    defaults_.setValidStrings("designer:separator", ListUtils::create<String>("tab,semi-colon,comma,whitespace"));

    defaults_.setSectionDescription("designer", "Additional options for quantitative experimental design");

    defaultsToParam_();
  }

  QuantitativeExperimentalDesign::~QuantitativeExperimentalDesign()
  {
  }

  void QuantitativeExperimentalDesign::applyDesign2Resolver(ProteinResolver& resolver, TextFile& file, StringList& file_paths)
  {

    //create mapping from experimental setting to all respective file names
    map<String, StringList> design2FileBaseName;
    mapFiles2Design_(design2FileBaseName, file);
    //filter out all non-existing files
    map<String, StringList> design2FilePath;
    findRelevantFilePaths_(design2FileBaseName, design2FilePath, file_paths);

    //determine wether we deal with idXML or featureXML
    FileTypes::Type in_type = FileHandler::getType(file_paths.front());

    if (in_type == FileTypes::IDXML)
    {
      vector<ProteinIdentification> proteins;
      vector<PeptideIdentification> peptides;

      for (map<String, StringList>::iterator iter =  design2FilePath.begin(); iter != design2FilePath.end(); ++iter)
      {
        // merge the respective files
        mergeIDFiles_(proteins, peptides, iter->first, iter->second);
      }

      resolver.resolveID(peptides);
    }
    else
    {
      ConsensusMap consensus;

      for (map<String, StringList>::iterator iter =  design2FilePath.begin(); iter != design2FilePath.end(); ++iter)
      {
        mergeConsensusMaps_(consensus, iter->first, iter->second);
      }

      resolver.resolveConsensus(consensus);
    }
  }

  void QuantitativeExperimentalDesign::mergeConsensusMaps_(ConsensusMap& out, const String& experiment, StringList& file_paths)
  {
    ConsensusMap map;

    OPENMS_LOG_INFO << "Merge consensus maps: " << endl;
    UInt counter = 1;
    for (StringList::iterator file_it = file_paths.begin(); file_it != file_paths.end(); ++file_it, ++counter)
    {
      //load should clear the map
      ConsensusXMLFile().load(*file_it, map);
      for (ConsensusMap::iterator it = map.begin(); it != map.end(); ++it)
      {
        it->setMetaValue("experiment", DataValue(experiment));
      }
      out.appendRows(map);
    }
    OPENMS_LOG_INFO << endl;
  }

  void QuantitativeExperimentalDesign::mergeIDFiles_(vector<ProteinIdentification>& proteins, vector<PeptideIdentification>& peptides, const String& experiment, StringList& file_paths)
  {
    set<String> used_ids;
    vector<ProteinIdentification> additional_proteins;
    vector<PeptideIdentification> additional_peptides;

    OPENMS_LOG_INFO << "Merge idXML-files:" << endl;
    for (StringList::iterator file_it = file_paths.begin(); file_it != file_paths.end(); ++file_it)
    {
      // load should clear the vectors
      IdXMLFile().load(*file_it, additional_proteins, additional_peptides);

      for (vector<ProteinIdentification>::iterator prot_it =
             additional_proteins.begin(); prot_it !=
           additional_proteins.end(); ++prot_it)
      {
        prot_it->setMetaValue("experiment", DataValue(experiment));
      }

      for (vector<PeptideIdentification>::iterator pep_it =
             additional_peptides.begin(); pep_it !=
           additional_peptides.end(); ++pep_it)
      {
        pep_it->setMetaValue("experiment", DataValue(experiment));
      }

      UInt counter = 1;
      for (vector<ProteinIdentification>::iterator prot_it = additional_proteins.begin(); prot_it != additional_proteins.end(); ++prot_it, ++counter)
      {
        String id = prot_it->getIdentifier();
        if (used_ids.find(id) != used_ids.end()) // ID used previously
        {
          OPENMS_LOG_INFO << "Warning: The identifier '" + id + "' was used before!" << endl;
          // generate a new ID:
          DateTime date_time = prot_it->getDateTime();
          String new_id;
          String search_engine = prot_it->getSearchEngine();
          
          do
          {
            date_time = date_time.addSecs(1);
            new_id = search_engine + "_" + date_time.toString(Qt::ISODate);
          } while (used_ids.find(new_id) != used_ids.end());

          OPENMS_LOG_INFO << "New identifier '" + new_id + "' generated as replacement." << endl;
          // update fields:
          prot_it->setIdentifier(new_id);
          prot_it->setDateTime(date_time);
          for (vector<PeptideIdentification>::iterator pep_it = additional_peptides.begin(); pep_it != additional_peptides.end(); ++pep_it)
          {
            if (pep_it->getIdentifier() == id)
              pep_it->setIdentifier(new_id);
          }
          used_ids.insert(new_id);
        }
        else
          used_ids.insert(id);
      }

      proteins.insert(proteins.end(), additional_proteins.begin(), additional_proteins.end());
      peptides.insert(peptides.end(), additional_peptides.begin(), additional_peptides.end());
    }
  }

  void QuantitativeExperimentalDesign::findRelevantFilePaths_(map<String, StringList>& design2FileBaseName, map<String, StringList>& design2FilePath, StringList& filePaths)
  {
    //find all file from the input file that belong to an experimental setting
    // files without a mapping are ignored

    // for every experimental setup
    for (map<String, StringList>::iterator iter =  design2FileBaseName.begin(); iter != design2FileBaseName.end(); ++iter)
    {
      StringList& files_base_name_design = iter->second;

      StringList existing_files_input;

      // for every base file name
      for (StringList::iterator it = files_base_name_design.begin(); it != files_base_name_design.end(); ++it)
      {
        // search against all files from the user input
        for (StringList::iterator it2 = filePaths.begin(); it2 != filePaths.end(); ++it2)
        {
          // QFileInfo fi("/tmp/archive.tar.gz");
          // QString name = fi.baseName(); --> name = "archive"
          const String file_ = QFileInfo(it2->toQString()).baseName().toStdString();
          // if given store file path in string list
          if (it->compare(file_) == 0)
          {
            existing_files_input.push_back(*it2);
          }
        }
      }
      // iff files are provided for an setup, create a map entry
      if (!existing_files_input.empty())
        design2FilePath.insert(make_pair(iter->first, existing_files_input));
    }
  }

  void QuantitativeExperimentalDesign::analyzeHeader_(UInt& expCol, UInt& fileCol, StringList& header)
  {
    // read parameter
    String experiment = param_.getValue("designer:experiment");
    String fileName = param_.getValue("designer:file");

    // iterate through header strings to look for matching identifier
    UInt col = 0;
    for (StringList::iterator iter = header.begin(); iter != header.end(); ++iter, ++col)
    {
      if (experiment.compare(*iter) == 0)
        expCol = col;
      if (fileName.compare(*iter) == 0)
        fileCol = col;
    }

    // in case one or all identifier could not be found throw an exception
    UInt invalid = -1;
    if (expCol == invalid || fileCol == invalid)
    {
      if (expCol == invalid && fileCol == invalid)
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Both identifier (experimental design and file name) are not correct");
      if (expCol == invalid)
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Identifier for experimental design is not correct");
      if (fileCol == invalid)
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Identifier for the file name is not correct");
    }
  }

  void QuantitativeExperimentalDesign::getSeparator_(String& separator)
  {
    // get separator from parameter setting
    String sep = param_.getValue("designer:separator");

    // assign
    if (sep.compare("tab") == 0)
      separator = "\t";
    else if (sep.compare("semi-colon") == 0)
      separator = ";";
    else if (sep.compare("comma") == 0)
      separator = ",";
    else if (sep.compare("whitespace") == 0)
      separator = " ";
  }

  void QuantitativeExperimentalDesign::mapFiles2Design_(map<String, StringList>& experiments, TextFile& file)
  {
    // get the defined separator from the parameter setting
    String separator;
    getSeparator_(separator);

    // read the header and split according separator
    StringList header;
    TextFile::ConstIterator titer = file.begin();
    titer->split(separator, header);
    ++titer;

    // define the column of file name and experimental setting
    UInt expCol = -1;
    UInt fileCol = -1;
    analyzeHeader_(expCol, fileCol, header);

    // read rest of the file, each row is already split according to separator
    vector<StringList> rows;
    for (; titer != file.end(); ++titer)
    {
      StringList column;
      titer->split(separator, column);
      rows.push_back(column);
    }

    // map all file names to the respective experimental setting
    map<String, StringList>::iterator it;

    for (vector<StringList>::iterator liter = rows.begin(); liter != rows.end(); ++liter)
    {
      // get experimental setting and file name
      String experiment = liter->at(expCol);
      String fileName = liter->at(fileCol);

      // search for experimental setting
      it = experiments.find(experiment);

      // if experimental setting is already present, add file name
      if (it != experiments.end())
      {
        StringList& list = it->second;
        list.push_back(fileName);
      }
      // otherwise create new list
      else
      {
        StringList newList;
        newList.push_back(fileName);
        experiments.insert(make_pair(experiment, newList));
      }
    }

    OPENMS_LOG_INFO << "\n Statistics: \n";
    for (it = experiments.begin(); it != experiments.end(); ++it)
    {
      OPENMS_LOG_INFO << "Experiment: " << it->first << ", number datasets: " << it->second.size() << endl;
    }
  }

} // namespace OpenMS
