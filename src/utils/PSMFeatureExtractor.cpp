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
// $Maintainer: Mathias Walzer $
// $Authors: Andreas Simon, Mathias Walzer, Matthew The $
// --------------------------------------------------------------------------
#include <OpenMS/config.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/ANALYSIS/ID/PercolatorFeatureSetHelper.h>
#include <QtCore/QFile>
#include <QtCore/QDir>
#include <QtCore/QProcess>

#include <iostream>
#include <cmath>
#include <string>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page UTILS_PSMFeatureExtractor PSMFeatureExtractor

  @brief PSMFeatureExtractor computes extra features for each input PSM

  @experimental Parts of this tool are still work in progress and usage and input requirements or output might change. (multiple_search_engine, Mascot support)

  <center>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ PSMFeatureExtractor \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeptideIndexer</td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PercolatorAdapter </td>
        </tr>
    </table>
  </center>

  <p>
PSMFeatureExtractor is search engine sensitive, i.e. it's extra features
vary, depending on the search engine. Thus, please make sure the input is
compliant with TOPP SearchengineAdapter output. Also, PeptideIndexer compliant
target/decoy annotation is mandatory.
Currently supported search engines are Comet, X!Tandem, MSGF+.
Mascot support is available but in beta development.
  </p>

  @note if you have extra features you want to pass to percolator, use the extra
 flag and list the MetaData entries containing the extra features.

  <B>The command line parameters of this tool are:</B>
  @verbinclude UTILS_PSMFeatureExtractor.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude UTILS_PSMFeatureExtractor.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class PSMFeatureExtractor :
  public TOPPBase
{
public:
  PSMFeatureExtractor() :
    TOPPBase("PSMFeatureExtractor", "Computes extra features for each input PSM.", false)
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<files>", StringList(), "Input file(s)", true);
    setValidFormats_("in", ListUtils::create<String>("mzid,idXML"));
    registerOutputFileList_("out", "<files>", StringList(), "Output file(s)", true);
    setValidFormats_("out", ListUtils::create<String>("mzid,idXML"));
    registerStringList_("extra", "<MetaData parameter>", vector<String>(), "List of the MetaData parameters to be included in a feature set for precolator.", false, false);
    // setValidStrings_("extra", ?);
    // TODO: add this MHC feature back in with TopPerc::hasMHCEnd_()
    //registerFlag_("MHC", "Add a feature for MHC ligand properties to the specific PSM.", true);
    registerFlag_("multiple_search_engines", "Combine PSMs from different search engines by merging on scan level.");
    registerFlag_("skip_db_check", "Manual override to skip the check if same settings for multiple search engines were applied. Only valid together with -multiple_search_engines flag.", true);
    registerFlag_("concat", "Naive merging of PSMs from different search engines: concatenate multiple search results instead of merging on scan level. Only valid together with -multiple_search_engines flag.", true);
    registerFlag_("impute", "Will instead of discarding all PSM not unanimously detected by all SE, impute missing values by their respective scores min/max observed. Only valid together with -multiple_search_engines flag.", true);
    registerFlag_("limit_imputation", "Will impute missing scores with the worst numerical limit (instead of min/max observed) of the respective score. Only valid together with -multiple_search_engines flag.", true);
    registerFlag_("replicates", "Add features across multiple replicate runs analysed using the same search engine. Cannot be combined with -multiple_search_engines.", false);
  }

  void writeOutputToFile_(const String out, vector<PeptideIdentification>& peptide_ids, vector<ProteinIdentification>& protein_ids, StringList& feature_set, StringList& extra_features)
  {
    ProteinIdentification::SearchParameters search_parameters = protein_ids.front().getSearchParameters();
    search_parameters.setMetaValue("feature_extractor", "TOPP_PSMFeatureExtractor");
    feature_set.insert(feature_set.end(), extra_features.begin(), extra_features.end());
    search_parameters.setMetaValue("extra_features", ListUtils::concatenate(feature_set, ","));
    protein_ids.front().setSearchParameters(search_parameters);

    // Storing the PeptideHits with calculated q-value, pep and svm score
    FileTypes::Type out_type = FileHandler::getType(out);

    LOG_INFO << "writing output file: " << out << endl;

    if (out_type == FileTypes::IDXML)
    {
      IdXMLFile().store(out, protein_ids, peptide_ids);
    }
    else if (out_type == FileTypes::MZIDENTML)
    {
      MzIdentMLFile().store(out, protein_ids, peptide_ids);
    }
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // general variables and data to perform PSMFeatureExtractor
    //-------------------------------------------------------------
    vector<PeptideIdentification> all_peptide_ids;
    vector<ProteinIdentification> all_protein_ids;
    
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    const StringList in_list = getStringList_("in");
    const StringList out_list = getStringList_("out");
    bool multiple_search_engines = getFlag_("multiple_search_engines");
    bool replicates = getFlag_("replicates");
    LOG_DEBUG << "Input file (of target?): " << ListUtils::concatenate(in_list, ",") << endl;

    if (in_list.size() > 1 && !multiple_search_engines && !replicates)
    {
      writeLog_("Fatal error: multiple input files given for -in, but -multiple_search_engines or -replicates flag not specified.");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    if (multiple_search_engines && replicates)
    {
      writeLog_("Fatal error: both -multiple_search_engines and -replicates flag specified.");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    if (replicates && in_list.size() != out_list.size())
    {
      writeLog_("Fatal error: If -replicates flag is specified, please provide one output file for each input file.");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    if (!replicates && out_list.size() > 1)
    {
      LOG_WARN << "Multiple output files given for -out, but -replicates flag not specified. Output files other than the first one will be ignored." << endl;
    }

    //-------------------------------------------------------------
    // read input
    //-------------------------------------------------------------
    bool skip_db_check = getFlag_("skip_db_check");
    bool concatenate = getFlag_("concat");
    StringList search_engines_used;
    vector<Int> peptide_ids_range {0};
    int file_count = 0;
    for (StringList::const_iterator fit = in_list.begin(); fit != in_list.end(); ++fit)
    {
      vector<PeptideIdentification> peptide_ids;
      vector<ProteinIdentification> protein_ids;
      String in = *fit;
      FileHandler fh;
      FileTypes::Type in_type = fh.getType(in);
      LOG_INFO << "Loading input file: " << in << endl;
      if (in_type == FileTypes::IDXML)
      {
        IdXMLFile().load(in, protein_ids, peptide_ids);
      }
      else if (in_type == FileTypes::MZIDENTML)
      {
        LOG_WARN << "Converting from mzid: possible loss of information depending on target format." << endl;
        MzIdentMLFile().load(in, protein_ids, peptide_ids);
      }
      //else caught by TOPPBase:registerInput being mandatory mzid or idxml

      // will check if all ProteinIdentifications have the same search db unless it is the first, in which case all_protein_ids is empty yet.
      if ((multiple_search_engines || replicates) && !skip_db_check && !all_protein_ids.empty())
      {
        ProteinIdentification::SearchParameters all_search_parameters = all_protein_ids.front().getSearchParameters();
        ProteinIdentification::SearchParameters search_parameters = protein_ids.front().getSearchParameters();
        if (search_parameters.db != all_search_parameters.db)
        {
          writeLog_("Input files are not searched with the same protein database, " + search_parameters.db + " vs. " + all_search_parameters.db + ". Set -skip_db_check flag to ignore this. Aborting!");
          return INCOMPATIBLE_INPUT_DATA;
        }
      }

      if (replicates)
      {
        // update file count (serves as identifier)
        file_count += 1;

        // for each peptide hit, remember the file it originated from (will be removed after processing)
        StringList sources;
        protein_ids.front().getPrimaryMSRunPath(sources);
        for (vector<PeptideIdentification>::iterator it = peptide_ids.begin(); it != peptide_ids.end(); ++it)
        {
          for (vector<PeptideHit>::iterator hit = it->getHits().begin(); hit != it->getHits().end(); ++hit)
          {
            hit->setMetaValue("TMP:sources", sources);

            if (sources.empty())
            {
              hit->setMetaValue("TMP:sources", "file_" + to_string(file_count));
            }
            else
            {
              hit->setMetaValue("TMP:sources", hit->getMetaValue("TMP:sources").toString());
            }
          }
        }

        all_protein_ids.insert(all_protein_ids.end(), protein_ids.begin(), protein_ids.end());
      }
      else
      {
        PercolatorFeatureSetHelper::mergeMULTISEProteinIds(all_protein_ids, protein_ids);
      }

      if (!multiple_search_engines)
      {
        all_peptide_ids.insert(all_peptide_ids.end(), peptide_ids.begin(), peptide_ids.end());
        peptide_ids_range.push_back(all_peptide_ids.size());
      }
      else
      {
        String search_engine = protein_ids.front().getSearchEngine();
        if (!ListUtils::contains(search_engines_used, search_engine))
        {
          search_engines_used.push_back(search_engine);
        }
        
        if (concatenate)
        {
          // will concatenate the list
          PercolatorFeatureSetHelper::concatMULTISEPeptideIds(all_peptide_ids, peptide_ids, search_engine);
        }
        else
        {
          // will collapse the list (based on spectrum_reference)
          PercolatorFeatureSetHelper::mergeMULTISEPeptideIds(all_peptide_ids, peptide_ids, search_engine);
        }
      }
    }
    
    if (all_protein_ids.empty())
    {
      writeLog_("No protein hits found in input file. Aborting!");
      printUsage_();
      return INPUT_FILE_EMPTY;
    }

    //-------------------------------------------------------------
    // extract search engine and prepare pin
    //-------------------------------------------------------------
    String search_engine = all_protein_ids.front().getSearchEngine();
    if (multiple_search_engines) search_engine = "multiple";
    LOG_DEBUG << "Registered search engine: " << search_engine << endl;
    
    StringList extra_features = getStringList_("extra");
    StringList feature_set;

    if (search_engine == "multiple")
    {
      if (getFlag_("concat"))
      {
        PercolatorFeatureSetHelper::addCONCATSEFeatures(all_peptide_ids, search_engines_used, feature_set);
      }
      else
      {
        bool impute = getFlag_("impute");
        bool limits = getFlag_("limit_imputation");
        PercolatorFeatureSetHelper::addMULTISEFeatures(all_peptide_ids, search_engines_used, feature_set, !impute, limits);
      }
    }
    else if (search_engine == "MS-GF+") PercolatorFeatureSetHelper::addMSGFFeatures(all_peptide_ids, feature_set);
    else if (search_engine == "Mascot") PercolatorFeatureSetHelper::addMASCOTFeatures(all_peptide_ids, feature_set);
    else if (search_engine == "XTandem") PercolatorFeatureSetHelper::addXTANDEMFeatures(all_peptide_ids, feature_set);
    else if (search_engine == "Comet") PercolatorFeatureSetHelper::addCOMETFeatures(all_peptide_ids, feature_set);
    else
    {
      LOG_ERROR << "No known input to create PSM features from. Aborting" << std::endl;
      return INCOMPATIBLE_INPUT_DATA;
    }

    if (replicates)
    {
      PercolatorFeatureSetHelper::addREPLICATEFeatures(all_peptide_ids, feature_set);

      // remove temporary meta value storing a peptide hit's source
      for (vector<PeptideIdentification>::iterator it = all_peptide_ids.begin(); it != all_peptide_ids.end(); ++it)
      {
        for (vector<PeptideHit>::iterator hit = it->getHits().begin(); hit != it->getHits().end(); ++hit)
        {
          hit->removeMetaValue("TMP:sources");
        }
      }

      for (string::size_type i = 0; i < in_list.size(); ++i)
      {
        ProteinIdentification p = all_protein_ids[i];
        vector<PeptideIdentification> peptide_ids(all_peptide_ids.begin() + peptide_ids_range[i], all_peptide_ids.begin() + peptide_ids_range[i + 1]);

        String run_identifier = p.getIdentifier();
        for (vector<PeptideIdentification>::iterator it = peptide_ids.begin(); it != peptide_ids.end(); ++it)
        {
          it->setIdentifier(run_identifier);
          PercolatorFeatureSetHelper::checkExtraFeatures(it->getHits(), extra_features);  // will remove inconsistently available features
        }
      }
    }
    else
    {
      String run_identifier = all_protein_ids.front().getIdentifier();
      for (vector<PeptideIdentification>::iterator it = all_peptide_ids.begin(); it != all_peptide_ids.end(); ++it)
      {
        it->setIdentifier(run_identifier);
        PercolatorFeatureSetHelper::checkExtraFeatures(it->getHits(), extra_features);  // will remove inconsistently available features
      }
    }

    if ((!replicates && all_protein_ids.size() > 1) || (replicates && all_protein_ids.size() != in_list.size()))
    {
      LOG_ERROR << "Multiple identifications in one file are not supported. Please resume with separate input files. Quitting." << std::endl;
      return INCOMPATIBLE_INPUT_DATA;
    }

    if (replicates)
    {
      for (string::size_type i = 0; i < out_list.size(); ++i)
      {
        vector<PeptideIdentification> peptide_ids(all_peptide_ids.begin() + peptide_ids_range[i], all_peptide_ids.begin() + peptide_ids_range[i + 1]);
        vector<ProteinIdentification> protein_ids(1, all_protein_ids[i]);

        writeOutputToFile_(out_list[i], peptide_ids, protein_ids, feature_set, extra_features);
      }
    }
    else
    {
      writeOutputToFile_(out_list[0], all_peptide_ids, all_protein_ids, feature_set, extra_features);
    }

    writeLog_("PSMFeatureExtractor finished successfully!");
    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  PSMFeatureExtractor tool;

  return tool.main(argc, argv);
}

/// @endcond
