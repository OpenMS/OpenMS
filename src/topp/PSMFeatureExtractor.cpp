// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Andreas Simon, Mathias Walzer, Matthew The $
// --------------------------------------------------------------------------
#include <OpenMS/config.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
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
@page TOPP_PSMFeatureExtractor PSMFeatureExtractor

@brief PSMFeatureExtractor computes extra features for each input PSM

@experimental Parts of this tool are still work in progress and usage and input requirements or output might change. (multiple_search_engine, Mascot support)

<center>
  <table>
      <tr>
          <th ALIGN = "center"> pot. predecessor tools </td>
          <td VALIGN="middle" ROWSPAN=2> &rarr; PSMFeatureExtractor &rarr;</td>
          <th ALIGN = "center"> pot. successor tools </td>
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
@verbinclude TOPP_PSMFeatureExtractor.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_PSMFeatureExtractor.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class PSMFeatureExtractor :
  public TOPPBase
{
public:
  PSMFeatureExtractor() :
    TOPPBase("PSMFeatureExtractor", "Computes extra features for each input PSM.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<files>", StringList(), "Input file(s)", true);
    setValidFormats_("in", ListUtils::create<String>("idXML,mzid"));
    registerOutputFile_("out", "<file>", "", "Output file in mzid or idXML format", true);
    setValidFormats_("out", ListUtils::create<String>("idXML,mzid"));    
    registerStringOption_("out_type", "<type>", "", "Output file type -- default: determined from file extension or content.", false);
    setValidStrings_("out_type", ListUtils::create<String>("idXML,mzid"));
    registerStringList_("extra", "<MetaData parameter>", vector<String>(), "List of the MetaData parameters to be included in a feature set for precolator.", false, false);
    registerFlag_("multiple_search_engines", "Combine PSMs from different search engines by merging on scan level.");
    registerFlag_("skip_db_check", "Manual override to skip the check if same settings for multiple search engines were applied. Only valid together with -multiple_search_engines flag.", true);
    registerFlag_("concat", "Naive merging of PSMs from different search engines: concatenate multiple search results instead of merging on scan level. Only valid together with -multiple_search_engines flag.", true);
    registerFlag_("impute", "Will instead of discarding all PSM not unanimously detected by all SE, impute missing values by their respective scores min/max observed. Only valid together with -multiple_search_engines flag.", true);
    registerFlag_("limit_imputation", "Will impute missing scores with the worst numerical limit (instead of min/max observed) of the respective score. Only valid together with -multiple_search_engines flag.", true);
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
    bool multiple_search_engines = getFlag_("multiple_search_engines");
    OPENMS_LOG_DEBUG << "Input file (of target?): " << ListUtils::concatenate(in_list, ",") << endl;
    if (in_list.size() > 1 && !multiple_search_engines)
    {
      writeLogError_("Error: multiple input files given for -in, but -multiple_search_engines flag not specified. If the same search engine was used, feed the input files into PSMFeatureExtractor one by one.");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }
    
    const String out(getStringOption_("out"));

    //-------------------------------------------------------------
    // read input
    //-------------------------------------------------------------
    bool skip_db_check = getFlag_("skip_db_check");
    bool concatenate = getFlag_("concat");
    StringList search_engines_used;
    for (StringList::const_iterator fit = in_list.begin(); fit != in_list.end(); ++fit)
    {
      vector<PeptideIdentification> peptide_ids;
      vector<ProteinIdentification> protein_ids;
      String in = *fit;
      FileHandler fh;
      FileTypes::Type in_type = fh.getType(in);
      OPENMS_LOG_INFO << "Loading input file: " << in << endl;
      if (in_type == FileTypes::IDXML || in_type == FileTypes::MZIDENTML)
      {
        FileHandler().loadIdentifications(in, protein_ids, peptide_ids, {FileTypes::IDXML, FileTypes::MZIDENTML});
      }
      if (in_type == FileTypes::MZIDENTML)
      {
        OPENMS_LOG_WARN << "Converting from mzid: possible loss of information depending on target format." << endl;
      }
      //else caught by TOPPBase:registerInput being mandatory mzid or idxml

      //check and warn if merged from multiple runs
      if (protein_ids.size() > 1)
      {
        throw Exception::InvalidValue(
            __FILE__,
            __LINE__,
            OPENMS_PRETTY_FUNCTION,
            "File '" + in + "' has more than one ProteinIDRun. This is currently not correctly handled."
            "Please use the merge_proteins_add_psms option if you used IDMerger. Alternatively, pass"
            " all original single-run idXML inputs as list to this tool.",
            "# runs: " + String(protein_ids.size()));
      }

      // will check if all ProteinIdentifications have the same search db unless it is the first, in which case all_protein_ids is empty yet.
      if (multiple_search_engines && !skip_db_check && !all_protein_ids.empty())
      {
        ProteinIdentification::SearchParameters all_search_parameters = all_protein_ids.front().getSearchParameters();
        ProteinIdentification::SearchParameters search_parameters = protein_ids.front().getSearchParameters();
        if (search_parameters.db != all_search_parameters.db)
        {
          writeLogError_("Error: Input files are not searched with the same protein database, " + search_parameters.db + " vs. " + all_search_parameters.db + ". Set -skip_db_check flag to ignore this. Aborting!");
          return INCOMPATIBLE_INPUT_DATA;
        }
      }
      
      if (!multiple_search_engines)
      {
        all_peptide_ids.insert(all_peptide_ids.end(), peptide_ids.begin(), peptide_ids.end());
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
      PercolatorFeatureSetHelper::mergeMULTISEProteinIds(all_protein_ids, protein_ids);
    }
    
    if (all_protein_ids.empty())
    {
      writeLogError_("Error: No protein hits found in input file. Aborting!");
      printUsage_();
      return INPUT_FILE_EMPTY;
    }

    //-------------------------------------------------------------
    // extract search engine and prepare pin
    //-------------------------------------------------------------
    String search_engine = all_protein_ids.front().getSearchEngine();
    if (multiple_search_engines) search_engine = "multiple";
    OPENMS_LOG_DEBUG << "Registered search engine: " << search_engine << endl;
    
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
    else if (search_engine == "MS-GF+") 
    {
      PercolatorFeatureSetHelper::addMSGFFeatures(all_peptide_ids, feature_set);
    }
    else if (search_engine == "Mascot") 
    {
      PercolatorFeatureSetHelper::addMASCOTFeatures(all_peptide_ids, feature_set);
    }
    else if (search_engine == "XTandem") 
    {
      PercolatorFeatureSetHelper::addXTANDEMFeatures(all_peptide_ids, feature_set);
    }
    else if (search_engine == "Comet") 
    {
      PercolatorFeatureSetHelper::addCOMETFeatures(all_peptide_ids, feature_set);
    }
    else if (search_engine == "MSFragger") 
    {
      PercolatorFeatureSetHelper::addMSFRAGGERFeatures(feature_set);
    }
    else
    {
      OPENMS_LOG_ERROR << "No known input to create PSM features from. Aborting" << std::endl;
      return INCOMPATIBLE_INPUT_DATA;
    }

    String run_identifier = all_protein_ids.front().getIdentifier();
    for (vector<PeptideIdentification>::iterator it = all_peptide_ids.begin(); it != all_peptide_ids.end(); ++it)
    {
      it->setIdentifier(run_identifier);
      PercolatorFeatureSetHelper::checkExtraFeatures(it->getHits(), extra_features);  // will remove inconsistently available features
    }
    
    if (all_protein_ids.size() > 1)
    {
      OPENMS_LOG_ERROR << "Multiple identifications in one file are not supported. Please resume with separate input files. Quitting." << std::endl;
      return INCOMPATIBLE_INPUT_DATA;
    }
    else
    {
      ProteinIdentification::SearchParameters search_parameters = all_protein_ids.front().getSearchParameters();
      
      search_parameters.setMetaValue("feature_extractor", "TOPP_PSMFeatureExtractor");
      feature_set.insert(feature_set.end(), extra_features.begin(), extra_features.end());
      search_parameters.setMetaValue("extra_features", ListUtils::concatenate(feature_set, ","));
      all_protein_ids.front().setSearchParameters(search_parameters);
    }
    
    // Storing the PeptideHits with calculated q-value, pep and svm score
    FileTypes::Type out_type = FileTypes::nameToType(getStringOption_("out_type"));

    if (out_type == FileTypes::UNKNOWN)
    {
      FileHandler fh;
      out_type = fh.getTypeByFileName(out);
    }

    if (out_type == FileTypes::UNKNOWN)
    {
      writeLogError_("Error: Could not determine output file type! Set 'out_type' parameter to desired file type.");
      return PARSE_ERROR;
    }
    OPENMS_LOG_INFO << "writing output file: " << out << endl;
    

    FileHandler().storeIdentifications(out, all_protein_ids, all_peptide_ids, {FileTypes::MZIDENTML, FileTypes::IDXML});


    writeLogInfo_("PSMFeatureExtractor finished successfully!");
    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  PSMFeatureExtractor tool;

  return tool.main(argc, argv);
}

/// @endcond
