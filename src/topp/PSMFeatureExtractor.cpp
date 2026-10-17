// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Andreas Simon, Mathias Walzer, Matthew The $
// --------------------------------------------------------------------------
#include <OpenMS/config.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/ANALYSIS/ID/PercolatorFeatureSetHelper.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
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

@experimental Parts of this tool are still work in progress and usage and input requirements or output might change. (Mascot support)

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
Currently supported search engines are Comet, X!Tandem, MSGF+, MSFragger, andes.
Mascot support is available but in beta development.

The tool processes exactly one search result at a time. To rescore results from
several search engines jointly, combine them with @ref TOPP_ConsensusID instead.
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
    registerInputFile_("in", "<file>", "", "Input file", true);
    setValidFormats_("in", ListUtils::create<std::string>("idXML,mzid,idparquet"));
    registerOutputFile_("out", "<file>", "", "Output file in idXML, mzid or idparquet format", true);
    setValidFormats_("out", ListUtils::create<std::string>("idXML,mzid,idparquet"));
    registerStringOption_("out_type", "<type>", "", "Output file type -- default: determined from file extension or content.", false);
    setValidStrings_("out_type", ListUtils::create<std::string>("idXML,mzid,idparquet"));
    registerStringList_("extra", "<MetaData parameter>", vector<std::string>(), "List of the MetaData parameters to be included in a feature set for precolator.", false, false);
  }
  
  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // general variables and data to perform PSMFeatureExtractor
    //-------------------------------------------------------------
    PeptideIdentificationList all_peptide_ids;
    vector<ProteinIdentification> all_protein_ids;
    
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    const std::string in(getStringOption_("in"));
    const std::string out(getStringOption_("out"));

    //-------------------------------------------------------------
    // read input
    //-------------------------------------------------------------
    {
      PeptideIdentificationList peptide_ids;
      vector<ProteinIdentification> protein_ids;
      FileHandler fh;
      FileTypes::Type in_type = fh.getType(in);
      OPENMS_LOG_INFO << "Loading input file: " << in << endl;
      if (in_type == FileTypes::IDXML || in_type == FileTypes::MZIDENTML || in_type == FileTypes::IDPARQUET)
      {
        FileHandler().loadIdentifications(in, protein_ids, peptide_ids, {FileTypes::IDXML, FileTypes::MZIDENTML, FileTypes::IDPARQUET});
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
            "Please use the merge_proteins_add_psms option if you used IDMerger. Alternatively, feed the"
            " original single-run idXML inputs into this tool one by one.",
            "# runs: " + StringUtils::toStr(protein_ids.size()));
      }

      all_peptide_ids.insert(all_peptide_ids.end(), peptide_ids.begin(), peptide_ids.end());
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
    const std::string search_engine = all_protein_ids.front().getSearchEngine();
    OPENMS_LOG_DEBUG << "Registered search engine: " << search_engine << endl;

    StringList extra_features = getStringList_("extra");
    StringList feature_set;

    if (search_engine == "MS-GF+")
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
    else if (search_engine == "andes")
    {
      PercolatorFeatureSetHelper::addANDESFeatures(all_peptide_ids, feature_set);
    }
    else
    {
      OPENMS_LOG_ERROR << "No known input to create PSM features from. Aborting" << std::endl;
      return INCOMPATIBLE_INPUT_DATA;
    }

    std::string run_identifier = all_protein_ids.front().getIdentifier();
    for (PeptideIdentificationList::iterator it = all_peptide_ids.begin(); it != all_peptide_ids.end(); ++it)
    {
      it->setIdentifier(run_identifier);
      PercolatorFeatureSetHelper::checkExtraFeatures(it->getHits(), extra_features);  // will remove inconsistently available features
    }
    
    ProteinIdentification::SearchParameters search_parameters = all_protein_ids.front().getSearchParameters();

    search_parameters.setMetaValue("feature_extractor", "TOPP_PSMFeatureExtractor");
    feature_set.insert(feature_set.end(), extra_features.begin(), extra_features.end());
    search_parameters.setMetaValue("extra_features", ListUtils::concatenate(feature_set, ","));
    all_protein_ids.front().setSearchParameters(search_parameters);


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
    

    FileHandler().storeIdentifications(out, all_protein_ids, all_peptide_ids, {out_type});


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
