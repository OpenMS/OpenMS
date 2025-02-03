// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/SearchEngineBase.h>

#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/SYSTEM/File.h>

using namespace std;

namespace OpenMS
{
  SearchEngineBase::SearchEngineBase(const String& tool_name, const String& tool_description, bool official, const std::vector<Citation>& citations, bool toolhandler_test) :
    TOPPBase(tool_name, tool_description, official, citations, toolhandler_test)
  {
  }

  SearchEngineBase::~SearchEngineBase() = default;

  
  String SearchEngineBase::getRawfileName(int ms_level) const
  {
    String inputfile_name = getStringOption_("in");
    FileHandler fh;
    auto type = fh.getType(inputfile_name);
    switch (type)
    {
      case FileTypes::MZML:
      {
        MzMLFile mzml;
        mzml.getOptions().setMSLevels({ ms_level }); // only query MS2 (or whatever ms_level is)
        const auto& centroid_info = mzml.getCentroidInfo(inputfile_name);
        const auto& lvl_info = centroid_info.find(ms_level);
        if (lvl_info == centroid_info.end())
        {
          throw Exception::FileEmpty(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Error: No MS" + String(ms_level) + " spectra in input file.");
        }

        if (lvl_info->second.count_profile > 0)
        {
          if (getFlag_("force"))
          {
            OPENMS_LOG_WARN << "Warning: Profile data found, but centroid MS spectra required. "
                               "Since '-force' flag is in effect, we will continue, but results are likely bogus." << std::endl;
          }
          else
          {
            throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "Error: Profile data provided but centroided MS" + String(ms_level) + " spectra required. To enforce processing (unwise!) of the data enable the -force flag (results will be bogus!).");
          }
        }
        if (lvl_info->second.count_centroided == 0)
        {
          if (getFlag_("force"))
          {
            OPENMS_LOG_WARN << "Warning: No centroided MS" + String(ms_level) + " were found, but are required. "
                               "Since '-force' flag is in effect, we will continue, but results might be bogus." << std::endl;
          }
          else
          {
            throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
              "Error: No centroided MS" + String(ms_level) + " spectra were found, but are required. To enforce processing of the data enable the -force flag (results will likely be bogus!).");
          }
        }
        // do no check for UNKNOWN, since it does not really tell much (UNKNOWN can only occur if meta data is missing and our peak type estimation fails (which only happens for (almost) empty spectra))
      }
      case FileTypes::MGF:
        // no warning required. MGF files should be centroided by definition
        break;
      default:
        OPENMS_LOG_WARN << "Warning: make sure that MS" << ms_level << " spectra in '" << inputfile_name << "' are centroided. Otherwise the results may be undefined!";
    }

      
    return inputfile_name;
  }

  String SearchEngineBase::getDBFilename(const String& db) const
  {
    String db_name(db.empty() ? getStringOption_("database") : db);
    if (!File::readable(db_name))
    {
      db_name = File::findDatabase(db_name);
    }
    return db_name;
  }

  SearchEngineBase::ExitCodes SearchEngineBase::reindex_(std::vector<ProteinIdentification>& protein_identifications, 
                                       std::vector<PeptideIdentification>& peptide_identifications) const
  {
    if (getStringOption_("reindex") == "true")
    {
      PeptideIndexing indexer;
      
      // extract parameter subtree
      Param param = getParam_().copy("PeptideIndexing:", true);
      
      Param param_pi = indexer.getParameters();
      // copy search engine specific default parameter for peptide indexing into param_pi
      param_pi.update(param, false, false, false, false, OpenMS_Log_debug); // suppress param. update message
      indexer.setParameters(param_pi);
      indexer.setLogType(this->log_type_);
      FASTAContainer<TFI_File> proteins(getDBFilename());
      PeptideIndexing::ExitCodes indexer_exit = indexer.run(proteins, protein_identifications, peptide_identifications);

      if ((indexer_exit != PeptideIndexing::EXECUTION_OK) &&
          (indexer_exit != PeptideIndexing::PEPTIDE_IDS_EMPTY))
      {
        if (indexer_exit == PeptideIndexing::DATABASE_EMPTY)
        {
          return INPUT_FILE_EMPTY;       
        }
        else if (indexer_exit == PeptideIndexing::UNEXPECTED_RESULT)
        {
          return UNEXPECTED_RESULT;
        }
        else
        {
          return UNKNOWN_ERROR;
        }
      } 
    }
    return EXECUTION_OK;
  }

  void SearchEngineBase::registerPeptideIndexingParameter_(Param peptide_indexing_parameter)
  {
    registerStringOption_("reindex", "<choice>", "true", "Recalculate peptide to protein association using OpenMS. Annotates target-decoy information.", false);
    setValidStrings_("reindex", { "true", "false" });
  
    peptide_indexing_parameter.setValue("missing_decoy_action", "warn");

    // hide entries
    for (const auto& s : {"decoy_string", "decoy_string_position", "missing_decoy_action", "enzyme:name", "enzyme:specificity",
                          "write_protein_sequence", "write_protein_description", "keep_unreferenced_proteins", "unmatched_action", 
                          "aaa_max","mismatches_max", "IL_equivalent"})
    {
      peptide_indexing_parameter.addTag(s, "advanced");
    }
    // move parameter to PeptideIndexing subtree so we don't accidently overwrite duplicate keys in the tool and indexer
    Param combined;
    combined.insert("PeptideIndexing:", peptide_indexing_parameter);
    registerFullParam_(combined);
  }

}
