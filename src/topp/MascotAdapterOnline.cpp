// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch, Daniel Jameson, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/PercolatorFeatureSetHelper.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/FORMAT/MascotXMLFile.h>
#include <OpenMS/FORMAT/MascotRemoteQuery.h>
#include <OpenMS/FORMAT/MascotGenericFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/SYSTEM/File.h>

#include <OpenMS/APPLICATIONS/SearchEngineBase.h>

#include <sstream>

#include <QtCore/QFile>
#include <QtCore/QCoreApplication>
#include <QtCore/QTimer>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**

@page TOPP_MascotAdapterOnline MascotAdapterOnline

@brief Identifies peptides in MS/MS spectra via Mascot.

<CENTER>
    <table>
        <tr>
            <th ALIGN = "center"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> &rarr; MascotAdapterOnline &rarr;</td>
            <th ALIGN = "center"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any signal-/preprocessing tool @n (that writes mzML format)</td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter or @n any protein/peptide processing tool</td>
        </tr>
    </table>
</CENTER>

This wrapper application generates peptide identifications for MS/MS
spectra using the search engine Mascot. It communicates with the Mascot
server over the network (i.e. it does not have to run on the server
itself).

The adapter supports Mascot security features as well as proxy connections.
Mascot versions 2.2.x up to 2.4.1 are supported and have been successfully
tested (to varying degrees).

@bug Running the adapter on Mascot 2.4 (possibly also other versions) produces the following error messages, which should be ignored:\n
MascotRemoteQuery: An error occurred (requestId=11): Request aborted (QT Error Code: 7)\n
MascotRemoteQuery: An error occurred (requestId=12): Request aborted (QT Error Code: 7)

@note Some Mascot server instances seem to fail without reporting back an
error message. In such cases, try to run the search on another Mascot
server or change/validate the search parameters (e.g. using modifications
that are known to Mascot and can thus be set in the INI file, but which are
unknown to Mascot, might pose a problem).

@note Mascot returns incomplete/incorrect protein assignments for most
identified peptides (due to protein-level grouping/filtering). Thus,
the protein associations are therefore not included in the output of this
adapter, only the peptide sequences. @ref TOPP_PeptideIndexer should be run
after this tool to get correct assignments.

@note Currently mzIdentML (mzid) is not directly supported as an
input/output format of this tool. Convert mzid files to/from idXML using
@ref TOPP_IDFileConverter if necessary.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_MascotAdapterOnline.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_MascotAdapterOnline.html

For the parameters of the algorithm section see the algorithms documentation: @n
@ref OpenMS::MascotRemoteQuery "Mascot_server" @n
@ref OpenMS::MascotGenericFile "Mascot_parameters" @n

*/

// We do not want this class to show up in the doc:
/// @cond TOPPCLASSES


class TOPPMascotAdapterOnline :
  public SearchEngineBase
{
public:
  TOPPMascotAdapterOnline() :
    SearchEngineBase("MascotAdapterOnline", "Annotates MS/MS spectra using Mascot.")
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input file in mzML format.\n");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "output file in idXML format.\n");
    setValidFormats_("out", ListUtils::create<String>("idXML"));

    registerSubsection_("Mascot_server", "Mascot server details");
    registerSubsection_("Mascot_parameters", "Mascot parameters used for searching");
  }

  Param getSubsectionDefaults_(const String& section) const override
  {
    if (section == "Mascot_server")
    {
      MascotRemoteQuery mascot_query;
      return mascot_query.getParameters();
    }

    if (section == "Mascot_parameters")
    {
      MascotGenericFile mgf_file;
      Param p = mgf_file.getParameters();
      p.remove("internal:");
      return p;
    }

    return Param();
  }

  void parseMascotResponse_(const PeakMap& exp, bool decoy, MascotRemoteQuery* mascot_query, ProteinIdentification& prot_id, vector<PeptideIdentification>& pep_ids)
  {   
    String mascot_tmp_file_name = decoy ? (File::getTempDirectory() + "/" + File::getUniqueName() + "_Mascot_decoy_response") : (File::getTempDirectory() + "/" + File::getUniqueName() + "_Mascot_response");
    QFile mascot_tmp_file(mascot_tmp_file_name.c_str());
    mascot_tmp_file.open(QIODevice::WriteOnly);
    if (decoy)
    {
      mascot_tmp_file.write(mascot_query->getMascotXMLDecoyResponse());
    }
    else
    {
      mascot_tmp_file.write(mascot_query->getMascotXMLResponse());
    }
    mascot_tmp_file.close();
  
    writeDebug_(String("\nMascot Server Response file saved to: '") + mascot_tmp_file_name + "'. If an error occurs, send this file to the OpenMS team.\n", 100);

    // set up helper object for looking up spectrum meta data:
    SpectrumMetaDataLookup lookup;
    MascotXMLFile::initializeLookup(lookup, exp);

    // read the response
    MascotXMLFile().load(mascot_tmp_file_name, prot_id, pep_ids, lookup);
    writeDebug_("Read " + String(pep_ids.size()) + " peptide ids and " + String(prot_id.getHits().size()) + " protein identifications from Mascot", 5);

    // for debugging errors relating to unexpected response files
    if (this->debug_level_ >= 100)
    {
      writeDebug_(String("\nMascot Server Response file saved to: '") + mascot_tmp_file_name + "'. If an error occurs, send this file to the OpenMS team.\n", 100);
    }
    else
    {
      mascot_tmp_file.remove(); // delete file
    }
  }

  // merge b into a
  void mergeIDs_(ProteinIdentification& p_a, const ProteinIdentification& p_b, vector<PeptideIdentification>& pep_a, const vector<PeptideIdentification>& pep_b)
  {
    // if p_a is empty use all meta values and hits from p_b to initialize p_a
    if (p_a.getHits().empty())
    {
      p_a = p_b;
    }
    else
    {
      // p_a already initialized? just add proteins of b to a
      for (const ProteinHit& p : p_b.getHits())
      {
        p_a.insertHit(p);
      }
    }
    
    map<String, size_t> native_id2id_index;
    size_t index{};
    String run_identifier;
    for (const PeptideIdentification& pep : pep_a)
    {
      const String& native_id = pep.getSpectrumReference();
      native_id2id_index[native_id] = index;
      ++index;
      if (run_identifier.empty()) run_identifier = pep.getIdentifier();
    }

    for (auto pep : pep_b) //OMS_CODING_TEST_EXCLUDE
    {
      auto it = native_id2id_index.find(pep.getSpectrumReference());
      if (it == native_id2id_index.end()) // spectrum not yet identified? add decoy id
      {
        pep.setIdentifier(run_identifier);
        pep_a.push_back(pep);
      }
      else
      { 
        const vector<PeptideHit>& hits = pep.getHits();
        if (hits.empty())
        {
          continue;
        }
        for (const PeptideHit& h : hits)
        {
          pep_a[it->second].insertHit(h);
        }
        pep_a[it->second].assignRanks();
      }
    }
  }
   
  ExitCodes main_(int argc, const char** argv) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    // input/output files
    String in = getRawfileName();
    String out(getStringOption_("out"));

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------

    PeakMap exp;
    // keep only MS2 spectra
    FileHandler fh;
    fh.getOptions().setMSLevels({2});
    fh.loadExperiment(in, exp, {FileTypes::Type::MZML}, log_type_, false, false);
    writeLogInfo_("Number of spectra loaded: " + String(exp.size()));


    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    Param mascot_param = getParam_().copy("Mascot_parameters:", true);

    // overwrite default search title with filename
    if (mascot_param.getValue("search_title") == "OpenMS_search")
    {
      mascot_param.setValue("search_title", FileHandler::stripExtension(File::basename(in)));
    }

    mascot_param.setValue("internal:HTTP_format", "true");

    SpectrumLookup lookup;        
    lookup.readSpectra(exp.getSpectra());

    Param mascot_query_param = getParam_().copy("Mascot_server:", true);
    size_t batch_size = (size_t)mascot_query_param.getValue("batch_size");
    size_t chunks = (exp.size() - 1) / batch_size + 1; // Note: safe as we have at least one spectrum

    vector<ProteinIdentification> all_prot_ids;
    ProteinIdentification all_prot_id;
    vector<PeptideIdentification> all_pep_ids;

    MSExperiment current_batch;
    for (size_t k = 0; k < chunks; ++k)
    {
      // get range for next set of n elements
      auto start_itr = std::next(exp.begin(), k*batch_size);
      auto end_itr = std::next(exp.begin(), k*batch_size + batch_size);

      // allocate memory for the current chunk
      current_batch.resize(batch_size);

      // code to handle the last sub-vector as it might
      // contain less elements
      if (k*batch_size + batch_size > exp.size()) 
      {
        end_itr = exp.end();
        current_batch.resize(exp.size() - k*batch_size);
      }

      // copy elements from the input range to the sub-vector        
      std::copy(start_itr, end_itr, current_batch.begin());

      // write mgf and run search
      MascotGenericFile mgf_file;
      mgf_file.setParameters(mascot_param);
      // get the spectra into string stream
      writeDebug_("Writing MGF file to stream", 1);
      stringstream ss;
      mgf_file.store(ss, in, current_batch, true); // write in compact format

      // Usage of a QCoreApplication is overkill here (and ugly too), but we just use the
      // QEventLoop to process the signals and slots and grab the results afterwards from
      // the MascotRemotQuery instance
      char** argv2 = const_cast<char**>(argv);
      QCoreApplication event_loop(argc, argv2);
      MascotRemoteQuery* mascot_query = new MascotRemoteQuery(&event_loop);
      writeDebug_("Setting parameters for Mascot query", 1);
      mascot_query->setParameters(mascot_query_param);
      
      bool internal_decoys = mascot_param.getValue("decoy") == "true";
      // We used internal decoy search. Set that we want to retrieve decoy search results during export.
      if (internal_decoys)
      {
        mascot_query->setExportDecoys(true);
      }

      writeDebug_("Setting spectra for Mascot query", 1);
      mascot_query->setQuerySpectra(ss.str());

      // remove unnecessary spectra
      ss.clear();

      QObject::connect(mascot_query, SIGNAL(done()), &event_loop, SLOT(quit()));
      QTimer::singleShot(1000, mascot_query, SLOT(run()));
      writeLogInfo_("Submitting Mascot query (now: " + DateTime::now().get() + ")...");
      event_loop.exec();
      writeLogInfo_("Mascot query finished");

      if (mascot_query->hasError())
      {
        writeLogError_("An error occurred during the query: " + mascot_query->getErrorMessage());
        delete mascot_query;
        return EXTERNAL_PROGRAM_ERROR;
      }

      vector<PeptideIdentification> pep_ids;
      ProteinIdentification prot_id;

      if (!mascot_query_param.exists("skip_export") ||
          !mascot_query_param.getValue("skip_export").toBool())
      {
        // write Mascot response to file
        parseMascotResponse_(current_batch, false, mascot_query, prot_id, pep_ids); // targets

        // reannotate proper spectrum native id if missing
        for (auto& pep : pep_ids)
        {
          // no need to reannotate
          if (pep.metaValueExists("spectrum_reference") 
            && !(static_cast<String>(pep.getMetaValue("spectrum_reference")).empty()))
          {
            continue;
          }
            
          try
          { 
            Size index = lookup.findByRT(pep.getRT());
            pep.setSpectrumReference( exp[index].getNativeID());
          }
          catch (Exception::ElementNotFound&)
          {
            OPENMS_LOG_ERROR << "Error: Failed to look up spectrum native ID for peptide identification with retention time '" + String(pep.getRT()) + "'." << endl;
          }
        }

        if (internal_decoys)
        {
          vector<PeptideIdentification> decoy_pep_ids;
          ProteinIdentification decoy_prot_id;
          parseMascotResponse_(current_batch, true, mascot_query, decoy_prot_id, decoy_pep_ids);  // decoys

          // reannotate proper spectrum native id if missing
          for (auto& pep : decoy_pep_ids)
          {
            // no need to reannotate
            if (pep.metaValueExists("spectrum_reference") 
              && !(static_cast<String>(pep.getMetaValue("spectrum_reference")).empty()))
            {
              continue;
            } 

            try
            { 
              Size index = lookup.findByRT(pep.getRT());
              pep.setSpectrumReference( exp[index].getNativeID());
            }
            catch (Exception::ElementNotFound&)
            {
              OPENMS_LOG_ERROR << "Error: Failed to look up spectrum native ID for peptide identification with retention time '" + String(pep.getRT()) + "'." << endl;
            }
          }
          mergeIDs_(prot_id, decoy_prot_id, pep_ids, decoy_pep_ids);
        }
      }

      String search_number = mascot_query->getSearchIdentifier();
      if (search_number.empty())
      {
        writeLogError_("Error: Failed to extract the Mascot search identifier (search number).");
        if (mascot_query_param.exists("skip_export") &&
            mascot_query_param.getValue("skip_export").toBool())
        {
          return PARSE_ERROR;
        }
      }
      else 
      {
        prot_id.setMetaValue("SearchNumber", search_number);
      }

      // clean up
      delete mascot_query;
      
      current_batch.clear(true); // clear meta data

      mergeIDs_(all_prot_id, prot_id, all_pep_ids, pep_ids);
    }

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    all_prot_id.setPrimaryMSRunPath({ in }, exp);

    DateTime now = DateTime::now();
    String date_string = now.get();
    String run_identifier("Mascot_" + date_string);

    // remove proteins as protein links seem are broken and reindexing is needed
    all_prot_id.getHits().clear(); 
    all_prot_id.setIdentifier(run_identifier);
    all_prot_ids.push_back(all_prot_id);

    // remove protein links from peptides as protein links seem are broken and reindexing is needed
    for (auto& pep : all_pep_ids)
    {
      pep.setIdentifier(run_identifier);
      for (auto& hit : pep.getHits())
      {
        hit.setPeptideEvidences({});
      }
    }        

    // write all (!) parameters as metavalues to the search parameters
    DefaultParamHandler::writeParametersToMetaValues(this->getParam_(), all_prot_ids[0].getSearchParameters(), this->getToolPrefix());

    // get feature set used in percolator
    StringList feature_set;
    PercolatorFeatureSetHelper::addMASCOTFeatures(all_pep_ids, feature_set);
    all_prot_ids.front().getSearchParameters().setMetaValue("extra_features", ListUtils::concatenate(feature_set, ","));
    
    FileHandler().storeIdentifications(out, all_prot_ids, all_pep_ids, {FileTypes::IDXML});
    
    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPMascotAdapterOnline tool;

  return tool.main(argc, argv);
}

/// @endcond
