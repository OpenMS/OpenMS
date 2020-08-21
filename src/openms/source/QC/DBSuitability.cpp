// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Tom Waschischeck $
// $Authors: Tom Waschischeck $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMDecoy.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FILTERING/ID/IDFilter.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/QC/DBSuitability.h>
#include <OpenMS/SYSTEM/ExternalProcess.h>
#include <OpenMS/SYSTEM/File.h>

using namespace std;

namespace OpenMS
{
  DBSuitability::DBSuitability()
    : DefaultParamHandler("DBSuitability"), results_{}
  {
    defaults_.setValue("no_re_rank", "false", "Enable/Disable re-ranking");
    defaults_.setValidStrings("no_re_rank", { "true", "false" });
    defaults_.setValue("cut_off_fract", 1., "Percentile to determine which decoy cut-off to use.");
    defaults_.setMinFloat("cut_off_fract", 0.);
    defaults_.setMaxFloat("cut_off_fract", 1.);
    defaults_.setValue("FDR", 0.01, "Filtering peptide hits based on this q-value");
    defaults_.setMinFloat("FDR", 0.);
    defaults_.setMaxFloat("FDR", 1.);
    defaultsToParam_();
  }
  
  void DBSuitability::compute(vector<PeptideIdentification> pep_ids, const MSExperiment& exp, const vector<FASTAFile::FASTAEntry>& original_fasta, const vector<FASTAFile::FASTAEntry>& novo_fasta, const ProteinIdentification::SearchParameters& search_params)
  {
    bool no_re_rank = param_.getValue("no_re_rank").toBool();
    double cut_off_fract = param_.getValue("cut_off_fract");
    double FDR = param_.getValue("FDR");

    SuitabilityData d;
    results_.push_back(d);
    SuitabilityData& data = results_.back();

    pair<String, Param> search_info = extractSearchAdapterInfoFromMetaValues_(search_params);

    // execute combined deNovo+database ID-search
    vector<FASTAFile::FASTAEntry> combined_db(original_fasta);
    combined_db.insert(combined_db.end(), novo_fasta.begin(), novo_fasta.end());
    vector<PeptideIdentification> combined_ids = runIdentificationSearch_(exp, combined_db, search_info.first, search_info.second);

    if (combined_ids.empty())
    {
      OPENMS_LOG_WARN << "No peptide identifications found with given mzML and database! No calculations performed." << endl;
      return;
    }

    if (!no_re_rank)
    {
      data.cut_off = getDecoyCutOff_(combined_ids, cut_off_fract);
    }

    for (PeptideIdentification& pep_id : combined_ids)
    {
      // sort hits by q-value
      pep_id.sort();

      vector<PeptideHit>& hits = pep_id.getHits();

      if (hits.empty()) continue;

      const PeptideHit& top_hit = hits[0];

      // skip if the top hit is a decoy hit
      if (!top_hit.metaValueExists("target_decoy"))
      {
        throw(Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No target/decoy information found! Make sure 'PeptideIndexer' is run beforehand."));
      }
      if (top_hit.getMetaValue("target_decoy") == "decoy") continue;

      // skip if top hit is out ouf FDR
      if (!passesFDR_(top_hit, FDR)) continue;

      // check if top hit is found in de novo protein
      if (!isNovoHit_(top_hit)) // top hit is db hit
      {
        ++data.num_top_db;
        continue;
      }

      // find the second target hit, skip all decoy or novo hits inbetween
      const PeptideHit* second_hit = nullptr;
      for (UInt i = 1; i < hits.size(); ++i)
      {
        // check for FDR
        if (!passesFDR_(hits[i], FDR)) break;

        // check if target, also check for "target+decoy" value
        String td_info(hits[i].getMetaValue("target_decoy"));
        if (td_info.find("target") != 0) continue;

        // check if hit is novo hit
        if (isNovoHit_(hits[i])) continue;

        second_hit = &hits[i];
        break;
      }
      if (second_hit == nullptr) // no second target hit with given FDR found
      {
        ++data.num_top_novo;
        continue;
      }

      // second hit is db hit
      ++data.num_interest;

      // check for re-ranking
      if (no_re_rank)
      {
        ++data.num_top_novo;
        continue;
      }

      // check for xcorr score
      if (!top_hit.metaValueExists("MS:1002252") || !second_hit->metaValueExists("MS:1002252"))
      {
        throw(Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No cross correlation score found at peptide hit. Only Comet search engine is supported right now."));
      }

      double top_xscore_mw = double(top_hit.getMetaValue("MS:1002252")) / top_hit.getSequence().getMonoWeight();
      double second_xscore_mw = double(second_hit->getMetaValue("MS:1002252")) / second_hit->getSequence().getMonoWeight();
      if (top_xscore_mw - second_xscore_mw <= data.cut_off)
      {
        ++data.num_top_db;
        ++data.num_re_ranked;
      }
      else
      {
        ++data.num_top_novo;
      }
    }

    if (data.num_top_db == 0 && data.num_top_novo == 0)
    {
      OPENMS_LOG_WARN << "Identifications could not be assigned to either the database or the deNovo protein. Probably your FDR threshold is too strict." << endl;
      data.suitability = DBL_MAX;
      return;
    }

    data.suitability = double(data.num_top_db) / (data.num_top_db + data.num_top_novo);

    // calculate correction of suitability
    if (pep_ids[0].getScoreType() == "q-value")
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "q-value found at PeptideIdentifications. That is not allowed! Please make sure FDR did not run previously.");
    }
    for (const auto& id : pep_ids)
    {
      if (id.getHits().empty()) continue;
      if (id.getHits()[0].metaValueExists("q-value"))
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "q-value found at PeptideIdentifications. That is not allowed! Please make sure FDR did not run previously.");
      }
    }

    Param p;
    p.setValue("use_all_hits", "true");
    p.setValue("add_decoy_peptides", "true");
    p.setValue("add_decoy_proteins", "true");

    FalseDiscoveryRate fdr;
    fdr.setParameters(p);
    fdr.apply(pep_ids);
    
    Size num_ids_only_db = countNumberOfIdentifications_(pep_ids);
    Size num_ids_only_novo = countNumberOfIdentifications_(runIdentificationSearch_(exp, novo_fasta, search_info.first, search_info.second));

    double factor = num_ids_only_novo / double(num_ids_only_db);
    data.corr_factor = factor;

    data.num_top_novo_corr = data.num_top_novo * factor;
    data.suitability_corr = double(data.num_top_db) / (data.num_top_db + (factor * data.num_top_novo));
  }

  const std::vector<DBSuitability::SuitabilityData>& DBSuitability::getResults() const
  {
    return results_;
  }

  double DBSuitability::getDecoyDiff_(const PeptideIdentification& pep_id)
  {
    double diff = DBL_MAX;

    // get the score of the first two decoy hits
    double decoy_1 = DBL_MAX;
    double decoy_2 = DBL_MAX;
    UInt curr_hit = 0;

    for (const auto& hit : pep_id.getHits())
    {
      if (curr_hit > 10) break;
      ++curr_hit;

      if (!hit.metaValueExists("target_decoy"))
      {
        throw(Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No target/decoy information found! Make sure 'PeptideIndexer' is run beforehand."));
      }

      if (!hit.metaValueExists("MS:1002252"))
      {
        throw(Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No cross correlation score found at peptide hit. Only Comet search engine is supported right now."));
      }

      if (decoy_1 == DBL_MAX && hit.getMetaValue("target_decoy") == "decoy")
      {
        decoy_1 = hit.getMetaValue("MS:1002252");
        continue;
      }
      if (decoy_1 < DBL_MAX && hit.getMetaValue("target_decoy") == "decoy")
      {
        decoy_2 = hit.getMetaValue("MS:1002252");
        break;
      }
    }

    if (decoy_2 < DBL_MAX) // if there are two decoy hits
    {
      diff = abs(decoy_1 - decoy_2) / pep_id.getHits()[0].getSequence().getMonoWeight(); // normalized by mw
    }

    // if there aren't two decoy hits DBL_MAX is returned
    return diff;
  }

  double DBSuitability::getDecoyCutOff_(const vector<PeptideIdentification>& pep_ids, double cut_off_fract)
  {
    if (cut_off_fract < 0 || cut_off_fract > 1)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "'cut_off_fract' is not within its allowed range [0,1]. Please select a valid value.");
    }

    // get all decoy diffs of peptide ids with at least two decoy hits
    vector<double> diffs;
    for (const auto& pep_id : pep_ids)
    {
      double diff = getDecoyDiff_(pep_id);
      if (diff < DBL_MAX)
      {
        diffs.push_back(diff);
      }
    }

    if (double(diffs.size()) / pep_ids.size() < 0.2)
    {
      throw(Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Under 20 % of peptide identifications have two decoy hits. This is not enough for re-ranking. Use the 'force_no_re_rank' flag to still compute a suitability score."));
    }

    // sort the diffs decreasing and get the (1-novo_fract)*N one
    UInt index = round((1 - cut_off_fract) * diffs.size());
    
    if (index >= diffs.size())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "'cut_off_fract' is set too low. Please set the parameter to a higher value.");
    }

    nth_element(diffs.begin(), diffs.begin() + index, diffs.end(), greater<double>());

    return diffs[index];
  }

  bool DBSuitability::isNovoHit_(const PeptideHit& hit)
  {
    const set<String>& accessions = hit.extractProteinAccessionsSet();
    for (const String& acc : accessions)
    {
      if (acc.find(Constants::UserParam::CONCAT_PEPTIDE) == String::npos)
      {
        return false;
      }
    }
    return true;
  }

  bool DBSuitability::passesFDR_(const PeptideHit& hit, double FDR) const
  {
    if (hit.getScore() > FDR) return false;
    return true;
  }

  pair<String, Param> DBSuitability::extractSearchAdapterInfoFromMetaValues_(const ProteinIdentification::SearchParameters& search_params) const
  {
    Param p;
    // list of all allowed adapters
    vector<String> working_adapters{ "CometAdapter", "CruxAdapter", "MSGFPlusAdapter", "MSFraggerAdapter", "MyriMatchAdapter", "OMSSAAdapter", "XTandemAdapter" };

    String adapter;
    vector<String> keys;
    search_params.getKeys(keys);
    for (const String& key : keys)
    {
      if (adapter.empty()) // trying to find which adapter was used
      {
        for (const String& a : working_adapters)
        {
          if (key.find(a) == 0)
          {
            // used adapter found
            adapter = a;
            break;
          }
        }
      }
      else // found adapter name
      {
        if (key.find(adapter) != 0) continue; // does adapter appear in meta value key?
        p.setValue(key, search_params.getMetaValue(key));
      }
    }

    if (p.empty()) // non of the allowed adapter names where found in the meta values
    {
      String message;
      message = "No parameters found for any of the allowed adapters in the given meta values. Allowed are:\n";
      for (const String& a : working_adapters)
      {
        message += a + " ";
      }
      message = "\n";
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message);
    }

    OPENMS_LOG_DEBUG << "Parameters for the following adapter were found: " << adapter << endl;

    return make_pair(adapter, p);
  }

  void DBSuitability::writeIniFile_(const Param& parameters, const String& filename) const
  {
    if (!File::writable(filename)) // check if file exists and is writable
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unable to create file: " + filename);
    }
    ParamXMLFile param_file;
    param_file.store(filename, parameters);
  }

  vector<PeptideIdentification> DBSuitability::runIdentificationSearch_(const MSExperiment& exp, const vector<FASTAFile::FASTAEntry>& fasta_data, const String& adapter_name, Param& parameters) const
  {
    if (adapter_name.empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No adapter name given. Aborting!");
    }

    // temporary folder for search in- und output files
    File::TempDir tmp_dir;
    String mzml_path = tmp_dir.getPath() + "spectra.mzML";
    String db_path = tmp_dir.getPath() + "database.FASTA";
    String out_path = tmp_dir.getPath() + "out.idXML";

    // override the in- and output files in the parameters
    parameters.setValue(adapter_name + ":1:in", mzml_path);
    parameters.setValue(adapter_name + ":1:database", db_path);
    parameters.setValue(adapter_name + ":1:out", out_path);

    // store data in temporary files
    MzMLFile spectra_file;
    spectra_file.store(mzml_path, exp);
    FASTAFile database;
    database.store(db_path, fasta_data);

    String ini_path = tmp_dir.getPath() + "parameters.INI";
    writeIniFile_(parameters, ini_path);

    // run identification search
    String proc_stdout;
    String proc_stderr;
    auto lam_out = [&](const String& out) { proc_stdout += out; };
    auto lam_err = [&](const String& out) { proc_stderr += out; };

    ExternalProcess ep(lam_out, lam_err);
    OPENMS_LOG_DEBUG << "Running " << adapter_name << "..." << endl << endl;
    const auto& rt = ep.run(adapter_name.toQString(), QStringList() << "-ini" << ini_path.toQString(), tmp_dir.getPath().toQString(), true);
    if (rt != ExternalProcess::RETURNSTATE::SUCCESS)
    { // error occured
      OPENMS_LOG_ERROR << "An error occured while running " << adapter_name << "." << endl;
      OPENMS_LOG_ERROR << "Standard output: " << proc_stdout << endl;
      OPENMS_LOG_ERROR << "Standard error: " << proc_stderr << endl;
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Return state was: " + static_cast<Int>(rt));
    }
    // search was successful

    // load result
    vector<ProteinIdentification> prot_ids;
    vector<PeptideIdentification> pep_ids;
    IdXMLFile comet_out;
    comet_out.load(out_path, prot_ids, pep_ids);

    // annotate target/decoy information
    PeptideIndexing indexer;
    FASTAContainer<TFI_Vector> proteins(fasta_data);
    OPENMS_LOG_DEBUG << "Running PeptideIndexer functionalities ..." << endl << endl;
    OPENMS_LOG_INFO.remove(cout); // prevent indexer from writing statistic
    PeptideIndexing::ExitCodes indexer_exit = indexer.run(proteins, prot_ids, pep_ids);
    OPENMS_LOG_INFO.insert(cout); // revert logging change
    if (indexer_exit != PeptideIndexing::ExitCodes::EXECUTION_OK)
    {
      OPENMS_LOG_ERROR << "An error occured while trying to index the search results." << endl;
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Return state was: " + static_cast<Int>(indexer_exit));
    }

    // calculate q-values
    Param p;
    p.setValue("use_all_hits", "true");
    p.setValue("add_decoy_peptides", "true");
    p.setValue("add_decoy_proteins", "true");

    FalseDiscoveryRate fdr;
    fdr.setParameters(p);
    OPENMS_LOG_DEBUG << "Calculating q-values ..." << endl << endl;
    fdr.apply(pep_ids);

    return pep_ids;
  }

  Size DBSuitability::countNumberOfIdentifications_(const vector<PeptideIdentification>& pep_ids) const
  {
    Size count{};
    for (const auto& pep_id : pep_ids)
    {
      if (pep_id.getHits().empty()) continue; // skip empty IDs without hits
      
      const PeptideHit& hit = pep_id.getHits()[0];

      if (!hit.metaValueExists("target_decoy"))
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No target/decoy annotation found. Make sure PeptideIndexer is run beforehand.");
      }
      if (hit.getMetaValue("target_decoy") == "decoy") continue; // skip decoy IDs

      if (!passesFDR_(hit, this->getParameters().getValue("FDR"))) continue;

      ++count; // count the rest
    }
    return count;
  }
}
