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

#include <random>

using namespace std;

namespace OpenMS
{
  DBSuitability::DBSuitability()
    : DefaultParamHandler("DBSuitability"), results_{}
  {
    defaults_.setValue("no_rerank", "false", "Use this flag if you want to disable re-ranking. Cases, where a de novo peptide scores just higher than the database peptide, are overlooked and counted as a de novo hit. This might underestimate the database quality.");
    defaults_.setValidStrings("no_rerank", { "true", "false" });
    defaults_.setValue("reranking_cutoff_percentile", 0.01, "Swap a top-scoring deNovo hit with a lower scoring DB hit if their xcorr score difference is in the given percentile of all score differences between the first two decoy hits of a PSM. The lower the value the lower the decoy cut-off will be. Therefore it will be harder for a lower scoring DB hit to be re-ranked to the top.");
    defaults_.setMinFloat("reranking_cutoff_percentile", 0.);
    defaults_.setMaxFloat("reranking_cutoff_percentile", 1.);
    defaults_.setValue("FDR", 0.01, "Filter peptide hits based on this q-value. (e.g., 0.05 = 5 % FDR)");
    defaults_.setMinFloat("FDR", 0.);
    defaults_.setMaxFloat("FDR", 1.);
    defaults_.setValue("force", "false", "Set this flag to enforce re-ranking when no cross correlation score is present. For re-ranking the default score found at each peptide hit is used. Use with care!");
    defaults_.setValidStrings("force", { "true", "false" });
    defaults_.setValue("keep_search_files", "false", "Set this flag if you wish to keep the files used by and produced by the internal ID search.");
    defaults_.setValidStrings("keep_search_files", { "true", "false" });
    defaultsToParam_();
  }
  
  void DBSuitability::compute(vector<PeptideIdentification> pep_ids, const MSExperiment& exp, vector<FASTAFile::FASTAEntry> original_fasta, std::vector<FASTAFile::FASTAEntry> novo_fasta, const ProteinIdentification::SearchParameters& search_params)
  {
    pair<String, Param> search_info = extractSearchAdapterInfoFromMetaValues_(search_params);

    if (pep_ids[0].getScoreType() == "q-value") // q-value as score?
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "q-value found at PeptideIdentifications. That is not allowed! Please make sure FDR did not run previously.");
    }

    for (const auto& id : pep_ids)
    {
      if (id.getHits().empty()) continue;
      if (id.getHits()[0].metaValueExists("q-value")) // q-value at meta values?
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "q-value found at PeptideIdentifications. That is not allowed! Please make sure FDR did not run previously.");
      }

      if (!id.getHits()[0].metaValueExists("MS:1002252")) // no xcorr at meta values?
      {
        if (!param_.getValue("force").toBool())
        {
          OPENMS_LOG_WARN << "No cross correlation score found. Make sure Comet is used for identification search. Re-ranking will be turned off. Set the 'force' flag to re-enable re-ranking. Use with care!" << endl;
          param_.setValue("no_rerank", "true");
        }
        else
        {
          OPENMS_LOG_WARN << "'force' flag is set. Re-ranking (if not disabled) will be done with the default score. Be aware that this may result in undefined behaviour." << endl;
        }
      }
      break;
    }

    Param p;
    p.setValue("use_all_hits", "true");
    p.setValue("add_decoy_peptides", "true");
    p.setValue("add_decoy_proteins", "true");

    FalseDiscoveryRate fdr;
    fdr.setParameters(p);
    fdr.apply(pep_ids);

    // calculate suitability
    results_.push_back(SuitabilityData());
    SuitabilityData& suitability_data_full = results_.back();
    calculateSuitability_(pep_ids, suitability_data_full);

    // calculate correction of suitability with extrapolation

    // sampled run
    // TODO: maybe multiple runs? could be controlled with a parameter
    double subsampling_rate = 0.5;
    vector<FASTAFile::FASTAEntry> sampled_db = getSubsampledFasta_(original_fasta, subsampling_rate);
    sampled_db.insert(sampled_db.end(), novo_fasta.begin(), novo_fasta.end());
    appendDecoys_(sampled_db);
    vector<PeptideIdentification> subsampled_ids = runIdentificationSearch_(exp, sampled_db, search_info.first, search_info.second);

    SuitabilityData suitability_data_sampled;
    calculateSuitability_(subsampled_ids, suitability_data_sampled);

    suitability_data_full.setCorrectionFactor(calculateCorrectionFactor_(suitability_data_full, suitability_data_sampled, subsampling_rate));

    // fill in theoretical suitability if re-ranking hadn't happen
    SuitabilityData no_rerank = simulateNoReRanking_(suitability_data_full);
    SuitabilityData no_rerank_sampled = simulateNoReRanking_(suitability_data_sampled);

    double factor_no_rerank = calculateCorrectionFactor_(no_rerank, no_rerank_sampled, subsampling_rate);

    suitability_data_full.suitability_corr_no_rerank = double(no_rerank.num_top_db) / (no_rerank.num_top_novo * factor_no_rerank + no_rerank.num_top_db);
  }

  const std::vector<DBSuitability::SuitabilityData>& DBSuitability::getResults() const
  {
    return results_;
  }

  double DBSuitability::getDecoyDiff_(const PeptideIdentification& pep_id) const
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

      if (decoy_1 == DBL_MAX && hit.getMetaValue("target_decoy") == "decoy")
      {
        decoy_1 = getRightScore_(hit);
        continue;
      }
      if (decoy_1 < DBL_MAX && hit.getMetaValue("target_decoy") == "decoy")
      {
        decoy_2 = getRightScore_(hit);
        break;
      }
    }

    if (decoy_2 < DBL_MAX) // if there are two decoy hits
    {
      diff = abs(decoy_1 - decoy_2);
    }

    // if there aren't two decoy hits DBL_MAX is returned
    return diff;
  }

  double DBSuitability::getDecoyCutOff_(const vector<PeptideIdentification>& pep_ids, double reranking_cutoff_percentile) const
  {
    if (reranking_cutoff_percentile < 0 || reranking_cutoff_percentile > 1)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "'reranking_cutoff_percentile' is not within its allowed range [0,1]. Please select a valid value.");
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
      throw(Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Under 20 % of peptide identifications have two decoy hits. This is not enough for re-ranking. Use the 'no_rerank' flag to still compute a suitability score."));
    }

    // sort the diffs decreasing and get the percentile one
    UInt index = round(reranking_cutoff_percentile * diffs.size());
    
    if (index >= diffs.size())
    {
      return *max_element(diffs.begin(), diffs.end());
    }

    nth_element(diffs.begin(), diffs.begin() + index, diffs.end());

    return diffs[index];
  }

  bool DBSuitability::isNovoHit_(const PeptideHit& hit) const
  {
    const boost::regex decoy_pattern(DecoyHelper::getPrefixRegex() + "|" + DecoyHelper::getSuffixRegex());
    const set<String>& accessions = hit.extractProteinAccessionsSet();
    for (const String& acc : accessions)
    {
      if (acc.find(Constants::UserParam::CONCAT_PEPTIDE) == String::npos && !boost::regex_search(acc, decoy_pattern))
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

    vector<String> keys;
    search_params.getKeys(keys);

    // find adapter name
    String adapter;
    for (const String& key : keys)
    {
      for (const String& a : working_adapters)
      {
        if (key.compare(0, a.size(), a) == 0)
        {
          // used adapter found
          adapter = a;
          break;
        }
      }
      if (!adapter.empty())
      {
        break;
      }
    }

    if (adapter.empty()) // non of the allowed adapter names where found in the meta values
    {
      String message;
      message = "No parameters found for any of the allowed adapters in the given meta values. Allowed are:\n";
      message += ListUtils::concatenate(working_adapters, ", ");
      message = "\n";
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message);
    }

    // extract parameters
    for (const String& key : keys)
    {
      if (key.compare(0, adapter.size(), adapter) != 0) continue; // does adapter appear in meta value key?
      p.setValue(key, search_params.getMetaValue(key));
    }

    OPENMS_LOG_DEBUG << "Parameters for the following adapter were found: " << adapter << endl;

    return make_pair(adapter, p);
  }

  void DBSuitability::writeIniFile_(const Param& parameters, const String& filename) const
  {
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
    bool keep_files = !param_.getValue("keep_search_files").toBool();
    File::TempDir tmp_dir(keep_files);
    String mzml_path = tmp_dir.getPath() + "spectra.mzML";
    String db_path = tmp_dir.getPath() + "database.FASTA";
    String out_path = tmp_dir.getPath() + "out.idXML";

    // override the in- and output files in the parameters
    if (!parameters.exists(adapter_name + ":1:in") || !parameters.exists(adapter_name + ":1:database") || !parameters.exists(adapter_name + ":1:out"))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "'in', 'out' or 'database' parameter not found! The search adapter is probably not supported anymore.");
    }
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
      throw Exception::InternalToolError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Return state was: " + static_cast<Int>(rt));
    }
    // search was successful

    // load result
    vector<ProteinIdentification> prot_ids;
    vector<PeptideIdentification> pep_ids;
    IdXMLFile id_file;
    id_file.load(out_path, prot_ids, pep_ids);

    if (keep_files)
    {
      id_file.store(tmp_dir.getPath() + "id_search_out.idXML", prot_ids, pep_ids);
    }

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
      throw Exception::InternalToolError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Return state was: " + static_cast<Int>(indexer_exit));
    }

    if (keep_files)
    {
      id_file.store(tmp_dir.getPath() + "indexed_pre_FDR.idXML", prot_ids, pep_ids);
    }

    // calculate q-values
    FalseDiscoveryRate fdr;
    Param p(fdr.getParameters());
    if (!p.exists("use_all_hits") || !p.exists("add_decoy_peptides") || !p.exists("add_decoy_proteins"))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "FDR parameters probably changed. 'use_all_hits', 'add_decoy_peptides' or 'add_decoy_proteins' not found.");
    }
    p.setValue("use_all_hits", "true");
    p.setValue("add_decoy_peptides", "true");
    p.setValue("add_decoy_proteins", "true");

    fdr.setParameters(p);
    OPENMS_LOG_DEBUG << "Calculating q-values ..." << endl << endl;
    fdr.apply(pep_ids);

    if (keep_files)
    {
      id_file.store(tmp_dir.getPath() + "indexed_post_FDR.idXML", prot_ids, pep_ids);
    }

    return pep_ids;
  }

  std::vector<FASTAFile::FASTAEntry> DBSuitability::getSubsampledFasta_(std::vector<FASTAFile::FASTAEntry> fasta_data, double subsampling_rate) const
  {
    if (subsampling_rate < 0 || subsampling_rate > 1)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Subsampling rate has to be between 0 and 1. Aborting!");
    }
    Size num_AA{};
    for (const auto& entry : fasta_data)
    {
      num_AA += entry.sequence.size();
    }
    double num_AS_written = num_AA * subsampling_rate;

    //random_device rd();
    //mt19937 g(rd());
    mt19937 g(0);
    shuffle(fasta_data.begin(), fasta_data.end(), g);

    Size curr_AA{};
    vector<FASTAFile::FASTAEntry> sampled_fasta;
    for (const auto& entry : fasta_data)
    {
      if (curr_AA >= num_AS_written) break;
      sampled_fasta.push_back(entry);
      curr_AA += entry.sequence.size();
    }
    return sampled_fasta;
  }

  void DBSuitability::calculateSuitability_(std::vector<PeptideIdentification> pep_ids, SuitabilityData& data) const
  {
    bool no_re_rank = param_.getValue("no_rerank").toBool();
    double cut_off_fract = param_.getValue("reranking_cutoff_percentile");
    double FDR = param_.getValue("FDR");

    if (pep_ids.empty())
    {
      OPENMS_LOG_WARN << "No peptide identifications found in given idXML! No calculations performed." << endl;
      return;
    }

    if (!no_re_rank)
    {
      data.cut_off = getDecoyCutOff_(pep_ids, cut_off_fract);
    }

    for (PeptideIdentification& pep_id : pep_ids)
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
      const PeptideHit* second_hit_ptr = nullptr;
      for (UInt i = 1; i < hits.size(); ++i)
      {
        // check for FDR
        if (!passesFDR_(hits[i], FDR)) break;

        // check if target, also check for "target+decoy" value
        String td_info(hits[i].getMetaValue("target_decoy"));
        if (td_info.find("target") != 0) continue;

        // check if hit is novo hit
        if (isNovoHit_(hits[i])) continue;

        second_hit_ptr = &hits[i];
        break;
      }
      if (second_hit_ptr == nullptr) // no second target hit with given FDR found
      {
        ++data.num_top_novo;
        continue;
      }

      const PeptideHit second_hit(*second_hit_ptr);

      // second hit is db hit
      ++data.num_interest;

      // check for re-ranking
      if (no_re_rank)
      {
        ++data.num_top_novo;
        continue;
      }

      // re-ranking
      if (getRightScore_(top_hit) - getRightScore_(second_hit) <= data.cut_off)
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

    data.suitability_no_rerank = double(data.num_top_db - data.num_re_ranked) / (data.num_top_db + data.num_top_novo);
  }

  void DBSuitability::appendDecoys_(std::vector<FASTAFile::FASTAEntry>& fasta) const
  {
    fasta.reserve(fasta.size() * 2);

    for (auto& entry : fasta)
    {
      ProteaseDigestion digestion;
      digestion.setEnzyme("Trypsin");
      std::vector<AASequence> peptides;
      digestion.digest(AASequence::fromString(entry.sequence), peptides);
      String new_sequence = "";
      for (auto const& peptide : peptides)
      {
        OpenMS::TargetedExperiment::Peptide p;
        p.sequence = peptide.toString();
        OpenMS::TargetedExperiment::Peptide decoy_p = MRMDecoy::reversePeptide(p, true, true, "");
        new_sequence += decoy_p.sequence;
      }
      FASTAFile::FASTAEntry decoy_entry;
      decoy_entry.sequence = new_sequence;
      decoy_entry.identifier = "DECOY_" + entry.identifier;
      fasta.push_back(decoy_entry);
    }
  }

  double DBSuitability::getRightScore_(const PeptideHit& pep_hit) const
  {
    if (pep_hit.metaValueExists("MS:1002252")) // use xcorr
    {
      return double(pep_hit.getMetaValue("MS:1002252")) / pep_hit.getSequence().getMonoWeight(); // normalized by mw
    }
    else
    {
      if (!param_.getValue("force").toBool()) // without xcorr, need to force re-ranking
      {
        throw(Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No cross correlation score found at peptide hit. Only Comet search engine is supported for re-ranking. Set 'force' flag to use the default score for this. This may result in undefined behaviour and is not advised."));
      }

      return pep_hit.getScore();
    }
  }

  DBSuitability::SuitabilityData DBSuitability::simulateNoReRanking_(const SuitabilityData& data) const
  {
    SuitabilityData simulated_data;
    simulated_data.num_top_db = data.num_top_db - data.num_re_ranked;
    simulated_data.num_top_novo = data.num_top_novo + data.num_re_ranked;    
    return simulated_data;
  }

  double DBSuitability::calculateCorrectionFactor_(const SuitabilityData& data_full, const SuitabilityData& data_sampled, double sampling_rate) const
  {
    if (sampling_rate >= 1 || sampling_rate < 0)
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The sampling rate has to be element of [0,1).");
    }

    double db_slope = (int(data_sampled.num_top_db) - int(data_full.num_top_db)) / (sampling_rate - 1);
    double deNovo_slope = (int(data_sampled.num_top_novo) - int(data_full.num_top_novo)) / (sampling_rate - 1);

    return -(db_slope) / (deNovo_slope);
  }

  void DBSuitability::SuitabilityData::setCorrectionFactor(double factor)
  {
    if (num_top_db == 0 || num_top_novo == 0)
    {
      throw(Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No suitability data found. Can't apply correction factor."));
    }
    corr_factor = factor;
    num_top_novo_corr = num_top_novo * factor;
    suitability_corr = num_top_db / (num_top_db + num_top_novo_corr);
  }

}// namespace OpenMS
