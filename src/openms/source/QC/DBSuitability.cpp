// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/MATH/StatisticFunctions.h>
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
    : DefaultParamHandler("DBSuitability"), results_{}, decoy_pattern_(string(DecoyHelper::regexstr_prefix + "|" + DecoyHelper::regexstr_suffix))
  {
    defaults_.setValue("no_rerank", "false", "Use this flag if you want to disable re-ranking. Cases, where a de novo peptide scores just higher than the database peptide, are overlooked and counted as a de novo hit. This might underestimate the database quality.");
    defaults_.setValidStrings("no_rerank", { "true", "false" });
    defaults_.setValue("reranking_cutoff_percentile", 0.01, "Swap a top-scoring deNovo hit with a lower scoring DB hit if their xcorr score difference is in the given percentile of all score differences between the first two decoy hits of a PSM. The lower the value the lower the decoy cut-off will be. Therefore it will be harder for a lower scoring DB hit to be re-ranked to the top.");
    defaults_.setMinFloat("reranking_cutoff_percentile", 0.);
    defaults_.setMaxFloat("reranking_cutoff_percentile", 1.);
    defaults_.setValue("FDR", 0.01, "Filter peptide hits based on this q-value. (e.g., 0.05 = 5 % FDR)");
    defaults_.setMinFloat("FDR", 0.);
    defaults_.setMaxFloat("FDR", 1.);
    defaults_.setValue("number_of_subsampled_runs", 1, "Controls how many runs should be done for calculating corrected suitability. (0 : number of runs will be estimated automaticly) ATTENTION: For each run a seperate ID-search is performed. This can result in some serious run time.");
    defaults_.setMinInt("number_of_subsampled_runs", 0);
    defaults_.setValue("keep_search_files", "false", "Set this flag if you wish to keep the files used by and produced by the internal ID search.");
    defaults_.setValidStrings("keep_search_files", { "true", "false" });
    defaults_.setValue("disable_correction", "false", "Set this flag to disable the calculation of the corrected suitability.");
    defaults_.setValidStrings("disable_correction", { "true", "false" });
    defaults_.setValue("force", "false", "Set this flag to enforce re-ranking when no cross correlation score is present. For re-ranking the default score found at each peptide hit is used. Use with care!");
    defaults_.setValidStrings("force", { "true", "false" });
    defaultsToParam_();
  }
  
  void DBSuitability::compute(vector<PeptideIdentification>&& pep_ids, const MSExperiment& exp, const vector<FASTAFile::FASTAEntry>& original_fasta, const std::vector<FASTAFile::FASTAEntry>& novo_fasta, const ProteinIdentification::SearchParameters& search_params)
  {
    for (const auto& id : pep_ids)
    {
      if (id.getScoreType() == "q-value") // q-value as score?
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "q-value found at PeptideIdentifications. That is not allowed! Please make sure FDR did not run previously.");
      }
      if (id.getHits().empty()) continue;
      if (id.getHits()[0].metaValueExists("q-value")) // q-value at meta values?
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "q-value found at PeptideIdentifications. That is not allowed! Please make sure FDR did not run previously.");
      }

      if (!id.getHits()[0].metaValueExists("MS:1002252")) // no xcorr at meta values?
      {
        if (!param_.getValue("force").toBool())
        {
          OPENMS_LOG_WARN << "No cross correlation score found. Comet is recommended for identification search. Re-ranking will be turned off. Set the 'force' flag to re-enable re-ranking. Use with care!" << endl;
          param_.setValue("no_rerank", "true");
        }
        else
        {
          OPENMS_LOG_WARN << "'force' flag is set. Re-ranking (if not disabled) will be done with the default score. Be aware that this may result in undefined behaviour." << endl;
        }
      }
      break;
    }

    // calculate suitability
    results_.emplace_back();
    SuitabilityData& suitability_data_full = results_.back();

    // make sure pep_ids are sorted
    for (auto& pep_id : pep_ids)
    {
      pep_id.sort();
    }
    calculateSuitability_(pep_ids, suitability_data_full);

    if(!param_.getValue("disable_correction").toBool())
    {
      pair<String, Param> search_info = extractSearchAdapterInfoFromMetaValues_(search_params);
      
      // calculate correction of suitability with extrapolation
      UInt number_of_runs = param_.getValue("number_of_subsampled_runs");
      if (number_of_runs == 0)
      {
        number_of_runs = ceil(original_fasta.size() / double(numberOfUniqueProteins_(pep_ids)));
      }

      // sampled run(s)
      // TODO: maybe multiple runs? could be controlled with a parameter
      double subsampling_rate = 0.5;
      vector<SuitabilityData> subsampled_results;

      UInt current_run = 0;
      while (current_run < number_of_runs)
      {
        ++current_run;
        vector<FASTAFile::FASTAEntry> sampled_db = getSubsampledFasta_(original_fasta, subsampling_rate);
        sampled_db.insert(sampled_db.end(), novo_fasta.begin(), novo_fasta.end());
        appendDecoys_(sampled_db);
        vector<PeptideIdentification> subsampled_ids = runIdentificationSearch_(exp, sampled_db, search_info.first, search_info.second);
        // make sure pep_ids are sorted
        for (auto& pep_id : subsampled_ids)
        {
          pep_id.sort();
        }

        SuitabilityData suitability_data_sampled;
        calculateSuitability_(subsampled_ids, suitability_data_sampled);
        subsampled_results.push_back(suitability_data_sampled);
      }

      SuitabilityData median_sampled_data = subsampled_results[getIndexWithMedianNovoHits_(subsampled_results)];

      suitability_data_full.setCorrectionFactor(calculateCorrectionFactor_(suitability_data_full, median_sampled_data, subsampling_rate));

      // fill in theoretical suitability if re-ranking hadn't happen
      SuitabilityData no_rerank = suitability_data_full.simulateNoReRanking();
      SuitabilityData no_rerank_sampled = median_sampled_data.simulateNoReRanking();

      double factor_no_rerank = calculateCorrectionFactor_(no_rerank, {no_rerank_sampled}, subsampling_rate);

      suitability_data_full.suitability_corr_no_rerank = double(no_rerank.num_top_db) / (no_rerank.num_top_novo * factor_no_rerank + no_rerank.num_top_db);
    }
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

    for (const auto& hit : pep_id.getHits())
    {
      if (!hit.metaValueExists("target_decoy"))
      {
        throw(Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No target/decoy information found! Make sure 'PeptideIndexer' is run beforehand."));
      }

      if (decoy_1 == DBL_MAX && hit.getMetaValue("target_decoy") == "decoy")
      {
        decoy_1 = extractScore_(hit);
        continue;
      }
      if (decoy_1 < DBL_MAX && hit.getMetaValue("target_decoy") == "decoy")
      {
        decoy_2 = extractScore_(hit);
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

    // sort the diffs decreasing and get the nth-percentile diff
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
    const set<String>& accessions = hit.extractProteinAccessionsSet();
    for (const String& acc : accessions)
    {
      if (acc.find(Constants::UserParam::CONCAT_PEPTIDE) == String::npos && !boost::regex_search(String(acc).toLower(), decoy_pattern_))
      {
        return false;
      }
    }
    return true;
  }

  bool DBSuitability::checkScoreBetterThanThreshold_(const PeptideHit& hit, double threshold, bool higher_score_better) const
  {
    if (higher_score_better)
    {
      if (hit.getScore() < threshold)
        return false;
      return true;
    }
    if (hit.getScore() > threshold)
      return false;
    return true;
  }

  pair<String, Param> DBSuitability::extractSearchAdapterInfoFromMetaValues_(const ProteinIdentification::SearchParameters& search_params) const
  {
    Param p;
    // list of all allowed adapters
    vector<String> working_adapters{ "CometAdapter", "MSGFPlusAdapter", "MSFraggerAdapter", "MyriMatchAdapter", "OMSSAAdapter", "XTandemAdapter" };

    vector<String> keys;
    search_params.getKeys(keys);

    // find adapter name
    String adapter;
    for (const String& key : keys)
    {
      for (const String& a : working_adapters)
      {
        // look for adapter name as a prefix of a search parameter
        if (key.compare(0, a.size() + 1, a + ":") == 0)
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
    bool keep_files = param_.getValue("keep_search_files").toBool();
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
    FileHandler spectra_file;
    spectra_file.storeExperiment(mzml_path, exp,{FileTypes::MZML});
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
      throw Exception::InternalToolError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Return state was: " + String(static_cast<Int>(rt)));
    }
    // search was successful

    // load result
    vector<ProteinIdentification> prot_ids;
    vector<PeptideIdentification> pep_ids;
    FileHandler id_file;
    id_file.loadIdentifications(out_path, prot_ids, pep_ids, {FileTypes::IDXML});

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
      throw Exception::InternalToolError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Return state was: " + String(static_cast<Int>(indexer_exit)));
    }

    if (keep_files)
    {
      id_file.storeIdentifications(tmp_dir.getPath() + "indexed_pre_FDR.idXML", prot_ids, pep_ids, {FileTypes::IDXML});
    }

    return pep_ids;
  }

  std::vector<FASTAFile::FASTAEntry> DBSuitability::getSubsampledFasta_(const std::vector<FASTAFile::FASTAEntry>& fasta_data, double subsampling_rate) const
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

    Math::RandomShuffler shuffler;
    shuffler.seed(UniqueIdGenerator::getUniqueId());
    std::vector<int> rnd_indices(fasta_data.size());
    std::iota(std::begin(rnd_indices), std::end(rnd_indices), 0);
    shuffler.portable_random_shuffle(rnd_indices.begin(), rnd_indices.end());

    Size curr_AA{};
    vector<FASTAFile::FASTAEntry> sampled_fasta;
    for (const int i : rnd_indices)
    {
      if (curr_AA >= num_AS_written) break;
      sampled_fasta.push_back(fasta_data[i]);
      curr_AA += fasta_data[i].sequence.size();
    }
    return sampled_fasta;
  }

  void DBSuitability::calculateSuitability_(const std::vector<PeptideIdentification>& pep_ids, SuitabilityData& data) const
  {
    // make sure no old data messes up the calculations
    data.clear();

    bool no_re_rank = param_.getValue("no_rerank").toBool();
    double cut_off_fract = param_.getValue("reranking_cutoff_percentile");

    if (pep_ids.empty())
    {
      OPENMS_LOG_WARN << "No peptide identifications found in given idXML! No calculations performed." << endl;
      return;
    }

    bool hsb = pep_ids[0].isHigherScoreBetter();

    // calculate score that corresponds to the FDR cut-off
    double score_cut_off;
    {
      vector<PeptideIdentification> ids_copy(pep_ids);

      FalseDiscoveryRate fdr;
      fdr.apply(ids_copy);

      score_cut_off = getScoreMatchingFDR_(ids_copy, param_.getValue("FDR"), pep_ids[0].getScoreType(), hsb);
    }

    if (!no_re_rank)
    {
      data.cut_off = getDecoyCutOff_(pep_ids, cut_off_fract);
    }

    for (const PeptideIdentification& pep_id : pep_ids)
    {
      const vector<PeptideHit>& hits = pep_id.getHits();

      if (hits.empty())
      {
        continue;
      }
      const PeptideHit& top_hit = hits[0];

      // skip if the top hit is a decoy hit
      if (!top_hit.metaValueExists("target_decoy"))
      {
        throw(Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No target/decoy information found! Make sure 'PeptideIndexer' is run beforehand."));
      }
      if (top_hit.getMetaValue("target_decoy") == "decoy") continue;

      // skip if top hit is out ouf FDR
      if (!checkScoreBetterThanThreshold_(top_hit, score_cut_off, hsb)) continue;

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
        if (!checkScoreBetterThanThreshold_(hits[i], score_cut_off, hsb)) break;

        // check if target, also check for "target+decoy" value
        String td_info(hits[i].getMetaValue("target_decoy"));
        if (td_info.find("target") != 0)
        {
          continue;
        }
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
      // score difference smaller than cut-off                              or       scores equal using floating point precision
      //                                                                             ('getMonoWeight()' can sometimes have a float error, was unable to fix it (maybe the error happens while reading the input idXML))
      if (extractScore_(top_hit) - extractScore_(second_hit) < data.cut_off || abs(extractScore_(top_hit) - extractScore_(second_hit)) < std::numeric_limits<float>::epsilon())
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

  double DBSuitability::extractScore_(const PeptideHit& pep_hit) const
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

      return pep_hit.getScore(); // uses q-value from FDR
    }
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

  UInt DBSuitability::numberOfUniqueProteins_(const std::vector<PeptideIdentification>& peps, UInt number_of_hits) const
  {
    set<String> proteins;

    for (const auto& pep : peps)
    {
      vector<PeptideHit> hits = pep.getHits();
      if (hits.empty())
        continue;

      for (Size i = 0; (i < number_of_hits && i < hits.size()); ++i)
      {
        const PeptideHit& hit = hits[i];

        if (!hit.metaValueExists("target_decoy"))
        {
          throw(Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No target/decoy information found! Make sure 'PeptideIndexer' is run beforehand."));
        }
        if (hit.getMetaValue("target_decoy") == "decoy")
          continue;// skip if the hit is a decoy hit

        // insert protein accessions
        const set<String> accessions = hit.extractProteinAccessionsSet();
        for (const String& acc : accessions)
        {
          if (acc.find(Constants::UserParam::CONCAT_PEPTIDE) != String::npos)// skip novo accessions
          {
            continue;
          }
          if (boost::regex_search(String(acc).toLower(), decoy_pattern_))// skip decoy accessions (this can happen if the hit is 'target+decoy'.)
          {
            continue;
          }

          proteins.insert(acc);// insert the rest
        }
      }
    }

    return proteins.size();
  }

  Size DBSuitability::getIndexWithMedianNovoHits_(const vector<SuitabilityData>& data) const
  {
    if (data.empty())
    {
      throw(Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No suitability data given!"));
    }

    vector<Size> novo_data;
    map<Size, Size> novo_hits_to_data;
    for (size_t i = 0; i < data.size(); ++i)
    {
      const Size num_top_novo = data[i].num_top_novo;
      novo_data.push_back(num_top_novo);
      novo_hits_to_data[num_top_novo] = i;
    }
    
    std::sort(novo_data.begin(), novo_data.end());

    return novo_hits_to_data.at(novo_data[ceil(novo_data.size() / 2)]);
  }

  double DBSuitability::getScoreMatchingFDR_(const std::vector<PeptideIdentification>& pep_ids, double FDR, const String& score_name, bool higher_score_better) const
  {
    double worst_score = DBL_MAX;
    if (!higher_score_better)
    {
      worst_score = -DBL_MAX;
    }

    for (const auto& id: pep_ids)
    {
      const vector<PeptideHit>& hits = id.getHits();

      if (hits.empty()) continue;

      if (id.getScoreType() != "q-value") // did FDR run?
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No q-value found at peptide identification.");
      }

      const PeptideHit& top_hit = hits[0]; // get first hit

      if (!checkScoreBetterThanThreshold_(top_hit, FDR, false)) continue; // does this hit make the FDR?

      // extract score from metavalues
      double score;
      bool score_found = false;

      vector<String> meta_keys;
      top_hit.getKeys(meta_keys);
      for (const String& key : meta_keys)
      {
        if (key.find(score_name) != String::npos)
        {
          score = top_hit.getMetaValue(key);
          score_found = true;
          break;
        }
      }

      if (!score_found)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "'" + score_name + "' not found. The given score name has to exist as a meta value.");
      }

      
      // is this score the worst score yet?
      if (higher_score_better)
      {
        if (score < worst_score)
        {
          worst_score = score;
        }
        continue;
      }
      if (score > worst_score)
      {
        worst_score = score;
      }
      continue;
    }
    return worst_score;
  }

  void DBSuitability::SuitabilityData::clear()
  {
    *this = SuitabilityData();
  }

  void DBSuitability::SuitabilityData::setCorrectionFactor(double factor)
  {
    if (num_top_db == 0 && num_top_novo == 0)
    {
      throw(Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No suitability data found. Can't apply correction factor."));
    }
    corr_factor = factor;
    num_top_novo_corr = num_top_novo * factor;
    suitability_corr = num_top_db / (num_top_db + num_top_novo_corr);
  }

  double DBSuitability::SuitabilityData::getCorrectionFactor() const
  {
    return this->corr_factor;
  }

  double DBSuitability::SuitabilityData::getCorrectedNovoHits() const
  {
    return this->num_top_novo_corr;
  }

  double DBSuitability::SuitabilityData::getCorrectedSuitability() const
  {
    return this->suitability_corr;
  }

  DBSuitability::SuitabilityData DBSuitability::SuitabilityData::simulateNoReRanking() const
  {
    SuitabilityData simulated_data;
    simulated_data.num_top_db = num_top_db - num_re_ranked;
    simulated_data.num_top_novo = num_top_novo + num_re_ranked;
    return simulated_data;
  }
}// namespace OpenMS
