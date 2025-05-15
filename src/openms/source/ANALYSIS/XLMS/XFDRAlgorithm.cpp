// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Lukas Zimmermann, Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/XLMS/XFDRAlgorithm.h>
#include <OpenMS/CONCEPT/Constants.h>

// using namespace std;
using namespace OpenMS;

  XFDRAlgorithm::XFDRAlgorithm()
    : DefaultParamHandler("XFDRAlgorithm")
  {
    defaults_.setValue(param_decoy_string_, "DECOY_", "Prefix of decoy protein ids. The correspondig target protein id should be retrievable by deleting this prefix.");
    defaults_.setValue(param_minborder_, -50.0, "Filter for minimum precursor mass error (ppm) before FDR estimation. Values outside of the tolerance window of the original search will effectively disable this filter.");
    defaults_.setValue(param_maxborder_, 50.0, "Filter for maximum precursor mass error (ppm) before FDR estimation. Values outside of the tolerance window of the original search will effectively disable this filter.");

    defaults_.setValue(param_mindeltas_, 0.0, "Filter for delta score, 0 disables the filter. Minimum delta score required, hits are rejected if larger or equal. The delta score is a ratio of the score of a hit and the score of the next best hit to the same spectrum, so the value range is between 0 and 1 with 1.0 meaning the scores are equal and 0.5 meaning the next best score is half as high as the current one.");
    defaults_.setMinFloat(param_mindeltas_, 0.0);
    defaults_.setMaxFloat(param_mindeltas_, 1.0);

    defaults_.setValue(param_minionsmatched_, 0, "Filter for minimum matched ions per peptide.");
    defaults_.setMinInt(param_minionsmatched_, 0);

    std::vector<std::string> bool_strings = {"true","false"};

    defaults_.setValue(param_uniquexl_, "false", "Calculate statistics based only on unique IDs. For a set of IDs from equal candidates (same pair of peptides, modifications and cross-linked positions), only the highest scoring hit will be considered. By default the score distribution will be estimated using all 1st ranked candidates.");
    defaults_.setValidStrings(param_uniquexl_, bool_strings);

    defaults_.setValue(param_no_qvalues_, "false", "Do not transform simple FDR to q-values");
    defaults_.setValidStrings(param_no_qvalues_, bool_strings);

    defaults_.setValue(param_minscore_, -10.0, "Minimum score to be considered for FDR calculation. A number lower than the lowest score will effectively disable this filter.");

    defaults_.setValue(param_binsize_, 0.0001, "Bin size for the cumulative histograms for score distributions. Should be about the same size as the smallest expected difference between scores. Smaller numbers will make XFDR more robust, but much slower. Negative numbers are not allowed. Should only be changed if the range of the main score changes or another score than the OpenPepXL score is used.");
    defaults_.setMinFloat(param_binsize_, 1e-15);

    defaultsToParam_();
  }

  XFDRAlgorithm::~XFDRAlgorithm() = default;

  void XFDRAlgorithm::updateMembers_()
  {
    decoy_string_ = static_cast<String>(param_.getValue(param_decoy_string_).toString());
    arg_mindeltas_ = static_cast<double>(param_.getValue(param_mindeltas_));
    arg_minborder_ = static_cast<double>(param_.getValue(param_minborder_));
    arg_maxborder_ = static_cast<double>(param_.getValue(param_maxborder_));
    arg_minionsmatched_ = static_cast<Int>(param_.getValue(param_minionsmatched_));
    arg_minscore_ = static_cast<double>(param_.getValue(param_minscore_));
    arg_uniquex_ = (param_.getValue(param_uniquexl_) == "true" ? true : false);
    arg_no_qvalues_ = (param_.getValue(param_no_qvalues_) == "true" ? true : false);
    arg_binsize_ = static_cast<double>(param_.getValue(param_binsize_));
    min_score_ = 0;
    max_score_ = arg_minscore_;
  }

  XFDRAlgorithm::ExitCodes XFDRAlgorithm::run(std::vector<PeptideIdentification>& peptide_ids, ProteinIdentification& protein_id)
  {
    writeArgumentsLog_();
    std::cout << "Initializing data structures..." << std::endl;
    // Initialize and validate data structures that are derived from the main peptide identification vector 'all_pep)ids'
    initDataStructures_(peptide_ids, protein_id);

    // Maps the cross link class to the encountered scores
    std::map<String, std::vector<double>> scores;
    UInt num_flagged(0);

    std::cout << "Collecting scores for each class..." << std::endl;
    // Loop through the peptides, apply filter, and assign cross-link types
    for (PeptideIdentification& pep_id : peptide_ids)
    {
      if (pep_id.getHits().empty())
      {
        continue;
      }

      // usually the 1st ranked hit should be the first in the list
      // but we make sure here and use the first hit with rank 1 that we find
      int rank_one_hit_index(0);
      for (Size i = 0; i < pep_id.getHits().size(); ++i)
      {
        PeptideHit& ph = pep_id.getHits()[i];
        if (int(ph.getMetaValue(Constants::UserParam::OPENPEPXL_XL_RANK)) == 1)
        {
          rank_one_hit_index = i;
          break;
        }
      }
      PeptideHit& ph = pep_id.getHits()[rank_one_hit_index];


      // if after the search above we don't have a rank 1 hit, skip this pep_id
      if (int(ph.getMetaValue(Constants::UserParam::OPENPEPXL_XL_RANK)) != 1)
      {
        continue;
      }

      // Attributes of peptide identification that can be used for filtering
      const double delta_score = ph.getMetaValue(Constants::UserParam::DELTA_SCORE);
      const double score = ph.getScore();

      double error_rel(0);
      if (ph.metaValueExists(Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM))
      {
        error_rel = ph.getMetaValue(Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM);
      }

      const Size min_ions_matched = getMinIonsMatched_(ph);
      const String id = ph.getMetaValue("OpenPepXL:id");

      // Only consider IDs which fullfill all filter criteria specified by the user
      if (   (arg_minborder_ <= error_rel)   // minborder fullfilled
          && (arg_maxborder_ >= error_rel)   // maxborder fullfilled
          && (arg_mindeltas_ == 0  || delta_score < arg_mindeltas_)
          && (min_ions_matched  >= (Size)arg_minionsmatched_)
          &&  score >= arg_minscore_
         )
      {
        // check for the unique ID criterion
        if (arg_uniquex_)
        {
          auto uid_it = std::find(this->unique_ids_.begin(), this->unique_ids_.end(), id);
          if (uid_it != this->unique_ids_.end())
          {
            int index = std::distance(this->unique_ids_.begin(), uid_it);
            if (this->unique_id_scores_[index] != ph.getScore())
            {
              // this is not the highest scoring ID for this candidate
              continue;
            }
          }
        }
        num_flagged++;

        ph.setMetaValue("XFDR:used_for_FDR", 1);

        for (const String &cross_link_class : this->cross_link_classes_[id])
        {
          scores[cross_link_class].push_back(score);
        }
      }
    }
    std::cout << "XFDR has used " + String(num_flagged) + " hits to calculate the FDR" << std::endl;

    // Log number of scores within each class
    std::cout << "Number of Scores for each class:" << std::endl;

    for (const auto &score : scores)
    {
      std::cout << score.first + ": " + score.second.size() << std::endl;
    }

    // Generate Histograms of the scores for each class
    // Use cumulative histograms to count the number of scores above consecutive thresholds
    std::map< String, Math::Histogram<> >  cum_histograms;
    for (const auto &class_scores: scores)
    {
      std::vector< double > current_scores = class_scores.second;

      Math::Histogram<> histogram(this->min_score_, this->max_score_, arg_binsize_);
      Math::Histogram<>::getCumulativeHistogram(current_scores.begin(), current_scores.end(), true, true, histogram);
      cum_histograms[class_scores.first] = histogram;
    }

    std::cout << "Calculating Score Distributions..." << std::endl;
    // Calculate FDR for interlinks
    std::vector< double > fdr_interlinks;
    this->fdr_xprophet_(cum_histograms, crosslink_class_interlinks_, crosslink_class_interdecoys_, crosslink_class_fulldecoysinterlinks_, fdr_interlinks, false);

    // Calculate FDR for intralinks
    std::vector< double > fdr_intralinks;
    this->fdr_xprophet_(cum_histograms, crosslink_class_intralinks_, crosslink_class_intradecoys_, crosslink_class_fulldecoysintralinks_, fdr_intralinks, false);

    // Calculate FDR for monolinks and looplinks
    std::vector< double > fdr_monolinks;
    this->fdr_xprophet_(cum_histograms, crosslink_class_monolinks_, crosslink_class_monodecoys_, "", fdr_monolinks, true);

    // Determine whether qTransform should be performed (and consequently the score type)
    // bool arg_no_qvalues = getFlag_(param_no_qvalues_);
    String score_type = arg_no_qvalues_ ? "FDR" : "q-value";

    if ( ! arg_no_qvalues_)
    {
      std::cout << "Performing qFDR transformation..." << std::endl;

      std::vector< double > qfdr_interlinks;
      this->calc_qfdr_(fdr_interlinks, qfdr_interlinks);

      std::vector< double > qfdr_intralinks;
      this->calc_qfdr_(fdr_intralinks, qfdr_intralinks);

      std::vector< double > qfdr_monolinks;
      this->calc_qfdr_(fdr_monolinks, qfdr_monolinks);

      fdr_interlinks = qfdr_interlinks;
      fdr_intralinks = qfdr_intralinks;
      fdr_monolinks = qfdr_monolinks;
    }

    std::cout << "Assigning FDRs..." << std::endl;
    // Assign FDR values to all identifications
    for (PeptideIdentification &pep_id : peptide_ids)
    {
      for (PeptideHit& ph : pep_id.getHits())
      {
        if ( ! ph.metaValueExists("XFDR:used_for_FDR"))
        {
          ph.setMetaValue("XFDR:used_for_FDR", 0);
        }
        double score = ph.getScore();

        StringList crosslink_types;
        assignTypes_(ph, crosslink_types);

        ph.setMetaValue("XFDR:fdr_type", score_type);

        // Assign FDR value as meta value and also set as score
        bool assigned = false;
        double fdr = 1;
        for (StringList::const_iterator crosslink_types_it = crosslink_types.begin();
            crosslink_types_it != crosslink_types.end(); ++crosslink_types_it)
        {
          String current_crosslink_type = *crosslink_types_it;
          Size idx = std::floor((score - this->min_score_) / arg_binsize_);
          if (   current_crosslink_type == crosslink_class_fulldecoysinterlinks_
              || current_crosslink_type == crosslink_class_hybriddecoysinterlinks_
              || current_crosslink_type == crosslink_class_interdecoys_
              || current_crosslink_type == crosslink_class_interlinks_)
          {
            fdr = fdr_interlinks[idx];
            assigned = true;
            break;
          }
          else if (   current_crosslink_type == crosslink_class_fulldecoysintralinks_
              || current_crosslink_type == crosslink_class_hybriddecoysintralinks_
              || current_crosslink_type == crosslink_class_intradecoys_
              || current_crosslink_type == crosslink_class_intralinks_)
          {
            fdr = fdr_intralinks[idx];
            assigned = true;
            break;
          }
          else if (   current_crosslink_type == crosslink_class_monodecoys_
              || current_crosslink_type == crosslink_class_monolinks_)
          {
            fdr = fdr_monolinks[idx];
            assigned = true;
            break;
          }
        }
        if ( assigned)
        {
          ph.setMetaValue("XFDR:FDR", fdr);
        }
        else
        {
          std::cout << "WARNING: A Crosslink could not be identified as either interlink, intralink, or monolink, so no FDR will be available for it." << std::endl;
        }
      }
    }
    return EXECUTION_OK;
  }

  void XFDRAlgorithm::initDataStructures_(std::vector<PeptideIdentification>& peptide_ids, ProteinIdentification& protein_id)
  {
    const String prot_identifier = protein_id.getIdentifier();

    // if the metaValue exists in search_params and the default value for XFDR was not changed, use the one in search_params
    ProteinIdentification::SearchParameters search_params = protein_id.getSearchParameters();
    if (search_params.metaValueExists("decoy_string") && decoy_string_ == "DECOY_")
    {
      decoy_string_ = search_params.getMetaValue("decoy_string");
    }

    // Preprocess all peptide identifications and construct derived data structures necessary for XFDR
    for (Size i = 0; i < peptide_ids.size(); ++i)
    {
      PeptideIdentification &pep_id = peptide_ids[i];
      if (pep_id.getHits().empty())
      {
        continue;
      }

      pep_id.setIdentifier(prot_identifier);

      std::vector< PeptideHit > &pep_hits = pep_id.getHits();

      for (PeptideHit& ph : pep_hits)
      {
        // Set the minScore and maxScore attribute depending on the input data
        const double score = ph.getScore();

        // Set score boundaries
        if (score < this->min_score_)
        {
          this->min_score_ = std::floor(score);
        }
        if (score > this->max_score_)
        {
          this->max_score_ = std::ceil(score);
        }
        assert(this->min_score_ <= this->max_score_);

        // figure out if crosslink is inter- or intra protein
        // for cases with multiple proteins, count as true, if any one possible combination of proteins fits the criteria
        // so both can be true at the same time (or false for mono-links)
        setIntraProtein_(ph, false);
        setInterProtein_(ph, false);

        if (ph.metaValueExists(Constants::UserParam::OPENPEPXL_XL_TYPE) && ph.getMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE) == "cross-link")
        {
          StringList alpha_prots;
          const std::vector<PeptideEvidence> pevs_alpha = ph.getPeptideEvidences();
          for (PeptideEvidence pev : pevs_alpha)
          {
            alpha_prots.push_back(pev.getProteinAccession());
          }
          StringList beta_prots = ListUtils::create<String>(ph.getMetaValue(Constants::UserParam::OPENPEPXL_BETA_ACCESSIONS).toString());

          for (String& alpha_prot : alpha_prots)
          {
            for (String& beta_prot : beta_prots)
            {
              if (isSameProtein_(alpha_prot, beta_prot, decoy_string_))
              {
                setIntraProtein_(ph, true);
              }
              else
              {
                setInterProtein_(ph, true);
              }
            }
          }
        }

        String id = getId_(ph);
        ph.setMetaValue("OpenPepXL:id", id);
        // candidates with the same ID will also have the same types
        if (this->cross_link_classes_.find(id) == this->cross_link_classes_.end())
        {
          assignTypes_(ph, this->cross_link_classes_[id]);
        }
      }
    }
    if (arg_uniquex_)
    {
      findTopUniqueHits_(peptide_ids);
    }
  }

  void XFDRAlgorithm::assignTypes_(PeptideHit &ph, StringList &types)
  {
    types.clear();
    bool xl_is_decoy = ph.getMetaValue(Constants::UserParam::TARGET_DECOY) == "decoy";

    // target or decoy
    if (xl_is_decoy)
    {
      types.push_back(crosslink_class_decoys_);
    }
    else
    {
      types.push_back(crosslink_class_targets_);
    }

    // intralinks
    if (ph.getMetaValue("XFDR:is_intraprotein").toBool() && (!xl_is_decoy))
    {
      types.push_back(crosslink_class_intralinks_);
    }

    // intradecoys
    if (ph.getMetaValue("XFDR:is_intraprotein").toBool() && xl_is_decoy)
    {
      types.push_back(crosslink_class_intradecoys_);
    }

    // interlinks
    if (ph.getMetaValue("XFDR:is_interprotein").toBool() && (!xl_is_decoy))
    {
      types.push_back(crosslink_class_interlinks_);
    }

    // interdecoys
    if (ph.getMetaValue("XFDR:is_interprotein").toBool() && xl_is_decoy)
    {
      types.push_back(crosslink_class_interdecoys_);
    }

    assert(ph.metaValueExists(Constants::UserParam::OPENPEPXL_XL_TYPE));
    String current_crosslink_type = ph.getMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE);

    // monolinks
    if ( (!xl_is_decoy) && (current_crosslink_type == "mono-link"
        ||  current_crosslink_type == "loop-link"))
    {
      types.push_back(crosslink_class_monolinks_);
    }

    // monodecoys
    if ( xl_is_decoy && (current_crosslink_type == "mono-link"
        ||  current_crosslink_type == "loop-link"))
    {
      types.push_back(crosslink_class_monodecoys_);
    }

    if (current_crosslink_type == "cross-link")
    {
      const bool alpha_is_decoy = ph.getMetaValue(Constants::UserParam::OPENPEPXL_TARGET_DECOY_ALPHA).toString() == "decoy";
      const bool beta_is_decoy = ph.getMetaValue(Constants::UserParam::OPENPEPXL_TARGET_DECOY_BETA).toString() == "decoy";

      // fulldecoysintralinks
      if (ph.getMetaValue("XFDR:is_intraprotein").toBool() && alpha_is_decoy && beta_is_decoy)
      {
        types.push_back(crosslink_class_fulldecoysintralinks_);
      }

      // fulldecoysinterlinks
      if (ph.getMetaValue("XFDR:is_interprotein").toBool() && alpha_is_decoy && beta_is_decoy)
      {
        types.push_back(crosslink_class_fulldecoysinterlinks_);
      }

      // hybriddecoysintralinks
      if (ph.getMetaValue("XFDR:is_intraprotein").toBool()
          && (( (!alpha_is_decoy) &&   beta_is_decoy)
          ||     (alpha_is_decoy  && (!beta_is_decoy))))
      {
        types.push_back(crosslink_class_hybriddecoysintralinks_);
      }

      // hybriddecoysinterlinks
      if (ph.getMetaValue("XFDR:is_interprotein").toBool()
          && (( (!alpha_is_decoy) &&   beta_is_decoy)
          ||     (alpha_is_decoy  && (!beta_is_decoy))))
      {
        types.push_back(crosslink_class_hybriddecoysinterlinks_);
      }
    }
  }

  void XFDRAlgorithm::fdr_xprophet_(std::map< String, Math::Histogram<> > & cum_histograms,
                    const String  & targetclass, const String & decoyclass, const String & fulldecoyclass,
                    std::vector< double > & fdr, bool mono) const
  {
    // Determine whether targetclass, decoyclass, and fulldecoyclass are present in the histogram map
    bool targetclass_present = cum_histograms.find(targetclass) != cum_histograms.end();
    bool decoyclass_present = cum_histograms.find(decoyclass) != cum_histograms.end();
    bool fulldecoyclass_present = cum_histograms.find(fulldecoyclass) != cum_histograms.end();

    for (double current_score = this->min_score_ +  (arg_binsize_/2);
        current_score <= this->max_score_ - (arg_binsize_/2);
        current_score += arg_binsize_)
    {
      double estimated_n_decoys = decoyclass_present ? cum_histograms[decoyclass].binValue(current_score) : 0;
      if ( ! mono)
      {
        estimated_n_decoys -= 2 * ( fulldecoyclass_present ? cum_histograms[fulldecoyclass].binValue(current_score) : 0);
      }
      double n_targets = targetclass_present ? cum_histograms[targetclass].binValue(current_score) : 0;
      fdr.push_back(n_targets > 0 ? estimated_n_decoys / (n_targets) : 0);
    }
  }

  void XFDRAlgorithm::calc_qfdr_(const std::vector< double > &fdr, std::vector< double > &qfdr)
  {
    qfdr.resize(fdr.size());
    for (Int i = fdr.size() - 1; i >= 0; --i)
    {
      double current_fdr = fdr[i];
      double smallest_fdr = current_fdr;
      for (Int j = i; j >= 0; j--)
      {
        double fdr_to_check = fdr[j];
        if (fdr_to_check < smallest_fdr)
        {
          smallest_fdr = fdr_to_check;
        }
      }
      qfdr[i] = smallest_fdr < current_fdr ? smallest_fdr : current_fdr;
    }
  }

  void XFDRAlgorithm::findTopUniqueHits_(std::vector<PeptideIdentification>& peptide_ids)
  {
    for (PeptideIdentification& pep_id : peptide_ids)
    {
      for (PeptideHit& ph : pep_id.getHits())
      {
        String id = ph.getMetaValue("OpenPepXL:id");
        auto uid_it = std::find(this->unique_ids_.begin(), this->unique_ids_.end(), id);
        // if an ID for this candidate already exists, check if the new score is higher than the last
        if (uid_it != this->unique_ids_.end())
        {
          int index = std::distance(this->unique_ids_.begin(), uid_it);
          if (this->unique_id_scores_[index] < ph.getScore())
          {
            this->unique_id_scores_[index] = ph.getScore();
          }
        }
        else
        {
          this->unique_ids_.push_back(id);
          this->unique_id_scores_.push_back(ph.getScore());
        }
      }
    }
  }

  void XFDRAlgorithm::writeArgumentsLog_() const
  {
    //-------------------------------------------------------------
    // Printing parameters to log
    //-------------------------------------------------------------
    std::cout << std::endl;
    std::cout << ((arg_minborder_ != -1) ? "Lower bound for precursor mass error for FDR calculation is " + String(arg_minborder_) + " ppm"
                                  : "No lower bound for precursor mass error for FDR calculation") << std::endl;
    std::cout << ((arg_maxborder_ != -1) ? "Upper bound for precursor mass error for FDR calculation is " + String(arg_maxborder_) + " ppm"
                                  : "No upper bound for precursor mass error for FDR calculation") << std::endl;
    std::cout << ((arg_mindeltas_ != 0)  ? "Filtering of hits by a deltascore of " + String(arg_mindeltas_) + " is used."
                                  : "No filtering of hits by deltascore") << std::endl;
    std::cout << ((arg_minionsmatched_ > 0) ? "Filtering of hits by minimum ions matched: " + String(arg_minionsmatched_) + " is used"
                                     : "No filtering of hits by minimum ions matched.") << std::endl;
    std::cout << ((arg_minscore_ > 0) ? "Filtering of hits by minimum score of " + String(arg_minscore_) + " is used."
                               : "No filtering of hits by minimum score.") << std::endl;
    std::cout << ((arg_uniquex_) ? "Error model is generated based on unique cross-links."
                          : "Error model is generated based on redundant cross-links.") << std::endl;
    std::cout << "Bin size for cumulative histograms is " + String(arg_binsize_) << std::endl;
  }

  XFDRAlgorithm::ExitCodes XFDRAlgorithm::validateClassArguments() const
  {
    if (arg_minborder_ >= arg_maxborder_)
    {
      std::cout << "Minborder cannot be larger or equal than Maxboder!" << std::endl;
      return ILLEGAL_PARAMETERS;
    }
    return EXECUTION_OK;
  }

  String XFDRAlgorithm::getId_(const PeptideHit& ph) const
  {
    if (ph.metaValueExists("OpenPepXL:id"))
    {
      return ph.getMetaValue("OpenPepXL:id").toString();
    }

    if (ph.getMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE) == "cross-link")
    {
      return   ph.getSequence().toUnmodifiedString()
               + "-" + AASequence::fromString(ph.getMetaValue(Constants::UserParam::OPENPEPXL_BETA_SEQUENCE)).toUnmodifiedString()
               + "-a" + String(ph.getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1))
               + "-b" + String(ph.getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2));

    }
    else if (ph.getMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE) == "loop-link")
    {
      return   ph.getSequence().toUnmodifiedString()
               + "-a" + String(ph.getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1))
               + "-b" + String(ph.getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2));

    }
    else if (ph.metaValueExists(Constants::UserParam::OPENPEPXL_XL_MASS))
    {
      return   ph.getSequence().toUnmodifiedString()
               + "-" + String(ph.getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1))
               + "-" + String(ph.getMetaValue(Constants::UserParam::OPENPEPXL_XL_MASS));
    }
    else
    {
      return   ph.getSequence().toUnmodifiedString()
               + "-" + String(ph.getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1));
    }
  }

  // Names of the class parameters
  const String XFDRAlgorithm::param_decoy_string_ = "decoy_string";
  const String XFDRAlgorithm::param_minborder_ = "minborder";
  const String XFDRAlgorithm::param_maxborder_ = "maxborder";
  const String XFDRAlgorithm::param_mindeltas_ = "mindeltas";
  const String XFDRAlgorithm::param_minionsmatched_ = "minionsmatched";
  const String XFDRAlgorithm::param_uniquexl_ = "uniquexl";
  const String XFDRAlgorithm::param_no_qvalues_ = "no_qvalues";
  const String XFDRAlgorithm::param_minscore_ = "minscore";
  const String XFDRAlgorithm::param_binsize_ = "binsize";

  // Names of cross-link classes
  const String XFDRAlgorithm::crosslink_class_intradecoys_ = "intradecoys";
  const String XFDRAlgorithm::crosslink_class_fulldecoysintralinks_ = "fulldecoysintralinks";
  const String XFDRAlgorithm::crosslink_class_interdecoys_ = "interdecoys";
  const String XFDRAlgorithm::crosslink_class_fulldecoysinterlinks_ = "fulldecoysinterlinks";
  const String XFDRAlgorithm::crosslink_class_monodecoys_ = "monodecoys";
  const String XFDRAlgorithm::crosslink_class_intralinks_ = "intralinks";
  const String XFDRAlgorithm::crosslink_class_interlinks_ = "interlinks";
  const String XFDRAlgorithm::crosslink_class_monolinks_  = "monolinks";
  const String XFDRAlgorithm::crosslink_class_decoys_ = "decoys";
  const String XFDRAlgorithm::crosslink_class_targets_ = "targets";
  const String XFDRAlgorithm::crosslink_class_hybriddecoysintralinks_ = "hybriddecoysintralinks";
  const String XFDRAlgorithm::crosslink_class_hybriddecoysinterlinks_ = "hybriddecoysinterlinks";
