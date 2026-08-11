// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer, Matthew The $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/PercolatorFeatureSetHelper.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Constants.h>

#include <set>

using namespace std;

namespace OpenMS
{    
    void PercolatorFeatureSetHelper::addMSGFFeatures(PeptideIdentificationList& peptide_ids, StringList& feature_set)
    {
      // MSGF+ does not always produce all scores so we focus on the main ones 
      // and make sure they are present and initalized
      feature_set.push_back("MS:1002049"); // MS-GF:RawScore
      feature_set.push_back("MS:1002050"); // MS-GF:DeNovoScore
      feature_set.push_back("MS:1002052"); // MS-GF:SpecEValue
      feature_set.push_back("MS:1002053"); // MS-GF:EValue
      feature_set.push_back(Constants::UserParam::ISOTOPE_ERROR);
      for (auto& p : peptide_ids)
      {
        for (auto& h : p.getHits())
        {
          if (!h.metaValueExists("MS:1002049")) h.setMetaValue("MS:1002049", 0.0);
          if (!h.metaValueExists("MS:1002050")) h.setMetaValue("MS:1002050", 0.0);
          if (!h.metaValueExists("MS:1002052")) h.setMetaValue("MS:1002052", 0.0);
          if (!h.metaValueExists("MS:1002053")) h.setMetaValue("MS:1002053", 0.0);
        }
      }
    }
    
    void PercolatorFeatureSetHelper::addXTANDEMFeatures(PeptideIdentificationList& peptide_ids, StringList& feature_set)
    {
      //TODO annotate isotope error in Adapter and add here as well.
      // Find out which ions are in XTandem-File and take only these as features
      StringList ion_types = ListUtils::create<std::string>("a,b,c,x,y,z");
      StringList ion_types_found;
      for (StringList::const_iterator ion = ion_types.begin(); ion != ion_types.end(); ++ion)
      {
        if (!peptide_ids.front().getHits().front().getMetaValue(*ion + "_score").toString().empty() &&
            !peptide_ids.front().getHits().front().getMetaValue(*ion + "_ions").toString().empty())
        {
          feature_set.push_back("XTANDEM:frac_ion_" + *ion);
          ion_types_found.push_back(*ion);
        }
      }

      feature_set.push_back("XTANDEM:deltascore");
      
      for (PeptideIdentificationList::iterator it = peptide_ids.begin(); it != peptide_ids.end(); ++it)
      {
        double hyper_score = it->getHits().front().getScore();
        double delta_score = hyper_score - StringUtils::toDouble(it->getHits().front().getMetaValue("nextscore").toString());
        it->getHits().front().setMetaValue("XTANDEM:deltascore", delta_score);
        
        std::string sequence = it->getHits().front().getSequence().toUnmodifiedString();
        int length = sequence.length();

        // Find out correct ion types and get its Values
        for (StringList::const_iterator ion = ion_types_found.begin(); ion != ion_types_found.end(); ++ion)
        {
          if (!peptide_ids.front().getHits().front().getMetaValue(*ion + "_score").toString().empty() &&
              !peptide_ids.front().getHits().front().getMetaValue(*ion + "_ions").toString().empty())
          {
            // recalculate ion score
            double ion_score = StringUtils::toDouble(it->getHits().front().getMetaValue(*ion + "_ions").toString()) / length;
            it->getHits().front().setMetaValue("XTANDEM:frac_ion_" + *ion, ion_score);
          }
        }
      }
    }

    void PercolatorFeatureSetHelper::addMSFRAGGERFeatures(StringList& feature_set)
    {
      feature_set.push_back("MS:1001330"); // expect_score
      feature_set.push_back("hyperscore");
      feature_set.push_back("nextscore");
      feature_set.push_back(Constants::UserParam::ISOTOPE_ERROR);
    }

    void PercolatorFeatureSetHelper::addANDESFeatures(PeptideIdentificationList& peptide_ids, StringList& feature_set)
    {
      // andes annotates each PSM with its own numeric features as MetaValues prefixed with "andes:".
      // We collect all numeric "andes:"-prefixed MetaValues present on the hits as Percolator features,
      // so the andes search engine stays the single source of truth for its feature set (currently 47).
      std::set<std::string> andes_features;
      for (const auto& p : peptide_ids)
      {
        for (const auto& h : p.getHits())
        {
          StringList keys;
          h.getKeys(keys);
          for (const std::string& key : keys)
          {
            if (key.compare(0, 6, "andes:") != 0)
            {
              continue;
            }
            const DataValue::DataType vt = h.getMetaValue(key).valueType();
            if (vt == DataValue::DOUBLE_VALUE || vt == DataValue::INT_VALUE)
            {
              andes_features.insert(key);
            }
          }
        }
      }
      feature_set.insert(feature_set.end(), andes_features.begin(), andes_features.end());
    }

    void PercolatorFeatureSetHelper::addCOMETFeatures(PeptideIdentificationList& peptide_ids, StringList& feature_set)
    {

      feature_set.push_back(Constants::UserParam::ISOTOPE_ERROR);
      feature_set.push_back("COMET:deltaCn"); // recalculated deltaCn = (current_XCorr - 2nd_best_XCorr) / max(current_XCorr, 1)
      feature_set.push_back("COMET:deltaLCn"); // deltaLCn = (current_XCorr - worst_XCorr) / max(current_XCorr, 1)
      feature_set.push_back("COMET:lnExpect"); // log(E-value)
      feature_set.push_back("MS:1002252"); // unchanged XCorr
      feature_set.push_back("MS:1002255"); // unchanged Sp = number of candidate peptides
      feature_set.push_back("COMET:lnNumSP"); // log(number of candidate peptides)
      feature_set.push_back("COMET:lnRankSP"); // log(rank based on Sp score)
      feature_set.push_back("COMET:IonFrac"); // matched_ions / total_ions
      
      for (PeptideIdentificationList::iterator it = peptide_ids.begin(); it != peptide_ids.end(); ++it)
      {
        double worst_xcorr = 0, second_xcorr = 0;
        Int cnt = 0;
        for (vector<PeptideHit>::iterator hit = it->getHits().begin(); hit != it->getHits().end(); ++hit)
        {
          double xcorr = StringUtils::toDouble(hit->getMetaValue("MS:1002252").toString());
          worst_xcorr = xcorr;
          if (cnt == 1) { second_xcorr = xcorr; }
          ++cnt;
        }
        
        for (vector<PeptideHit>::iterator hit = it->getHits().begin(); hit != it->getHits().end(); ++hit)
        {

          double xcorr = StringUtils::toDouble(hit->getMetaValue("MS:1002252").toString());

          if (!hit->metaValueExists("COMET:deltaCn"))
          {
            double delta_cn = (xcorr - second_xcorr) / max(1.0, xcorr);
            hit->setMetaValue("COMET:deltaCn", delta_cn);
          }

          if (!hit->metaValueExists("COMET:deltaLCn"))
          {
            double delta_last_cn = (xcorr - worst_xcorr) / max(1.0, xcorr);
            hit->setMetaValue("COMET:deltaLCn", delta_last_cn);
          }
          
          double ln_expect = log(StringUtils::toDouble(hit->getMetaValue("MS:1002257").toString()));
          hit->setMetaValue("COMET:lnExpect", ln_expect);

          if (!hit->metaValueExists("COMET:lnNumSP"))
          {
            double ln_num_sp;   
            if (hit->metaValueExists("num_matched_peptides"))
            {
              double num_sp = StringUtils::toDouble(hit->getMetaValue("num_matched_peptides").toString());
              ln_num_sp = log(max(1.0, num_sp));  // if recorded, one can be safely assumed
            }
            else // fallback TODO: remove?
            {
              ln_num_sp = StringUtils::toDouble(hit->getMetaValue("MS:1002255").toString());
            }  
            hit->setMetaValue("COMET:lnNumSP", ln_num_sp);
          }

          if (!hit->metaValueExists("COMET:lnRankSP"))
          {          
            double ln_rank_sp = log(max(1.0, StringUtils::toDouble(hit->getMetaValue("MS:1002256").toString())));
            hit->setMetaValue("COMET:lnRankSP", ln_rank_sp);
          }

          if (!hit->metaValueExists("COMET:IonFrac"))
          {
            double num_matched_ions = StringUtils::toDouble(hit->getMetaValue("MS:1002258").toString());
            double num_total_ions = StringUtils::toDouble(hit->getMetaValue("MS:1002259").toString());
            double ion_frac = num_matched_ions / num_total_ions;
            hit->setMetaValue("COMET:IonFrac", ion_frac);
          }
        }
      }
    }

    /**
    Features 1-9 Represent the Basic Feature Set

    feature abbreviation  feature description
    1. mass        Calculated monoisotopic mass of the identified peptide. Present as generic feature.
    2. charge      Precursor ion charge. Present as generic feature.
    3. mScore      Mascot score. Added in this function.
    4. dScore      Mascot score minus Mascot score of next best non isobaric peptide hit. Added in this function.
    5. deltaM      Calculated minus observed peptide mass (in Dalton and ppm). Present as generic feature.
    6. absDeltaM   Absolute value of calculated minus observed peptide mass (in Dalton and ppm). Present as generic feature.
    7. isoDeltaM   Calculated minus observed peptide mass, isotope error corrected (in Dalton and ppm)
    8. uniquePeps  None (0), one (1), two or more (2) distinct peptide sequences match same protein. Added in this function.
    9. mc          Missed tryptic cleavages. Present as generic feature.

    Features 10-18 Represent the Extended Feature Set As Used in Mascot Percolator

    feature abbreviation  feature description
    10. totInt            Total ion intensity (log). Not available in mascot adapter.
    11. intMatchedTot     Total matched ion intensity (log). Not available in mascot adapter.
    12. relIntMatchedTot  Total matched ion intensity divided by total ion intensity. Not available in mascot adapter.
    13. binom             Peptide Score as described in ref 28. Not available in mascot adapter.
    14. fragMassError     Mean fragment mass error (in Dalton and ppm). Not available in mascot adapter.
    15. absFragMassError  Mean absolute fragment mass error (in Dalton and ppm). Not available in mascot adapter.
    16. fracIonsMatched   Fraction of calculated ions matched (per ion series). Not available in mascot adapter.
    17. seqCov            Sequence coverage of matched ions (per ion series). Not available in mascot adapter.
    18. intMatched        Matched ion intensity (per ion series). Not available in mascot adapter.
    */
    void PercolatorFeatureSetHelper::addMASCOTFeatures(PeptideIdentificationList& peptide_ids, StringList& feature_set)
    {      
      feature_set.push_back("MS:1001171"); // unchanged mScore
      feature_set.push_back("MASCOT:delta_score"); // delta score based on mScore
      feature_set.push_back("MASCOT:hasMod"); // bool: has post translational modification
      
      for (PeptideIdentificationList::iterator it = peptide_ids.begin(); it != peptide_ids.end(); ++it)
      {
        it->sort();
        std::vector<PeptideHit> hits = it->getHits();
        assignDeltaScore_(hits, "MS:1001171", "MASCOT:delta_score");
        
        for (vector<PeptideHit>::iterator hit = hits.begin(); hit != hits.end(); ++hit)
        {
          bool has_mod = hit->getSequence().isModified();
          hit->setMetaValue("MASCOT:hasMod", has_mod);
        }
      }
    }

    // references from PeptideHits to ProteinHits work with the protein accessions, so no need to update the PeptideHits
    void PercolatorFeatureSetHelper::mergeMULTISEProteinIds(vector<ProteinIdentification>& all_protein_ids, vector<ProteinIdentification>& new_protein_ids)
    {      
      OPENMS_LOG_DEBUG << "merging search parameters" << endl;
      
      std::string SE = new_protein_ids.front().getSearchEngine();  
      if (all_protein_ids.empty())
      {
        all_protein_ids.emplace_back();
        DateTime now = DateTime::now();
        std::string date_string = now.getDate();
        std::string identifier = "TopPerc_" + date_string;
        all_protein_ids.front().setDateTime(now);
        all_protein_ids.front().setIdentifier(identifier);
        all_protein_ids.front().setSearchEngine(SE);
        OPENMS_LOG_DEBUG << "Setting search engine to " << SE << endl;
        all_protein_ids.front().setSearchParameters(new_protein_ids.front().getSearchParameters());
      }
      else if (all_protein_ids.front().getSearchEngine() != SE)
      {
        all_protein_ids.front().setSearchEngine("multiple");
      }
      std::vector<ProteinHit>& all_protein_hits = all_protein_ids.front().getHits();
      std::vector<ProteinHit>& new_protein_hits = new_protein_ids.front().getHits();
      
      OPENMS_LOG_DEBUG << "Sorting " << new_protein_hits.size() << " new ProteinHits." << endl;
      std::sort(new_protein_hits.begin(), new_protein_hits.end(), PercolatorFeatureSetHelper::lq_ProteinHit());
      
      OPENMS_LOG_DEBUG << "Melting with " << all_protein_hits.size() << " previous ProteinHits." << endl;
      if (all_protein_hits.empty())
      {
        all_protein_hits.swap(new_protein_hits);
      }
      else
      {
        std::vector<ProteinHit> tmp_protein_hits(new_protein_hits.size() + all_protein_hits.size());
        std::vector<ProteinHit>::iterator uni = set_union(
            all_protein_hits.begin(), all_protein_hits.end(),
            new_protein_hits.begin(), new_protein_hits.end(), tmp_protein_hits.begin(),
            PercolatorFeatureSetHelper::lq_ProteinHit() );
        tmp_protein_hits.resize(uni - tmp_protein_hits.begin());
        all_protein_hits.swap(tmp_protein_hits);
      }
      OPENMS_LOG_DEBUG << "Done with next ProteinHits." << endl;
    
      StringList keys;
      all_protein_ids.front().getSearchParameters().getKeys(keys);      
      if (!ListUtils::contains(keys, "SE:" + SE)) 
      {
        OPENMS_LOG_DEBUG << "Melting Parameters from " << SE << " into MetaInfo." << endl;
        
        //insert into MetaInfo as SE:param
        ProteinIdentification::SearchParameters sp = new_protein_ids.front().getSearchParameters();
        ProteinIdentification::SearchParameters all_sp = all_protein_ids.front().getSearchParameters();
        all_sp.setMetaValue("SE:"+SE,new_protein_ids.front().getSearchEngineVersion());
        all_sp.setMetaValue(SE+":db",sp.db);
        all_sp.setMetaValue(SE+":db_version",sp.db_version);
        all_sp.setMetaValue(SE+":taxonomy",sp.taxonomy);
        all_sp.setMetaValue(SE+":charges",sp.charges);
        all_sp.setMetaValue(SE+":fixed_modifications",ListUtils::concatenate(sp.fixed_modifications, ","));
        all_sp.setMetaValue(SE+":variable_modifications",ListUtils::concatenate(sp.variable_modifications, ","));
        all_sp.setMetaValue(SE+":missed_cleavages",sp.missed_cleavages);
        all_sp.setMetaValue(SE+":fragment_mass_tolerance",sp.fragment_mass_tolerance);
        all_sp.setMetaValue(SE+":fragment_mass_tolerance_unit", sp.fragment_mass_tolerance_ppm ? "ppm" : "Da");
        all_sp.setMetaValue(SE+":precursor_mass_tolerance",sp.precursor_mass_tolerance);
        all_sp.setMetaValue(SE+":precursor_mass_tolerance_unit", sp.precursor_mass_tolerance_ppm ? "ppm" : "Da");
        all_sp.setMetaValue(SE+":digestion_enzyme",sp.digestion_enzyme.getName());
        all_sp.setMetaValue(SE+":enzyme_term_specificity",sp.enzyme_term_specificity);
        //TODO maybe add all the files in file origin that were searched with this SE. then we can do a lookup later
        // for every PepID based on its file_origin, with which SEs and settings it was identified.
        
        OPENMS_LOG_DEBUG << "Done with next Parameters." << endl;
        all_protein_ids.front().setSearchParameters(all_sp);
      }
      OPENMS_LOG_DEBUG << "Merging primaryMSRunPaths." << endl;
      try
      {
        StringList all_primary_ms_run_path;
        all_protein_ids.front().getPrimaryMSRunPath(all_primary_ms_run_path);
        StringList new_primary_ms_run_path;
        new_protein_ids.front().getPrimaryMSRunPath(new_primary_ms_run_path);

        all_primary_ms_run_path.insert(all_primary_ms_run_path.end(), new_primary_ms_run_path.begin(), new_primary_ms_run_path.end());
        all_protein_ids.front().setPrimaryMSRunPath(all_primary_ms_run_path);

        OPENMS_LOG_DEBUG << "New primary run paths: " << ListUtils::concatenate(new_primary_ms_run_path,",") << endl;
        OPENMS_LOG_DEBUG << "All primary run paths: " << ListUtils::concatenate(all_primary_ms_run_path,",") << endl;
      }
      catch (Exception::BaseException& e)
      {
        OPENMS_LOG_DEBUG << "faulty primary MS run path: " << endl;
        Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, e.what(), "");
      }
      OPENMS_LOG_DEBUG << "Merging for this file finished." << endl;
    }
    
    void PercolatorFeatureSetHelper::checkExtraFeatures(const vector<PeptideHit>& psms, StringList& extra_features)
    {
      set<StringList::iterator> unavail;
      for (vector<PeptideHit>::const_iterator hit = psms.begin(); hit != psms.end(); ++hit)
      {
        for (StringList::iterator ef = extra_features.begin(); ef != extra_features.end(); ++ef)
        {
          if (!hit->metaValueExists(*ef))
          {
            unavail.insert(ef);
          }
        }
      }
      for (set<StringList::iterator>::reverse_iterator rit = unavail.rbegin(); rit != unavail.rend(); ++rit)
      {
        OPENMS_LOG_WARN << "A extra_feature requested (" << *(*rit) << ") was not available - removed." << endl;
        extra_features.erase(*rit);
      }
    }

    
    // Function adapted from MSGFPlusReader in Percolator converter
    double PercolatorFeatureSetHelper::rescaleFragmentFeature_(double featureValue, int NumMatchedMainIons)
    {
      // Rescale the fragment features to penalize features calculated by few ions
      int numMatchedIonLimit = 7;
      int numerator = (1 + numMatchedIonLimit) * (1 + numMatchedIonLimit);
      int denominator = (1 + (min)(NumMatchedMainIons, numMatchedIonLimit)) * (1 + (min)(NumMatchedMainIons, numMatchedIonLimit));
      return featureValue * ((double)numerator / denominator);
    }
    
    void PercolatorFeatureSetHelper::assignDeltaScore_(vector<PeptideHit>& hits, const std::string& score_ref, const std::string& output_ref)
    {
      if (!hits.empty())
      {
        vector<PeptideHit>::iterator prev = hits.begin();
        double prev_score = double(prev->getMetaValue(score_ref));
        for (vector<PeptideHit>::iterator hit = hits.begin()+1; hit != hits.end(); ++hit)
        {
          double cur_score = double(hit->getMetaValue(score_ref));
          double value = prev_score - cur_score;
          prev->setMetaValue(output_ref, value);
          prev = hit;
        }
        (hits.end()-1)->setMetaValue(output_ref, 0.0); //if last hit or only one hit
      }
    }
}
