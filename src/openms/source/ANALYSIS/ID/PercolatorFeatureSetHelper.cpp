// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer, Matthew The $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/PercolatorFeatureSetHelper.h>

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Constants.h>

#include <boost/lexical_cast.hpp>

using namespace std;

namespace OpenMS
{    
    void PercolatorFeatureSetHelper::addMSGFFeatures(vector<PeptideIdentification>& peptide_ids, StringList& feature_set)
    {
      feature_set.push_back("MS:1002049"); // unchanged RawScore
      feature_set.push_back("MS:1002050"); // unchanged DeNovoScore
      feature_set.push_back("MSGF:ScoreRatio");
      feature_set.push_back("MSGF:Energy");
      feature_set.push_back("MSGF:lnEValue");
      feature_set.push_back(Constants::UserParam::ISOTOPE_ERROR);
      feature_set.push_back("MSGF:lnExplainedIonCurrentRatio");
      feature_set.push_back("MSGF:lnNTermIonCurrentRatio");
      feature_set.push_back("MSGF:lnCTermIonCurrentRatio");
      feature_set.push_back("MSGF:lnMS2IonCurrent");
      feature_set.push_back("MSGF:MeanErrorTop7");
      feature_set.push_back("MSGF:sqMeanErrorTop7");
      feature_set.push_back("MSGF:StdevErrorTop7");
      
      for (vector<PeptideIdentification>::iterator it = peptide_ids.begin(); it != peptide_ids.end(); ++it)
      {
        for (vector<PeptideHit>::iterator hit = it->getHits().begin(); hit != it->getHits().end(); ++hit)
        {
          // Some Hits have no NumMatchedMainIons, and MeanError, etc. values. Have to ignore them!
          if (hit->metaValueExists("NumMatchedMainIons"))
          {
            // only take features from first ranked entries and only with meanerrortop7 != 0.0
            if (hit->getMetaValue("MeanErrorTop7").toString().toDouble() != 0.0)
            {
              double raw_score = hit->getMetaValue("MS:1002049").toString().toDouble();
              double denovo_score = hit->getMetaValue("MS:1002050").toString().toDouble();
              
              double energy = denovo_score - raw_score;
              double score_ratio = raw_score * 10000;
              if (denovo_score > 0)
              {
                score_ratio = (raw_score / denovo_score);
              }
              hit->setMetaValue("MSGF:ScoreRatio", score_ratio);
              hit->setMetaValue("MSGF:Energy", energy);
              
              double ln_eval = -log(hit->getMetaValue("MS:1002053").toString().toDouble());
              hit->setMetaValue("MSGF:lnEValue", ln_eval);
              
              double ln_explained_ion_current_ratio = log(hit->getMetaValue("ExplainedIonCurrentRatio").toString().toDouble() + 0.0001);
              double ln_NTerm_ion_current_ratio = log(hit->getMetaValue("NTermIonCurrentRatio").toString().toDouble() + 0.0001);
              double ln_CTerm_ion_current_ratio = log(hit->getMetaValue("CTermIonCurrentRatio").toString().toDouble() + 0.0001);
              hit->setMetaValue("MSGF:lnExplainedIonCurrentRatio", ln_explained_ion_current_ratio);
              hit->setMetaValue("MSGF:lnNTermIonCurrentRatio", ln_NTerm_ion_current_ratio);
              hit->setMetaValue("MSGF:lnCTermIonCurrentRatio", ln_CTerm_ion_current_ratio);
              
              double ln_MS2_ion_current = log(hit->getMetaValue("MS2IonCurrent").toString().toDouble());
              hit->setMetaValue("MSGF:lnMS2IonCurrent", ln_MS2_ion_current);
              
              double mean_error_top7 = hit->getMetaValue("MeanErrorTop7").toString().toDouble();
              int num_matched_main_ions =  hit->getMetaValue("NumMatchedMainIons").toString().toInt();

              double stdev_error_top7 = 0.0;
              if (hit->getMetaValue("StdevErrorTop7").toString() != "NaN")
              {
                stdev_error_top7 = hit->getMetaValue("StdevErrorTop7").toString().toDouble();
                if (stdev_error_top7 == 0.0)
                {
                  stdev_error_top7 = mean_error_top7;
                }
              }
              else
              {
                stdev_error_top7 = mean_error_top7;
                OPENMS_LOG_WARN << "StdevErrorTop7 is NaN, setting as MeanErrorTop7 instead." << endl;
              }
              
              mean_error_top7 = rescaleFragmentFeature_(mean_error_top7, num_matched_main_ions);
              double sq_mean_error_top7 = rescaleFragmentFeature_(mean_error_top7 * mean_error_top7, num_matched_main_ions);
              stdev_error_top7 = rescaleFragmentFeature_(stdev_error_top7, num_matched_main_ions);
              hit->setMetaValue("MSGF:MeanErrorTop7", mean_error_top7);
              hit->setMetaValue("MSGF:sqMeanErrorTop7", sq_mean_error_top7);
              hit->setMetaValue("MSGF:StdevErrorTop7", stdev_error_top7);
            }
          }
          else OPENMS_LOG_WARN << "MS-GF+ PSM with missing NumMatchedMainIons skipped." << endl;
        }
      }
    }
    
    void PercolatorFeatureSetHelper::addXTANDEMFeatures(vector<PeptideIdentification>& peptide_ids, StringList& feature_set)
    {
      //TODO annotate isotope error in Adapter and add here as well.
      // Find out which ions are in XTandem-File and take only these as features
      StringList ion_types = ListUtils::create<String>("a,b,c,x,y,z");
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
      
      for (vector<PeptideIdentification>::iterator it = peptide_ids.begin(); it != peptide_ids.end(); ++it)
      {
        double hyper_score = it->getHits().front().getScore();
        double delta_score = hyper_score - it->getHits().front().getMetaValue("nextscore").toString().toDouble();
        it->getHits().front().setMetaValue("XTANDEM:deltascore", delta_score);
        
        String sequence = it->getHits().front().getSequence().toUnmodifiedString();
        int length = sequence.length();

        // Find out correct ion types and get its Values
        for (StringList::const_iterator ion = ion_types_found.begin(); ion != ion_types_found.end(); ++ion)
        {
          if (!peptide_ids.front().getHits().front().getMetaValue(*ion + "_score").toString().empty() &&
              !peptide_ids.front().getHits().front().getMetaValue(*ion + "_ions").toString().empty())
          {
            // recalculate ion score
            double ion_score = it->getHits().front().getMetaValue(*ion + "_ions").toString().toDouble() / length;
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

    void PercolatorFeatureSetHelper::addCOMETFeatures(vector<PeptideIdentification>& peptide_ids, StringList& feature_set)
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
      
      for (vector<PeptideIdentification>::iterator it = peptide_ids.begin(); it != peptide_ids.end(); ++it)
      {
        double worst_xcorr = 0, second_xcorr = 0;
        Int cnt = 0;
        for (vector<PeptideHit>::iterator hit = it->getHits().begin(); hit != it->getHits().end(); ++hit)
        {
          double xcorr = hit->getMetaValue("MS:1002252").toString().toDouble();
          worst_xcorr = xcorr;
          if (cnt == 1) { second_xcorr = xcorr; }
          ++cnt;
        }
        
        for (vector<PeptideHit>::iterator hit = it->getHits().begin(); hit != it->getHits().end(); ++hit)
        {

          double xcorr = hit->getMetaValue("MS:1002252").toString().toDouble();

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
          
          double ln_expect = log(hit->getMetaValue("MS:1002257").toString().toDouble());
          hit->setMetaValue("COMET:lnExpect", ln_expect);

          if (!hit->metaValueExists("COMET:lnNumSP"))
          {
            double ln_num_sp;   
            if (hit->metaValueExists("num_matched_peptides"))
            {
              double num_sp = hit->getMetaValue("num_matched_peptides").toString().toDouble();
              ln_num_sp = log(max(1.0, num_sp));  // if recorded, one can be safely assumed
            }
            else // fallback TODO: remove?
            {
              ln_num_sp = hit->getMetaValue("MS:1002255").toString().toDouble();
            }  
            hit->setMetaValue("COMET:lnNumSP", ln_num_sp);
          }

          if (!hit->metaValueExists("COMET:lnRankSP"))
          {          
            double ln_rank_sp = log(max(1.0, hit->getMetaValue("MS:1002256").toString().toDouble()));
            hit->setMetaValue("COMET:lnRankSP", ln_rank_sp);
          }

          if (!hit->metaValueExists("COMET:IonFrac"))
          {
            double num_matched_ions = hit->getMetaValue("MS:1002258").toString().toDouble();
            double num_total_ions = hit->getMetaValue("MS:1002259").toString().toDouble();
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
    void PercolatorFeatureSetHelper::addMASCOTFeatures(vector<PeptideIdentification>& peptide_ids, StringList& feature_set)
    {      
      feature_set.push_back("MS:1001171"); // unchanged mScore
      feature_set.push_back("MASCOT:delta_score"); // delta score based on mScore
      feature_set.push_back("MASCOT:hasMod"); // bool: has post translational modification
      
      for (vector<PeptideIdentification>::iterator it = peptide_ids.begin(); it != peptide_ids.end(); ++it)
      {
        it->sort();
        it->assignRanks();
        std::vector<PeptideHit> hits = it->getHits();
        assignDeltaScore_(hits, "MS:1001171", "MASCOT:delta_score");
        
        for (vector<PeptideHit>::iterator hit = hits.begin(); hit != hits.end(); ++hit)
        {
          bool has_mod = hit->getSequence().isModified();
          hit->setMetaValue("MASCOT:hasMod", has_mod);
        }
      }
    }

    void PercolatorFeatureSetHelper::addCONCATSEFeatures(vector<PeptideIdentification>& peptide_ids, StringList& search_engines_used, StringList& feature_set)
    {     
      for (StringList::iterator it = search_engines_used.begin(); it != search_engines_used.end(); ++it) {
        feature_set.push_back("CONCAT:" + *it);
      }
      OPENMS_LOG_INFO << "Using " << ListUtils::concatenate(search_engines_used, ", ") << " as source for search engine specific features." << endl;
      feature_set.push_back("CONCAT:lnEvalue");
      feature_set.push_back("CONCAT:deltaLnEvalue");
      
      // feature values have been set in concatMULTISEids
      for (vector<PeptideIdentification>::iterator it = peptide_ids.begin(); it != peptide_ids.end(); ++it)
      {
        it->sort();
        it->assignRanks();
        assignDeltaScore_(it->getHits(), "CONCAT:lnEvalue", "CONCAT:deltaLnEvalue");
      }
    }

    void PercolatorFeatureSetHelper::mergeMULTISEPeptideIds(vector<PeptideIdentification>& all_peptide_ids, vector<PeptideIdentification>& new_peptide_ids, const String& search_engine)
    {
      OPENMS_LOG_DEBUG << "creating spectrum map" << endl;
      
      std::map<String,PeptideIdentification> unified;
      //setup map of merge characteristics per spectrum
      for (vector<PeptideIdentification>::iterator pit = all_peptide_ids.begin(); pit != all_peptide_ids.end(); ++pit)
      {
        PeptideIdentification ins = *pit;
        ins.setScoreType("multiple");
        ins.setIdentifier("TopPerc_multiple_SE_input");
        String spectrum_reference = getScanMergeKey_(pit, all_peptide_ids.begin());
        unified[spectrum_reference] = ins;
      }
      OPENMS_LOG_DEBUG << "filled spectrum map with previously observed PSM: " << unified.size() << endl;

      int nc = 0;
      int mc = 0;
      OPENMS_LOG_DEBUG << "About to merge in:" << new_peptide_ids.size() << "PSMs." << endl;
      for (vector<PeptideIdentification>::iterator pit = new_peptide_ids.begin(); pit != new_peptide_ids.end(); ++pit)
      {
        PeptideIdentification ins = *pit;
        String st = pit->getScoreType();
        //prepare for merge
        for (vector<PeptideHit>::iterator hit = ins.getHits().begin(); hit != ins.getHits().end(); ++hit)
        {
          // keep the hit score as meta value
          if (st == "MS-GF:RawScore")
          {
            st = "MS:1002049";
          }
          else if (st == "XTandem")
          {
            st = "MS:1001331";
          }
          else if (st == "Mascot")
          {
            st = "MS:1001171";
          }
          else if ((st == "expect" && search_engine == "Comet" )|| st == "Comet")
          {
            st = "MS:1002257";
          }

          if (!hit->metaValueExists(st))
          {
            hit->setMetaValue(st, hit->getScore());
          }

          hit->setScore(1);  // new 'multiple' score is just the number of times identified by different SE

          // rename ambiguous meta value names to PSI cv term ids
          if (search_engine == "MS-GF+" && hit->metaValueExists("EValue"))  // MS-GF should have all values as PSI cv terms available anyway
          {
            hit->setMetaValue("MS:1002053", hit->getMetaValue("EValue"));
          }
          if (search_engine == "Mascot" && hit->metaValueExists("EValue"))
          {
            hit->setMetaValue("MS:1001172", hit->getMetaValue("EValue"));
          }
          if (search_engine == "Comet" && hit->metaValueExists("xcorr"))
          {
            hit->setMetaValue("MS:1002252", hit->getMetaValue("xcorr"));
          }
          if (search_engine == "XTandem" && hit->metaValueExists("E-Value"))
          {
            hit->setMetaValue("MS:1001330", hit->getMetaValue("E-Value"));
          }
        }
        ins.setScoreType("multiple");
        ins.setIdentifier("TopPerc_multiple_SE_input");
        String spectrum_reference = getScanMergeKey_(pit, new_peptide_ids.begin());
        //merge in unified map
        // insert newly identified spectra (PeptideIdentifications) or ..
        if (unified.find(spectrum_reference) == unified.end())
        {
          unified[spectrum_reference] = ins;
          ++nc;
        }
        // .. add PSMs to already identified spectrum
        else
        {
          //find corresponding hit (i.e. sequences must match)
          for (vector<PeptideHit>::iterator hit = ins.getHits().begin(); hit != ins.getHits().end(); ++hit)
          {
            for (vector<PeptideHit>::iterator merger = unified[spectrum_reference].getHits().begin(); merger != unified[spectrum_reference].getHits().end(); ++merger)
            {
              if (hit->getSequence()==merger->getSequence())
              {
                //care for peptide evidences!! set would be okay if checked for same search db in parameters,
//                  vector<PeptideEvidence> pev;
//                  pev.reserve(max(hit->getPeptideEvidences().size(),merger->getPeptideEvidences().size()));
//                  std::vector<ProteinHit>::iterator uni;
//                  std::sort(merger->getPeptideEvidences().begin(),merger->getPeptideEvidences().end(), TopPerc::lq_PeptideEvidence);
//                  std::sort(hit->getPeptideEvidences().begin(),hit->getPeptideEvidences().end(), TopPerc::lq_PeptideEvidence);
//                  uni = std::set_union(swop.front().getHits().begin(), swop.front().getHits().end(),
//                                       it->front().getHits().begin(),it->front().getHits().end(), pev.begin(),
//                                       TopPerc::lq_PeptideEvidence);
//                  pev.resize(uni-pev.begin());
//                  merger->setPeptideEvidences(pev);
                //There is no mutable getPeptideEvidences() accessor in PeptideHit - above will not work, but so long:
                //Implying PeptideIndexer was applied (with the same search db each) will care for that all PeptideEvidences from two hits with equal AASequence are the same

                //merge meta values
                StringList keys;
                hit->getKeys(keys);
                for (StringList::const_iterator kt = keys.begin(); kt != keys.end(); ++kt)
                {
                  if (!merger->metaValueExists(*kt))
                  {
                    merger->setMetaValue(*kt, hit->getMetaValue(*kt));
                  }
                }

                // adds up the number of hits, as the score of each separate (new) hit is 1
                merger->setScore(merger->getScore() + hit->getScore());
                ++mc;
                break;
              }
            }
          }
        }
      }
      
      OPENMS_LOG_DEBUG << "filled spectrum map" << endl;
      std::vector<PeptideIdentification> swip;
      swip.reserve(unified.size());
      OPENMS_LOG_DEBUG << "merging spectrum map" << endl;
      for (std::map<String,PeptideIdentification>::iterator it = unified.begin(); it != unified.end(); ++it)
      {
        swip.push_back(it->second);
      }
      all_peptide_ids.swap(swip);
      OPENMS_LOG_DEBUG << "Now containing " << all_peptide_ids.size() << " spectra identifications, new: " << nc << " merged in: " << mc << endl;
    }
    
    // references from PeptideHits to ProteinHits work with the protein accessions, so no need to update the PeptideHits
    void PercolatorFeatureSetHelper::mergeMULTISEProteinIds(vector<ProteinIdentification>& all_protein_ids, vector<ProteinIdentification>& new_protein_ids)
    {      
      OPENMS_LOG_DEBUG << "merging search parameters" << endl;
      
      String SE = new_protein_ids.front().getSearchEngine();  
      if (all_protein_ids.empty())
      {
        all_protein_ids.emplace_back();
        DateTime now = DateTime::now();
        String date_string = now.getDate();
        String identifier = "TopPerc_" + date_string;
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
    
    void PercolatorFeatureSetHelper::concatMULTISEPeptideIds(vector<PeptideIdentification>& all_peptide_ids, vector<PeptideIdentification>& new_peptide_ids, const String& search_engine)
    {      
      for (vector<PeptideIdentification>::iterator pit = new_peptide_ids.begin(); pit != new_peptide_ids.end(); ++pit)
      {
        for (vector<PeptideHit>::iterator hit = pit->getHits().begin(); hit != pit->getHits().end(); ++hit)
        {
          double evalue = 1000.0;
          if (search_engine == "MS-GF+")
          {
            hit->setMetaValue("CONCAT:" + search_engine, hit->getMetaValue("MS:1002049"));  // rawscore
            evalue = hit->getMetaValue("MS:1002049").toString().toDouble();  // evalue
          }
          if (search_engine == "Mascot")
          {
            hit->setMetaValue("CONCAT:" + search_engine, hit->getMetaValue("MS:1001171")); // mscore
            evalue = hit->getMetaValue("EValue").toString().toDouble();
          }
          if (search_engine == "Comet")
          {
            hit->setMetaValue("CONCAT:" + search_engine, hit->getMetaValue("MS:1002252"));  // xcorr
            evalue = hit->getMetaValue("MS:1002257").toString().toDouble();
          }
          if (search_engine == "XTandem")
          {
            hit->setMetaValue("CONCAT:" + search_engine, hit->getMetaValue("XTandem_score"));  // xtandem score
            evalue = hit->getMetaValue("E-Value").toString().toDouble();
          }
          hit->setMetaValue("CONCAT:lnEvalue", log(evalue));  // log(evalue)
        }
      }
      all_peptide_ids.insert(all_peptide_ids.end(), new_peptide_ids.begin(), new_peptide_ids.end());
    }

    void PercolatorFeatureSetHelper::addMULTISEFeatures(vector<PeptideIdentification>& peptide_ids, StringList& search_engines_used, StringList& feature_set, bool complete_only, bool limits_imputation)
    {
      map<String,vector<double> > extremals;  // will have as keys the below SE cv terms
      vector<String> max_better, min_better;
      // This is the minimum set for each SE that should be available in all openms id files in one way or another
      if (ListUtils::contains(search_engines_used, "MS-GF+"))
      {
        feature_set.push_back("MS:1002049");  // rawscore
        feature_set.push_back("MS:1002053");  // evalue
        max_better.emplace_back("MS:1002049");  // higher is better - start high, get min
        min_better.emplace_back("MS:1002053");  // lower is better - start low, get max
      }
      if (ListUtils::contains(search_engines_used, "Mascot"))
      {
        feature_set.push_back("MS:1001171");  // score aka Mascot
        feature_set.push_back("MS:1001172");  // evalue aka EValue
        max_better.emplace_back("MS:1001171");  // higher is better - start high, get min
        min_better.emplace_back("MS:1001172");  // lower is better - start low, get max
      }
      if (ListUtils::contains(search_engines_used, "Comet"))
      {
        feature_set.push_back("MS:1002252");  // xcorr
        feature_set.push_back("MS:1002257");  // evalue
        max_better.emplace_back("MS:1002252");  // higher is better - start high, get min
        min_better.emplace_back("MS:1002257");  // lower is better - start low, get max
      }
      if (ListUtils::contains(search_engines_used, "XTandem"))
      {
        feature_set.push_back("MS:1001331");  // hyperscore aka XTandem
        feature_set.push_back("MS:1001330");  // evalue aka E-Value
        max_better.emplace_back("MS:1001331");  // higher is better - start high, get min
        min_better.emplace_back("MS:1001330");  // lower is better - start low, get max
      }
      //feature_set.push_back("MULTI:ionFrac");
      //feature_set.push_back("MULTI:numHits"); // this is not informative if we only keep PSMs with hits for all search engines
      
      OPENMS_LOG_INFO << "Using " << ListUtils::concatenate(search_engines_used, ", ") << " as source for search engine specific features." << endl;

      // get all the feature values
      if (!complete_only)
      {
        for (vector<PeptideIdentification>::iterator it = peptide_ids.begin(); it != peptide_ids.end(); ++it)
        {
          for (vector<PeptideHit>::iterator hit = it->getHits().begin(); hit != it->getHits().end(); ++hit)
          {
            for (StringList::iterator feats = feature_set.begin(); feats != feature_set.end(); ++feats)
            {
              if (hit->metaValueExists(*feats))
              {
                // TODO raise issue: MS-GF raw score values are sometimes registered as string DataValues and henceforth casted defectively
                if (hit->getMetaValue(*feats).valueType() == DataValue::STRING_VALUE)
                {
                  String recast = hit->getMetaValue(*feats);
                  double d = boost::lexical_cast<double>(recast);
                  OPENMS_LOG_DEBUG << "recast: "
                            << recast << " "
                            << double(hit->getMetaValue(*feats)) << "* ";
                  hit->setMetaValue(*feats,d);
                  OPENMS_LOG_DEBUG << hit->getMetaValue(*feats).valueType() << " "
                            << hit->getMetaValue(*feats)
                            << endl;
                }
                extremals[*feats].push_back(hit->getMetaValue(*feats));
              }
            }
          }
        }
        // TODO : add optional manual extremal values settings for 'data imputation' instead of min/max or numeric_limits value
        for (vector<String>::iterator maxbt = max_better.begin(); maxbt != max_better.end(); ++maxbt)
        {
          map<String,vector<double> >::iterator fi = extremals.find(*maxbt);
          if (fi != extremals.end())
          {
            vector<double>::iterator mymax = min_element(fi->second.begin(), fi->second.end());
            iter_swap(fi->second.begin(), mymax);
            if (limits_imputation)
            {
              fi->second.front() = -std::numeric_limits<float>::max();
            }
          }
        }
        for (vector<String>::iterator minbt = min_better.begin(); minbt != min_better.end(); ++minbt)
        {
          map<String,vector<double> >::iterator fi = extremals.find(*minbt);
          if (fi != extremals.end())
          {
            vector<double>::iterator mymin = max_element(fi->second.begin(), fi->second.end());
            iter_swap(fi->second.begin(), mymin);
            if (limits_imputation)
            {
              fi->second.front() = std::numeric_limits<float>::max();
            }
          }
        }
      }

      size_t sum_removed = 0;
      size_t imputed_values = 0;
      size_t observed_values = 0;
      size_t complete_spectra = 0;
      size_t incomplete_spectra = 0;

      OPENMS_LOG_DEBUG << "Looking for minimum feature set:" << ListUtils::concatenate(feature_set, ", ") << "." << endl;

      for (vector<PeptideIdentification>::iterator pi = peptide_ids.begin(); pi != peptide_ids.end(); ++pi)
      {
        pi->sort();
        pi->assignRanks();
        vector<vector<PeptideHit>::iterator> incompletes;

        size_t imputed_back = imputed_values;
        for (vector<PeptideHit>::iterator hit = pi->getHits().begin(); hit != pi->getHits().end(); ++hit)
        {
          //double ion_frac = hit->getMetaValue("matched_intensity").toString().toDouble() / hit->getMetaValue("sum_intensity").toString().toDouble();  // also consider "matched_ion_number"/"peak_number"
          //hit->setMetaValue("MULTI:ionFrac", ion_frac);
          
          for (StringList::iterator feats = feature_set.begin(); feats != feature_set.end(); ++feats)
          {
            if (complete_only && !hit->metaValueExists(*feats))
            {
              incompletes.push_back(hit);  // mark for removal
              break;
            }
            else if (!hit->metaValueExists(*feats))
            {
              hit->setMetaValue(*feats, extremals[*feats].front());  // imputation
              ++imputed_values;
            }
            else
            {
              ++observed_values;
            }
          }
          int num_hits = hit->getScore();
          hit->setMetaValue("MULTI:numHits", num_hits);
        }
        if (complete_only)
        {
          // remove incompletes
          for (vector<vector<PeptideHit>::iterator>::reverse_iterator rit = incompletes.rbegin(); rit != incompletes.rend(); ++rit)
          {
            pi->getHits().erase(*rit);
          }
          sum_removed += incompletes.size();
        }
        if (!incompletes.empty() || imputed_back < imputed_values)
          ++incomplete_spectra;
        else
          ++complete_spectra;
      }
      if (sum_removed > 0)
      {
        OPENMS_LOG_WARN << "Removed " << sum_removed << " incomplete cases of PSMs." << endl;
      }
      if (imputed_values > 0)
      {
        OPENMS_LOG_WARN << "Imputed " << imputed_values << " of " << observed_values+imputed_values
                 << " missing values. ("
                 << imputed_values*100.0/(imputed_values+observed_values)
                 << "%)" << endl;
        OPENMS_LOG_WARN << "Affected " << incomplete_spectra << " of " << incomplete_spectra+complete_spectra
                 << " spectra. ("
                 << incomplete_spectra*100.0/(incomplete_spectra+complete_spectra)
                 << "%)" << endl;
      }
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
    
    void PercolatorFeatureSetHelper::assignDeltaScore_(vector<PeptideHit>& hits, const String& score_ref, const String& output_ref)
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
        
    // TODO: this is code redundancy to PercolatorAdapter
    // TODO: in case of merged idXML files from fractions and/or replicates make sure that you also consider the file origin
    //  this is usually stored in the map_index MetaValue of a PeptideIdentification (PSM) object.
    String PercolatorFeatureSetHelper::getScanMergeKey_(vector<PeptideIdentification>::iterator it, vector<PeptideIdentification>::iterator start)
    {
      // MSGF+ uses this field, is empty if not specified
      String scan_identifier = it->getSpectrumReference();
      if (scan_identifier.empty())
      {
        // XTandem uses this (integer) field
        // these ids are 1-based in contrast to the index which is 0-based, so subtract 1.
        if (it->metaValueExists("spectrum_id") && !it->getMetaValue("spectrum_id").toString().empty())
        {
          scan_identifier = "index=" + String(it->getMetaValue("spectrum_id").toString().toInt() - 1);
        }
        else
        {
          scan_identifier = "index=" + String(it - start + 1);
          OPENMS_LOG_WARN << "no known spectrum identifiers, using index [1,n] - use at own risk." << endl;
        }
      }
      
      Int scan = 0;
      StringList fields = ListUtils::create<String>(scan_identifier);
      for (StringList::const_iterator it = fields.begin(); it != fields.end(); ++it)
      {
        Size idx = 0;
        if ((idx = it->find("scan=")) != String::npos)
        {
          scan = it->substr(idx + 5).toInt();
          break;
        }  // only if scan number is not available, use the scan index
        else if ((idx = it->find("index=")) != String::npos)
        {
          scan = it->substr(idx + 6).toInt();
        }
      }
      return String(scan);
    }
}
