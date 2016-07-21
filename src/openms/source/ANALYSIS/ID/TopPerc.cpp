// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer, Matthew The $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>
#include <OpenMS/ANALYSIS/ID/TopPerc.h>

using namespace std;

namespace OpenMS
{
    /*
    void TopPerc::prepareCUSTOMpin(vector<PeptideIdentification>& peptide_ids, vector<String>& user_param_features)
    {
      // Create header for the features
      string min_featureset = "SpecId, Label, ScanNr";
      StringList txt_header = ListUtils::create<String>(min_featureset);
      txt_header.insert(txt_header.end(), user_param_features.begin(), user_param_features.end() );
      txt.addLine(ListUtils::concatenate(txt_header, out_sep));

      for (vector<PeptideIdentification>::iterator it = peptide_ids.begin(); it != peptide_ids.end(); ++it)
      {
        for (vector<PeptideHit>::const_iterator hit = it->getHits().begin(); hit != it->getHits().end(); ++hit)
        {
          String scan_identifier = getScanIdentifier(it, peptide_ids.begin());
          Int scan_number = getScanNumber(scan_identifier);
          int label = 1;
          if (hit->metaValueExists("target_decoy") && String(hit->getMetaValue("target_decoy")).hasSubstring("decoy"))
          {
            label = -1;
          }

          StringList collected_feats;
          collected_feats.push_back(scan_identifier);
          collected_feats.push_back(String(label));
          collected_feats.push_back(String(scan_number));

          for (vector<String>::const_iterator feat = user_param_features.begin(); feat != user_param_features.end(); ++feat)
          {
          // Some Hits have no NumMatchedMainIons, and MeanError, etc. values. Have to ignore them!
            if (hit->metaValueExists(*feat))
            {
              collected_feats.push_back(hit->getMetaValue(*feat).toString());
            }
          }
          if (collected_feats.size() == user_param_features.size())
          { // only if all feats were present add
            txt.addLine(ListUtils::concatenate(collected_feats, out_sep));
          }
        }
      }
    }
    */
    
    void TopPerc::addMSGFFeatures(vector<PeptideIdentification>& peptide_ids, StringList& feature_set)
    {
      feature_set.push_back("MS:1002049"); // unchanged RawScore
      feature_set.push_back("MS:1002050"); // unchanged DeNovoScore
      feature_set.push_back("MSGF:ScoreRatio");
      feature_set.push_back("MSGF:Energy");
      feature_set.push_back("MSGF:lnEValue");
      feature_set.push_back("IsotopeError"); // unchanged IsotopeError
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
              
              double ln_explained_ion_current_ratio = log(hit->getMetaValue("ExplainedIonCurrentRatio").toString().toDouble() + 0.0001);           // @andsi: wtf?!
              double ln_NTerm_ion_current_ratio = log(hit->getMetaValue("NTermIonCurrentRatio").toString().toDouble() + 0.0001);           // @andsi: wtf?!
              double ln_CTerm_ion_current_ratio = log(hit->getMetaValue("CTermIonCurrentRatio").toString().toDouble() + 0.0001);           // @andsi: wtf?!
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
                LOG_WARN << "StdevErrorTop7 is NaN, setting as MeanErrorTop7 instead." << endl;
              }
              
              mean_error_top7 = rescaleFragmentFeature_(mean_error_top7, num_matched_main_ions);
              double sq_mean_error_top7 = rescaleFragmentFeature_(mean_error_top7 * mean_error_top7, num_matched_main_ions);
              stdev_error_top7 = rescaleFragmentFeature_(stdev_error_top7, num_matched_main_ions);
              hit->setMetaValue("MSGF:MeanErrorTop7", mean_error_top7);
              hit->setMetaValue("MSGF:sqMeanErrorTop7", sq_mean_error_top7);
              hit->setMetaValue("MSGF:StdevErrorTop7", stdev_error_top7);
            }
          }
        }
      }
    }
    
    void TopPerc::addXTANDEMFeatures(vector<PeptideIdentification>& peptide_ids, StringList& feature_set)
    {
      // Find out which ions are in XTandem-File and take only these as features
      StringList ion_types = ListUtils::create<String>("a,b,c,x,y,z");
      StringList ion_types_found;
      for (StringList::const_iterator ion = ion_types.begin(); ion != ion_types.end(); ++ion)
      {
        if (peptide_ids.front().getHits().front().getMetaValue(*ion + "_score").toString() != "" &&
            peptide_ids.front().getHits().front().getMetaValue(*ion + "_ions").toString() != "")
        {
          feature_set.push_back("XTANDEM:frac_ion_" + *ion);
          ion_types_found.push_back(*ion);
        }
      }
      feature_set.push_back("XTANDEM:hyperscore");
      feature_set.push_back("XTANDEM:deltascore");
      
      for (vector<PeptideIdentification>::iterator it = peptide_ids.begin(); it != peptide_ids.end(); ++it)
      {
        double hyper_score = it->getHits().front().getScore();
        double delta_score = hyper_score - it->getHits().front().getMetaValue("nextscore").toString().toDouble();
        it->getHits().front().setMetaValue("XTANDEM:hyperscore", hyper_score);
        it->getHits().front().setMetaValue("XTANDEM:deltascore", delta_score);
        
        String sequence = it->getHits().front().getSequence().toUnmodifiedString();
        int length = sequence.length();

        // Find out correct ion types and get its Values
        for (StringList::const_iterator ion = ion_types_found.begin(); ion != ion_types_found.end(); ++ion)
        {
          if (peptide_ids.front().getHits().front().getMetaValue(*ion + "_score").toString() != "" &&
              peptide_ids.front().getHits().front().getMetaValue(*ion + "_ions").toString() != "")
          {
            // recalculate ion score
            double ion_score = it->getHits().front().getMetaValue(*ion + "_ions").toString().toDouble() / length;
            it->getHits().front().setMetaValue("XTANDEM:frac_ion_" + *ion, ion_score);
          }
        }
      }
    }

    void TopPerc::addCOMETFeatures(vector<PeptideIdentification>& peptide_ids, StringList& feature_set)
    {
      feature_set.push_back("COMET:deltCn"); // recalculated deltCn = (current_XCorr - 2nd_best_XCorr) / max(current_XCorr, 1)
      feature_set.push_back("COMET:deltLCn"); // deltLCn = (current_XCorr - worst_XCorr) / max(current_XCorr, 1)
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
          if (cnt == 1) second_xcorr = xcorr;
          ++cnt;
        }
        
        for (vector<PeptideHit>::iterator hit = it->getHits().begin(); hit != it->getHits().end(); ++hit)
        {
          double xcorr = hit->getMetaValue("MS:1002252").toString().toDouble();
          double delta_cn = (xcorr - second_xcorr) / max(1.0, xcorr);
          double delta_last_cn = (xcorr - worst_xcorr) / max(1.0, xcorr);
          hit->setMetaValue("COMET:deltCn", delta_cn);
          hit->setMetaValue("COMET:deltLCn", delta_last_cn);
          
          double ln_expect = log(hit->getMetaValue("MS:1002257").toString().toDouble());
          hit->setMetaValue("COMET:lnExpect", ln_expect);
          
          double ln_num_sp = log(hit->getMetaValue("MS:1002255").toString().toDouble());
          double ln_rank_sp = log(max(1.0, hit->getMetaValue("MS:1002256").toString().toDouble()));
          hit->setMetaValue("COMET:lnNumSP", ln_num_sp);
          hit->setMetaValue("COMET:lnRankSP", ln_rank_sp);
          
          double num_matched_ions = hit->getMetaValue("MS:1002258").toString().toDouble();
          double num_total_ions = hit->getMetaValue("MS:1002259").toString().toDouble();
          double ion_frac = num_matched_ions / num_total_ions;
          hit->setMetaValue("COMET:IonFrac", ion_frac);
        }
      }
    }

    /**
    Features 1-9 Represent the Basic Feature Set

    feature abbreviation	feature description
    1. mass	Calculated monoisotopic mass of the identified peptide. Present as generic feature.
    2. charge	Precursor ion charge. Present as generic feature.
    3. mScore	Mascot score. Added in this function.
    4. dScore	Mascot score minus Mascot score of next best nonisobaric peptide hit. Added in this function.
    5. deltaM	Calculated minus observed peptide mass (in Dalton and ppm). Present as generic feature.
    6. absDeltaM	Absolute value of calculated minus observed peptide mass (in Dalton and ppm). Present as generic feature.
    7. isoDeltaM	Calculated minus observed peptide mass, isotope error corrected (in Dalton and ppm)
    8. uniquePeps	None (0), one (1), two or more (2) distinct peptide sequences match same protein. Added in this function.
    9. mc	Missed tryptic cleavages. Present as generic feature.

    Features 10-18 Represent the Extended Feature Set As Used in Mascot Percolator

    feature abbreviation	feature description
    10. totInt	Total ion intensity (log). Not available in mascot adapter.
    11. intMatchedTot	Total matched ion intensity (log). Not available in mascot adapter.
    12. relIntMatchedTot	Total matched ion intensity divided by total ion intensity. Not available in mascot adapter.
    13. binom	Peptide Score as described in ref 28. Not available in mascot adapter.
    14. fragMassError	Mean fragment mass error (in Dalton and ppm). Not available in mascot adapter.
    15. absFragMassError	Mean absolute fragment mass error (in Dalton and ppm). Not available in mascot adapter.
    16. fracIonsMatched	Fraction of calculated ions matched (per ion series). Not available in mascot adapter.
    17. seqCov	Sequence coverage of matched ions (per ion series). Not available in mascot adapter.
    18. intMatched	Matched ion intensity (per ion series). Not available in mascot adapter.
    */
    void TopPerc::addMASCOTFeatures(vector<PeptideIdentification>& peptide_ids, StringList& feature_set)
    {      
      feature_set.push_back("MS:1001171"); // unchanged mScore
      feature_set.push_back("MASCOT:delta_score"); // delta score based on mScore
      feature_set.push_back("MASCOT:uniqueToProt"); // bool: peptide unique to protein
      feature_set.push_back("MASCOT:hasMod"); // bool: has post translational modification
      
      for (vector<PeptideIdentification>::iterator it = peptide_ids.begin(); it != peptide_ids.end(); ++it)
      {
        it->sort();
        it->assignRanks();
        std::vector<PeptideHit> hits = it->getHits();
        assignDeltaScore_(hits, "MS:1001171", "MASCOT:delta_score");
        
        for (vector<PeptideHit>::iterator hit = hits.begin(); hit != hits.end(); ++hit)
        {
          bool unique_to_protein = (String(hit->getMetaValue("protein_references")) == "unique");
          bool has_mod = hit->getSequence().isModified();
          hit->setMetaValue("COMET:uniqueToProt", unique_to_protein);
          hit->setMetaValue("COMET:hasMod", has_mod);
        }
      }
    }

    void TopPerc::addCONCATSEFeatures(vector<PeptideIdentification>& peptide_ids, StringList& search_engines_used, StringList& feature_set)
    {     
      for (StringList::iterator it = search_engines_used.begin(); it != search_engines_used.end(); ++it) {
        feature_set.push_back("CONCAT:" + *it);
      }
      LOG_INFO << "Using " << ListUtils::concatenate(search_engines_used, ", ") << " as source for search engine specific features." << endl;
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

    void TopPerc::mergeMULTISEPeptideIds(vector<PeptideIdentification>& all_peptide_ids, vector<PeptideIdentification>& new_peptide_ids)
    {
      LOG_DEBUG << "creating spectrum map" << endl;
      
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
      
      for (vector<PeptideIdentification>::iterator pit = new_peptide_ids.begin(); pit != new_peptide_ids.end(); ++pit)
      {
        PeptideIdentification ins = *pit;
        //prepare for merge
        for (vector<PeptideHit>::iterator hit = ins.getHits().begin(); hit != ins.getHits().end(); ++hit)
        {
          hit->setScore(1);
        }
        ins.setScoreType("multiple");
        ins.setIdentifier("TopPerc_multiple_SE_input");
        String spectrum_reference = getScanMergeKey_(pit, new_peptide_ids.begin());
        //merge in unified map
        if (unified.find(spectrum_reference) == unified.end())
        {
          unified[spectrum_reference] = ins;
        }
        else
        {
          //find corresponding hit
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
                //There is no mutable getPeptideEvidences() accessor in PeptideHit - above will not werk, but so long:
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
                // adds up the number of hits, as the score of each separate hit is 1
                merger->setScore(merger->getScore() + hit->getScore());
                break;
              }
            }
          }
        }
      }
      
      LOG_DEBUG << "filled spectrum map" << endl;
      std::vector<PeptideIdentification> swip;
      swip.reserve(unified.size());
      LOG_DEBUG << "merging spectrum map" << endl;
      for (std::map<String,PeptideIdentification>::iterator it = unified.begin(); it != unified.end(); ++it)
      {
        swip.push_back(it->second);
      }
      all_peptide_ids.swap(swip);
      LOG_DEBUG << "Now containing " << all_peptide_ids.size() << " spectra identifications."<< endl;        
    }
    
    // references from PeptideHits to ProteinHits work with the protein accessions, so no need to update the PeptideHits
    void TopPerc::mergeMULTISEProteinIds(vector<ProteinIdentification>& all_protein_ids, vector<ProteinIdentification>& new_protein_ids)
    {      
      LOG_DEBUG << "merging search parameters" << endl;
      //care for search parameters!!
      
      String SE = new_protein_ids.front().getSearchEngine();  
      if (all_protein_ids.empty())
      {
        all_protein_ids.push_back(ProteinIdentification());
        DateTime now = DateTime::now();
        String date_string = now.getDate();
        String identifier = "TopPerc_" + date_string;
        all_protein_ids.front().setDateTime(now);
        all_protein_ids.front().setIdentifier(identifier);
        all_protein_ids.front().setSearchEngine(SE);
        LOG_DEBUG << "Setting search engine to " << SE << endl;
        all_protein_ids.front().setSearchParameters(new_protein_ids.front().getSearchParameters());
      }
      else if (all_protein_ids.front().getSearchEngine() != SE)
      {
        all_protein_ids.front().setSearchEngine("multiple");
      }
      std::vector<ProteinHit>& all_protein_hits = all_protein_ids.front().getHits();
      std::vector<ProteinHit>& new_protein_hits = new_protein_ids.front().getHits();
      
      LOG_DEBUG << "Sorting " << new_protein_hits.size() << " new ProteinHits." << endl;
      std::sort(new_protein_hits.begin(), new_protein_hits.end(), TopPerc::lq_ProteinHit());
      
      LOG_DEBUG << "Melting with " << all_protein_hits.size() << " previous ProteinHits." << endl;
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
            TopPerc::lq_ProteinHit() );
        tmp_protein_hits.resize(uni - tmp_protein_hits.begin());
        all_protein_hits.swap(tmp_protein_hits);
      }
      LOG_DEBUG << "Done with next ProteinHits." << endl;
    
      StringList keys;
      all_protein_ids.front().getSearchParameters().getKeys(keys);      
      if (!ListUtils::contains(keys, "SE:" + SE)) 
      {
        LOG_DEBUG << "Melting Parameters from " << SE << " into MetaInfo." << endl;
        
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
        all_sp.setMetaValue(SE+":fragment_mass_tolerance_ppm",sp.fragment_mass_tolerance_ppm);
        all_sp.setMetaValue(SE+":precursor_tolerance",sp.precursor_tolerance);
        all_sp.setMetaValue(SE+":precursor_mass_tolerance_ppm",sp.precursor_mass_tolerance_ppm);
        all_sp.setMetaValue(SE+":digestion_enzyme",sp.digestion_enzyme.getName());
        
        LOG_DEBUG << "Done with next Parameters." << endl;
        all_protein_ids.front().setSearchParameters(all_sp);
      }
      
      StringList all_primary_ms_run_path = all_protein_ids.front().getPrimaryMSRunPath();
      StringList new_primary_ms_run_path = new_protein_ids.front().getPrimaryMSRunPath();
      all_primary_ms_run_path.insert(all_primary_ms_run_path.end(), new_primary_ms_run_path.begin(), new_primary_ms_run_path.end());
      all_protein_ids.front().setPrimaryMSRunPath(all_primary_ms_run_path);
      LOG_DEBUG << "New primary run paths: " << ListUtils::concatenate(new_primary_ms_run_path,",") << endl;
      LOG_DEBUG << "All primary run paths: " << ListUtils::concatenate(all_primary_ms_run_path,",") << endl;
      
      LOG_DEBUG << "All merging finished." << endl;
    }
    
    void TopPerc::concatMULTISEPeptideIds(vector<PeptideIdentification>& all_peptide_ids, vector<PeptideIdentification>& new_peptide_ids, String search_engine)
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

    void TopPerc::addMULTISEFeatures(vector<PeptideIdentification>& peptide_ids, StringList& search_engines_used, StringList& feature_set)
    {
      if (ListUtils::contains(search_engines_used, "MS-GF+"))
      {
        feature_set.push_back("MS:1002049");  // rawscore
        feature_set.push_back("MS:1002053");  // evalue
      }
      if (ListUtils::contains(search_engines_used, "Mascot"))
      {
        feature_set.push_back("MS:1001171");
        feature_set.push_back("EValue");
      }
      if (ListUtils::contains(search_engines_used, "Comet"))
      {
        feature_set.push_back("MS:1002252");  //xcorr
        feature_set.push_back("MS:1002257");  //evalue
      }
      if (ListUtils::contains(search_engines_used, "XTandem"))
      {
        //TODO: create XTandem score
        //feature_set.push_back("XTandem_score");
        feature_set.push_back("E-Value");
      }
      //feature_set.push_back("MULTI:ionFrac");
      //feature_set.push_back("MULTI:numHits"); // this is not informative if we only keep PSMs with hits for all search engines
      
      LOG_INFO << "Using " << ListUtils::concatenate(search_engines_used, ", ") << " as source for search engine specific features." << endl;

      // get all the feature values
      for (vector<PeptideIdentification>::iterator it = peptide_ids.begin(); it != peptide_ids.end(); ++it)
      {
        it->sort();
        it->assignRanks();
        for (vector<PeptideHit>::iterator hit = it->getHits().begin(); hit != it->getHits().end(); ++hit)
        {
          //double ion_frac = hit->getMetaValue("matched_intensity").toString().toDouble() / hit->getMetaValue("sum_intensity").toString().toDouble();  // also consider "matched_ion_number"/"peak_number"
          //hit->setMetaValue("MULTI:ionFrac", ion_frac);
          
          int num_hits = hit->getScore();
          hit->setMetaValue("MULTI:numHits", num_hits);
        }
      }
    }
    
    // Function adapted from MsgfplusReader in Percolator converter
    double TopPerc::rescaleFragmentFeature_(double featureValue, int NumMatchedMainIons)
    {
      // Rescale the fragment features to penalize features calculated by few ions
      int numMatchedIonLimit = 7;
      int numerator = (1 + numMatchedIonLimit) * (1 + numMatchedIonLimit);
      int denominator = (1 + (min)(NumMatchedMainIons, numMatchedIonLimit)) * (1 + (min)(NumMatchedMainIons, numMatchedIonLimit));
      return featureValue * ((double)numerator / denominator);
    }
    
    void TopPerc::assignDeltaScore_(vector<PeptideHit>& hits, String score_ref, String output_ref)
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
    
    bool TopPerc::hasMHCEnd_(String peptide)
    {
      bool suf = false;
      static const string arr[] = {"A", "F", "I", "K", "M", "L", "R", "W", "V"};
      vector<string> mhcends (arr, arr + sizeof(arr) / sizeof(arr[0]) );
      for (std::vector<string>::iterator eit = mhcends.begin(); eit != mhcends.end(); ++eit)
      {
        if (peptide.hasSuffix(string(*eit)))
        {
          suf = true;
          break;
        }
      }
      return suf;
    }
    
    // TODO: check if this is consistent for all search engines. MSGF+ and X!Tandem have been checked.
    String TopPerc::getScanMergeKey_(vector<PeptideIdentification>::iterator it, vector<PeptideIdentification>::iterator start)
    {
      // MSGF+ uses this field, is empty if not specified
      String scan_identifier = it->getMetaValue("spectrum_reference");
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
          LOG_WARN << "no known spectrum identifiers, using index [1,n] - use at own risk." << endl;
        }
      }
      
      Int scan_number = 0;
      StringList fields = ListUtils::create<String>(scan_identifier);
      for (StringList::const_iterator it = fields.begin(); it != fields.end(); ++it)
      {
        // if scan number is not available, use the scan index
        Size idx = 0;
        if ((idx = it->find("index=")) != string::npos) 
        {
          scan_number = it->substr(idx + 6).toInt();
        }
      }
      return String(scan_number);
    }
    

}
