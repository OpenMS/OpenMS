// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include "NuXLFDR.h"
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/CONCEPT/Exception.h>
using namespace std;

namespace OpenMS
{
  NuXLFDR::NuXLFDR(size_t report_top_hits, XLGroups xl_groups):
    report_top_hits_(report_top_hits),
    xl_groups_(xl_groups)
  {
  }

  NuXLFDR::XLGroups NuXLFDR::xlGroupsFromString(const std::string& s)
  {
    if (s == "auto") return XLGroups::AUTO;
    if (s == "on") return XLGroups::ON;
    if (s == "off") return XLGroups::OFF;
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Expected 'auto', 'on' or 'off'.", s);
  }

  bool NuXLFDR::isShortXL(const PeptideHit& hit)
  {
    return static_cast<int>(hit.getMetaValue("NuXL:isXL", 0)) != 0
      && hit.getSequence().size() <= SHORT_XL_MAX_LENGTH
      && static_cast<int>(hit.getMetaValue("NuXL:NA_length", 0)) >= static_cast<int>(SHORT_XL_MIN_NA_LENGTH);
  }

  void NuXLFDR::setUpFDR_(FalseDiscoveryRate& fdr) const
  {
    Param p = fdr.getParameters();
    p.setValue("add_decoy_proteins", "true"); // we still want decoys in the result (e.g., to run percolator)
    p.setValue("add_decoy_peptides", "true");
    if (report_top_hits_ >= 2)
    {
      p.setValue("use_all_hits", "true");
    }
    fdr.setParameters(p);
  }

  bool NuXLFDR::useShortXLGroup(const PeptideIdentificationList& xl_pi, Size& short_decoys) const
  {
    short_decoys = 0;
    if (xl_groups_ == XLGroups::OFF) return false;
    const bool has_short = std::any_of(xl_pi.begin(), xl_pi.end(), [](const PeptideIdentification& pi)
      { return !pi.getHits().empty() && isShortXL(pi.getHits()[0]); });
    if (!has_short) return false;
    if (xl_groups_ == XLGroups::ON) return true;

    // AUTO: short decoys among the cross-links that pass AUTO_FDR with one FDR over all cross-links
    PeptideIdentificationList all = xl_pi;
    FalseDiscoveryRate fdr;
    setUpFDR_(fdr);
    fdr.apply(all, false);
    for (const auto& pi : all)
    {
      if (pi.getHits().empty()) continue;
      const PeptideHit& h = pi.getHits()[0];
      if (h.isDecoy() && h.getScore() <= AUTO_FDR && isShortXL(h)) { ++short_decoys; }
    }
    return short_decoys >= AUTO_MIN_DECOYS;
  }

  void NuXLFDR::QValueAtPSMLevel(PeptideIdentificationList& peptide_ids) const
  {
    FalseDiscoveryRate fdr;
    setUpFDR_(fdr);
    fdr.apply(peptide_ids, true); // also calculate Constants::UserParam::PEPTIDE_Q_VALUE
  }

  void NuXLFDR::splitIntoPeptidesAndXLs(const PeptideIdentificationList& peptide_ids, PeptideIdentificationList& pep_pi, PeptideIdentificationList& xl_pi) const
  {
    pep_pi.clear();
    xl_pi.clear();
    for (const auto & pi : peptide_ids)
    {
      vector<PeptideHit> pep_ph, xl_ph;
      for (const auto & ph : pi.getHits())
      {
         if (static_cast<int>(ph.getMetaValue("NuXL:isXL")) == 0)
         { // only add best hit
           if (pep_ph.empty() && xl_ph.empty()) pep_ph.push_back(ph); 
         }
         else
         {
           if (pep_ph.empty() && xl_ph.empty()) xl_ph.push_back(ph); 
         }
      }
      if (!pep_ph.empty()) { pep_pi.push_back(pi); pep_pi.back().setHits(pep_ph); }
      if (!xl_ph.empty()) { xl_pi.push_back(pi); xl_pi.back().setHits(xl_ph); }
    }
  }

  void NuXLFDR::mergePeptidesAndXLs(const PeptideIdentificationList& pep_pi, const PeptideIdentificationList& xl_pi, PeptideIdentificationList& peptide_ids) const
  {
    peptide_ids.clear();
    map<std::string, size_t> native_id_to_id_index;
    size_t id_index = 0;
    for (const auto & pi : pep_pi)
    {
      peptide_ids.push_back(pi);
      native_id_to_id_index[pi.getMetaValue("spectrum_reference")] = id_index;
      ++id_index;
    }

    for (const auto & pi : xl_pi)
    {
      if (!native_id_to_id_index.contains(pi.getMetaValue("spectrum_reference")))
      {
        peptide_ids.push_back(pi);
      }
      else
      {
        // spectrum already identified
        size_t index = native_id_to_id_index.at(pi.getMetaValue("spectrum_reference"));
        auto& hits = peptide_ids[index].getHits();
        for (const auto& h : pi.getHits())
        {
          hits.push_back(h);
        }
        peptide_ids[index].sort();
      }
    }
  }

  void NuXLFDR::calculatePeptideAndXLQValueAtPSMLevel(const PeptideIdentificationList& peptide_ids, PeptideIdentificationList& pep_pi, PeptideIdentificationList& xl_pi) const
  {
    FalseDiscoveryRate fdr;
    setUpFDR_(fdr);

    splitIntoPeptidesAndXLs(peptide_ids, pep_pi, xl_pi);

    // calculate PSM and peptide FDRs separately
    fdr.apply(pep_pi, true);

    Size short_decoys = 0;
    const bool short_group = useShortXLGroup(xl_pi, short_decoys);
    if (xl_groups_ == XLGroups::AUTO)
    {
      OPENMS_LOG_INFO << "Cross-link FDR: " << short_decoys << " decoys of short cross-links (at most " << SHORT_XL_MAX_LENGTH
                      << " residues, at least " << SHORT_XL_MIN_NA_LENGTH << " nucleotides) at " << AUTO_FDR * 100.0
                      << "% FDR over all cross-links: "
                      << (short_group ? "short cross-links form their own FDR group." : "one FDR over all cross-links (fewer than " + std::to_string(AUTO_MIN_DECOYS) + ").")
                      << endl;
    }
    if (!short_group)
    {
      fdr.apply(xl_pi, true);
      return;
    }

    // separate FDR groups: short cross-links and the others, merged back in the input order
    std::vector<char> is_short;
    is_short.reserve(xl_pi.size());
    PeptideIdentificationList groups[2]; // 0: other cross-links, 1: short cross-links
    for (auto& pi : xl_pi)
    {
      const bool s = !pi.getHits().empty() && isShortXL(pi.getHits()[0]);
      is_short.push_back(s);
      groups[s ? 1 : 0].push_back(std::move(pi));
    }
    xl_pi.clear();
    Size next[2] = {0, 0};
    for (int g = 0; g < 2; ++g)
    {
      fdr.apply(groups[g], true);
      for (auto& pi : groups[g])
      {
        for (auto& h : pi.getHits()) { h.setMetaValue("NuXL:fdr_group", g == 1 ? "short_XL" : "XL"); }
      }
    }
    for (char s : is_short) { xl_pi.push_back(std::move(groups[s ? 1 : 0][next[s ? 1 : 0]++])); }
  }

  void NuXLFDR::calculatePeptideAndXLQValueAndFilterAtPSMLevel(
      const vector<ProteinIdentification>& protein_ids,
      const PeptideIdentificationList& peptide_ids, 
      PeptideIdentificationList& pep_pi,
      double peptide_PSM_qvalue_threshold,
      double peptide_peptide_qvalue_threshold,
      PeptideIdentificationList& xl_pi,
      vector<double> xl_PSM_qvalue_thresholds,
      vector<double> xl_peptidelevel_qvalue_thresholds,
      const std::string& out_idxml,
      int decoy_factor) const
  {
    // calculate separate FDRs
    calculatePeptideAndXLQValueAtPSMLevel(peptide_ids, pep_pi, xl_pi);

    // correct for decoy factor
    if (decoy_factor != 1)
    {
      for (auto & pi : xl_pi)
      {
        for (auto & p : pi.getHits())
        {
          p.setScore(p.getScore() / (double)decoy_factor);
        }      
      }
    }

    // q-values might not be unique so we choose a second score to make them slightly different    
    bool svm_score_exists = false;
    for (auto & pi : xl_pi)
    {
      if (svm_score_exists) break;
      for (auto & p : pi.getHits())
      {
        if (p.metaValueExists("svm_score")) 
        {
          svm_score_exists = true;
          break;
        }
      }      
    }

    // maximum value of the "main" score besides the q-value
    double max_score = -1e10;    
    double min_score = 1e10; 
    for (auto & pi : xl_pi)
    {
      for (auto & p : pi.getHits())
      {
        if (svm_score_exists)
        {
          max_score = std::max(max_score, (double)p.getMetaValue("svm_score"));
          min_score = std::min(min_score, (double)p.getMetaValue("svm_score"));
        }
        else
        {
          max_score = std::max(max_score, (double)p.getMetaValue("NuXL:score"));
          min_score = std::min(min_score, (double)p.getMetaValue("NuXL:score"));
        }
      }
    }

    // add a very small value to q-value to break ties between same q-value but different main score
    double score_range = max_score - min_score;
    for (auto & pi : xl_pi)
    {
      for (auto & p : pi.getHits())
      {
        if (svm_score_exists)
        {
          double small_value = (1.0 - ((double)p.getMetaValue("svm_score") - min_score) / score_range) * 1e-5; // a high score will not or only slightly increase the q-value, lower scores will increase it more
          p.setScore(p.getScore() + small_value);
        }
        else 
        {
          double small_value = (1.0 - ((double)p.getMetaValue("NuXL:score") - min_score) / score_range) * 1e-5; // a high score will not or only slightly increase the q-value, lower scores will increase it more          
          p.setScore(p.getScore() + small_value);
        }
      }
    }

    // filter
    IDFilter::removeDecoyHits(xl_pi);
    IDFilter::removeDecoyHits(pep_pi);

    // filter on peptide-level q-value
    if (peptide_peptide_qvalue_threshold > 0.0 && peptide_peptide_qvalue_threshold < 1.0)
    {
      auto chechBadPeptideQValue = [&peptide_peptide_qvalue_threshold](PeptideHit& ph)->bool
      {
        return (double)ph.getMetaValue(Constants::UserParam::PEPTIDE_Q_VALUE) >= peptide_peptide_qvalue_threshold; 
      }; // of lambda

      for (auto & pid : pep_pi)
      {
        vector<PeptideHit>& phs = pid.getHits();
        phs.erase(remove_if(phs.begin(), phs.end(), chechBadPeptideQValue), phs.end());
      }
    }

    // filter on PSM-level q-value
    if (peptide_PSM_qvalue_threshold > 0.0 && peptide_PSM_qvalue_threshold < 1.0)
    {
      IDFilter::filterHitsByScore(pep_pi, peptide_PSM_qvalue_threshold); 
    }


    // store peptide PSM result
    {
      vector<ProteinIdentification> tmp_prots = protein_ids;
      IDFilter::removeUnreferencedProteins(tmp_prots, pep_pi);
      IdXMLFile().store(out_idxml + StringUtils::number(peptide_PSM_qvalue_threshold, 4) + "_peptides.idXML", tmp_prots, pep_pi);
    }

    // treat disabled filtering as 100% FDR
    std::replace(xl_PSM_qvalue_thresholds.begin(), xl_PSM_qvalue_thresholds.end(), 0.0, 1.0);
    std::sort(xl_PSM_qvalue_thresholds.begin(), xl_PSM_qvalue_thresholds.end(), greater<double>()); // important: sort by threshold (descending) to generate results by applying increasingly stringent q-value filters

    if (xl_PSM_qvalue_thresholds.size() != xl_peptidelevel_qvalue_thresholds.size())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
          "q-value list for PSMs and peptides differ in size.", 
          StringUtils::toStr(xl_PSM_qvalue_thresholds.size()) + "!=" + StringUtils::toStr(xl_peptidelevel_qvalue_thresholds.size()));
    }

    for (Size i = 0; i != xl_PSM_qvalue_thresholds.size(); ++i)
    { 
      double xlFDR = xl_PSM_qvalue_thresholds[i];
      double xl_peptidelevel_FDR = xl_peptidelevel_qvalue_thresholds[i];
      OPENMS_LOG_INFO << "Writing XL results at xl-FDR: " << xlFDR << endl;

      // filter cross-links on peptide-level q-value
      if (xl_peptidelevel_FDR > 0.0 && xl_peptidelevel_FDR < 1.0)
      {
        auto chechBadPeptideQValue = [&xl_peptidelevel_FDR](PeptideHit& ph)->bool
        {
          return (double)ph.getMetaValue(Constants::UserParam::PEPTIDE_Q_VALUE) >= xl_peptidelevel_FDR;
        }; // of lambda

        for (auto & pid : xl_pi)
        {
          vector<PeptideHit>& phs = pid.getHits();
          phs.erase(remove_if(phs.begin(), phs.end(), chechBadPeptideQValue), phs.end());
        }
      }

      // filter cross-links on PSM-level q-value
      if (xlFDR > 0.0 && xlFDR < 1.0)
      {
        IDFilter::filterHitsByScore(xl_pi, xlFDR);
      }

      vector<ProteinIdentification> tmp_prots = protein_ids;
      IDFilter::removeUnreferencedProteins(tmp_prots, xl_pi);

      // compute coverage by cross-linked peptides
      tmp_prots[0].computeCoverage(xl_pi);

      // write out XL PSM result
      IdXMLFile().store(out_idxml + StringUtils::number(xlFDR, 4) + "_XLs.idXML", tmp_prots, xl_pi);

      // write out XL protein result only for results with FDR < 10% otherwise we get to many protein associations and large memory consumption
      if (xlFDR <= 0.1)
      {
        OPENMS_LOG_INFO << "Writing XL protein results at xl-FDR: " << xlFDR << endl;
        TextFile tsv_file;
        NuXLProteinReport::annotateProteinModificationForTopHits(tmp_prots, xl_pi, tsv_file);
        tsv_file.store(out_idxml + "proteins" + StringUtils::number(xlFDR, 4) + "_XLs.tsv");
      }      
   }
  }

} // namespace

