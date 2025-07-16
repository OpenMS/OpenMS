// Copyright (c) 2002-2024, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/AhoCorasickAmbiguous.h>
#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHExtenderAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTaggerAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTnTAlgorithm.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

namespace OpenMS
{
inline const int min_tag_count = 2;
inline const int max_hit_count = 20;
FLASHTnTAlgorithm::FLASHTnTAlgorithm(): DefaultParamHandler("FLASHTnTAlgorithm"), ProgressLogger()
{
  setDefaultParams_();
}

FLASHTnTAlgorithm::FLASHTnTAlgorithm(const FLASHTnTAlgorithm& other): DefaultParamHandler(other), ProgressLogger(other)
{
}

FLASHTnTAlgorithm& FLASHTnTAlgorithm::operator=(const FLASHTnTAlgorithm& rhs)
{
  if (this == &rhs) return *this;

  DefaultParamHandler::operator=(rhs);
  return *this;
}

void FLASHTnTAlgorithm::setDefaultParams_()
{
  defaults_.setValue("prsm_fdr", 1.0, "Specifies the PrSM-level FDR.");
  defaults_.setMinFloat("prsm_fdr", 0.0);

  defaults_.setValue("pro_fdr", 1.0, "Specifies the proteoform-level FDR.");
  defaults_.setMinFloat("pro_fdr", 0.0);

  defaults_.setValue("only_single_hit", "false", "Allows only a single hit per spectrum.");
  defaults_.setValidStrings("only_single_hit", {"true", "false"});

  defaults_.setValue("discard_underdetermined", "false",
                     "Discards underdetermined proteoform IDs (e.g., those without exact precursor masses or start/end positions).");
  defaults_.setValidStrings("discard_underdetermined", {"true", "false"});

  defaults_.setValue("keep_decoy", "false", "Retains decoy hits in the results.");
  defaults_.setValidStrings("keep_decoy", {"true", "false"});

  defaults_.setValue("ion_type", std::vector<std::string> {"b", "y"}, "Specifies ion types to consider.");
  defaults_.setValidStrings("ion_type", {"b", "c", "a", "y", "z", "x", "zp1", "zp2"});

  auto tparam = FLASHTaggerAlgorithm().getDefaults();
  tparam.remove("ion_type");
  defaults_.insert("tag:", tparam);
  auto eparam = FLASHExtenderAlgorithm().getDefaults();
  eparam.remove("ion_type");

  defaults_.insert("ex:", eparam);
  defaultsToParam_();
}

void FLASHTnTAlgorithm::updateMembers_()
{
  tagger_param_ = param_.copy("tag:", true);
  tagger_param_.setValue("ion_type", param_.getValue("ion_type"));
  extender_param_ = param_.copy("ex:", true);
  extender_param_.setValue("ion_type", param_.getValue("ion_type"));
  prsm_fdr_ = param_.getValue("prsm_fdr");
  keep_decoy_ = param_.getValue("keep_decoy").toString() == "true";
  keep_underdetermined_ = param_.getValue("discard_underdetermined").toString() == "false";
  multiple_hits_per_spec_ = param_.getValue("only_single_hit").toString() == "false";
}

bool FLASHTnTAlgorithm::areConsistent_(const ProteinHit& a, const ProteinHit& b, double tol) const
{
  double mass1 = a.getMetaValue("Mass");
  double mass2 = b.getMetaValue("Mass");
  if (mass1 * mass2 < 0) return false;
  if (mass1 > 0 && mass2 > 0 && std::abs(mass1 - mass2) > std::max(mass1, mass2) * tol / 1e6 * 2) return false;

  int sp1 = a.getMetaValue("StartPosition");
  int ep1 = a.getMetaValue("EndPosition");

  int sp2 = b.getMetaValue("StartPosition");
  int ep2 = b.getMetaValue("EndPosition");

  // double rt1 = a.getMetaValue("RT");
  // double rt2 = b.getMetaValue("RT");
  // if (std::abs(rt1 - rt2) < 30.0) return true; // if rts are within 30 sec, true

  if (a.metaValueExists("mod_masses") && b.metaValueExists("mod_masses"))
  {
    std::vector<double> mod_masses1 = a.getMetaValue("mod_masses");
    std::vector<double> mod_masses2 = b.getMetaValue("mod_masses");
    if (mod_masses1.size() != mod_masses2.size()) return false;
    else
    {
      std::vector<int> mod_starts1 = a.getMetaValue("ModificationStarts");
      std::vector<int> mod_starts2 = b.getMetaValue("ModificationStarts");
      std::vector<int> mod_ends1 = a.getMetaValue("ModificationEnds");
      std::vector<int> mod_ends2 = b.getMetaValue("ModificationEnds");

      for (Size i = 0; i < mod_masses1.size(); i++)
      {
        if (std::abs(mod_masses1[i] - mod_masses2[i]) > std::max(mod_masses1[i], mod_masses2[i]) * tol / 1e6 * 2) { return false; }
        else if (mod_starts1[i] > mod_ends2[i] || mod_starts2[i] > mod_ends1[i]) { return false; }
      }
    }
    return true;
  }
  else if (! a.metaValueExists("mod_masses") && ! a.metaValueExists("mod_masses"))
  {
    if (sp1 == sp2 && ep1 == ep2) return true;
  }
  return false;
}

void FLASHTnTAlgorithm::markRepresentativeProteoformHits_(double tol)
{
  std::sort(proteoform_hits_.begin(), proteoform_hits_.end(),
            [](const ProteinHit& left, const ProteinHit& right) { return left.getScore() > right.getScore(); });
  std::map<String, std::vector<ProteinHit>> proteoform_map;
  for (auto& hit : proteoform_hits_)
  {
    const auto& acc = hit.getAccession();

    if (proteoform_map.find(acc) != proteoform_map.end())
    {
      bool skip = false;
      for (const auto& hit2 : proteoform_map[acc])
      {
        if (areConsistent_(hit, hit2, tol))
        {
          skip = true;
          break;
        }
      }
      if (skip) continue;
    }
    hit.setMetaValue("Representative", true);

    proteoform_map[acc].push_back(hit);
  }
}


void FLASHTnTAlgorithm::run(const MSExperiment& map, const std::vector<FASTAFile::FASTAEntry>& original_fasta_entry)
{
  setLogType(CMD);

  int max_mod_cntr = extender_param_.getValue("max_mod_count");
  double max_mod_mass = max_mod_cntr * (double)extender_param_.getValue("max_mod_mass") + 1.0;
  std::map<double, std::vector<ResidueModification>> mod_map;
  const auto inst = ModificationsDB::getInstance();             // give this from outside ...
  std::map<String, std::vector<Size>> tag_to_protein_indices;   // tag to protein index in fasta
  std::map<String, std::vector<Size>> tag_to_protein_positions; // tag to protein position
  std::map<int, std::vector<Size>> scan_to_tag_indices;         // scan number to tag indices
  std::vector<std::string> tag_seqs, protein_seqs;
  std::vector<DeconvolvedSpectrum> dspecs;

  std::vector<FASTAFile::FASTAEntry> fasta_entry;
  fasta_entry.reserve(original_fasta_entry.size());
  for (const auto& fe : original_fasta_entry)
  {
    auto seq = fe.sequence;
    std::replace(seq.begin(), seq.end(), 'I', 'L');
    fasta_entry.emplace_back(fe.identifier, fe.description, seq);
  }

  {
    protein_seqs.reserve(fasta_entry.size());
    double taget_count = 0;
    double decoy_count = 0;
    for (const auto& fe : fasta_entry)
    {
      if (fe.identifier.hasPrefix("DECOY")) { decoy_count++; }
      else { taget_count++; }
      protein_seqs.push_back(fe.sequence);
    }

    decoy_factor_ = decoy_count / taget_count;
  }

  std::vector<String> mod_strs;
  inst->getAllSearchModifications(mod_strs);
  for (int i = 0; i < mod_strs.size(); i++)
  {
    const auto mod = *inst->getModification(mod_strs[i]);
    if (std::abs(mod.getDiffMonoMass()) > max_mod_mass) continue;
    mod_map[mod.getDiffMonoMass()].push_back(mod);
  }
  double precursor_tol = -1, tol;

  startProgress(0, (SignedSize)map.size(), "Finding sequence tags ...");

  for (int index = 0; index < map.size(); index++)
  {
    auto spec = map[index];
    nextProgress();

    if (spec.size() < 5) continue;

    int scan = FLASHDeconvAlgorithm::getScanNumber(map, index);

    DeconvolvedSpectrum dspec(scan);
    dspec.setOriginalSpectrum(spec);
    String deconv_meta_str = spec.getMetaValue("DeconvMassInfo").toString();

    int tol_loc_s = deconv_meta_str.find("tol=") + 4;
    int tol_loc_e = deconv_meta_str.find(";", tol_loc_s);

    tol = stod(deconv_meta_str.substr(tol_loc_s, tol_loc_e - tol_loc_s));
    if (spec.getMSLevel() == 1 && precursor_tol < 0)
    {
      precursor_tol = tol;
      // continue;
    }

    int q_loc_s = deconv_meta_str.find("qscore=") + 7;
    int q_loc_e = deconv_meta_str.find(";", q_loc_s);
    auto q_str = deconv_meta_str.substr(q_loc_s, q_loc_e - q_loc_s);
    Size pos = 0;
    std::vector<double> qscores;
    while (true)
    {
      Size pos_t = q_str.find(",", pos);
      if (pos_t == String::npos) break;
      auto token = q_str.substr(pos, pos_t - pos);
      qscores.push_back(stod(token));
      pos = pos_t + 1;
    }

    int s_loc_s = deconv_meta_str.find("snr=") + 4;
    int s_loc_e = deconv_meta_str.find(";", s_loc_s);
    auto s_str = deconv_meta_str.substr(s_loc_s, s_loc_e - s_loc_s);
    pos = 0;
    std::vector<float> snrs;
    while (true)
    {
      Size pos_t = s_str.find(",", pos);
      if (pos_t == String::npos) break;
      auto token = s_str.substr(pos, pos_t - pos);
      snrs.push_back(stof(token));
      pos = pos_t + 1;
    }

    int s_loc_pre_s = deconv_meta_str.find("precursorscan=") + 14;
    int s_loc_pre_e = deconv_meta_str.find(";", s_loc_pre_s);
    int precursor_scan = stoi(deconv_meta_str.substr(s_loc_pre_s, s_loc_pre_e - s_loc_pre_s));

    if (precursor_scan > 0)
    {
      int s_loc_prem_s = deconv_meta_str.find("precursormass=") + 14;
      int s_loc_prem_e = deconv_meta_str.find(";", s_loc_prem_s);
      double precursor_mass = stod(deconv_meta_str.substr(s_loc_prem_s, s_loc_prem_e - s_loc_prem_s));
      PeakGroup pg;
      pg.setMonoisotopicMass(precursor_mass);
      if (deconv_meta_str.hasSubstring("precursorscore="))
      {
        int s_loc_preq_s = deconv_meta_str.find("precursorscore=") + 15;
        int s_loc_preq_e = deconv_meta_str.find(";", s_loc_preq_s);
        double precursor_qscore = stod(deconv_meta_str.substr(s_loc_preq_s, s_loc_preq_e - s_loc_preq_s));
        pg.setQscore2D(precursor_qscore);
      }
      if (deconv_meta_str.hasSubstring("precursorSNR="))
      {
        int s_loc_pres_s = deconv_meta_str.find("precursorSNR=") + 13;
        int s_loc_pres_e = deconv_meta_str.find(";", s_loc_pres_s);
        double precursor_snr = stod(deconv_meta_str.substr(s_loc_pres_s, s_loc_pres_e - s_loc_pres_s));
        pg.setSNR(precursor_snr);
      }

      dspec.setPrecursorPeakGroup(pg);
    }
    for (Size i = 0; i < spec.size(); i++)
    {
      PeakGroup peak;
      peak.setQscore(qscores[i]);
      peak.setSNR(snrs[i]);
      peak.setMonoisotopicMass(spec[i].getMZ());
      peak.setScanNumber(scan);
      dspec.push_back(peak);
    }

    if (dspec.size() > max_node_cntr_)
    {
      dspec.sortByQscore();
      std::vector<PeakGroup> filtered_peaks;
      filtered_peaks.reserve(max_node_cntr_);
      for (const auto& pg : dspec)
      {
        if (filtered_peaks.size() >= max_node_cntr_) break;
        filtered_peaks.push_back(pg);
      }
      dspec.setPeakGroups(filtered_peaks);
    }

    dspec.sort();
    dspecs.push_back(dspec);
    FLASHTaggerAlgorithm tagger;
    // Run tagger
    tagger.setParameters(tagger_param_);
    tagger.run(dspec, tol);
    tagger.fillTags(tags_);
  }

  endProgress();

  for (Size i = 0; i < tags_.size(); i++)
  {
    const auto& seq = tags_[i].getUppercaseSequence();
    scan_to_tag_indices[tags_[i].getScan()].push_back(i);
    tags_[i].setIndex((int)i);
    if (tag_to_protein_indices.find(seq) != tag_to_protein_indices.end()) continue;
    tag_to_protein_indices[seq] = std::vector<Size>();
    tag_to_protein_positions[seq] = std::vector<Size>();
    tag_seqs.push_back(seq);
  } // now we have all the tags

  // tag_seqs protein_seqs

  ACTrie ac_trie(0, 0); // no ambiguities, no mismatches
  ac_trie.addNeedlesAndCompress(tag_seqs);
  OpenMS::ACTrieState ac_state;

  startProgress(0, (SignedSize)protein_seqs.size(), "Searching tags against database ...");
  // #pragma omp parallel for default(none) shared(protein_seqs, ac_state, ac_trie)
  for (int i = 0; i < protein_seqs.size(); i++) // // Natural19WithoutI
  {
    nextProgress();
    ac_state.setQuery(protein_seqs[i]);
    ac_trie.getAllHits(ac_state);

    for (const auto& h : ac_state.hits)
    {
      String seq = tag_seqs[h.needle_index];
      tag_to_protein_indices[seq].push_back(i);
      tag_to_protein_positions[seq].push_back(h.query_pos);
    }
  }

  endProgress();
  startProgress(0, (SignedSize)dspecs.size(), "Running candidate protein filtration and extension algorithm ...");

  for (const auto& dspec : dspecs)
  {
    nextProgress();
    //if (dspec.getScanNumber() > 993)break;
    if (scan_to_tag_indices.find(dspec.getScanNumber()) == scan_to_tag_indices.end()) continue;
    const auto& tag_indices = scan_to_tag_indices[dspec.getScanNumber()];
    std::vector<ProteinHit> hits;
    std::map<Size, std::map<Size, int>> pi_to_pos;         // protein index to covered positions
    std::map<Size, std::set<Size>> pi_to_tag_indices; // protein index to tag indices

    for (const auto& tag_index : tag_indices)
    {
      const auto& tag = tags_[tag_index];
      const auto& seq = tag.getUppercaseSequence();

      const auto& pis = tag_to_protein_indices[seq];
      const auto& pps = tag_to_protein_positions[seq];
      for (Size i = 0; i < pis.size(); i++)
      {
        pi_to_tag_indices[pis[i]].insert(tag_index);
        for (int j = 0; j <= seq.length(); j++) //
        {
          pi_to_pos[pis[i]][pps[i] + j] = 1;//tag.getScore(j);
        }
      }
    }

    for (const auto& [pi, pp] : pi_to_pos)
    {
      const auto& fe = fasta_entry[pi];
      double covered_pos = (double)pp.size();
      double coverage = (double)covered_pos / (double)fe.sequence.length();

      std::vector<int> positions;
      std::transform(pp.begin(), pp.end(), std::back_inserter(positions),
                     [](const std::pair<const Size, int>& p) {
                       return p.first;
                     });
      //if (positions.size() < min_tag_count) continue;
      ProteinHit hit(coverage, 0, fe.identifier, fe.sequence); //
      hit.setDescription(fe.description);
      hit.setMetaValue("Scan", dspec.getScanNumber());
      hit.setMetaValue("MatchedAA", covered_pos);
      hit.setCoverage(coverage);
      hit.setMetaValue("TagIndices", std::vector<int>(pi_to_tag_indices[pi].begin(), pi_to_tag_indices[pi].end()));
      hit.setMetaValue("TagPositions", positions);
      hit.setMetaValue("FastaIndex", pi);
      hits.push_back(hit);
    }

    if (hits.empty()) continue;

        std::sort(hits.begin(), hits.end(), [](const ProteinHit& left, const ProteinHit& right) {
          return left.getScore() == right.getScore() ? (left.getCoverage() == right.getCoverage() ? (left.getDescription() > right.getDescription())
                                                                                                  : (left.getCoverage() > right.getCoverage()))
                                                     : (left.getScore() > right.getScore());
        });
   // std::cout<<hits.size()<<std::endl; // 예전 runmathcing 으로 좀 돌아갈 듯.. tag 매칭도 필요함. 그래도 빠르긴 하다만....
    //if (hits.size() > 5 * max_hit_count) { hits.resize(5 * max_hit_count); }
        if (hits.size() > 100*max_hit_count) { hits.resize(100*max_hit_count); }

    FLASHTaggerAlgorithm::runMatching(hits, dspec, max_mod_mass);

    std::sort(hits.begin(), hits.end(), [](const ProteinHit& left, const ProteinHit& right) {
      return left.getScore() == right.getScore() ? (left.getCoverage() == right.getCoverage() ? (left.getDescription() > right.getDescription())
                                                                                              : (left.getCoverage() > right.getCoverage()))
                                                 : (left.getScore() > right.getScore());
    });

    if (hits.size() > max_hit_count) { hits.resize(max_hit_count); }

    FLASHExtenderAlgorithm extender;
    extender.setParameters(extender_param_);
    extender.setModificationMap(mod_map);
    extender.run(hits, tags_, dspec, tol, multiple_hits_per_spec_);
    extender.fillProteoforms(proteoform_hits_);
  }

  endProgress();

  markRepresentativeProteoformHits_(precursor_tol);
  std::sort(proteoform_hits_.begin(), proteoform_hits_.end(), [](const ProteinHit& left, const ProteinHit& right) {
    return left.getScore() == right.getScore() ? (left.getCoverage() == right.getCoverage() ? (left.getMetaValue("Scan") > right.getMetaValue("Scan"))
                                                                                            : (left.getCoverage() > right.getCoverage()))
                                               : (left.getScore() > right.getScore());
  });

  if (decoy_factor_ > 0 || ! keep_underdetermined_)
  {
    std::vector<ProteinHit> filtered_proteoform_hits;
    filtered_proteoform_hits.reserve(proteoform_hits_.size());

    for (int k = 0; k < (keep_underdetermined_ ? 2 : 1); k++)
    {
      for (Size mod = 0; mod <= max_mod_cntr; mod++)
      {
        double taget_count = 0;
        double decoy_count = 0;

        double taget_count_pro = 0;
        double decoy_count_pro = 0;

        std::map<double, double> map_qvalue;
        std::map<double, double> map_qvalue_pro;

        for (auto& hit : proteoform_hits_)
        {
          std::vector<double> mod_masses = hit.getMetaValue("Modifications");
          if (mod_masses.size() != mod) continue;
          if (k == 0
              && ((double)hit.getMetaValue("Mass") < 0 || (int)hit.getMetaValue("StartPosition") < 0 || (int)hit.getMetaValue("EndPosition") < 0))
          {
            continue;
          }
          else if (k == 1
                   && ((double)hit.getMetaValue("Mass") > 0 && (int)hit.getMetaValue("StartPosition") > 0
                       && (int)hit.getMetaValue("EndPosition") > 0))
            continue;

          bool is_decoy = hit.getAccession().hasPrefix("DECOY");
          bool is_rep = hit.metaValueExists("Representative");
          if (is_decoy)
          {
            decoy_count += 1.0 / decoy_factor_;
            if (is_rep) decoy_count_pro += 1.0 / decoy_factor_;
          }
          else
          {
            taget_count++;
            if (is_rep) taget_count_pro++;
          }

          double tmp_qvalue = decoy_count / (decoy_count + taget_count);
          map_qvalue[hit.getScore()] = std::min(1.0, tmp_qvalue);

          if (! is_rep) continue;
          double tmp_qvalue_pro = decoy_count_pro / (decoy_count_pro + taget_count_pro);
          map_qvalue_pro[hit.getScore()] = std::min(1.0, tmp_qvalue_pro);
        }

        double cummin = 1.0;
        for (auto&& rit = map_qvalue.begin(); rit != map_qvalue.end(); ++rit)
        {
          cummin = std::min(rit->second, cummin);
          rit->second = cummin;
        }

        cummin = 1.0;
        for (auto&& rit = map_qvalue_pro.begin(); rit != map_qvalue_pro.end(); ++rit)
        {
          cummin = std::min(rit->second, cummin);
          rit->second = cummin;
        }

        for (auto& hit : proteoform_hits_)
        {
          std::vector<double> mod_masses = hit.getMetaValue("Modifications");
          if (mod_masses.size() != mod) continue;
          if (k == 0
              && ((double)hit.getMetaValue("Mass") < 0 || (int)hit.getMetaValue("StartPosition") < 0 || (int)hit.getMetaValue("EndPosition") < 0))
            continue;
          if (k == 1
              && ((double)hit.getMetaValue("Mass") > 0 && (int)hit.getMetaValue("StartPosition") > 0 && (int)hit.getMetaValue("EndPosition") > 0))
            continue;

          bool is_decoy = hit.getAccession().hasPrefix("DECOY");
          double qvalue = map_qvalue[hit.getScore()];
          hit.setMetaValue("qvalue", qvalue);
          if (! keep_decoy_ && is_decoy) continue;

          if (prsm_fdr_ < 1 && qvalue > prsm_fdr_) continue;

          auto iter = map_qvalue_pro.lower_bound(hit.getScore());
          if (iter != map_qvalue_pro.end())
          {
            double qvalue_pro = iter->second;
            hit.setMetaValue("proqvalue", qvalue_pro);
          }
          else
            hit.setMetaValue("proqvalue", 1.0);
          // if (pro_fdr_ < 1 && qvalue > pro_fdr_) continue;

          filtered_proteoform_hits.push_back(hit);
        }
      }
    }
    proteoform_hits_.swap(filtered_proteoform_hits);
  }
  std::sort(proteoform_hits_.begin(), proteoform_hits_.end(),
            [](const ProteinHit& left, const ProteinHit& right) { return left.getMetaValue("RT") < right.getMetaValue("RT"); });

  //  define proteoform index and tag to proteoform indices.
  int proteoform_index = 0;
  for (auto& hit : proteoform_hits_)
  {
    int fasta_index = hit.getMetaValue("FastaIndex");
    hit.setSequence(original_fasta_entry[fasta_index].sequence);
    hit.setMetaValue("Index", proteoform_index);
    for (int tag_index : (std::vector<int>)hit.getMetaValue("TagIndices").toIntList())
    {
      matching_hits_indices_[tag_index].push_back(proteoform_index);
    }
    proteoform_index++;
  }

  // TODO
  // per scan, get the scan number and the precursor mass - maybe given or maybe deconvolved.
  // get all protoeforms - mass and score.
  // get all combinations with a certain tolerance.
  // record the best ranking combination(s) along with summed score (can score be non overlapping?) and mass difference
}

void FLASHTnTAlgorithm::getProteoformHitsMatchedBy(const FLASHHelperClasses::Tag& tag, std::vector<ProteinHit>& hits) const
{
  int index = tag.getIndex();

  if (index < 0 || matching_hits_indices_.find(index) == matching_hits_indices_.end()) return;

  for (auto i : matching_hits_indices_.at(index))
  {
    hits.push_back(proteoform_hits_[i]);
  }
}

void FLASHTnTAlgorithm::getTags(std::vector<FLASHHelperClasses::Tag>& tags) const
{
  for (const auto& tag : tags_)
  {
    tags.push_back(tag);
  }
}

} // namespace OpenMS