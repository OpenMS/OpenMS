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

bool FLASHTnTAlgorithm::areConsistent_(const ProteinHit& a, const ProteinHit& b, double tol)
{
  double mass1 = a.getMetaValue("Mass");
  double mass2 = b.getMetaValue("Mass");
  if (mass1 * mass2 < 0) return false;

  int sp1 = a.getMetaValue("StartPosition");
  int ep1 = a.getMetaValue("EndPosition");

  int sp2 = b.getMetaValue("StartPosition");
  int ep2 = b.getMetaValue("EndPosition");

  std::vector<String> mod_accs1, mod_accs2;

  if (a.metaValueExists("Modifications"))
  {
    mod_accs1 = a.getMetaValue("ModificationACCs");
    std::sort(mod_accs1.begin(), mod_accs1.end());
  }

  if (b.metaValueExists("Modifications"))
  {
    mod_accs2 = b.getMetaValue("ModificationACCs");
    std::sort(mod_accs2.begin(), mod_accs2.end());
  }

  if (mod_accs1.empty() && mod_accs2.empty()) return sp1 == sp2 && ep1 == ep2; // unmodified ones

  // at least one is modified
  bool mass_matched = std::abs(mass1 - mass2) < std::max(mass1, mass2) * tol / 1e6 * 2;
  bool mod_matched = std::equal(mod_accs1.begin(), mod_accs1.end(), mod_accs2.begin(), mod_accs2.end());
  bool mod_loc_matched = mod_matched;

  if (mod_matched && mod_accs1.size() == mod_accs2.size())
  {
    std::vector<int> mod_starts1 = a.getMetaValue("ModificationStarts");
    std::vector<int> mod_starts2 = b.getMetaValue("ModificationStarts");
    std::vector<int> mod_ends1 = a.getMetaValue("ModificationEnds");
    std::vector<int> mod_ends2 = b.getMetaValue("ModificationEnds");

    for (Size i = 0; i < mod_starts1.size(); i++)
    {
      if (mod_starts1[i] > mod_ends2[i] || mod_starts2[i] > mod_ends1[i]) { mod_loc_matched = false; break; }
    }
  }

  if (mass1 > 0 && mass_matched)
  {
    //false if mods are different or mods are the same but mod locs are different
    if (!mod_matched || !mod_loc_matched) return false;
    return true;
  }

  // masses are underdetermined
  if (sp1 == sp2 && ep1 == ep2 && mod_matched && mod_loc_matched)  return true;
  // true if sp ep are the same and mod accs are the same and locs are consistent.
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


void FLASHTnTAlgorithm::vectorizeProteinSequence_(const std::vector<std::string>& cleaned_protein_seqs,
                                                     std::vector<std::unordered_set<int>>& vec_pro,
                                                     std::vector<std::unordered_set<int>>& rev_vec_pro)
{
  vec_pro.reserve(cleaned_protein_seqs.size());
  rev_vec_pro.reserve(cleaned_protein_seqs.size());
  for (const auto& seq : cleaned_protein_seqs)
  { //ignore X
    std::unordered_set<int> vec, rev_vec;
    double mass = 0;
    vec.insert(0);
    rev_vec.insert(0);

    for (Size j = 0; j < seq.size(); j++)
    {
      Size index = j;
      mass += ResidueDB::getInstance()->getResidue(seq[index])->getMonoWeight(Residue::Internal);
      int vindex = int(round(mass));
      vec.insert(vindex);
    }
    mass = 0;
    for (Size j = 0; j < seq.size(); j++)
    {
      Size index = seq.size() - j - 1;
      mass += ResidueDB::getInstance()->getResidue(seq[index])->getMonoWeight(Residue::Internal);
      int vindex = int(round(mass));
      rev_vec.insert(vindex);
    }
    vec_pro.push_back(vec);
    rev_vec_pro.push_back(rev_vec);
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
  std::map<String, std::vector<double>> tag_to_n_flanking_masses, tag_to_c_flanking_masses; // tag to protein flanking masses

  std::map<int, std::vector<Size>> scan_to_tag_indices;         // scan number to tag indices
  std::vector<std::string> tag_seqs, protein_seqs, cleaned_protein_seqs;
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
    cleaned_protein_seqs.reserve(fasta_entry.size());
    double taget_count = 0;
    double decoy_count = 0;
    for (const auto& fe : fasta_entry)
    {
      if (fe.identifier.hasPrefix("DECOY")) { decoy_count++; }
      else { taget_count++; }
      protein_seqs.push_back(fe.sequence);
      String cleaned_seq;
      std::remove_copy(fe.sequence.begin(), fe.sequence.end(),
                       std::back_inserter(cleaned_seq), 'X');
      cleaned_protein_seqs.push_back(cleaned_seq);
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
    tag_to_n_flanking_masses[seq] = std::vector<double>();
    tag_to_c_flanking_masses[seq] = std::vector<double>();

    tag_seqs.push_back(seq);
  } // now we have all the tags

  // tag_seqs protein_seqs

  ACTrie ac_trie(0, 0); // no ambiguities, no mismatches
  ac_trie.addNeedlesAndCompress(tag_seqs);
  OpenMS::ACTrieState ac_state;

  startProgress(0, (SignedSize)protein_seqs.size(), "Searching tags against database ...");

  for (int i = 0; i < protein_seqs.size(); i++) // // Natural19WithoutI
  {
    nextProgress();
    const String& pseq = protein_seqs[i];
    const String& cleanseq = cleaned_protein_seqs[i];
    ac_state.setQuery(pseq);
    ac_trie.getAllHits(ac_state);

    for (const auto& h : ac_state.hits)
    {
      String seq = tag_seqs[h.needle_index];
      tag_to_protein_indices[seq].push_back(i);
      tag_to_protein_positions[seq].push_back(h.query_pos);

      double nmass = h.query_pos <= 0 ? .0 : AASequence::fromString(cleanseq.substr(0, h.query_pos)).getMonoWeight(Residue::Internal);
      double cmass = cleanseq.length() - seq.length() <= h.query_pos ? .0 :
         AASequence::fromString(cleanseq.substr(h.query_pos + seq.length(), cleanseq.length() - h.query_pos - seq.length())).getMonoWeight(Residue::Internal);

      tag_to_n_flanking_masses[seq].push_back(nmass);
      tag_to_c_flanking_masses[seq].push_back(cmass);
    }
  }

  endProgress();
  startProgress(0, (SignedSize)dspecs.size(), "Running candidate protein filtration and extension algorithm ...");

  std::vector<std::unordered_set<int>> vec_pro, rev_vec_pro;
  vectorizeProteinSequence_(cleaned_protein_seqs, vec_pro, rev_vec_pro);

  for (const auto& dspec : dspecs)
  {
    nextProgress();
    if (scan_to_tag_indices.find(dspec.getScanNumber()) == scan_to_tag_indices.end()) continue;
    const auto& tag_indices = scan_to_tag_indices[dspec.getScanNumber()];
    std::vector<ProteinHit> hits;
    std::map<Size, std::set<Size>> pi_to_pos;         // protein index to covered positions
    std::map<Size, std::set<Size>> pi_to_tag_indices; // protein index to tag indices
    std::map<Size, std::set<int>> pi_to_n_flanking_masses, pi_to_c_flanking_masses;
    for (const auto& tag_index : tag_indices)
    {
      const auto& tag = tags_[tag_index];
      const auto& seq = tag.getUppercaseSequence();

      const auto& pis = tag_to_protein_indices[seq];
      const auto& pps = tag_to_protein_positions[seq];

      for (Size i = 0; i < pis.size(); i++)
      {
        const auto& pseq = cleaned_protein_seqs[pis[i]];

        if (tag.getNtermMass() >= 0 && pps[i] > 0 && pps[i] <= pseq.length() + seq.length())
        {
          double nmass = tag_to_n_flanking_masses[seq][i];
          if (tag.getNtermMass() > nmass + max_mod_mass + 200) continue;
          pi_to_n_flanking_masses[pis[i]].insert((int)round(nmass - tag.getNtermMass()));
        }
        if (tag.getCtermMass() >= 0 && pps[i] + seq.length() >= 0 && pseq.length() > pps[i] + seq.length())
        {
          double cmass = tag_to_c_flanking_masses[seq][i];
          if (tag.getCtermMass() > cmass + max_mod_mass + 200) continue;
          pi_to_c_flanking_masses[pis[i]].insert((int)round(cmass - tag.getCtermMass()));
        }

        pi_to_tag_indices[pis[i]].insert(tag_index);
        pi_to_pos[pis[i]].insert(pps[i]);
      }
    }

    for (const auto& [pi, pp] : pi_to_pos)
    {
      const auto& fe = fasta_entry[pi];
      Size covered_pos = pp.size();
      if (covered_pos < min_tag_count) continue;

      // check flanking masses ...
      double coverage = (double)covered_pos / (double)fe.sequence.length();
      std::vector<int> positions(pp.begin(), pp.end());

      //if (positions.size() < min_tag_count) continue;
      ProteinHit hit(coverage, 0, fe.identifier, fe.sequence); //
      hit.setDescription(fe.description);
      hit.setMetaValue("Scan", dspec.getScanNumber());
      hit.setMetaValue("MatchedAA", covered_pos);
      hit.setCoverage(coverage);
      hit.setMetaValue("NtermFlankingMasses", std::vector<int>(pi_to_n_flanking_masses[pi].begin(), pi_to_n_flanking_masses[pi].end()));
      hit.setMetaValue("CtermFlankingMasses", std::vector<int>(pi_to_c_flanking_masses[pi].begin(), pi_to_c_flanking_masses[pi].end()));

      hit.setMetaValue("TagIndices", std::vector<int>(pi_to_tag_indices[pi].begin(), pi_to_tag_indices[pi].end()));
      hit.setMetaValue("TagPositions", positions);
      hit.setMetaValue("FastaIndex", pi);
      hits.push_back(hit);
    }

    if (hits.empty()) continue;

    FLASHTaggerAlgorithm::runMatching(hits, vec_pro, rev_vec_pro, dspec, max_mod_mass);

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