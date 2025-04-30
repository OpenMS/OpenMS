// Copyright (c) 2002-2024, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHExtenderAlgorithm.h>
#include <queue>
#include <utility>

namespace OpenMS
{
inline const int multi_ion_score = 1;
inline const bool debug = false;
inline const double i2f_mass = Residue::getInternalToFull().getMonoWeight();
FLASHExtenderAlgorithm::FLASHExtenderAlgorithm(): DefaultParamHandler("FLASHExtenderAlgorithm"), ProgressLogger()
{
  setDefaultParams_();
}

FLASHExtenderAlgorithm& FLASHExtenderAlgorithm::operator=(const FLASHExtenderAlgorithm& rhs)
{
  if (this == &rhs) return *this;

  DefaultParamHandler::operator=(rhs);
  return *this;
}

void FLASHExtenderAlgorithm::setDefaultParams_()
{
  defaults_.setValue("max_mod_mass", 500.0, "Maximum mass shift allowed for modifications.");
  defaults_.setValue("max_mod_count", 2, "Maximum number of blind modifications.");

  defaults_.setValue("ion_type", std::vector<std::string> {"b", "y"}, "Specifies ion types to consider");
  defaults_.setValidStrings("ion_type", {"b", "c", "a", "y", "z", "x", "zp1", "zp2"});

  defaultsToParam_();
}

void FLASHExtenderAlgorithm::updateMembers_()
{
  max_mod_cntr_ = param_.getValue("max_mod_count");
  max_mod_mass_ = param_.getValue("max_mod_mass");
}

inline Size FLASHExtenderAlgorithm::getVertex_(int node_index, int pro_index, int score, int num_mod, Size pro_mass_size) const
{
  return ((node_index * pro_mass_size + pro_index) * (max_mod_cntr_ + 1) + num_mod) * (max_path_score_ - min_path_score_ + 1)
         + (std::min(max_path_score_, std::max(min_path_score_, score)) - min_path_score_);
}

inline int FLASHExtenderAlgorithm::getNodeIndex_(Size vertex, Size pro_mass_size) const
{
  return (vertex / (max_path_score_ - min_path_score_ + 1) / ((Size)max_mod_cntr_ + 1)) / pro_mass_size;
}

inline int FLASHExtenderAlgorithm::getProIndex_(Size vertex, Size pro_mass_size) const
{
  return ((vertex / (max_path_score_ - min_path_score_ + 1) / (max_mod_cntr_ + 1))) % pro_mass_size;
}

inline int FLASHExtenderAlgorithm::getScore_(Size vertex) const
{
  return vertex % (max_path_score_ - min_path_score_ + 1) + min_path_score_;
}

inline int FLASHExtenderAlgorithm::getModNumber_(Size vertex) const
{
  return (vertex / (max_path_score_ - min_path_score_ + 1)) % (max_mod_cntr_ + 1);
}

// take the hits. Just calculate the mass of truncated protein. Then add modification masses if they are disjoint. If they overlap and the same mass,
// we have a single one. If they are overlapping but different, we add all of them. the max mod count is also adjusted.
void FLASHExtenderAlgorithm::calculatePrecursorMass_(const ProteinHit& hit,
                                                     const std::map<int, std::vector<Size>>& best_path_map,
                                                     HitInformation& hi)
{
  hi.calculated_precursor_mass_ = -1;
  if (hi.protein_start_position_ < 0 || hi.protein_end_position_ < 0) return;

  const auto& bp0 = best_path_map.at(0);
  const auto& bp1 = best_path_map.at(1);

  //PeakGroup pg;
  //pg.setQscore(.8);
  const int terminal_score_threshold = 1;//(int)round(FLASHTaggerAlgorithm::getNodeScore(pg) * 3); // n term or c term should have at least this score to be considered = three high quality masses
  int max_score = terminal_score_threshold; // need to exceed 20 for precursor correction
  double min_tol = 1000;
  Size pro_size = hi.pro_mass_map_[0].size();
  int min_excessive_aa = pro_size;
  for(auto iter0 = bp0.rbegin(); iter0 != bp0.rend(); iter0++)
  {
    int proi = getProIndex_(*iter0, pro_size);
    int proj = 0;
    auto iter1 = bp1.rbegin();
    for(; iter1 != bp1.rend(); iter1++)
    {
      proj = getProIndex_(*iter1, pro_size);
      if (proi + proj >= (int)hit.getSequence().size()) break;
    }
    if (getScore_(*iter0) < terminal_score_threshold || getScore_(*iter1) < terminal_score_threshold) continue; //
    int excessive_aa = (proi + proj) - (int)hit.getSequence().size();

    if (excessive_aa < 0) continue;
    if (excessive_aa > min_excessive_aa) continue;
    if (excessive_aa < min_excessive_aa)
    {
      max_score = terminal_score_threshold;
      min_tol = 1000;
      min_excessive_aa = excessive_aa;
    }

    int score = getScore_(*iter0) + getScore_(*iter1);
    if (max_score > score) continue;

    if (max_score == score)
    {
      double tol = std::max(hi.tol_spec_map_.at(0)[getNodeIndex_(*iter0, pro_size)].getIntensity(),
        hi.tol_spec_map_.at(1)[getNodeIndex_(*iter1, pro_size)].getIntensity());
      if (tol > min_tol) continue;
      min_tol = tol;
    }
    max_score = score;
    hi.calculated_precursor_mass_ = hi.node_spec_map_.at(0)[getNodeIndex_(*iter0, pro_size)].getMZ()
                                    + hi.node_spec_map_.at(1)[getNodeIndex_(*iter1, pro_size)].getMZ();

    if (proi + proj > (int)hit.getSequence().size())
    {
      hi.calculated_precursor_mass_ += hi.pro_mass_map_[0].back() - (hi.pro_mass_map_[0][proi] + hi.pro_mass_map_[1][proj]);
    }
  }
  if (hi.calculated_precursor_mass_ > 0) hi.calculated_precursor_mass_ += i2f_mass;
}

void FLASHExtenderAlgorithm::getProMasses_(const ProteinHit& hit, std::vector<double>& pro_masses, int mode)
{
  pro_masses.reserve(hit.getSequence().size() + 1);
  pro_masses.push_back(0);

  auto seq = hit.getSequence();
  pro_masses.reserve(seq.size());

  if (mode == 0) seq = seq.reverse();
  for (const auto& aa : seq)
  {
    if (aa == 'X') pro_masses.push_back(pro_masses.back()); // repeat the previous mass
    else
      pro_masses.push_back(pro_masses.back() + AASequence::fromString(aa, true).getMonoWeight(Residue::Internal));
  }
}

double FLASHExtenderAlgorithm::getSpecMassSpan_(const std::vector<Size>& path, const MSSpectrum& node_spec, int pro_mass_size) const
{
  return node_spec[(getNodeIndex_(path[0], pro_mass_size))].getMZ();
}

double FLASHExtenderAlgorithm::getProteinMassSpan_(const std::vector<Size>& path, const std::vector<double>& pro_masses) const
{
  double minmass = 0, maxmass = 0;
  maxmass = pro_masses[getProIndex_(path[0], pro_masses.size())];

  for (int i = path.size() - 1; i >= 0; i--)
  {
    if (getNodeIndex_(path[i], pro_masses.size()) != 0) break;
    minmass = pro_masses[getProIndex_(path[i], pro_masses.size())];
  }
  return std::abs(maxmass - minmass);
}

int FLASHExtenderAlgorithm::getModifiedAACount_(const std::vector<Size>& path) const
{
  int cntr = 0;
  for (auto vertex : path)
  {
    if (getModNumber_(vertex) == 0) continue;
    cntr ++;
  }
  return cntr;
}

int FLASHExtenderAlgorithm::getProteinLength_(const std::vector<Size>& path, const std::vector<double>& pro_masses) const
{
  int minindex = 0, maxindex = 0;
  maxindex = getProIndex_(path[0], pro_masses.size());

  for (int i = path.size() - 1; i >= 0; i--)
  {
    if (getNodeIndex_(path[i], pro_masses.size()) != 0) break;
    minindex = getProIndex_(path[i], pro_masses.size());
  }
  return std::abs(maxindex - minindex);
}

void FLASHExtenderAlgorithm::defineNodes_(const DeconvolvedSpectrum& dspec, HitInformation& hi, double max_mass)
{
  // 0 for suffix 1 for prefix 2 for suffix and prefix if precursor mass is available
  double max_proteoform_mass = -1;
  float max_score = -(float)max_path_score_;

  std::vector<int> isos{0};
  if (hi.mode_ == 2 && hi.calculated_precursor_mass_ > 0)
  {
    for (int i = 1; i <= allowed_isotope_error_; ++i)
    {
      isos.push_back(i);
      isos.push_back(-i);
    }
  }
  for (int iso : isos)
  {
    MSSpectrum node_spec, tol_spec;
    node_spec.reserve(dspec.size() * ion_types_str_.size() + 1);
    tol_spec.reserve(dspec.size() * ion_types_str_.size() + 1);
    double proteoform_mass = hi.calculated_precursor_mass_ + iso * Constants::C13C12_MASSDIFF_U;

    for (const auto& pg : dspec)
    {
      const auto& shifts = (hi.mode_ == 0) ? suffix_shifts_ :
                           (hi.mode_ == 1) ? prefix_shifts_ :
                           (hi.mode_ == 2) ? prefix_shifts_ :
                                           std::vector<double>{};

      for (const auto& shift : shifts)
      {
        double mass = pg.getMonoMass() - shift;
        if (mass <= 0 || mass > max_mass + 1) continue;
        node_spec.emplace_back(mass, std::max(1, FLASHTaggerAlgorithm::getNodeScore(pg)));
        tol_spec.emplace_back(mass, tol_ * pg.getMonoMass());
      }

      if (hi.mode_ == 2 && hi.calculated_precursor_mass_ > 0)
      {
        for (const auto& shift : suffix_shifts_)
        {
          double mass = pg.getMonoMass() - shift;
          if (mass <= 0 || mass >= proteoform_mass - i2f_mass) continue;
          double prefix_mass = proteoform_mass - i2f_mass - mass;
          node_spec.emplace_back(prefix_mass, std::max(1, FLASHTaggerAlgorithm::getNodeScore(pg)));
          tol_spec.emplace_back(prefix_mass, tol_ * pg.getMonoMass());
        }
      }
    }

    node_spec.sortByPosition();
    tol_spec.sortByPosition();

    MSSpectrum t_node_spec, t_tol_spec;
    t_node_spec.reserve(node_spec.size() + 1);
    t_tol_spec.reserve(node_spec.size() + 1);

    t_node_spec.emplace_back(0, 0);
    t_tol_spec.emplace_back(0, 0);
    float overlapped_score = 0;

    for (Size k = 0; k < node_spec.size(); k++)
    {
      const auto& p = node_spec[k];
      double mass = p.getMZ();
      float score = p.getIntensity();
      double prev_margin = k == 0 ? .0 : tol_spec[k - 1].getIntensity();
      double margin = tol_spec[k].getIntensity();

      if (mass - margin < t_node_spec.back().getMZ() + prev_margin) // they are the same
      {
        float prev_score = t_node_spec.back().getIntensity();

        margin = (mass + margin - t_node_spec.back().getMZ() + prev_margin) / 2.0;
        mass = (mass + margin + t_node_spec.back().getMZ() - prev_margin) / 2.0;
        if (t_node_spec.size() > 1)
        {
          t_node_spec.pop_back();
          t_tol_spec.pop_back();
        }
        score = multi_ion_score + std::max(prev_score, score);
        //score = suffix_prefix == prev_suffix_prefix ? multi_ion_score + std::max(prev_score, score) : multi_ion_score + prev_score + score;
        overlapped_score += score;
      }

      if (mass <= 0) continue;

      t_node_spec.emplace_back(mass, score);
      t_tol_spec.emplace_back(mass, margin);
    }

    if (max_score < overlapped_score - multi_ion_score)
    {
      hi.node_spec_map_[hi.mode_] = t_node_spec;
      hi.tol_spec_map_[hi.mode_] = t_tol_spec;
      max_score = overlapped_score;
      max_proteoform_mass = proteoform_mass;
    }
  }

  hi.node_spec_map_[hi.mode_].sortByPosition();
  hi.tol_spec_map_[hi.mode_].sortByPosition();
  if (hi.calculated_precursor_mass_ > 0)
  {
    hi.calculated_precursor_mass_ = max_proteoform_mass;
    hi.node_spec_map_[hi.mode_].emplace_back(hi.calculated_precursor_mass_ - i2f_mass, 1);
    hi.tol_spec_map_[hi.mode_].emplace_back(hi.calculated_precursor_mass_ - i2f_mass, tol_ * hi.calculated_precursor_mass_);
  }
}

void FLASHExtenderAlgorithm::run_(const ProteinHit& hit,
                                  HitInformation& hi,
                                  const std::vector<FLASHHelperClasses::Tag>& matched_tags,
                                  std::map<int, std::vector<Size>>& all_paths_per_mode,
                                  int max_mod_cntr_for_last_mode) // per hit
{
  std::vector<std::vector<int>> tag_edges(4); // pro start end indices and node start end indices
  std::set<Size> sinks;
  std::set<String> upper_case_seqs;

  for (const auto& tag : matched_tags)
  {
    if (std::any_of(tag.getSequence().begin(), tag.getSequence().end(), ::islower)) continue;
    upper_case_seqs.insert(tag.getSequence());
  }

  if (! matched_tags.empty())
  {
    for (auto& edge : tag_edges)
    {
      edge.reserve(matched_tags.size() * 2);
    }
  }
  const auto& node_spec = hi.node_spec_map_[hi.mode_];
  const auto& pro_masses = hi.pro_mass_map_[hi.mode_];
  hi.dag_ = FLASHHelperClasses::DAG((1 + node_spec.size()) * (1 + pro_masses.size()) * (1 + max_mod_cntr_) * (1 + max_path_score_ - min_path_score_));

  bool tag_found = false;
  for (const auto& tag : matched_tags)
  {
    if (hi.mode_ == 2 && hi.calculated_precursor_mass_ <= 0) continue;
    const auto& upper_seq = tag.getUppercaseSequence();
    bool has_lower_case = std::any_of(tag.getSequence().begin(), tag.getSequence().end(), ::islower) ;
    if (has_lower_case && upper_case_seqs.find(upper_seq) != upper_case_seqs.end()) continue;

    tag_found = true;
    std::vector<int> positions;
    std::vector<double> masses;

    FLASHTaggerAlgorithm::fillMatchedPositionsAndFlankingMassDiffs(positions, masses, -1, hit.getSequence(), tag);
    auto tag_masses = tag.getMzs();
    std::sort(tag_masses.begin(), tag_masses.end());

    double seq_mass = 0;
    if (has_lower_case)
    {
      seq_mass = AASequence::fromString(upper_seq, true).getMonoWeight(Residue::Internal);
    }
    std::vector<double> start_masses, end_masses, start_tols, end_tols;
    if (tag.getCtermMass() >= 0) // suffix
    {
      for (const auto& shift : suffix_shifts_)
      {
        double start_mass = tag_masses[0] - shift;
        double end_mass = tag_masses.back() - shift;

        if (seq_mass > 0)
        {
          if (std::abs(start_mass + seq_mass - end_mass) > end_mass * tol_) end_mass = start_mass + seq_mass;
        }

        start_tols.push_back(start_mass * tol_);
        end_tols.push_back(end_mass * tol_);

        if (hi.mode_ == 2)
        {
          start_mass = hi.calculated_precursor_mass_ - i2f_mass - start_mass;
          end_mass = hi.calculated_precursor_mass_ - i2f_mass - end_mass;
          start_masses.push_back(end_mass);
          end_masses.push_back(start_mass);
        }
        else
        {
          start_masses.push_back(start_mass);
          end_masses.push_back(end_mass);
        }
      }
    }
    else // prefix
    {
      for (const auto& shift : prefix_shifts_)
      {
        double start_mass = tag_masses[0] - shift;
        double end_mass = tag_masses.back() - shift;

        if (seq_mass > 0)
        {
          if (std::abs(start_mass + seq_mass - end_mass) > end_mass * tol_) end_mass = start_mass + seq_mass;
        }

        start_tols.push_back(start_mass * tol_);
        end_tols.push_back(end_mass * tol_);
        start_masses.push_back(start_mass);
        end_masses.push_back(end_mass);
      }
    }

    bool same_terminal_tag = hi.mode_ == 2 || (hi.mode_ == 0 && tag.getCtermMass() >= 0) || (hi.mode_ == 1 && tag.getNtermMass() >= 0);

    for (Size l = 0; l < start_masses.size(); l++)
    {
      int highest_score_start = -1, highest_score_end = -1;
      if (same_terminal_tag)
      {
        double delta_start = start_tols[l];
        double delta_end = end_tols[l];

        highest_score_start = node_spec.findHighestInWindow(start_masses[l], delta_start, delta_start);
        highest_score_end = node_spec.findHighestInWindow(end_masses[l], delta_end, delta_end);

        if (highest_score_start < 0 || highest_score_end < 0 || highest_score_start >= (int)node_spec.size()
            || highest_score_end >= (int)node_spec.size())
          continue;
      }
      for (int pos : positions)
      {
        if (hi.mode_ == 0) // suffix inverted
        {
          pos = (int)pro_masses.size() - 1 - pos; // invert pos
          if (pos - tag.getLength() >= 0 && pos < (int)pro_masses.size())
          {
            tag_edges[0].emplace_back(pos - tag.getLength());
            tag_edges[1].emplace_back(pos);
            tag_edges[2].emplace_back(highest_score_start);
            tag_edges[3].emplace_back(highest_score_end); // check...
          }
        }
        else
        {
          if (pos >= 0 && pos + tag.getLength() < pro_masses.size())
          {
            tag_edges[0].emplace_back(pos);
            tag_edges[1].emplace_back(pos + tag.getLength()); // this can be much faster...
            tag_edges[2].emplace_back(highest_score_start);
            tag_edges[3].emplace_back(highest_score_end);
          }
        }
      }
    }
  }

  constructDAG_(sinks, hi, tag_edges, max_mod_cntr_for_last_mode, tag_found);

  Size src = getVertex_(0, 0, 0, 0, pro_masses.size());
  std::vector<int> max_scores(max_mod_cntr_ + 1, 0);
  for (Size sink : sinks)
  {
    int num_mod = getModNumber_(sink);
    if (sink == src || getScore_(sink) < max_scores[num_mod]) continue; // getModNumber_(sink) == 0 ||
    if (hi.calculated_precursor_mass_ > 0 && getNodeIndex_(sink, pro_masses.size()) < (int)node_spec.size() - 1) continue;
    max_scores[num_mod] = getScore_(sink);
  }
  std::vector<std::vector<std::vector<Size>>> paths(max_mod_cntr_ + 1, std::vector<std::vector<Size>>());
  for (Size sink : sinks)
  {
    int num_mod = getModNumber_(sink);
    if (sink == src || getScore_(sink) < max_scores[num_mod]) continue; //
    if (hi.calculated_precursor_mass_ > 0 && getNodeIndex_(sink, pro_masses.size()) < (int)node_spec.size() - 1) continue;
    //std::vector<std::vector<Size>> sub_paths;
    hi.dag_.findAllPaths(sink, src, paths[num_mod], 0);
  }
  for (int num_mod = 0; num_mod <= max_mod_cntr_; num_mod++)
  {
    for (const auto& path : paths[num_mod])
    {
      double mass = getSpecMassSpan_(path, node_spec, pro_masses.size());
      double pro_mass = getProteinMassSpan_(path, pro_masses);
      int pro_len = getProteinLength_(path, pro_masses);

//      if (hi.mode_ == 2 && hi.calculated_precursor_mass_ > 0 // TODO is this necessary?
//          && std::abs(hi.calculated_precursor_mass_ - i2f_mass - mass) > 1.1)
//      {
//        continue;
//      }

      auto iter = all_paths_per_mode.find(num_mod);

      if (iter == all_paths_per_mode.end())
      {
        all_paths_per_mode[num_mod] = path;
      }
      else
      {
        double mod_mass_ = num_mod == 0 ? 0 : std::abs(getSpecMassSpan_(iter->second, node_spec, pro_masses.size()) - getProteinMassSpan_(iter->second, pro_masses));
        double mod_mass = num_mod == 0 ? 0 : std::abs(mass - pro_mass);
        if (mod_mass_ < mod_mass) continue;

        // Prefer the path with smaller mod mass
        if (mod_mass_ > mod_mass ||
            (mod_mass_ == mod_mass &&
             (getProteinLength_(iter->second, pro_masses) < pro_len ||
              (getProteinLength_(iter->second, pro_masses) == pro_len &&
               getModifiedAACount_(path) < getModifiedAACount_(iter->second)))))
        {
          all_paths_per_mode[num_mod] = path;
        }
      }
    }
  }
}

void FLASHExtenderAlgorithm::run(std::vector<ProteinHit>& hits,
                                 const std::vector<FLASHHelperClasses::Tag>& tags,
                                 const DeconvolvedSpectrum& dspec,
                                 double ppm,
                                 bool multiple_hits_per_spec)
{
  if (hits.empty()) return;
  // setLogType(CMD);

  ion_types_str_ = param_.getValue("ion_type").toStringVector();
  std::sort(ion_types_str_.begin(), ion_types_str_.end());
  for (const auto& ion_str : ion_types_str_)
  {
    if (ion_str == "a") { prefix_shifts_.push_back(Residue::getInternalToAIon().getMonoWeight()); }
    else if (ion_str == "b") { prefix_shifts_.push_back(Residue::getInternalToBIon().getMonoWeight()); }
    else if (ion_str == "c") { prefix_shifts_.push_back(Residue::getInternalToCIon().getMonoWeight()); }
    else if (ion_str == "x") { suffix_shifts_.push_back(Residue::getInternalToXIon().getMonoWeight()); }
    else if (ion_str == "y") { suffix_shifts_.push_back(Residue::getInternalToYIon().getMonoWeight()); }
    else if (ion_str == "z") { suffix_shifts_.push_back(Residue::getInternalToZIon().getMonoWeight()); }
    else if (ion_str == "zp1") { suffix_shifts_.push_back(Residue::getInternalToZp1Ion().getMonoWeight()); }
    else if (ion_str == "zp2") { suffix_shifts_.push_back(Residue::getInternalToZp2Ion().getMonoWeight()); }
    else { continue; }
  }

  tol_ = ppm / 1e6;
  proteoform_hits_.clear();

  tags_ = tags;
  std::vector<double> mzs;
  std::vector<int> scores;

  proteoform_hits_.reserve(hits.size());
  if (dspec.getPrecursorPeakGroup().getMonoMass() > 0) { given_precursor_mass_ = dspec.getPrecursorPeakGroup().getMonoMass(); }

  startProgress(0, (int)hits.size(), "running FLASHExtender ...");

#pragma omp parallel for default(none) shared(hits, dspec, multiple_hits_per_spec, i2f_mass, std::cout)
  for (int i = 0; i < ((int)hits.size()); i++)
  {
    nextProgress();
    auto& hit = hits[i];
    HitInformation hi;
    int total_score = 0;
    std::vector<int> mod_starts, mod_ends;
    std::vector<double> mod_masses, mod_tols;
    int max_nterm_index = 0, max_cterm_rindex = 0;

    std::map<int, std::map<int, std::vector<Size>>> all_path_map; // mode, num_mod, path
    std::map<int, std::vector<Size>> best_path_map;               // mode, best paths

    std::vector<int> used_mode;
    std::vector<FLASHHelperClasses::Tag> matched_tags;

    if (! hit.metaValueExists("TagIndices")) continue;

    const std::vector<int>& tag_indices = hit.getMetaValue("TagIndices");

    matched_tags.reserve(tag_indices.size());
    for (const auto& index : tag_indices)
    {
      matched_tags.push_back(tags_[index]);
    }

    std::map<int, std::set<int>> matched_position_map;
    bool precursor_by_fragment = false;

    for (hi.mode_ = 0; hi.mode_ <= 2; hi.mode_++)
    {
      int max_mod_cntr_for_last_mode = -1;
      if (hi.mode_ == 2 && hi.calculated_precursor_mass_ <= 0)
      {//const ProteinHit& hit,
        if (max_nterm_index + max_cterm_rindex >= (int)hit.getSequence().size()) calculatePrecursorMass_(hit, best_path_map, hi);
        max_mod_cntr_for_last_mode = std::min(max_mod_cntr_, (int)mod_starts.size() + 1);

        if (hi.calculated_precursor_mass_ <= 0) hi.calculated_precursor_mass_ = given_precursor_mass_;
        else precursor_by_fragment = true;
        if (hi.calculated_precursor_mass_ <= 0) break;
      }

      auto& pro_masses = hi.pro_mass_map_[hi.mode_] = std::vector<double>();
     // auto& node_spec = hi.node_spec_map_[hi.mode_] = MSSpectrum();
      //auto& tol_spec = hi.tol_spec_map_[hi.mode_] = MSSpectrum();

      getProMasses_(hit, pro_masses, hi.mode_);
      defineNodes_(dspec, hi, pro_masses.back());

      if (hi.visited_.empty())
        hi.visited_ = boost::dynamic_bitset<>((3 + dspec.size() * ion_types_str_.size()) * (1 + pro_masses.size()) * (1 + max_mod_cntr_)
                                              * (1 + max_path_score_ - min_path_score_));
      run_(hit, hi, matched_tags, all_path_map[hi.mode_], max_mod_cntr_for_last_mode);

      if (hi.mode_ < 2)
      {
        const auto paths_c = all_path_map.find(0);
        const auto paths_n = all_path_map.find(1);

        std::map<int, int> nscores, cscores; // mod and score

        if (paths_n != all_path_map.end())
        {
          for (const auto& [mod, path] : paths_n->second)
          {
            if (path.empty()) continue;
            nscores[mod] = getScore_(path[0]);
          }
        }
        if (paths_c != all_path_map.end())
        {
          for (const auto& [mod, path] : paths_c->second)
          {
            if (path.empty()) continue;
            cscores[mod] = getScore_(path[0]);
          }
        }
        int max_score = 0;

        if (nscores.empty()) // only c term
        {
          for (const auto& [mod, path] : paths_c->second)
          {
            if (path.empty()) continue;
            if (max_score >= getScore_(path[0])) continue;
            max_score = getScore_(path[0]);
            best_path_map[0] = path;
          }
        }
        else if (cscores.empty()) // only n term
        {
          for (const auto& [mod, path] : paths_n->second)
          {
            if (path.empty()) continue;
            if (max_score >= getScore_(path[0])) continue;
            max_score = getScore_(path[0]);
            best_path_map[1] = path;
          }
        }
        else // both terms
        {
          for (int mc = 0; mc <= max_mod_cntr_; mc++)
          {
            if (cscores.find(mc) == cscores.end()) continue;
            const auto& cpath = paths_c->second[mc];
            for (int mn = 0; mc + mn <= max_mod_cntr_; mn++)
            {
              if (nscores.find(mn) == nscores.end()) continue;
              int sum_score = nscores[mn] + cscores[mc];
              if (max_score >= sum_score) continue;
              max_score = sum_score;
              best_path_map[0] = cpath;
              best_path_map[1] = paths_n->second[mn];
            }
          }
        }
      }
      else if (hi.mode_ == 2)
      {
        int max_score = 0;
        const auto paths = all_path_map.find(2);

        for (const auto& [mod, path] : paths->second)
        {
          if (path.empty()) continue;
          if (max_score >= getScore_(path[0])) continue;
          max_score = getScore_(path[0]);
          best_path_map[2] = path;
        }
        mod_starts.clear();
        mod_ends.clear();
        mod_masses.clear();
        mod_tols.clear();
//
//        int mode0_score = best_path_map.find(0) == best_path_map.end() ||  best_path_map[0].empty() ? 0 :  getScore_(best_path_map[0][0]);
//        int mode1_score = best_path_map.find(1) == best_path_map.end() ||  best_path_map[1].empty() ? 0 :  getScore_(best_path_map[1][0]);
//
//        if (max_score < mode0_score + mode1_score)
//        {
//          precursor_by_fragment = false;
//          hi.calculated_precursor_mass_ = -1;
//          break;
//        }
      }
      if (hi.mode_ == 0) continue;
      total_score = 0;
      // find the best paths per mode. Mode 0 and 1 should be considered together (since the modification counts for N C term paths should be summed up).

      for (int m = hi.mode_ == 2 ? 2 : 0; m <= hi.mode_; m++)
      {
        if (best_path_map.empty() || best_path_map.find(m) == best_path_map.end() || best_path_map[m].empty()) continue;
        auto& best_path = best_path_map[m];
        auto& t_pro_masses = hi.pro_mass_map_[m];
        auto& t_node_spec = hi.node_spec_map_[m];

        double prev_mass_shift = 0;
        int prev_mod_count = 0;
        int pre_pro_index = 0;
        int pre_node_index = 0;
        double mod_mass = 0;//, total_mod_mass = 0;
        //int total_mod_count = getModNumber_(*best_path.begin());

        for (auto iter = best_path.rbegin(); iter != best_path.rend(); iter++)
        {
          auto pro_index = getProIndex_(*iter, t_pro_masses.size());
          auto node_index = getNodeIndex_(*iter, t_pro_masses.size());
          auto mass_shift = t_node_spec[node_index].getMZ() - t_pro_masses[pro_index];
          auto mod_count = getModNumber_(*iter);

          if (node_index == 0)
          {
            if (m > 0) hi.protein_start_position_ = pro_index;
            if (m == 0) hi.protein_end_position_ = (int)hit.getSequence().size() - pro_index; //
            if (mod_count == 0) prev_mass_shift = mass_shift;
          }

          int pro_seq_index = m > 0 ? pro_index : ((int)hit.getSequence().size() - pro_index);
          if (node_index > 0 && t_node_spec[node_index].getMZ() != hi.calculated_precursor_mass_ - i2f_mass)
          {
            matched_position_map[m].insert(pro_seq_index);
          }

          if (m == 0) max_cterm_rindex = std::max(max_cterm_rindex, pro_index);
          if (m == 1) max_nterm_index = std::max(max_nterm_index, pro_index);
          if (m == 2) hi.protein_end_position_ = pro_index;
          if (mod_count != prev_mod_count)
          {
            mod_mass = mass_shift - prev_mass_shift;
           // total_mod_mass += mod_mass;
            int end = m > 0 ? (pro_index - 1) : ((int)hit.getSequence().size() - 1 - pre_pro_index);
            int start = m > 0 ? (pre_pro_index - 1) : ((int)hit.getSequence().size() - 1 - pro_index);

            for (int pi = pre_pro_index + 1; pi < pro_index; pi++)
            {
              double pm = mass_shift + t_pro_masses[pi];
              for (int ni = pre_node_index + 1; ni < node_index; ni++)
              {
                double nm = t_node_spec[ni].getMZ();
                if (std::abs(pm - nm) > hi.tol_spec_map_[m][ni].getIntensity()) continue;
                int t_end = m > 0 ? (pi - 1) : ((int)hit.getSequence().size() - 1 - pre_pro_index);
                int t_start = m > 0 ? (pre_pro_index - 1) : ((int)hit.getSequence().size() - 1 - pi);
                end = std::min(t_end, end);
                start = std::max(start, t_start);
              }
            }

            mod_starts.push_back(start + 1);
            mod_ends.push_back(end);
            mod_masses.push_back(mod_mass);
            mod_tols.push_back(hi.tol_spec_map_[m][node_index].getIntensity());
          }

          if (debug)
          {
            std::cout << hit.getAccession() << "\tmode\t" << m << "\tinput pre\t" << given_precursor_mass_ << "\tcal pre\t"
                      << std::to_string(hi.calculated_precursor_mass_) << "\tscore\t" << getScore_(*iter) << "\t" << node_index << "\t" << pro_index
                      << "\tin\t" << t_node_spec.size() << "\t" << t_pro_masses.size() << "\tmasses\t" << t_pro_masses.back() << "\t"
                      << std::to_string(t_pro_masses[pro_index]) << "\t" << std::to_string(t_pro_masses.back() - t_pro_masses[pro_index]) << "\t"
                      << std::to_string(t_node_spec[node_index].getMZ()) << " node score " << t_node_spec[node_index].getIntensity()
                      << "\t tolspec: " << hi.tol_spec_map_[m][node_index].getIntensity() << " tol: " << t_node_spec[node_index].getMZ() / 1e5 << "\t"
                      << std::to_string(mass_shift) << "\tmod mass: " << std::to_string(mod_mass) << "\t" << mod_count << std::endl;
          }

          if (mod_count > 0 && prev_mod_count != mod_count) prev_mass_shift = mass_shift;
          prev_mod_count = mod_count;
          pre_pro_index = pro_index;
          pre_node_index = node_index;
        }
        if (debug)
        {
          std::cout<<std::endl;
        }
        int mode_score = getScore_(best_path[0]);
        if (m == 1 && hi.protein_start_position_ >= 0 && hi.protein_end_position_ >= 0 && hi.protein_start_position_ >= hi.protein_end_position_)
        {
          if (total_score > mode_score) // mode 0 wins
          {
            hi.protein_start_position_ = -1;
            break;
          }
          else // mode 1 wins
          {
            hi.protein_end_position_ = -1;
            total_score = mode_score;
            used_mode.pop_back();
            used_mode.push_back(m);
            break;
          }
        }
        else
        {
          total_score += mode_score;
          used_mode.push_back(m);
        }
      }
    }
    if (hi.protein_start_position_ >= 0 && hi.protein_end_position_ >= 0 && hi.protein_start_position_ >= hi.protein_end_position_) { continue; }
    if (used_mode.empty()) continue;

    const auto t_mod_masses = mod_masses, t_mod_tols = mod_tols;
    const auto t_mod_starts = mod_starts, t_mod_ends = mod_ends;

    mod_masses.clear();
    mod_starts.clear();
    mod_ends.clear();
    mod_tols.clear();

    for (int k = 0; k < t_mod_masses.size(); k++)
    {
      if (hi.protein_start_position_ >= 0 && t_mod_starts[k] < hi.protein_start_position_) continue;
      if (hi.protein_end_position_ >= 0 && t_mod_ends[k] > hi.protein_end_position_) continue;
      mod_masses.push_back(t_mod_masses[k]);
      mod_starts.push_back(t_mod_starts[k]);
      mod_ends.push_back(t_mod_ends[k]);
      mod_tols.push_back(t_mod_tols[k]);
    }
    std::vector<String> mod_ids, mod_accs;

    for (int k = 0; k < mod_masses.size(); k++)
    {
      auto mod_mass = mod_masses[k];
      auto iter = mod_map_.lower_bound(mod_mass - mod_tols[k]);
      String mod_id = "";
      String mod_acc = "";
      std::set<int> mod_int_acc;
      while (iter != mod_map_.end())
      {
        double diff = iter->first - mod_mass;
        if (diff > mod_tols[k]) break;
        if (diff > -mod_tols[k])
        {
          for (const auto& mod : iter->second)
          {
            if (mod_int_acc.find(mod.getUniModRecordId()) != mod_int_acc.end()) continue;
            mod_int_acc.insert(mod.getUniModRecordId());
            mod_acc += std::to_string(mod.getUniModRecordId()) + ",";
            mod_id += mod.getId() + ",";
          }
        }
        iter++;
      }
      mod_ids.push_back(mod_id);
      mod_accs.push_back(mod_acc);
    }

    // remove unmatched tags.
    std::set<int> to_exclude_tag_indices, matched_positions;

    for (int m = (used_mode.back() == 2 ? 2 : 0); m <= (used_mode.back() == 2 ? 2 : 1); m++)
    {
      std::vector<Size> best_path;
      const auto& t_pro_masses = hi.pro_mass_map_[m];
      if (std::find(used_mode.begin(), used_mode.end(), m) == used_mode.end() || best_path_map.empty() || best_path_map.find(m) == best_path_map.end()
          || best_path_map[m].empty())
        ;
      else
        best_path = best_path_map[m];

      if (std::find(used_mode.begin(), used_mode.end(), m) != used_mode.end())
      {
        for (int pos : matched_position_map[m])
          matched_positions.insert(pos);
      }
      for (int j = 0; j < matched_tags.size(); j++) // for each tag
      {
        auto tag = matched_tags[j];

        if ((tag.getNtermMass() > 0 && m == 0) || (tag.getCtermMass() > 0 && m == 1)) { continue; }
        bool tag_matched = false;
        for (auto iter = best_path.rbegin(); iter != best_path.rend(); iter++) // compare against each path
        {
          auto node_index = getNodeIndex_(*iter, t_pro_masses.size());
          double node_mz = hi.node_spec_map_[m][node_index].getMZ();

          if (tag.getNtermMass() > 0)
          {
            for (const auto& shift : prefix_shifts_)
            {
              double t_mass = tag.getNtermMass() - shift;
              if (std::abs(t_mass - node_mz) > 1.1) continue;
              tag_matched = true;
              break;
            }
          }
          else
          {
            for (const auto& shift : suffix_shifts_)
            {
              double t_mass = tag.getCtermMass() - shift;
              if (m == 2 && hi.calculated_precursor_mass_ > 0)
                t_mass = hi.calculated_precursor_mass_ - i2f_mass - t_mass;
              if (std::abs(t_mass - node_mz) > 1.1) continue;
              tag_matched = true;
              break;
            }
          }

          if (tag_matched)
          {
            std::vector<int> positions;
            std::vector<double> masses;
            String seq = hit.getSequence();

            if (hi.protein_end_position_ >= 0) seq = seq.substr(0, hi.protein_end_position_);
            if (hi.protein_start_position_ >= 0) seq = seq.substr(hi.protein_start_position_);

            FLASHTaggerAlgorithm::fillMatchedPositionsAndFlankingMassDiffs(positions, masses, max_mod_mass_ * max_mod_cntr_ + 1, seq, tag);
            tag_matched = ! positions.empty();
            break;
          }
        }
        if (! tag_matched) to_exclude_tag_indices.insert(tag_indices[j]);
      }
    }

    std::vector<int> refined_tag_indices;
    for (auto index : tag_indices)
    {
      if (to_exclude_tag_indices.find(index) != to_exclude_tag_indices.end()) continue;
      refined_tag_indices.push_back(index);
    }

    if (refined_tag_indices.empty() || total_score <= 0) continue;

    hi.protein_start_position_ += hi.protein_start_position_ >= 0 ? 1 : 0;
    hit.setMetaValue("ModificationIDs", mod_ids); // TODO matching masses vs. all masses?
    hit.setMetaValue("ModificationACCs", mod_accs);
    hit.setMetaValue("Modifications", mod_masses);
    hit.setMetaValue("ModificationStarts", mod_starts);
    hit.setMetaValue("ModificationEnds", mod_ends);
    hit.setMetaValue("MatchedAA", matched_positions.size());
    hit.setMetaValue("TagIndices", refined_tag_indices);

    double protein_len = hit.getSequence().size();
    if (hi.protein_end_position_ > 0) { protein_len -= (protein_len - hi.protein_end_position_); }
    if (hi.protein_start_position_ > 0) { protein_len -= hi.protein_start_position_ - 1; }

    hit.setCoverage((double)matched_positions.size() / protein_len);
    hit.setScore(total_score);
    hit.setMetaValue("StartPosition", hi.protein_start_position_);
    hit.setMetaValue("EndPosition", hi.protein_end_position_);
    hit.setMetaValue("GivenMass", given_precursor_mass_);
    hit.setMetaValue("Mass", hi.calculated_precursor_mass_);
    hit.setMetaValue("RT", dspec.getOriginalSpectrum().getRT());
    hit.setMetaValue("NumMass", dspec.size());
    hit.setMetaValue("PrecursorScore", dspec.getPrecursorPeakGroup().getQscore2D());//dspec.getPrecursorPeakGroup().getChargeSNR(dspec.getPrecursor().getCharge())
    hit.setMetaValue("PrecursorSNR", dspec.getPrecursorPeakGroup().getSNR()); // it should be charge SNR in theory. But mass SNR was written for ease of coding.
    hit.setMetaValue("ProteoformMassByFragmentMass", precursor_by_fragment? 1 : 0);
    // hit.setMetaValue("Proforma", string)
#pragma omp critical
    {
      bool insert = true;

      if (! multiple_hits_per_spec && ! proteoform_hits_.empty()) // when multiple hits are not allowed
      {
        if (proteoform_hits_.back().getScore() >= hit.getScore()) insert = false;
        else
          proteoform_hits_.pop_back();
      }
      if (insert) { proteoform_hits_.push_back(hit); }
    }
  }
  endProgress();
}

void FLASHExtenderAlgorithm::constructDAG_(std::set<Size>& sinks,
                                           HitInformation& hi,
                                           const std::vector<std::vector<int>>& tag_edges,
                                           int max_mod_cntr_for_last_mode,
                                           bool use_tags)
{
  Size src = getVertex_(0, 0, 0, 0, hi.pro_mass_map_[hi.mode_].size());
  hi.visited_[src] = true;
  std::set<Size> visited_tag_edges;
  std::map<Size, std::set<std::pair<double, double>>> sink_map;
  std::map<Size, std::map<Size, int>> node_max_score_map; // node, cumulative mass, score

  connectBetweenTags_(visited_tag_edges, hi, sink_map, src, 0, 0, node_max_score_map, tag_edges, max_mod_cntr_for_last_mode, use_tags);

  for (const auto& sink : sink_map)
  {
    sinks.insert(sink.first);
    //    if (hi.mode_ == 2)
    //    {
    //      const auto & [t, c] = sink.second;
    //      std::cout<< t << " " << c << std::endl;
    //      std::cout<< std::to_string(hi.node_spec_map_[2][getNodeIndex_(sink.first, hi.pro_mass_map_[2].size())].getMZ() -
    //      hi.pro_mass_map_[2][getProIndex_(sink.first, hi.pro_mass_map_[2].size())]) << std::endl;
    //    }
  }
}

void FLASHExtenderAlgorithm::connectBetweenTags_(std::set<Size>& visited_tag_edges,
                                                 HitInformation& hi,
                                                 std::map<Size, std::set<std::pair<double, double>>>& sinks,
                                                 Size vertex,
                                                 const double truncation_mass,
                                                 const double cumulative_mod_mass,
                                                 std::map<Size, std::map<Size, int>>& node_max_score_map,
                                                 const std::vector<std::vector<int>>& tag_edges,
                                                 int max_mod_cntr_for_last_mode,
                                                 bool use_tags)
{
  const auto& pro_masses = hi.pro_mass_map_[hi.mode_];
  int node_index = getNodeIndex_(vertex, pro_masses.size());
  int pro_index = getProIndex_(vertex, pro_masses.size());

  int tag_start_index = -1;
  int tag_end_index = -1;

  const auto& tag_pro_starts = tag_edges[0];
  const auto& tag_pro_ends = tag_edges[1];
  const auto& tag_node_starts = tag_edges[2];
  const auto& tag_node_ends = tag_edges[3];

  int max_node_start = -1;
  for (int tag_node_start : tag_node_starts) // for all reachable tag starting point, run extension
  {
    max_node_start = std::max(max_node_start, tag_node_start);
  }

  for (int i = 0; i < (int)tag_pro_starts.size(); i++)
  {
    if (tag_start_index < 0 && tag_node_starts[i] == node_index && tag_pro_starts[i] == pro_index) { tag_start_index = i; }
    if (tag_end_index < 0 && tag_node_ends[i] == node_index && tag_pro_ends[i] == pro_index) { tag_end_index = i; }
  }

  Size src = getVertex_(0, 0, 0, 0, pro_masses.size());

  if (tag_start_index >= 0) // within tag
  {
    int node_end = -1;
    Size i = tag_start_index;
    int pro_end = -1;
    while ((i < tag_node_starts.size()) && (tag_node_starts[i] == node_index) && (tag_pro_starts[i] == pro_index))
    {
      if (pro_end < tag_pro_ends[i]) pro_end = tag_pro_ends[i];
      i++;
    }
    i = tag_start_index;
    while ((i < tag_node_starts.size()) && (tag_node_starts[i] == node_index) && (tag_pro_starts[i] == pro_index))
    {
      if (pro_end == tag_pro_ends[i]) { node_end = std::max(tag_node_ends[i], node_end); }
      i++;
    }

    std::map<Size, std::set<std::pair<double, double>>> next_vertices;
    extendBetweenTags_(next_vertices, hi, vertex, node_end, pro_end, 0, truncation_mass, cumulative_mod_mass, node_max_score_map,
                       max_mod_cntr_for_last_mode);

    std::map<std::set<std::pair<double, double>>, Size> mass_sink;
    for (const auto& [next_vertex, next_cumulative_shift] : next_vertices)
    {
      if (next_vertex == src) continue;
      if (mass_sink.find(next_cumulative_shift) == mass_sink.end() || getScore_(mass_sink[next_cumulative_shift]) < getScore_(next_vertex))
        mass_sink[next_cumulative_shift] = next_vertex;
    }

    for (const auto& [mass_set, next_vertex] : mass_sink)
    {
      for (const auto& masses : mass_set)
      {
        const auto& [t, c] = masses;
        connectBetweenTags_(visited_tag_edges, hi, sinks, next_vertex, t, c, node_max_score_map, tag_edges, max_mod_cntr_for_last_mode, use_tags);
      }
    }
  }

  if (vertex == src || tag_end_index >= 0 || max_node_start < 0) // between tag.
  {
    std::set<Size> reachable_vertices;

    for (Size tag_index = 0; tag_index < tag_node_starts.size(); tag_index++) // for all reachable tag starting point, run extension
    {
      int node_start = tag_node_starts[tag_index];
      int pro_start = tag_pro_starts[tag_index];
      int mod_num = getModNumber_(vertex);
      if (pro_index > pro_start) continue;
      if (max_node_start >= 0)
      {
        if (node_index > node_start) continue;
      }
      else if (node_start < 0)
        continue;

      bool is_visited_start = false;

      for (int mn = 0; mn <= mod_num; mn++)
      {
        if (visited_tag_edges.find(getVertex_(max_node_start >= 0 ? node_start : 0, pro_start, 0, mn, pro_masses.size())) != visited_tag_edges.end())
        {
          is_visited_start = true;
          break;
        }
      }

      if (is_visited_start) continue;
      visited_tag_edges.insert(getVertex_(max_node_start >= 0 ? node_start : 0, pro_start, 0, mod_num, pro_masses.size()));

      std::map<Size, std::set<std::pair<double, double>>> next_vertices;
      extendBetweenTags_(next_vertices, hi, vertex, node_start, pro_start, 0, truncation_mass, cumulative_mod_mass, node_max_score_map,
                         max_mod_cntr_for_last_mode);

      if (node_start < 0)
      {
        sinks = next_vertices;
        return;
      }

      std::map<std::set<std::pair<double, double>>, Size> mass_sink;

      for (const auto& [next_vertex, next_cumulative_shift] : next_vertices)
      {
        if (next_vertex == src) continue;
        if (mass_sink.find(next_cumulative_shift) == mass_sink.end() || getScore_(mass_sink[next_cumulative_shift]) < getScore_(next_vertex))
          mass_sink[next_cumulative_shift] = next_vertex;
      }

      for (const auto& [mass_set, next_vertex] : mass_sink)
      {
        for (const auto& masses : mass_set)
        {
          const auto& [t, c] = masses;
          connectBetweenTags_(visited_tag_edges, hi, sinks, next_vertex, t, c, node_max_score_map, tag_edges, max_mod_cntr_for_last_mode, use_tags);
        }
        reachable_vertices.insert(next_vertex);
      }
    }
    if ((vertex != src || ! use_tags) && reachable_vertices.empty())
    {
      if (! use_tags)
      {
        std::map<int, int> diff_count;
        for (const auto& p : hi.node_spec_map_.at(hi.mode_))
        {
          for (const auto& m : pro_masses)
          {
            int diff = (int)round((m - p.getMZ())); // can do log transform to reflect ppm error later.
            if (diff_count.find(diff) == diff_count.end()) { diff_count[diff] = 0; }
            diff_count[diff] += (int)p.getIntensity();
          }
        }

        std::priority_queue<std::pair<int, int>> maxHeap;

        for (const auto& entry : diff_count)
        {
          // Push pairs into the heap with count as the key (max heap based on count)
          maxHeap.push({entry.second, entry.first});
        }

        // Collect the top 3 most frequent differences
        std::vector<int> top_diffs;
        for (int i = 0; i < 1 && ! maxHeap.empty(); ++i)
        {
          top_diffs.push_back(maxHeap.top().second); // Get the difference
          maxHeap.pop();
        }
      }

      if (hi.mode_ != 2)
      {
        extendBetweenTags_(sinks, hi, vertex, -1, -1, use_tags ? 1e5 : 0, truncation_mass, cumulative_mod_mass, node_max_score_map,
                           max_mod_cntr_for_last_mode);
      }
      else
      {
        for (int j = 0; j < pro_masses.size(); j++)
        {
          if (std::abs(hi.calculated_precursor_mass_ + truncation_mass - pro_masses[j]) > max_mod_mass_ * (max_mod_cntr_ - getModNumber_(vertex)))
            continue;

          extendBetweenTags_(sinks, hi, vertex, hi.node_spec_map_[2].size() - 1, j, 0, truncation_mass, cumulative_mod_mass, node_max_score_map,
                             max_mod_cntr_for_last_mode);
        }
      }
    }
  }
}

void FLASHExtenderAlgorithm::extendBetweenTags_(std::map<Size, std::set<std::pair<double, double>>>& sinks,
                                                HitInformation& hi,
                                                Size start_vertex,
                                                const int end_node_index,
                                                int end_pro_index,
                                                int diagonal_counter,
                                                const double truncation_mass,
                                                const double cumulative_mod_mass,
                                                std::map<Size, std::map<Size, int>>& node_max_score_map,
                                                const int max_mod_cntr_for_last_mode)
{
  // TODO N term mod vs. 1st amino acid mod distinction
  if (! hi.visited_[start_vertex]) return;
  const auto& pro_masses = hi.pro_mass_map_[hi.mode_];
  const auto pro_mass_size = pro_masses.size();
  int max_mod_cntr = max_mod_cntr_for_last_mode >= 0 ? max_mod_cntr_for_last_mode : max_mod_cntr_;
  int start_node_index = getNodeIndex_(start_vertex, pro_mass_size);
  int start_pro_index = getProIndex_(start_vertex, pro_mass_size);
  int start_score = getScore_(start_vertex);
  int start_num_mod = getModNumber_(start_vertex);
  if (start_num_mod == max_mod_cntr) diagonal_counter = 1e5;

  auto src = getVertex_(0, 0, 0, 0, pro_mass_size);
  const auto& node_spec = hi.node_spec_map_.at(hi.mode_);
  const auto& tol_spec = hi.tol_spec_map_.at(hi.mode_);

  if (end_node_index < 0) //
  {
    end_pro_index = 0;
    if (hi.protein_end_position_ >= 0)
      while (end_pro_index < pro_mass_size
             && pro_masses[end_pro_index] - pro_masses[hi.protein_end_position_ - 1] < max_mod_mass_ * (max_mod_cntr - start_num_mod) - 1)
        end_pro_index++;
    else
      end_pro_index = ((start_num_mod < max_mod_cntr) && diagonal_counter == 0)
                        ? start_pro_index + max_extension_stretch_
                        : ((int)pro_mass_size - 1); // if sink is not specified, stretch up to max_extension_stretch_ amino acids.
    end_pro_index = std::min(end_pro_index, (int)pro_mass_size - 1);
  }
  // make the range of truncation well...  make use of the positional information
  if (start_vertex == src)
  {
    {
      for (int pro_i = start_pro_index + 1; pro_i <= end_pro_index; pro_i++) // change later
      {
        if (hi.protein_start_position_ >= 0
            && pro_masses[pro_i] - pro_masses[hi.protein_start_position_] > max_mod_mass_ * (max_mod_cntr - start_num_mod) + 1.1)
          break;
        if (hi.protein_start_position_ >= 0
            && pro_masses[hi.protein_start_position_] - pro_masses[pro_i] > max_mod_mass_ * (max_mod_cntr - start_num_mod) + 1.1)
          continue;

        if (end_node_index >= 0
            && std::abs(node_spec[end_node_index].getMZ() - pro_masses[end_pro_index] + pro_masses[pro_i])
                 > max_mod_mass_ * (max_mod_cntr - start_num_mod) + tol_spec[end_node_index].getIntensity() + 1.1)
        {
          continue;
        }
        Size vertex2 = getVertex_(0, pro_i, 0, 0, pro_mass_size); //
        bool connected = hi.dag_.addEdge(vertex2, start_vertex, hi.visited_);

        if (vertex2 >= hi.dag_.size() || ! connected) continue;
        extendBetweenTags_(sinks, hi, vertex2, end_node_index, end_pro_index, diagonal_counter, pro_masses[pro_i], cumulative_mod_mass,
                           node_max_score_map, max_mod_cntr_for_last_mode);
      }
    }
  }

  // double start_delta_mass = start_node_mass - pro_masses[start_pro_index];
  double end_node_mass = node_spec[end_node_index].getMZ();
  double end_delta_mass = end_node_mass - pro_masses[end_pro_index];
  if (end_node_index >= 0)
  {
    double margin = tol_spec[end_node_index].getIntensity(); // tol_spec[end_node_index].getIntensity();
    if (std::abs(end_delta_mass - cumulative_mod_mass + truncation_mass) > max_mod_mass_ * (max_mod_cntr - start_num_mod) + margin) { return; }
    if (std::abs(end_delta_mass - cumulative_mod_mass + truncation_mass) > margin)
    {
      if (diagonal_counter > 0) return; //
    }
    else
      diagonal_counter = 1e5; // if the start and end points make a diagonal line, go through the diagonal line.
  }

  bool same_score_node = false;
  Size key_to_score_map = start_pro_index * pro_mass_size + end_pro_index; // - truncation_mass; //
  for (int nm = 0; nm <= start_num_mod; nm++)
  {
    Size zero_score_vertex = getVertex_(start_node_index, start_pro_index, 0, nm, pro_mass_size);
    if (node_max_score_map.find(zero_score_vertex) != node_max_score_map.end())
    {
      auto& kts = node_max_score_map[zero_score_vertex];
      if (kts.find(key_to_score_map) != kts.end())
      {
        if (kts[key_to_score_map] > start_score) { return; }
        else if (key_to_score_map != 0 && kts[key_to_score_map] == start_score)
          same_score_node = true;
      }
    }
  }

  if (start_node_index == end_node_index && start_pro_index == end_pro_index && start_vertex != src)
  {
    sinks[start_vertex].insert(std::pair<double, double> {truncation_mass, cumulative_mod_mass});
    return;
  }

  if (end_node_index < 0)
  {
    if (start_vertex != src)
    {
      sinks[start_vertex].insert(std::pair<double, double> {truncation_mass, cumulative_mod_mass});
    }
  }
  else if (start_node_index == node_spec.size() - 1)
  {
    sinks[start_vertex].insert(std::pair<double, double> {truncation_mass, cumulative_mod_mass});
    return;
  }
  else if (start_node_index > end_node_index || start_pro_index > end_pro_index)
    return;

  if (start_num_mod > 0 && same_score_node) return; // exclude truncation

  node_max_score_map[getVertex_(start_node_index, start_pro_index, 0, start_num_mod, pro_mass_size)][key_to_score_map] = start_score;

  for (int node_i = start_node_index + 1; node_i <= (end_node_index < 0 ? ((int)node_spec.size() - 1) : end_node_index); node_i++)
  {
    int score = start_score + (int)node_spec[node_i].getIntensity();

    double t_node_mass = node_spec[node_i].getMZ();
    double t_margin = tol_spec[node_i].getIntensity();

    for (int pro_i = start_pro_index; pro_i <= end_pro_index; pro_i++) //
    {
      double t_delta_mass = t_node_mass - pro_masses[pro_i];
      double delta_delta = t_delta_mass - cumulative_mod_mass + truncation_mass;
      if (delta_delta > max_mod_mass_ + t_margin) continue;
      if (-delta_delta > max_mod_mass_ + t_margin) break;

      int num_mod = start_num_mod;
      int next_score = score;
      double next_cumulative_mod_mass = cumulative_mod_mass;

      if (std::abs(delta_delta) > t_margin)
      {
        if (pro_i == start_pro_index) continue;
        if (std::abs(delta_delta) < 0.036386 - t_margin) continue;
        if (std::abs(delta_delta) > 0.036386 + t_margin && std::abs(delta_delta) < 0.947630 - t_margin) continue;
        num_mod++;
        if (num_mod > max_mod_cntr) continue;
        if (diagonal_counter > 0) continue; //

        next_cumulative_mod_mass = t_delta_mass + truncation_mass;

        if (hi.mode_ == 2)
        {
          for (int pi = std::max(pro_i, hi.protein_end_position_ - 5); pi < std::min(hi.protein_end_position_ + 5, (int)pro_mass_size); pi++)
          {
            double corrected_mod_mass
              = hi.calculated_precursor_mass_ - i2f_mass - pro_masses[pi] + truncation_mass;
            double delta = corrected_mod_mass - next_cumulative_mod_mass;
            if (delta > 1.1) continue;
            if (delta < -1.1) break;

            next_cumulative_mod_mass = corrected_mod_mass;
            break;
          }
          if (std::abs(t_delta_mass - next_cumulative_mod_mass + truncation_mass) > t_margin) continue;
        }

        next_score -= 1 + 2 * (multi_ion_score + FLASHTaggerAlgorithm::max_node_score); // at least two best scoring peaks should exist
        auto iter = mod_map_.lower_bound(delta_delta - t_margin);
        if (std::abs(delta_delta - iter->first) > t_margin)
        {
          next_score -= multi_ion_score + FLASHTaggerAlgorithm::max_node_score;
        } // if it was not found in unimod, subtract more
      }

      next_score = std::min(next_score, max_path_score_);
      if (next_score < min_path_score_) continue;
      auto next_vertex = getVertex_(node_i, pro_i, next_score, num_mod, pro_mass_size);
      if (next_vertex >= hi.dag_.size()) continue;

      if (next_score >= max_path_score_ && ! sinks.empty() && hi.visited_[next_vertex]) continue;
      if (!hi.dag_.addEdge(next_vertex, start_vertex, hi.visited_)) continue;

      int next_diagonal_counter = diagonal_counter;
      if (diagonal_counter > 0) next_diagonal_counter--;
      else if (num_mod != start_num_mod)
        next_diagonal_counter = 1; // at least one stretch after modification

      extendBetweenTags_(sinks, hi, next_vertex, end_node_index, end_pro_index, next_diagonal_counter, truncation_mass, next_cumulative_mod_mass,
                         node_max_score_map, max_mod_cntr_for_last_mode);
    }
  }
}
} // namespace OpenMS