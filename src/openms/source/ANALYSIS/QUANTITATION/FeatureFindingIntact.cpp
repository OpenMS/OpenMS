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
// $Maintainer: Jihyung Kim $
// $Authors: Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/QUANTITATION/FeatureFindingIntact.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>

#include <queue>

// test purpose
#include <iostream>
#include <fstream>
#include <sstream>

namespace OpenMS
{
  FeatureFindingIntact::FeatureFindingIntact():
      ProgressLogger()
  {
    this->setLogType(CMD);
  }

  FeatureFindingIntact::~FeatureFindingIntact(){}


  void getCoordinatesOfFeaturesForPython(std::vector<FeatureFindingIntact::FeatureHypothesis> hypo)
  {
    String out_path = "/Users/jeek/Documents/A4B_UKE/filg/rp/centroid/190904_Proswift50mm_7min_Filg_500ng.pphr.featHypo.tsv";
    ofstream out;
    out.open(out_path, ios::out);

    // header
    out << "feature_label\tscore\tbounding_box_pos\tbounding_box_width\tbounding_box_height\t"
           "iso_position\tmasstrace_centroid_rts\tmasstrace_centroid_mzs\n";

    for (const auto& feat : hypo)
    {
      double mz_upper_limit = 0.0;
      double rt_upper_limit  = 0.0;
      double mz_lower_limit = std::numeric_limits<double>::max();
      double rt_lower_limit = std::numeric_limits<double>::max();
      stringstream rts;
      stringstream mzs;
      for (const auto& mt : feat.getMassTraces())
      {
        const auto& boundingbox = mt->getConvexhull().getBoundingBox();

        if (boundingbox.minX() < rt_lower_limit)
          rt_lower_limit = boundingbox.minX();
        if (boundingbox.maxX() > rt_upper_limit)
          rt_upper_limit = boundingbox.maxX();
        if (boundingbox.minY() < mz_lower_limit)
          mz_lower_limit = boundingbox.minY();
        if (boundingbox.maxY() > mz_upper_limit)
          mz_upper_limit = boundingbox.maxY();

        mzs << mt->getCentroidMZ() << ",";
        rts << mt->getCentroidRT() << ",";
      }
      std::string centroids = rts.str();
      centroids.pop_back();
      centroids = centroids + "\t" + mzs.str();
      centroids.pop_back();

      stringstream isos;
      for (const auto& pair : feat.getIndicesOfMassTraces())
      {
        isos << pair.first << ",";
      }
      std::string iso_str = isos.str();
      iso_str.pop_back();

      String label = std::to_string(feat.getFeatureMass()) + "(cs" + std::to_string(feat.getCharge()) + ")";
      out << label << "\t" << to_string(feat.getScore()) << "\t"
          << to_string(rt_lower_limit) << "," << to_string(mz_lower_limit) << "\t"
          << to_string((rt_upper_limit-rt_lower_limit)) << "\t" << to_string((mz_upper_limit-mz_lower_limit)) << "\t"
          << iso_str << "\t" << centroids + "\n";
    }
    out.close();
  }

  void FeatureFindingIntact::run(std::vector<MassTrace> &input_mtraces, FeatureMap &output_featmap)
  {
    // *********************************************************** //
    // Step 1 building FeatureHypotheses
    // *********************************************************** //
    // output feature hypotheses
    std::vector<FeatureHypothesis> feat_hypos;

    // for tracking shared masstraces between features with mass traces indices
    // link masstrace_index -> hypothesis index
    std::vector<std::vector<Size>> shared_m_traces_indices(input_mtraces.size(), std::vector<Size>());

    // prepare isotope model for scoring method
    setAveragineModel();

    buildFeatureHypotheses_(input_mtraces, feat_hypos, shared_m_traces_indices);

    // *********************************************************** //
    // Step 2 clustering FeatureHypotheses
    // *********************************************************** //
//    getCoordinatesOfFeaturesForPython(feat_hypos);
    clusterFeatureHypotheses_(feat_hypos, shared_m_traces_indices);

    // *********************************************************** //
    // Step 3 resolving conflicts
    // *********************************************************** //
  }

  void FeatureFindingIntact::buildFeatureHypotheses_(std::vector<MassTrace>& input_mtraces,
                                                     std::vector<FeatureHypothesis>& output_hypotheses,
                                                     std::vector<std::vector<Size>>& shared_m_traces_indices) const
  {
    if (input_mtraces.empty())
    {
      return;
    }
    this->startProgress(0, input_mtraces.size(), "assembling mass traces to candidate features");

    // *********************************************************** //
    // Step 1 Preparation
    // *********************************************************** //
    // mass traces must be sorted by their centroid MZ
    std::sort(input_mtraces.begin(), input_mtraces.end(), CmpMassTraceByMZ());

    // total_intensity is needed for calculating feature hypothesis scores
    double total_intensity(0.0);
    for (Size i = 0; i < input_mtraces.size(); ++i)
    {
      total_intensity += input_mtraces[i].getIntensity(use_smoothed_intensities_);
    }

    // temporary vector
    std::vector<FeatureHypothesis> candidate_hypotheses;

    // *********************************************************** //
    // Step 2 Iterate through all mass traces to find likely matches
    // and generate isotopic / charge hypotheses
    // *********************************************************** //
    Size progress(0);
//    Size m_trace_size = input_mtraces.size();
//    Size progress_chunk = std::floor(m_trace_size/100);
    for (SignedSize i = 0; i < (SignedSize)input_mtraces.size(); ++i)
    {
        this->setProgress(progress);

      ++progress;

      std::vector<std::pair<const MassTrace*, Size>> local_traces;
      double ref_trace_mz(input_mtraces[i].getCentroidMZ());
      double ref_trace_rt(input_mtraces[i].getCentroidRT());

      local_traces.push_back(std::make_pair(&input_mtraces[i], i));

      for (Size ext_idx = i + 1; ext_idx < input_mtraces.size(); ++ext_idx)
      {
        // traces are sorted by m/z, so we can break when we leave the allowed window
        double diff_mz = std::fabs(input_mtraces[ext_idx].getCentroidMZ() - ref_trace_mz);
        if (diff_mz > local_mz_range_) break;

        double diff_rt = std::fabs(input_mtraces[ext_idx].getCentroidRT() - ref_trace_rt);
        if (diff_rt <= local_rt_range_)
        {
          // std::cout << " accepted!" << std::endl;
          local_traces.push_back(std::make_pair(&input_mtraces[ext_idx], ext_idx));
        }
      }
      findLocalFeatures_(local_traces, total_intensity, candidate_hypotheses, shared_m_traces_indices);
    }
    this->endProgress();
    OPENMS_LOG_INFO << "feature hypotheses size:" << candidate_hypotheses.size() << std::endl;

    // *********************************************************** //
    // Step 3 add charge score
    // *********************************************************** //
    // TODO: pregress logger?

    Size candi_size = candidate_hypotheses.size();
//    std::vector<double> feature_masses; // keeping track of deconvoluted feature masses
//    feature_masses.reserve(candi_size);
    std::vector<std::pair<double, Size>> feat_and_charges; // first : feature mass, second : charge state
    feat_and_charges.reserve(candi_size);

    // collect candidate masses
    for (const auto& hypo : candidate_hypotheses)
    {
//      feature_masses.push_back(hypo.getFeatureMass());
      feat_and_charges.push_back(std::make_pair(hypo.getFeatureMass(), hypo.getCharge()));
    }
//    std::sort(feature_masses.begin(), feature_masses.end());
    std::sort(feat_and_charges.begin(), feat_and_charges.end());

    // cluster mean_masses within tolerance
    std::vector<std::pair<double, std::set<Size>>> collected_masses; // first : average mass of current cluster, second : vector of charges
    collected_masses.reserve(candi_size);
    double tolerance = 1.5; // Da
    Size max_nr_charges = 0; // get the maximum size of charge set

    for(auto f_iter=feat_and_charges.begin(); f_iter != feat_and_charges.end(); f_iter++ )
    {
      auto f_iter_end = f_iter;

      double average = 0.0;
      int counter = 1;
      for(; f_iter_end!=feat_and_charges.end(); f_iter_end++)
      {
        average = 0.0;
        for (auto i=f_iter; i!=(f_iter_end+1); i++) average += i->first;
        average /= counter;
        if ((f_iter_end+1)->first-average > tolerance)
        {
          break;
        }
        counter++;
      }

      // save found charge states
      std::set<Size> tmp_charges;
      for(auto it=f_iter; it!=(f_iter_end+1); it++)
      {
        tmp_charges.insert(it->second);
      }
      collected_masses.push_back(std::make_pair(average, tmp_charges));
      if (tmp_charges.size() > max_nr_charges)
      {
        max_nr_charges = tmp_charges.size();
      }
      f_iter = f_iter_end;
//      OPENMS_LOG_INFO << std::distance(feat_and_charges.end(), f_iter) << endl;
    }
    collected_masses.shrink_to_fit();

    // give charge score for each feature
    for (auto& feature : candidate_hypotheses)
    {
      double cs_score=0;
      for (auto f_it = collected_masses.begin(); f_it != collected_masses.end(); f_it++)
      {
        if (std::fabs(f_it->first-feature.getFeatureMass()) > tolerance)
          continue;
        cs_score = f_it->second.size() / max_nr_charges;
        break;
      }
      feature.setChargeScore(cs_score);
    }

//    {
////      bool is_artifact = false;
//
//      bool found = false;
//      double matched_mass;
//      for (const auto& feat: feat_and_charges)
//      {
//        if ((feat.first-curr_feat.getFeatureMass()) > tolerance) break; // to large
//        if ((curr_feat.getFeatureMass()-feat.first) > tolerance) continue;  // to small
//
//        found = true;
//        if (matched_mass)
//      }
//
//      if (!found){
//        std::vector<Size> tmp_vec;
//        tmp_vec.push_back(curr_feat.getCharge());
//        feat_and_charges.push_back(std::make_pair(curr_feat.getFeatureMass(), tmp_vec));
//        std::sort(feat_and_charges.begin(), feat_and_charges.end());
//      }


//      for (Size j = i+1; j < candi_size; j++)
//      {
//        // check if low harmonic artifact
//        // check if high harmonic artifact
//        // check if charge-off-by-one
//      }
//      if (!is_artifact)
//      {
//        // add to collect charge states of features
//        feature_masses.push_back(curr_feat.getFeatureMass());
//        feature_masses.push_back(curr_feat.getFeatureMass());
//
//        // add to output
//        output_hypotheses.push_back(curr_feat);
//      }
//    }

    output_hypotheses = candidate_hypotheses;
  }

  void FeatureFindingIntact::findLocalFeatures_(const std::vector<std::pair<const MassTrace*, Size>>& candidates,
                                                const double total_intensity,
                                                std::vector<FeatureHypothesis>& output_hypotheses,
                                                std::vector<std::vector<Size>>& shared_m_traces_indices) const
  {
    // not storing hypothesis with only one mass trace (with only mono), while FeatureFindingMetabo does

    // compute maximum m/z window size

    for (Size charge = charge_lower_bound_; charge <= charge_upper_bound_; ++charge)
    {
      FeatureHypothesis fh_tmp;
      fh_tmp.addMassTrace(*candidates[0].first); // ref_mtrace (which is mono here)
      fh_tmp.setScore((candidates[0].first->getIntensity(use_smoothed_intensities_)) / total_intensity);
      fh_tmp.setCharge(charge);
      double mol_weight = (candidates[0].first->getCentroidMZ()-Constants::PROTON_MASS_U) * charge;
      if (mol_weight > mass_upper_bound)
        break;
      fh_tmp.setFeatureMass(mol_weight);

      // for shared mass traces tracker
      std::vector<std::pair<Size, Size>> used_mass_trace_indices; // first: iso index of current feature, second: mt index
      used_mass_trace_indices.push_back(std::make_pair(0, candidates[0].second));

      std::vector<double> feature_iso_intensities; // collection of selected isos (if nothing's in such iso index, intensity of it is 0)
      feature_iso_intensities.push_back(candidates[0].first->getIntensity(use_smoothed_intensities_));

      // expected m/z window for iso_pos -> 13C isotope peak position
      double mz_window = Constants::C13C12_MASSDIFF_U * max_nr_traces_ / charge;

      // calculate averagine isotope distribution here (based on mono trace)
      auto iso_dist = iso_model_.get(mol_weight);
      Size iso_size = iso_dist.size();
      double iso_norm = iso_model_.getNorm(mol_weight); // l2 norm
      int iso_index_in_avg_model = 0;

      Size last_iso_idx(0); // largest index of found iso index
      for (Size iso_pos = 1; iso_pos < iso_size; ++iso_pos)
      {
        // Find mass trace that best agrees with current hypothesis of charge & isotopic position
        double best_so_far(0.0);
        Size best_idx(0);
        for (Size mt_idx = last_iso_idx + 1; mt_idx < candidates.size(); ++mt_idx)
        {
          // if out of mz_window, pass this mass trace
          if(std::fabs(candidates[0].first->getCentroidMZ() - candidates[mt_idx].first->getCentroidMZ()) >= mz_window)
            break; // candidates is sorted. So if reached here, next candidate definitely exceed mz_window

          // Score current mass trace candidates against hypothesis
          double rt_score(scoreRT_(*candidates[0].first, *candidates[mt_idx].first));
          double mz_score(scoreMZ_(*candidates[0].first, *candidates[mt_idx].first, iso_pos, charge));

          if (rt_score <= 0.0 | mz_score <= 0)
          {
            continue;
          }

          double int_score(1.0);
          int offset; // getting isotope index against iso_dist
          std::vector<double> tmp_ints(feature_iso_intensities); // intensities up to the last isotope
          tmp_ints.push_back(candidates[mt_idx].first->getIntensity(use_smoothed_intensities_));
          int_score = computeAveragineCosineSimScore_(tmp_ints, iso_dist, iso_size, iso_norm, offset);

          double total_pair_score(0.0);
          if (int_score > 0.0)
          {
            total_pair_score = std::exp(std::log(rt_score) + log(mz_score) + log(int_score));
          }
          if (total_pair_score > best_so_far)
          {
            iso_index_in_avg_model = -offset;
            best_so_far = total_pair_score;
            best_idx = mt_idx;
          }
        } // end of mt_idx

        // Store mass trace that best agrees with current hypothesis of charge and isotopic position
        if (best_so_far > 0.0)
        {
          double weighted_score(((candidates[best_idx].first->getIntensity(use_smoothed_intensities_)) * best_so_far) / total_intensity);
          fh_tmp.setScore(fh_tmp.getScore() + weighted_score);
          fh_tmp.addMassTrace(*candidates[best_idx].first);

          // save up information for trackers
          used_mass_trace_indices.push_back(std::make_pair(iso_pos, candidates[best_idx].second));
          feature_iso_intensities.push_back(candidates[best_idx].first->getIntensity(use_smoothed_intensities_));

          // modify isotope index for this feature
          auto tmp_indices(used_mass_trace_indices);
          for (Size i = 0; i< used_mass_trace_indices.size(); i++)
          {
            tmp_indices[i] = std::make_pair(tmp_indices[i].first + iso_index_in_avg_model, tmp_indices[i].second);
          }
          fh_tmp.setIndicesOfMassTraces(tmp_indices);

          // update feature mass with modified isotope index
          if (iso_index_in_avg_model != 0)
          {
            fh_tmp.updateFeatureMass();
          }

          // update shared mass trace tracker
          Size h_index = output_hypotheses.size();// current hypothesis's index
          for (const auto& mt : used_mass_trace_indices )
          {
            shared_m_traces_indices[mt.second].push_back(h_index);
          }

          output_hypotheses.push_back(fh_tmp);
          last_iso_idx = best_idx;
        }
        else
        {
          feature_iso_intensities.push_back(0.);
          // if too many missing iso positions, break; TODO: find an elegant way...
          if (feature_iso_intensities.size() > 3 && feature_iso_intensities.rbegin()[0] == 0
              && feature_iso_intensities.rbegin()[1] == 0 && feature_iso_intensities.rbegin()[2] == 0)
          {
            break;
          }
        }
      } // end for iso_pos
    } // end of charge
  } // end of findLocalFeatures_(...)

  double FeatureFindingIntact::scoreRT_(const MassTrace& tr1, const MassTrace& tr2) const
  {
    std::map<double, std::vector<double> > coinciding_rts;

    std::pair<Size, Size> tr1_fwhm_idx(tr1.getFWHMborders());
    std::pair<Size, Size> tr2_fwhm_idx(tr2.getFWHMborders());

    double tr1_length(tr1.getFWHM());
    double tr2_length(tr2.getFWHM());
    double max_length = (tr1_length > tr2_length) ? tr1_length : tr2_length;

    // Extract peak shape between FWHM borders for both peaks
    for (Size i = tr1_fwhm_idx.first; i <= tr1_fwhm_idx.second; ++i)
    {
      coinciding_rts[tr1[i].getRT()].push_back(tr1[i].getIntensity());
    }
    for (Size i = tr2_fwhm_idx.first; i <= tr2_fwhm_idx.second; ++i)
    {
      coinciding_rts[tr2[i].getRT()].push_back(tr2[i].getIntensity());
    }

    // Look at peaks at the same RT
    std::vector<double> x, y, overlap_rts;
    for (std::map<double, std::vector<double> >::const_iterator m_it = coinciding_rts.begin(); m_it != coinciding_rts.end(); ++m_it)
    {
      if (m_it->second.size() == 2)
      {
        x.push_back(m_it->second[0]);
        y.push_back(m_it->second[1]);
        overlap_rts.push_back(m_it->first);
      }
    }

    double overlap(0.0);
    if (overlap_rts.size() > 0)
    {
      double start_rt(*(overlap_rts.begin())), end_rt(*(overlap_rts.rbegin()));
      overlap = std::fabs(end_rt - start_rt);
    }

    double proportion(overlap / max_length);
    if (proportion < 0.7)
    {
      return 0.0;
    }
    return computeCosineSim_(x, y);
  }

  double FeatureFindingIntact::scoreMZ_(const MassTrace& tr1, const MassTrace& tr2, Size iso_pos, Size charge) const
  {
    double diff_mz(std::fabs(tr2.getCentroidMZ() - tr1.getCentroidMZ()));

    double mt_sigma1(tr1.getCentroidSD());
    double mt_sigma2(tr2.getCentroidSD());
    double mt_variances(std::exp(2 * std::log(mt_sigma1)) + std::exp(2 * std::log(mt_sigma2)));

    double mz_score(0.0);
    /// mz scoring by expected mean w/ C13
    double mu = (Constants::C13C12_MASSDIFF_U * iso_pos) / charge; // using '1.0033548378'
    double sd = (0.0016633 * iso_pos - 0.0004751) / charge;
    double sigma_mult(3.0);

    //standard deviation including the estimated isotope deviation
    double score_sigma(std::sqrt(std::exp(2 * std::log(sd)) + mt_variances));

    if ((diff_mz < mu + sigma_mult * score_sigma) && (diff_mz > mu - sigma_mult * score_sigma))
    {
      double tmp_exponent((diff_mz - mu) / score_sigma);
      mz_score = std::exp(-0.5 * tmp_exponent * tmp_exponent);
    }

    return mz_score;
  }

  // from FLASHDeconvAlgorithm::getCosine
  double FeatureFindingIntact::computeCosineSimOfDiffSizedVector_(const std::vector<double>& a,
                                         const IsotopeDistribution& b,
                                         const int& b_size,
                                         const double& b_norm,
                                         const int offset) const
  {
    double n = .0, a_norm= .0;
    for (Size j = 0; j < a.size(); j++)
    {
      a_norm += a[j] * a[j];
      int i = j - offset;
      if (i < 0 || i >= b_size)
      {
        continue;
      }
      n += a[j] * b[i].getIntensity(); //
    }
    double d = (a_norm * b_norm);
    if (d <= 0)
    {
      return 0;
    }
    return n / sqrt(d);
  }

  // modified based on FLASHDeconvAlgorithm::getIsotopeCosineAndDetermineIsotopeIndex
  // isotopic intensity distribution comparison using the cosine similarity
  double FeatureFindingIntact::computeAveragineCosineSimScore_(const std::vector<double>& hypo_ints,
                                                               const IsotopeDistribution& iso_dist,
                                                               const Size& iso_size,
                                                               const double& iso_norm,
                                                               int& offset) const
  {
    // Unlike FLASHDeconv, candidate should be subset of averagine.
    // iso_size is always larger than hypo_size
    offset = 0; // initialize
    double max_cosine = -1;
    int hypo_size = hypo_ints.size();

    // find largest intensity in hypo_ints for normalization
    double max_inty_of_hypo = 0;
    for (int i = 0; i < hypo_size; ++i)
    {
      if (hypo_ints[i] > max_inty_of_hypo)
      {
        max_inty_of_hypo = hypo_ints[i];
      }
    }

    // normalize hypo_ints
    std::vector<double> hypo_norm_ints;
    for (int i = 0; i < hypo_size; ++i)
    {
      hypo_norm_ints.push_back(hypo_ints[i] / max_inty_of_hypo);
    }

    // determining isotope index within averagine guaranteeing highest score
    int max_cntr = 0;
    vector<double> tmp_cos_vec;
    for (int tmp_offset = -(iso_size-hypo_size); tmp_offset <= 0; tmp_offset++)
    {
      double tmp_cos = computeCosineSimOfDiffSizedVector_(hypo_norm_ints,
                                                          iso_dist,
                                                          iso_size,
                                                          iso_norm,
                                                          tmp_offset);
      tmp_cos_vec.push_back(tmp_cos);

      if (max_cosine <= tmp_cos)
      {
        if (max_cosine == tmp_cos)
        {
          max_cntr++;
          offset += tmp_offset;
        }
        else
        {
          max_cosine = tmp_cos;
          max_cntr = 1;
          offset = tmp_offset;
        }
      }
    }

    offset /= max_cntr;
    return max_cosine;
  }

  double FeatureFindingIntact::computeCosineSim_(const std::vector<double>& x, const std::vector<double>& y) const
  {
    if (x.size() != y.size())
    {
      return 0.0;
    }

    double mixed_sum(0.0);
    double x_squared_sum(0.0);
    double y_squared_sum(0.0);

    for (Size i = 0; i < x.size(); ++i)
    {
      mixed_sum += x[i] * y[i];
      x_squared_sum += x[i] * x[i];
      y_squared_sum += y[i] * y[i];
    }

    double denom(std::sqrt(x_squared_sum) * std::sqrt(y_squared_sum));
    return (denom > 0.0) ? mixed_sum / denom : 0.0;
  }

  void FeatureFindingIntact::setAveragineModel()
  {
    auto generator = new CoarseIsotopePatternGenerator();
    auto maxIso = generator->estimateFromPeptideWeight(mass_upper_bound);
    maxIso.trimRight(0.01 * maxIso.getMostAbundant().getIntensity());

    generator->setMaxIsotope(maxIso.size());
    iso_model_ = PrecalculatedAveragine(50, mass_upper_bound, 25, generator);
    iso_model_.setMaxIsotopeIndex(maxIso.size() - 1);

    max_nr_traces_ = iso_model_.getMaxIsotopeIndex();
  }

  void FeatureFindingIntact::clusterFeatureHypotheses_(vector<FeatureHypothesis>& hypotheses,
                                                       const std::vector<std::vector<Size>>& shared_m_traces_indices) const
  {
    // *********************************************************** //
    // Step 1 constructing hypergraph from featurehypotheses
    //        node = mass traces
    //        hyperedge = hypotheses
    // *********************************************************** //
    Size num_nodes = shared_m_traces_indices.size();
    std::vector<bool> bfs_visited;
    bfs_visited.resize(num_nodes, false);
    std::queue<Size> bfs_queue;
    Size search_pos = 0; // keeping track of mass trace index to look for seed

    std::vector<FeatureHypothesis> out_features;

    // BFS
    // TODO : progress logger
//    this->startProgress(0, shared_m_traces_indices.size(), "clustering features based on the shared mass traces");
//    Size progress = 0;
    while (true)
    {
//      this->setProgress(progress);
//      ++progress;
      // finding a seed 'shared_mass_trace' to start with (for constructing a cluster)
      bool finished = true;
      for (Size i = search_pos; i < num_nodes; ++i)
      {
        if (!bfs_visited[i])
        {
          bfs_queue.push(i);
          bfs_visited[i] = true;
          finished = false;
          search_pos = i + 1;
          break;
        }
      }
      if (finished) // if no possible seed is left
        break;

      set<Size> hypo_indices_in_current_cluster;

      while (!bfs_queue.empty())
      {
        Size i = bfs_queue.front(); // index of seed
        bfs_queue.pop();

        // get feature indices sharing this seed
        for (vector<Size>::const_iterator it = shared_m_traces_indices[i].begin();
             it != shared_m_traces_indices[i].end();
             ++it)
        {
          hypo_indices_in_current_cluster.insert(*it);
          FeatureHypothesis &current_hypo = hypotheses[*it];

          for (const auto &iso_mt_pair : current_hypo.getIndicesOfMassTraces())
          {
            Size j = iso_mt_pair.second;
            if (!bfs_visited[j])
            {
              bfs_queue.push(j);
              bfs_visited[j] = true;
            }
          }
        }
      }

      // current cluster
      if (hypo_indices_in_current_cluster.empty()) continue; // no cluster out of current mt
      if (hypo_indices_in_current_cluster.size() == 1){
        // no conflict, but cannot happen.
        out_features.push_back(hypotheses[*(hypo_indices_in_current_cluster.begin())]);
        break;
      }
      resolveConflictInCluster_(hypotheses, shared_m_traces_indices, hypo_indices_in_current_cluster, out_features);
    }
//    this->endProgress();
  }

  void FeatureFindingIntact::resolveConflictInCluster_(const std::vector<FeatureHypothesis>& feat_hypo,
                                 const std::vector<std::vector<Size> >& shared_m_traces_indices,
                                 const std::set<Size>& hypo_indices,
                                 std::vector<FeatureHypothesis>& out_features) const
  {
    // compute scores for each features
    vector<double> feat_scores(hypo_indices.size());
    for (const auto& h_index : hypo_indices)
    {
      feat_scores.push_back( feat_hypo[h_index].getScore() );
    }
  }

}