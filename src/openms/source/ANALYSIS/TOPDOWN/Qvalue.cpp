// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong$
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/Qvalue.h>

namespace OpenMS
{
  void Qvalue::removeOutliers(std::vector<double> qscores, uint bin_number)
  {
    std::vector<double> tmp_qscores(qscores);
    std::map<uint, int> bin_cntr;
    int max_bin_cntr = 0;
    qscores.clear();
    for (const auto s : tmp_qscores)
    {
      uint bin = getBinNumber(s, bin_number);
      if (bin_cntr.find(bin) == bin_cntr.end())
        bin_cntr[bin] = 0;
      bin_cntr[bin]++;
      max_bin_cntr = std::max(max_bin_cntr, bin_cntr[bin]);
    }
    for (const auto s : tmp_qscores)
    {
      uint bin = getBinNumber(s, bin_number);
      if (bin_cntr[bin] < max_bin_cntr / 20)
        continue;
      qscores.push_back(s);
    }
  }

  void Qvalue::updatePeakGroupQvalues(std::vector<DeconvolvedSpectrum>& deconvolved_spectra) // per ms level + precursor update as well.
  {
    uint bin_number;
    const uint min_bin_number = 100;

    std::map<uint, std::vector<double>> weights_map;
    std::map<uint, std::vector<double>> tscore_map; // per ms level

    std::map<uint, std::vector<double>> dscore_iso_decoy_map;
    std::map<uint, std::vector<double>> dscore_noise_decoy_map;
    std::map<uint, std::vector<double>> dscore_charge_decoy_map;

    std::map<uint, std::map<double, double>> qscore_qvalue_map; //

    // to calculate qvalues per ms level, store Qscores per ms level
    std::set<uint> used_feature_indices;

    for (auto& deconvolved_spectrum : deconvolved_spectra)
    {
      if (deconvolved_spectrum.empty() || deconvolved_spectrum.isDecoy())
        continue;

      uint ms_level = deconvolved_spectrum.getOriginalSpectrum().getMSLevel();
      for (auto& pg : deconvolved_spectrum)
      {
        if (pg.getFeatureIndex() > 0 && used_feature_indices.find(pg.getFeatureIndex()) != used_feature_indices.end())
          continue;
        used_feature_indices.insert(pg.getFeatureIndex());
        tscore_map[ms_level].push_back(pg.getQscore2D());
      }
    }

    for (auto& decoy_deconvolved_spectrum : deconvolved_spectra)
    {
      if (decoy_deconvolved_spectrum.empty() || !decoy_deconvolved_spectrum.isDecoy())
        continue;

      uint ms_level = decoy_deconvolved_spectrum.getOriginalSpectrum().getMSLevel();
      for (auto& pg : decoy_deconvolved_spectrum)
      {
        if (pg.getFeatureIndex() > 0 && used_feature_indices.find(pg.getFeatureIndex()) != used_feature_indices.end())
          continue;
        used_feature_indices.insert(pg.getFeatureIndex());
        if (pg.getTargetDecoyType() == PeakGroup::TargetDecoyType::isotope_decoy)
        {
          dscore_iso_decoy_map[ms_level].push_back(pg.getQscore2D());
        }
        else if (pg.getTargetDecoyType() == PeakGroup::TargetDecoyType::noise_decoy)
        {
          dscore_noise_decoy_map[ms_level].push_back(pg.getQscore2D());
        }
        else if (pg.getTargetDecoyType() == PeakGroup::TargetDecoyType::charge_decoy)
        {
          dscore_charge_decoy_map[ms_level].push_back(pg.getQscore2D());
        }
      }
    }

    for (auto& [ms_level, qscores] : tscore_map)
    {
      auto& dscore_iso = dscore_iso_decoy_map[ms_level];
      auto& dscore_charge = dscore_charge_decoy_map[ms_level];
      auto& dscore_noise = dscore_noise_decoy_map[ms_level];
      bin_number = min_bin_number; // std::max(min_bin_number, (uint)(qscores.size()/50));

      //removeOutliers(dscore_charge, bin_number);
      //removeOutliers(dscore_noise, bin_number);

      auto mixed_dist = getDistribution(qscores, bin_number);
      const auto charge_dist = getDistribution(dscore_charge, bin_number);
      const auto noise_dist = getDistribution(dscore_noise, bin_number);
      const auto iso_dist = getDistribution(dscore_iso, bin_number);

      std::vector<double> true_positive_dist(bin_number);
      std::vector<std::vector<double>> comp_dists {};

      for (int i = 0; i < mixed_dist.size(); i++)
      {
        mixed_dist[i] -= iso_dist[i] / (dscore_iso.empty() ? .0 : (((double)(qscores.size())) / dscore_iso.size()));
        mixed_dist[i] = mixed_dist[i] < .0 ? .0 : mixed_dist[i];
      }

      comp_dists.clear();
      comp_dists.push_back(charge_dist);
      comp_dists.push_back(noise_dist);
      auto weights = getDistributionWeights(mixed_dist, comp_dists, bin_number / 3);

      weights[0] *= dscore_charge.empty() ? .0 : (((double)(qscores.size())) / dscore_charge.size());
      weights[1] *= dscore_noise.empty() ? .0 : (((double)(qscores.size())) / dscore_noise.size());
      weights_map[ms_level] = weights;

      std::sort(qscores.begin(), qscores.end());
      std::sort(dscore_iso.begin(), dscore_iso.end());
      std::sort(dscore_noise.begin(), dscore_noise.end());
      std::sort(dscore_charge.begin(), dscore_charge.end());

      auto& map_qvalue = qscore_qvalue_map[ms_level];

      // calculate q values using targets and decoys
      for (size_t i = 0; i < qscores.size(); i++)
      {
        double ts = qscores[i];
        if (map_qvalue.find(ts) != map_qvalue.end())
          continue;

        size_t tindex = qscores.size() - i;
        double nom = 0;
        size_t dindex = dscore_charge.size() == 0 ? 0 : std::distance(std::lower_bound(dscore_charge.begin(), dscore_charge.end(), ts), dscore_charge.end()); // very inefficient... fix later.
        nom += weights[0] * (double)dindex;
        dindex = dscore_noise.size() == 0 ? 0 : std::distance(std::lower_bound(dscore_noise.begin(), dscore_noise.end(), ts), dscore_noise.end());
        nom += weights[1] * (double)dindex;
        dindex = dscore_iso.size() == 0 ? 0 : std::distance(std::lower_bound(dscore_iso.begin(), dscore_iso.end(), ts), dscore_iso.end());
        nom += (double)dindex;

        double tmp_q = (nom / (double)tindex);
        map_qvalue[ts] = std::min(1.0, tmp_q);
      }
    }

    for (auto& titem : tscore_map)
    {
      uint ms_level = titem.first;
      auto& map_qvalue = qscore_qvalue_map[ms_level];

      double cummin = 1.0;
      {
        for (auto&& rit = map_qvalue.begin(); rit != map_qvalue.end(); ++rit)
        {
          cummin = std::min(rit->second, cummin);
          rit->second = cummin;
        }
      }

      for (auto& deconvolved_spectrum : deconvolved_spectra)
      {
        if (deconvolved_spectrum.empty() || deconvolved_spectrum.isDecoy())
          continue;

        if (deconvolved_spectrum.getOriginalSpectrum().getMSLevel() == ms_level + 1 && !deconvolved_spectrum.getPrecursorPeakGroup().empty())
        {
          auto precursor_pg = deconvolved_spectrum.getPrecursorPeakGroup();
          double qs = precursor_pg.getQscore2D();

          precursor_pg.setQvalue(map_qvalue[qs]);
          deconvolved_spectrum.setPrecursorPeakGroup(precursor_pg);
        }

        if (deconvolved_spectrum.getOriginalSpectrum().getMSLevel() != ms_level)
        {
          continue;
        }

        for (auto& pg : deconvolved_spectrum)
        {
          pg.setQvalue(map_qvalue[pg.getQscore2D()]);
        }
      }
    }
  }

  uint Qvalue::getBinNumber(double qscore, uint total_bin_number)
  {
    return (uint)round(pow(qscore, 1) * (total_bin_number - 1.0));
  }


  float Qvalue::getBinValue(uint bin_number, uint total_bin_number)
  {
    return (float)pow((double)(bin_number) / (total_bin_number - 1.0), 1.0);
  }

  std::vector<double> Qvalue::getDistribution(const std::vector<double>& qscores, uint bin_number)
  {
    std::vector<double> ret(bin_number, .0);

    for (double qscore : qscores)
    {
      if (qscore < 0 || qscore > 1.0f)
      {
        continue;
      }
      uint bin = getBinNumber(qscore, bin_number);
      ret[bin]++;
    }

    double csum = 0;

    for (uint i = 0; i < bin_number; i++)
    {
      csum += ret[i];
    }

    if (csum > 0)
    {
      for (auto& r : ret)
        r /= csum;
    }
    return ret;
  }

  std::vector<double> Qvalue::getDistributionWeights(
    const std::vector<double>& mixed_dist, const std::vector<std::vector<double>>& comp_dists, int bin_threshold,
    uint num_iterations) // Richardson-Lucy algorithm https://stats.stackexchange.com/questions/501288/estimating-weights-of-known-component-distributions-in-a-mixture-distribution
  {
    uint weight_cntr = comp_dists.size();
    std::vector<double> weights(weight_cntr, 1.0 / weight_cntr);
    std::vector<double> c_sums(weight_cntr, .0);
    for (uint k = 0; k < bin_threshold; k++) //
    {
      for (uint i = 0; i < weight_cntr; i++)
        c_sums[i] += comp_dists[i][k];
    }

    for (uint n = 0; n < num_iterations; n++)
    {
      std::vector<double> tmp_weights(weights);

      for (uint i = 0; i < weight_cntr; i++)
      {
        double t = .0;
        for (uint k = 0; k < bin_threshold; k++) //
        {
          double denom = .0;
          for (uint j = 0; j < weight_cntr; j++)
          {
            denom += weights[j] * comp_dists[j][k];
          }
          if (denom > 0)
          {
            t += comp_dists[i][k] * mixed_dist[k] / denom;
          }
        }
        tmp_weights[i] *= t;
      }

      if (weights == tmp_weights)
      {
        break;
      }
      weights = tmp_weights;
    }
    return weights;
  }
} // namespace OpenMS