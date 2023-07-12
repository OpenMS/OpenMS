//--------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Maintainer: Kyowon Jeong$
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/Qvalue.h>

namespace OpenMS
{
  void Qvalue::updatePeakGroupQvalues(std::vector<DeconvolvedSpectrum>& deconvolved_spectra,
                                      std::vector<DeconvolvedSpectrum>& deconvolved_decoy_spectra) // per ms level + precursor update as well.
  {
    uint bin_number = 100;                         // 100 is enough resolution for qvalue calculation. In most cases FDR 5% will be used.
    std::map<uint, std::vector<float>> tscore_map; // per ms level

    std::map<uint, std::vector<float>> dscore_iso_decoy_map;
    std::map<uint, std::map<float, float>> qscore_iso_decoy_map; // maps for isotope decoy only qvalues

    std::map<uint, std::vector<float>> dscore_noise_decoy_map;
    std::map<uint, std::map<float, float>> qscore_noise_decoy_map; // maps for noise decoy only qvalues

    std::map<uint, std::vector<float>> dscore_charge_decoy_map;
    std::map<uint, std::map<float, float>> qscore_charge_decoy_map; // maps for charge decoy only qvalues

    // to calculate qvalues per ms level, store Qscores per ms level
    for (auto& deconvolved_spectrum : deconvolved_spectra)
    {
      uint ms_level = deconvolved_spectrum.getOriginalSpectrum().getMSLevel();
      for (auto& pg : deconvolved_spectrum)
      {
        tscore_map[ms_level].push_back(pg.getQscore());
      }
    }

    for (auto& decoy_deconvolved_spectrum : deconvolved_decoy_spectra)
    {
      uint ms_level = decoy_deconvolved_spectrum.getOriginalSpectrum().getMSLevel();
      for (auto& pg : decoy_deconvolved_spectrum)
      {
        if (pg.getTargetDummyType() == PeakGroup::TargetDummyType::isotope_dummy)
        {
          dscore_iso_decoy_map[ms_level].push_back(pg.getQscore());
        }
        else if (pg.getTargetDummyType() == PeakGroup::TargetDummyType::noise_dummy)
        {
          dscore_noise_decoy_map[ms_level].push_back(pg.getQscore());
        }
        else if (pg.getTargetDummyType() == PeakGroup::TargetDummyType::charge_dummy)
        {
          dscore_charge_decoy_map[ms_level].push_back(pg.getQscore());
        }
      }
    }

    for (auto& [ms_level, qscores] : tscore_map)
    {
      auto& dscore_iso = dscore_iso_decoy_map[ms_level];
      auto& dscore_charge = dscore_charge_decoy_map[ms_level];
      auto& dscore_noise = dscore_noise_decoy_map[ms_level];

      auto mixed_dist = getDistribution(qscores, bin_number);
      auto charge_dist = getDistribution(dscore_charge, bin_number);
      auto noise_dist = getDistribution(dscore_noise, bin_number);
      auto iso_dist = getDistribution(dscore_iso, bin_number);
      float iso_threshold = 1.0f;

      for (int i = (int)bin_number - 1; i >= 0; i--)
      {
        if ((iso_dist[i] > 0 && iso_dist[i] < charge_dist[i]))
        {
          iso_threshold = getBinValue(i, bin_number);
          break;
        }
      }
      std::vector<float> tmp_dscore_iso(dscore_iso);
      dscore_iso.clear();
      for (auto q : tmp_dscore_iso)
      {
        if (q < iso_threshold)
          continue;
        dscore_iso.push_back(q);
      }

      iso_dist = getDistribution(dscore_iso, bin_number);

      std::vector<float> weights(4, .25f);
      std::vector<float> target_dist(bin_number, .0f);

      for (uint iteration = 0; iteration < 100; iteration++)
      {
        std::fill(target_dist.begin(), target_dist.end(), .0f);
        float max_t = .0;

        for (int i = (int)bin_number - 1; i >= 0; i--)
        {
          float fp = (charge_dist[i] * weights[0] + noise_dist[i] * weights[1] + iso_dist[i] * weights[2]);
          if (mixed_dist[i] > 0 && mixed_dist[i] < fp)
          {
            break;
          }
          target_dist[i] = mixed_dist[i] - fp;
          if (target_dist[i] > 0 && target_dist[i] < max_t * .1f)
          {
            break;
          }
          max_t = std::max(max_t, target_dist[i]);
        }

        float csum = .0f;

        for (auto r : target_dist)
          csum += r;

        if (csum > 0)
        {
          for (auto& r : target_dist)
            r /= csum;
        }

        std::vector<std::vector<float>> comp_dists {};
        comp_dists.push_back(charge_dist);
        comp_dists.push_back(noise_dist);
        comp_dists.push_back(iso_dist);
        comp_dists.push_back(target_dist);

        auto tmp_weight = getDistributionWeights(mixed_dist, comp_dists);

        if (tmp_weight == weights)
        {
          break;
        }
        weights = tmp_weight;
      }

      weights[0] *= dscore_charge.empty() ? .0 : ((double)qscores.size() / dscore_charge.size());
      weights[1] *= dscore_noise.empty() ? .0 : ((double)qscores.size() / dscore_noise.size());
      weights[2] *= dscore_iso.empty() ? .0 : ((double)qscores.size() / dscore_iso.size());

      std::sort(qscores.begin(), qscores.end());
      std::sort(dscore_iso.begin(), dscore_iso.end());
      std::sort(dscore_noise.begin(), dscore_noise.end());
      std::sort(dscore_charge.begin(), dscore_charge.end());

      auto& map_charge = qscore_charge_decoy_map[ms_level];
      float tmp_q_charge = 1;

      // calculate q values using targets and charge dummies
      for (size_t i = 0; i < qscores.size(); i++)
      {
        float ts = qscores[i];
        size_t dindex = dscore_charge.size() == 0 ? 0 : std::distance(std::lower_bound(dscore_charge.begin(), dscore_charge.end(), ts), dscore_charge.end());
        size_t tindex = qscores.size() - i;
        float nom = weights[0] * (float)dindex;
        float denom = (float)(tindex);
        tmp_q_charge = std::min(tmp_q_charge, (nom / denom));
        map_charge[ts] = tmp_q_charge;
      }

      auto& map_iso = qscore_iso_decoy_map[ms_level];
      float tmp_q_iso = 1;

      // calculate q values using targets and isotope dummies
      for (size_t i = 0; i < qscores.size(); i++)
      {
        float ts = qscores[i];
        size_t dindex = dscore_iso.size() == 0 ? 0 : std::distance(std::lower_bound(dscore_iso.begin(), dscore_iso.end(), ts), dscore_iso.end());
        size_t tindex = qscores.size() - i;
        float nom = weights[2] * (float)dindex;
        float denom = (float)(tindex);
        tmp_q_iso = std::min(tmp_q_iso, (nom / denom));
        map_iso[ts] = tmp_q_iso;
      }

      auto& map_noise = qscore_noise_decoy_map[ms_level];
      float tmp_q_noise = 1;

      // calculate q values using targets and noise dummies
      for (size_t i = 0; i < qscores.size(); i++)
      {
        float ts = qscores[i];
        size_t dindex = dscore_noise.size() == 0 ? 0 : std::distance(std::lower_bound(dscore_noise.begin(), dscore_noise.end(), ts), dscore_noise.end());
        size_t tindex = qscores.size() - i;
        float nom = weights[1] * (float)dindex;
        float denom = (float)(tindex);
        tmp_q_noise = std::min(tmp_q_noise, (nom / denom));
        map_noise[ts] = tmp_q_noise;
      }
    }

    for (auto& titem : tscore_map)
    {
      uint ms_level = titem.first;
      for (auto& deconvolved_spectrum : deconvolved_spectra)
      {
        if (deconvolved_spectrum.getOriginalSpectrum().getMSLevel() != ms_level)
        {
          continue;
        }
        auto& map_iso = qscore_iso_decoy_map[ms_level];
        auto& map_noise = qscore_noise_decoy_map[ms_level];
        auto& map_charge = qscore_charge_decoy_map[ms_level];

        for (auto& pg : deconvolved_spectrum)
        {
          pg.setQvalue(map_charge[pg.getQscore()], PeakGroup::TargetDummyType::charge_dummy);
          pg.setQvalue(map_noise[pg.getQscore()], PeakGroup::TargetDummyType::noise_dummy);
          pg.setQvalue(map_iso[pg.getQscore()], PeakGroup::TargetDummyType::isotope_dummy);

          if (deconvolved_spectrum.getOriginalSpectrum().getMSLevel() > 1 && !deconvolved_spectrum.getPrecursorPeakGroup().empty())
          {
            float qs = deconvolved_spectrum.getPrecursorPeakGroup().getQscore();
            auto& pmap_iso = qscore_iso_decoy_map[ms_level - 1];
            auto& pmap_charge = qscore_charge_decoy_map[ms_level - 1];
            auto& pmap_noise = qscore_noise_decoy_map[ms_level - 1];
            deconvolved_spectrum.getPrecursorPeakGroup().setQvalue(pmap_iso[qs], PeakGroup::TargetDummyType::isotope_dummy);
            deconvolved_spectrum.getPrecursorPeakGroup().setQvalue(pmap_noise[qs], PeakGroup::TargetDummyType::noise_dummy);
            deconvolved_spectrum.getPrecursorPeakGroup().setQvalue(pmap_charge[qs], PeakGroup::TargetDummyType::charge_dummy);
          }
        }
      }
    }
  }

  uint Qvalue::getBinNumber(float qscore, uint total_bin_number)
  {
    return (uint)(pow(qscore, 3.0) * (total_bin_number - 1.0) + .5f);
  }

  float Qvalue::getBinValue(uint bin_number, uint total_bin_number)
  {
    return pow((double)(bin_number) / (total_bin_number - 1.0), 1.0 / 3.0);
  }

  std::vector<float> Qvalue::getDistribution(const std::vector<float>& qscores, uint bin_number)
  {
    std::vector<float> ret(bin_number, .0f);

    for (float qscore : qscores)
    {
      if (qscore < 0.0f || qscore > 1.0f)
      {
        continue;
      }
      uint bin = getBinNumber(qscore, bin_number);
      ret[bin]++;
    }
    if (qscores.size() > 0)
    {
      for (uint i = 0; i < bin_number; i++)
      {
        ret[i] /= qscores.size();
      }
    }
    float csum = .0f;
    for (auto r : ret)
      csum += r;
    if (csum > 0)
    {
      for (auto& r : ret)
        r /= csum;
    }
    // ret = smoothByMovingAvg(ret);
    return ret;
  }

  std::vector<float> Qvalue::getDistributionWeights(const std::vector<float>& mixed_dist, const std::vector<std::vector<float>>& comp_dists, uint num_iterations)
  {
    uint weight_cntr = comp_dists.size();
    uint bin_number = mixed_dist.size();
    std::vector<float> weights(weight_cntr, 1.0f / weight_cntr);

    for (uint n = 0; n < num_iterations; n++)
    {
      std::vector<float> tmp_weights(weights);
      float tmp_weight_sum = .0f;

      for (uint i = 0; i < weight_cntr; i++)
      {
        float t = .0f;
        for (uint k = 0; k < bin_number; k++)
        {
          float denom = .0f;
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
        tmp_weights[i] = std::max(.0f, tmp_weights[i]);
        tmp_weight_sum += tmp_weights[i];
      }

      if (tmp_weight_sum > 0)
      {
        for (float& tmp_weight : tmp_weights)
        {
          tmp_weight /= tmp_weight_sum;
        }
      }
      if (weights == tmp_weights)
      {
        break;
      }
      weights.swap(tmp_weights);
    }
    return weights;
  }
} // namespace OpenMS